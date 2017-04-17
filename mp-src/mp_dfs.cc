// Copyright (c) 2016, Kazuki Yoshizoe
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
// may be used to endorse or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// AREDISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>

#include <boost/array.hpp>
#include <boost/random.hpp>

#include "gflags/gflags.h"

#include "mp_dfs.h"

#ifdef __CDT_PARSER__
#undef DBG
#define DBG(a)  a
#endif

#ifdef __CDT_PARSER__
#undef LOG
#define LOG(a)  ;
#endif

DECLARE_bool(lcm); // false, "item file is lcm style", mp-lamp.cc

DECLARE_string(item);// "", "filename of item set"
DECLARE_string(pos);// "", "filename of positive / negative file"

DECLARE_double(a);// significance level alpha

DECLARE_bool(show_progress);// "show progress at each iteration"

DECLARE_bool(second_phase);// true, "do second phase"
DECLARE_bool(third_phase);// true, "do third phase"

DECLARE_int32(stack_size);// 1024*1024*64, used as int[stack_size], lamp.cc
DEFINE_int32(give_size_max, 1024 * 1024 * 4, "maximum size of one give");
DECLARE_int32(freq_max); // 1024*1024*64, "stack size for holding freq sets", lamp.cc
DEFINE_int32(sig_max, 1024 * 1024 * 64,
		"stack size for holding significant sets");

DEFINE_int32(bsend_buffer_size, 1024 * 1024 * 64, "size of bsend buffer");

DEFINE_int32(d, 10, "debug level. 0: none, higher level produce more log");
DEFINE_string(debuglogfile, "d", "base filename for debug log");
DEFINE_bool(log, false, "show log");

DEFINE_int32(probe_period, 128, "probe period during process node");
DEFINE_bool(probe_period_is_ms, false,
		"true: probe period is milli sec, false: num loops");

namespace lamp_search {

const int MP_LAMP::k_echo_tree_branch = 3;

const int MP_LAMP::k_int_max = std::numeric_limits<int>::max();
// assuming (digits in long long int) > (bits of double mantissa)
const long long int MP_LAMP::k_cs_max = 1ll
		<< (std::numeric_limits<double>::digits - 1);

const int MP_LAMP::k_probe_period = 128;

MP_LAMP::MP_LAMP(int rank, int nu_proc, int n, bool n_is_ms, int w, int l,
		int m) :
		mpi_data_(FLAGS_bsend_buffer_size, rank, nu_proc, n, n_is_ms, w, l, m,
				k_echo_tree_branch), dtd_(k_echo_tree_branch), d_(
		NULL), g_(NULL), bsh_(NULL), timer_(Timer::GetInstance()), cs_thr_(
		NULL), dtd_accum_array_base_(
		NULL), accum_array_(NULL), dtd_accum_recv_base_(NULL), accum_recv_(
		NULL), give_stack_(NULL), waiting_(false), stealer_(
				mpi_data_.nRandStealTrials_, mpi_data_.hypercubeDimension_), phase_(
				0), sup_buf_(
		NULL), child_sup_buf_(NULL), freq_stack_(NULL), significant_stack_(
		NULL), total_expand_num_(0ll), expand_num_(0ll), closed_set_num_(0ll), final_closed_set_num_(
				0ll), final_support_(0), final_sig_level_(0.0), last_bcast_was_dtd_(
				false) {
	printf("initializing MP_LAMP\n");
	if (FLAGS_d > 0) {
		SetLogFileName();
		DBG(lfs_.open(log_file_name_.c_str(), std::ios::out)
		;);
	}

//	mpi_data_.bsend_buffer_ = new int[FLAGS_bsend_buffer_size];
//
//	int ret = MPI_Buffer_attach(mpi_data_.bsend_buffer_,
//			FLAGS_bsend_buffer_size * sizeof(int));
//	if (ret != MPI_SUCCESS) {
//		throw std::bad_alloc();
//	}
//	printf("bsend_buffer\n");

// note: should remove victim if in lifeline ???



	dtd_.Init();

	printf("dtd\n");

	InitTreeRequest();

	Init();

	DBG(D(1) << "lifelines:"
	;);
	DBG(for (int i = 0; i < mpi_data_.hypercubeDimension_; i++)
		D(1, false) << "\t" << mpi_data_.lifelines_[i]
		;);
	DBG(D(1, false) << std::endl
	;);

	DBG(D(1) << "thieves_ size=" << mpi_data_.thieves_->Size() << std::endl
	;);
	DBG(
			D(1) << "mpi_data_.lifeline_thieves_ size="
					<< mpi_data_.lifeline_thieves_->Size() << std::endl
			;);

	DBG(
			D(1) << "bcast_targets:" << std::endl; for (std::size_t i = 0; i < k_echo_tree_branch; i++) { D(1) << "\t[" << i << "]=" << mpi_data_.bcast_targets_[i] << std::endl; }
			// << "\t[0]=" << bcast_targets_[0]
			// << "\t[1]=" << bcast_targets_[1]
			// << "\t[2]=" << bcast_targets_[2]
			// D(1) << std::endl;
			);

	printf("initialized MP_LAMP\n");
}

void MP_LAMP::InitTreeRequest() {
	// pushing tree to lifeline for the 1st wave
	int num = 0;

	for (std::size_t i = 1; i <= k_echo_tree_branch; i++) {
		if (k_echo_tree_branch * mpi_data_.mpiRank_ + i
				< mpi_data_.nTotalProc_) {
			mpi_data_.lifeline_thieves_->Push(
					k_echo_tree_branch * mpi_data_.mpiRank_ + i);
			mpi_data_.bcast_targets_[num++] = k_echo_tree_branch
					* mpi_data_.mpiRank_ + i;
		}
	}
	// if (3*h_+1 < p_) {
	//   mpi_data_.lifeline_thieves_->Push(3*h_+1);
	//   bcast_targets_[num++] = 3*h_+1;
	// }
	// if (3*h_+2 < p_) {
	//   mpi_data_.lifeline_thieves_->Push(3*h_+2);
	//   bcast_targets_[num++] = 3*h_+2;
	// }
	// if (3*h_+3 < p_) {
	//   mpi_data_.lifeline_thieves_->Push(3*h_+3);
	//   bcast_targets_[num++] = 3*h_+3;
	// }
	if (mpi_data_.mpiRank_ > 0) {
		mpi_data_.lifelines_activated_[(mpi_data_.mpiRank_ - 1)
				/ k_echo_tree_branch] = true;
		mpi_data_.bcast_source_ = (mpi_data_.mpiRank_ - 1) / k_echo_tree_branch;
	}

	for (int i = 0; i < k_echo_tree_branch; ++i) {
		printf("bcast[%d] = %d\n", i, mpi_data_.bcast_targets_[i]);
	}
}

void MP_LAMP::SetTreeRequest() {
	for (std::size_t i = 1; i <= k_echo_tree_branch; i++) {
		if (k_echo_tree_branch * mpi_data_.mpiRank_ + i < mpi_data_.nTotalProc_)
			mpi_data_.lifeline_thieves_->Push(
					k_echo_tree_branch * mpi_data_.mpiRank_ + i);
	}
	// if (3*h_+1 < p_) mpi_data_.lifeline_thieves_->Push(3*h_+1);
	// if (3*h_+2 < p_) mpi_data_.lifeline_thieves_->Push(3*h_+2);
	// if (3*h_+3 < p_) mpi_data_.lifeline_thieves_->Push(3*h_+3);
	if (mpi_data_.mpiRank_ > 0)
		mpi_data_.lifelines_activated_[(mpi_data_.mpiRank_ - 1)
				/ k_echo_tree_branch] = true;
}

MP_LAMP::~MP_LAMP() {
	DBG(D(2) << "MP_LAMP destructor begin" << std::endl
	;);

//	if (pmin_thr_)
//		delete[] pmin_thr_;
	if (cs_thr_)
		delete[] cs_thr_;
	if (dtd_accum_array_base_)
		delete[] dtd_accum_array_base_;
	if (dtd_accum_recv_base_)
		delete[] dtd_accum_recv_base_;
	if (node_stack_)
		delete node_stack_;
	if (give_stack_)
		delete give_stack_;
	if (sup_buf_)
		bsh_->Delete(sup_buf_);
	if (child_sup_buf_)
		bsh_->Delete(child_sup_buf_);
	if (freq_stack_)
		delete freq_stack_;
	if (significant_stack_)
		delete significant_stack_;

	if (d_)
		delete d_;
	if (g_)
		delete g_;

	if (bsh_)
		delete bsh_;
	int size;

//	delete mpi_data_;
//	delete treesearch_data_;

	DBG(D(2) << "MP_LAMP destructor end" << std::endl
	;);
	if (FLAGS_d > 0) {
		DBG(lfs_.close()
		;);
	}

}

void MP_LAMP::InitDatabaseRoot(std::istream & is1, std::istream & is2) {
	uint64 * data = NULL;
	uint64 * positive = NULL;
	boost::array<int, 5> counters; // nu_bits, nu_trans, nu_items, max_item_in_transaction
	counters.assign(-1);

	int nu_trans = 0;
	int nu_items = 0;
	int nu_pos_total = 0;
	int max_item_in_transaction = 0;

	std::vector<std::string> * item_names = NULL;
	std::vector<std::string> * transaction_names = NULL;

	DatabaseReader<uint64> reader;

	assert(mpi_data_.mpiRank_ == 0);
	{
		item_names = new std::vector<std::string>;
		transaction_names = new std::vector<std::string>;

		if (FLAGS_lcm) {
			reader.ReadFilesLCM(&bsh_, is1, &data, &nu_trans, &nu_items, is2,
					&positive, &nu_pos_total, item_names,
					&max_item_in_transaction);
		} else {
			reader.ReadFiles(&bsh_, is1, &data, &nu_trans, &nu_items, is2,
					&positive, &nu_pos_total, item_names, transaction_names,
					&max_item_in_transaction);
		}

		counters[0] = (int) (bsh_->nu_bits);
		counters[1] = nu_trans;
		counters[2] = nu_items;
		counters[3] = nu_pos_total;
		counters[4] = max_item_in_transaction;

		CallBcast(&counters, 5, MPI_INT);
	}

	CallBcast(data, bsh_->NewArraySize(nu_items), MPI_UNSIGNED_LONG_LONG);
	CallBcast(positive, bsh_->NuBlocks(), MPI_UNSIGNED_LONG_LONG);

	long long int start_time = timer_->Elapsed();
	d_ = new Database<uint64>(bsh_, data, nu_trans, nu_items, positive,
			nu_pos_total, max_item_in_transaction, item_names,
			transaction_names);
	log_.d_.pval_table_time_ = timer_->Elapsed() - start_time;

	g_ = new LampGraph<uint64>(*d_);

	// Lambda max is initially the number of total transactions in the database.
	lambda_max_ = d_->MaxX(); // used for getMinSup
	// D() << "max_x=lambda_max=" << lambda_max_ << std::endl;
	// D() << "pos_total=" << d_->PosTotal() << std::endl;
	// d_->DumpItems(D(false));
	// d_->DumpPMinTable(D(false));

	double *pmin_thr_ = new double[lambda_max_ + 1]; // used for getMinSup
	cs_thr_ = new long long int[lambda_max_ + 1]; // used for getMinSup

	dtd_accum_array_base_ = new long long int[lambda_max_ + 4];
	dtd_accum_recv_base_ = new long long int[lambda_max_ + 4];
	accum_array_ = &(dtd_accum_array_base_[3]); // TODO: ???
	accum_recv_ = &(dtd_accum_recv_base_[3]);

	node_stack_ = new VariableLengthItemsetStack(FLAGS_stack_size);
	give_stack_ = new VariableLengthItemsetStack(FLAGS_give_size_max);
	freq_stack_ = new VariableLengthItemsetStack(FLAGS_freq_max);
	// node_stack_ =  new VariableLengthItemsetStack(FLAGS_stack_size, lambda_max_);
	// give_stack_ =  new VariableLengthItemsetStack(FLAGS_give_size_max, lambda_max_);
	// freq_stack_ = new VariableLengthItemsetStack(FLAGS_freq_max, lambda_max_);

	for (int i = 0; i <= lambda_max_; i++)
		pmin_thr_[i] = d_->PMin(i);

	cs_thr_[0] = 0ll; // should not be used ???
	for (int i = 1; i <= lambda_max_; i++) {
		cs_thr_[i] = (long long int) (std::min(
				std::floor(FLAGS_a / pmin_thr_[i - 1]), (double) (k_cs_max)));
		DBG(
				D(2) << "cs_thr[" << i << "]=\t" << cs_thr_[i] << "\tpmin_thr="
						<< pmin_thr_[i - 1] << std::endl
				;);
	}

	for (int i = 0; i <= lambda_max_; i++)
		accum_array_[i] = 0ll;

	sup_buf_ = bsh_->New();
	child_sup_buf_ = bsh_->New();

	delete pmin_thr_;
}

void MP_LAMP::InitDatabaseRoot(std::istream & is1, int posnum) {
	uint64 * data = NULL;
	uint64 * positive = NULL;
	boost::array<int, 5> counters; // nu_bits, nu_trans, nu_items, max_item_in_transaction
	counters.assign(-1);

	int nu_trans = 0;
	int nu_items = 0;
	int nu_pos_total = posnum;
	int max_item_in_transaction = 0;

	std::vector<std::string> * item_names = NULL;
	std::vector<std::string> * transaction_names = NULL;

	DatabaseReader<uint64> reader;

	assert(mpi_data_.mpiRank_ == 0);
	{
		item_names = new std::vector<std::string>;
		transaction_names = new std::vector<std::string>;

		if (FLAGS_lcm) {
			reader.ReadFilesLCM(&bsh_, is1, &data, &nu_trans, &nu_items,
					item_names, &max_item_in_transaction);
		} else {
			reader.ReadFiles(&bsh_, is1, &data, &nu_trans, &nu_items,
					item_names, transaction_names, &max_item_in_transaction);
		}

		counters[0] = (int) (bsh_->nu_bits);
		counters[1] = nu_trans;
		counters[2] = nu_items;
		counters[3] = nu_pos_total;
		counters[4] = max_item_in_transaction;

		CallBcast(&counters, 5, MPI_INT);
	}

	CallBcast(data, bsh_->NewArraySize(nu_items), MPI_UNSIGNED_LONG_LONG);

	d_ = new Database<uint64>(bsh_, data, nu_trans, nu_items, positive,
			nu_pos_total, max_item_in_transaction, item_names,
			transaction_names);
	g_ = new LampGraph<uint64>(*d_);

	lambda_max_ = d_->MaxX();
	// D() << "max_x=lambda_max=" << lambda_max_ << std::endl;
	// D() << "pos_total=" << d_->PosTotal() << std::endl;
	// d_->DumpItems(D(false));
	// d_->DumpPMinTable(D(false));

	double* pmin_thr_ = new double[lambda_max_ + 1];
	cs_thr_ = new long long int[lambda_max_ + 1];
	dtd_accum_array_base_ = new long long int[lambda_max_ + 4];
	dtd_accum_recv_base_ = new long long int[lambda_max_ + 4];
	accum_array_ = &(dtd_accum_array_base_[3]);
	accum_recv_ = &(dtd_accum_recv_base_[3]);

	node_stack_ = new VariableLengthItemsetStack(FLAGS_stack_size);
	give_stack_ = new VariableLengthItemsetStack(FLAGS_give_size_max);
	freq_stack_ = new VariableLengthItemsetStack(FLAGS_freq_max);
	// node_stack_ =  new VariableLengthItemsetStack(FLAGS_stack_size, lambda_max_);
	// give_stack_ =  new VariableLengthItemsetStack(FLAGS_give_size_max, lambda_max_);
	// freq_stack_ = new VariableLengthItemsetStack(FLAGS_freq_max, lambda_max_);

	for (int i = 0; i <= lambda_max_; i++)
		pmin_thr_[i] = d_->PMin(i);

	cs_thr_[0] = 0ll; // should not be used ???
	for (int i = 1; i <= lambda_max_; i++) {
		cs_thr_[i] = (long long int) (std::min(
				std::floor(FLAGS_a / pmin_thr_[i - 1]), (double) (k_cs_max)));
		DBG(
				D(2) << "cs_thr[" << i << "]=\t" << cs_thr_[i] << "\tpmin_thr="
						<< pmin_thr_[i - 1] << std::endl
				;);
	}

	for (int i = 0; i <= lambda_max_; i++)
		accum_array_[i] = 0ll;

	sup_buf_ = bsh_->New();
	child_sup_buf_ = bsh_->New();

	delete pmin_thr_;
}

void MP_LAMP::InitDatabaseSub(bool pos) {
	uint64 * data = NULL;
	uint64 * positive = NULL;
	boost::array<int, 5> counters; // nu_bits, nu_trans, nu_items, max_item_in_transaction
	counters.assign(-1);

	int nu_trans = 0;
	int nu_items = 0;
	int nu_pos_total = 0;
	int max_item_in_transaction = 0;

	assert(mpi_data_.mpiRank_ != 0);
	{
		CallBcast(&counters, 5, MPI_INT);

		int nu_bits = counters[0];
		nu_trans = counters[1];
		nu_items = counters[2];
		nu_pos_total = counters[3];
		max_item_in_transaction = counters[4];

		bsh_ = new VariableBitsetHelper<uint64>(nu_bits);
		data = bsh_->NewArray(nu_items);
		if (pos)
			positive = bsh_->New();
	}

	CallBcast(data, bsh_->NewArraySize(nu_items), MPI_UNSIGNED_LONG_LONG);
	if (pos)
		CallBcast(positive, bsh_->NuBlocks(), MPI_UNSIGNED_LONG_LONG);

	d_ = new Database<uint64>(bsh_, data, nu_trans, nu_items, positive,
			nu_pos_total, max_item_in_transaction,
			NULL, NULL);
	g_ = new LampGraph<uint64>(*d_);

	lambda_max_ = d_->MaxX();
// D() << "max_x=lambda_max=" << lambda_max_ << std::endl;
// D() << "pos_total=" << d_->PosTotal() << std::endl;
// d_->DumpItems(D(false));
// d_->DumpPMinTable(D(false));

	double* pmin_thr_ = new double[lambda_max_ + 1];
	cs_thr_ = new long long int[lambda_max_ + 1];
	dtd_accum_array_base_ = new long long int[lambda_max_ + 4];
	dtd_accum_recv_base_ = new long long int[lambda_max_ + 4];
	accum_array_ = &(dtd_accum_array_base_[3]);
	accum_recv_ = &(dtd_accum_recv_base_[3]);

	node_stack_ = new VariableLengthItemsetStack(FLAGS_stack_size);
	give_stack_ = new VariableLengthItemsetStack(FLAGS_give_size_max);
	freq_stack_ = new VariableLengthItemsetStack(FLAGS_freq_max);
// node_stack_ =  new VariableLengthItemsetStack(FLAGS_stack_size, lambda_max_);
// give_stack_ =  new VariableLengthItemsetStack(FLAGS_give_size_max, lambda_max_);
// freq_stack_ = new VariableLengthItemsetStack(FLAGS_freq_max, lambda_max_);

	for (int i = 0; i <= lambda_max_; i++)
		pmin_thr_[i] = d_->PMin(i);

	cs_thr_[0] = 0ll; // should not be used ???
	for (int i = 1; i <= lambda_max_; i++) {
		cs_thr_[i] = (long long int) (std::min(
				std::floor(FLAGS_a / pmin_thr_[i - 1]), (double) (k_cs_max)));
		DBG(
				D(2) << "cs_thr[" << i << "]=\t" << cs_thr_[i] << "\tpmin_thr="
						<< pmin_thr_[i - 1] << std::endl
				;);
	}

	for (int i = 0; i <= lambda_max_; i++)
		accum_array_[i] = 0ll;

	sup_buf_ = bsh_->New();
	child_sup_buf_ = bsh_->New();

	delete pmin_thr_;
}

//int MP_LAMP::ComputeZ(int p, int l) {
//	int z0 = 1;
//	int zz = l;
//	while (zz < p) {
//		z0++;
//		zz *= l;
//	}
//	return z0;
//}

void MP_LAMP::Init() {
	// todo: move accum cs count variable to inner class
	mpi_data_.echo_waiting_ = false;
	for (int i = 0; i < k_echo_tree_branch; i++)
		mpi_data_.accum_flag_[i] = false;

	dtd_.Init();

	log_.Init();

	stealer_.Init();
	waiting_ = false;

	last_bcast_was_dtd_ = false;
}

void MP_LAMP::CheckPoint() {
	DBG(D(4) << "CheckPoint Started" << std::endl
	;);
	assert(mpi_data_.echo_waiting_ == false);

	// does this hold?
	for (int i = 0; i < k_echo_tree_branch; i++)
		assert(mpi_data_.accum_flag_[i] == false);

	dtd_.CheckPoint();
	stealer_.Init();
	waiting_ = false;

	mpi_data_.thieves_->Clear();
	mpi_data_.lifeline_thieves_->Clear();
	for (int pi = 0; pi < mpi_data_.nTotalProc_; pi++)
		mpi_data_.lifelines_activated_[pi] = false;

	SetTreeRequest();
	DBG(D(4) << "CheckPoint Finished" << std::endl
	;);
}

void MP_LAMP::ClearTasks() {
	MPI_Status probe_status, recv_status;
	int data_count, src, tag;
	int flag;
	int error;

	while (true) {
		error = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag,
				&probe_status);
		if (error != MPI_SUCCESS) {
			DBG(
					D(1) << "error in MPI_Iprobe in ClearTasks: " << error
							<< std::endl
					;);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		if (!flag)
			break;

		log_.d_.cleared_tasks_++;

		error = MPI_Get_count(&probe_status, MPI_INT, &data_count);
		if (error != MPI_SUCCESS) {
			DBG(
					D(1) << "error in MPI_Get_count in ClearTasks: " << error
							<< std::endl
					;);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		src = probe_status.MPI_SOURCE;
		tag = probe_status.MPI_TAG;

		error = MPI_Recv(give_stack_->Stack(), data_count, MPI_INT, src, tag,
		MPI_COMM_WORLD, &recv_status);
		if (error != MPI_SUCCESS) {
			DBG(
					D(1) << "error in MPI_Recv in ClearTasks: " << error
							<< std::endl
					;);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		give_stack_->Clear(); // drop messages
	}
}

//==============================================================================

void MP_LAMP::Probe(MPI_Data& mpi_data) {
	DBG(D(3) << "Probe" << std::endl
	;);
	MPI_Status probe_status;
	int probe_src, probe_tag;

	long long int start_time;
	start_time = timer_->Elapsed();

	log_.d_.probe_num_++;

	while (CallIprobe(&probe_status, &probe_src, &probe_tag)) {
		DBG(
				D(4) << "CallIprobe returned src=" << probe_src << "\ttag="
						<< probe_tag << std::endl
				;);
		switch (probe_tag) {
		// control tasks
		case Tag::DTD_REQUEST:
			RecvDTDRequest(mpi_data, probe_src);
			break;
		case Tag::DTD_REPLY:
			RecvDTDReply(mpi_data, probe_src);
			break;

		case Tag::DTD_ACCUM_REQUEST:
			assert(phase_ == 1);
			RecvDTDAccumRequest(mpi_data, probe_src);
			break;
		case Tag::DTD_ACCUM_REPLY:
			assert(phase_ == 1);
			RecvDTDAccumReply(mpi_data, probe_src);
			break;

		case Tag::BCAST_FINISH:
			RecvBcastFinish(mpi_data, probe_src);
			break;

// basic tasks
		case Tag::LAMBDA:
			assert(phase_ == 1);
			RecvLambda(mpi_data, probe_src);
			break;

		case Tag::REQUEST:
			RecvRequest(mpi_data, probe_src);
			break;
		case Tag::REJECT:
			RecvReject(mpi_data, probe_src);
			break;
		case Tag::GIVE:
			RecvGive(mpi_data, probe_src, probe_status);
			break;

// third phase tasks
		case Tag::RESULT_REQUEST:
			RecvResultRequest(mpi_data, probe_src);
			break;
		case Tag::RESULT_REPLY:
			RecvResultReply(mpi_data, probe_src, probe_status);
			break;
		default:
			DBG(
					D(1) << "unknown Tag=" << probe_tag
							<< " received in Probe: " << std::endl
					;
			)
			;
			MPI_Abort(MPI_COMM_WORLD, 1);
			break;
		}
	}

	// capacity, lambda, phase
	if (phase_ == 1 || phase_ == 2) {
		assert(node_stack_);
		log_.TakePeriodicLog(node_stack_->NuItemset(), lambda_, phase_);
	}

	if (mpi_data.mpiRank_ == 0) {
		// initiate termination detection

		// note: for phase_ 1, accum request and dtd request are unified
		if (!mpi_data.echo_waiting_ && !dtd_.terminated_) {
			if (phase_ == 1) {
				SendDTDAccumRequest(mpi_data);
				// log_.d_.dtd_accum_request_num_++;
			} else if (phase_ == 2) {
				if (node_stack_->Empty()) {
					SendDTDRequest(mpi_data);
					// log_.d_.dtd_request_num_++;
				}
			} else {
				// unknown phase
				assert(0);
			}
		}
	}

	long long int elapsed_time = timer_->Elapsed() - start_time;
	log_.d_.probe_time_ += elapsed_time;
	log_.d_.probe_time_max_ = std::max(elapsed_time, log_.d_.probe_time_max_);
}

// call this from root rank to start DTD
void MP_LAMP::SendDTDRequest(MPI_Data& mpi_data) {
	int message[1];
	message[0] = 1; // dummy

	mpi_data.echo_waiting_ = true;

	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data.bcast_targets_[i] < 0)
			break;
		assert(
				mpi_data.bcast_targets_[i] < mpi_data.nTotalProc_
						&& "SendDTDRequest");
		CallBsend(message, 1, MPI_INT, mpi_data.bcast_targets_[i],
				Tag::DTD_REQUEST);
		DBG(
				D(3) << "SendDTDRequest: dst=" << mpi_data.bcast_targets_[i]
						<< "\ttimezone=" << dtd_.time_zone_ << std::endl
				;);
	}
}

void MP_LAMP::RecvDTDRequest(MPI_Data& mpi_data, int src) {
	DBG(
			D(3) << "RecvDTDRequest: src=" << src << "\ttimezone="
					<< dtd_.time_zone_ << std::endl
			;);
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::DTD_REQUEST, &recv_status);

	if (IsLeaf(mpi_data))
		SendDTDReply(mpi_data);
	else
		SendDTDRequest(mpi_data);
}

bool MP_LAMP::DTDReplyReady(MPI_Data& mpi_data) const {
	for (int i = 0; i < k_echo_tree_branch; i++)
		if (mpi_data.bcast_targets_[i] >= 0 && !(dtd_.accum_flag_[i]))
			return false;
	return true;
	// check only valid bcast_targets_
	// always true if leaf
}

void MP_LAMP::SendDTDReply(MPI_Data& mpi_data) {
	int message[3];
	// using reduced vars
	message[0] = dtd_.count_ + dtd_.reduce_count_;
	bool tw_flag = dtd_.time_warp_ || dtd_.reduce_time_warp_;
	message[1] = (tw_flag ? 1 : 0);
	// for Steal
	dtd_.not_empty_ = !(node_stack_->Empty()) || (mpi_data.thieves_->Size() > 0)
			|| stealer_.StealStarted() || mpi_data.processing_node_; // thieves_ and stealer state check
	// for Steal2
	// dtd_.not_empty_ =
	//     !(node_stack_->Empty()) || (thieves_->Size() > 0) ||
	//     waiting_ || mpi_data.processing_node_;
	bool em_flag = dtd_.not_empty_ || dtd_.reduce_not_empty_;
	message[2] = (em_flag ? 1 : 0);
	DBG(
			D(3) << "SendDTDReply: dst = " << mpi_data.bcast_source_
					<< "\tcount=" << message[0] << "\ttw=" << tw_flag << "\tem="
					<< em_flag << std::endl
			;);

	assert(mpi_data.bcast_source_ < mpi_data.nTotalProc_ && "SendDTDReply");
	CallBsend(message, 3, MPI_INT, mpi_data.bcast_source_, Tag::DTD_REPLY);

	dtd_.time_warp_ = false;
	dtd_.not_empty_ = false;
	dtd_.IncTimeZone();

	mpi_data.echo_waiting_ = false;
	dtd_.ClearAccumFlags();
	dtd_.ClearReduceVars();
}

void MP_LAMP::RecvDTDReply(MPI_Data& mpi_data, int src) {
	MPI_Status recv_status;
	int message[3];
	CallRecv(&message, 3, MPI_INT, src, Tag::DTD_REPLY, &recv_status);
	assert(src == recv_status.MPI_SOURCE);

	// reduce reply (count, time_warp, not_empty)
	dtd_.Reduce(message[0], (message[1] != 0), (message[2] != 0));

	DBG(
			D(3) << "RecvDTDReply: src=" << src << "\tcount=" << message[0]
					<< "\ttw=" << message[1] << "\tem=" << message[2]
					<< "\treduced_count=" << dtd_.reduce_count_
					<< "\treduced_tw=" << dtd_.reduce_time_warp_
					<< "\treduced_em=" << dtd_.reduce_not_empty_ << std::endl
			;);

	bool flag = false;
	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data.bcast_targets_[i] == src) {
			flag = true;
			dtd_.accum_flag_[i] = true;
			break;
		}
	}
	assert(flag);

	if (DTDReplyReady(mpi_data)) {
		if (mpi_data.mpiRank_ == 0) {
			DTDCheck(mpi_data); // at root
			log_.d_.dtd_phase_num_++;
		} else
			SendDTDReply(mpi_data);
	}
}

void MP_LAMP::SendDTDAccumRequest(MPI_Data& mpi_data) {
	int message[1];
	message[0] = 1; // dummy

	mpi_data.echo_waiting_ = true;

	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data.bcast_targets_[i] < 0)
			break;
		assert(
				mpi_data.bcast_targets_[i] < mpi_data.nTotalProc_
						&& "SendDTDAccumRequest");
		CallBsend(message, 1, MPI_INT, mpi_data.bcast_targets_[i],
				Tag::DTD_ACCUM_REQUEST);

		DBG(
				D(3) << "SendDTDAccumRequest: dst="
						<< mpi_data.bcast_targets_[i] << "\ttimezone="
						<< dtd_.time_zone_ << std::endl
				;);
	}
}

void MP_LAMP::RecvDTDAccumRequest(MPI_Data& mpi_data, int src) {
	DBG(
			D(3) << "RecvDTDAccumRequest: src=" << src << "\ttimezone="
					<< dtd_.time_zone_ << std::endl
			;);
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::DTD_ACCUM_REQUEST, &recv_status);

	if (IsLeaf(mpi_data))
		SendDTDAccumReply(mpi_data);
	else
		SendDTDAccumRequest(mpi_data);
}

void MP_LAMP::SendDTDAccumReply(MPI_Data& mpi_data) {
	dtd_accum_array_base_[0] = dtd_.count_ + dtd_.reduce_count_;
	bool tw_flag = dtd_.time_warp_ || dtd_.reduce_time_warp_;
	dtd_accum_array_base_[1] = (tw_flag ? 1 : 0);
	// for Steal
	dtd_.not_empty_ = !(node_stack_->Empty()) || (mpi_data.thieves_->Size() > 0)
			|| stealer_.StealStarted() || mpi_data.processing_node_; // thieves_ and stealer state check
	// for Steal2
	// dtd_.not_empty_ =
	//     !(node_stack_->Empty()) || (thieves_->Size() > 0) ||
	//     waiting_ || mpi_data.processing_node_;
	bool em_flag = dtd_.not_empty_ || dtd_.reduce_not_empty_;
	dtd_accum_array_base_[2] = (em_flag ? 1 : 0);

	DBG(
			D(3) << "SendDTDAccumReply: dst = " << mpi_data.bcast_source_
					<< "\tcount=" << dtd_accum_array_base_[0] << "\ttw="
					<< tw_flag << "\tem=" << em_flag << std::endl
			;);

	assert(
			mpi_data.bcast_source_ < mpi_data.nTotalProc_
					&& "SendDTDAccumReply");
	CallBsend(dtd_accum_array_base_, lambda_max_ + 4, MPI_LONG_LONG_INT,
			mpi_data.bcast_source_, Tag::DTD_ACCUM_REPLY);

	dtd_.time_warp_ = false;
	dtd_.not_empty_ = false;
	dtd_.IncTimeZone();

	mpi_data.echo_waiting_ = false;
	dtd_.ClearAccumFlags();
	dtd_.ClearReduceVars();

	for (int l = 0; l <= lambda_max_; l++)
		accum_array_[l] = 0ll;
}

void MP_LAMP::RecvDTDAccumReply(MPI_Data& mpi_data, int src) {
	MPI_Status recv_status;

	CallRecv(dtd_accum_recv_base_, lambda_max_ + 4, MPI_LONG_LONG_INT, src,
			Tag::DTD_ACCUM_REPLY, &recv_status);
	assert(src == recv_status.MPI_SOURCE);

	int count = (int) (dtd_accum_recv_base_[0]);
	bool time_warp = (dtd_accum_recv_base_[1] != 0);
	bool not_empty = (dtd_accum_recv_base_[2] != 0);

	dtd_.Reduce(count, time_warp, not_empty);

	DBG(
			D(3) << "RecvDTDAccumReply: src=" << src << "\tcount=" << count
					<< "\ttw=" << time_warp << "\tem=" << not_empty
					<< "\treduced_count=" << dtd_.reduce_count_
					<< "\treduced_tw=" << dtd_.reduce_time_warp_
					<< "\treduced_em=" << dtd_.reduce_not_empty_ << std::endl
			;);

	for (int l = lambda_ - 1; l <= lambda_max_; l++)
		accum_array_[l] += accum_recv_[l];

	bool flag = false;
	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data.bcast_targets_[i] == src) {
			flag = true;
			dtd_.accum_flag_[i] = true;
			break;
		}
	}
	assert(flag);

	if (mpi_data.mpiRank_ == 0) {
		if (ExceedCsThr()) {
			int new_lambda = NextLambdaThr();
			SendLambda(mpi_data, new_lambda);
			lambda_ = new_lambda;
		}
		// if SendLambda is called, dtd_.count_ is incremented and DTDCheck will always fail
		if (DTDReplyReady(mpi_data)) {
			DTDCheck(mpi_data);
			log_.d_.dtd_accum_phase_num_++;
		}
	} else {  // not root
		if (DTDReplyReady(mpi_data))
			SendDTDAccumReply(mpi_data);
	}
}

void MP_LAMP::DTDCheck(MPI_Data& mpi_data) {
	assert(mpi_data.mpiRank_ == 0);
	// (count, time_warp, not_empty)
	dtd_.Reduce(dtd_.count_, dtd_.time_warp_, dtd_.not_empty_);

	if (dtd_.reduce_count_ == 0 && dtd_.reduce_time_warp_ == false
			&& dtd_.reduce_not_empty_ == false) {
		// termination
		SendBcastFinish(mpi_data);
		dtd_.terminated_ = true;
		DBG(D(1) << "terminated" << std::endl
		;);
	}
	// doing same thing as SendDTDReply
	dtd_.time_warp_ = false;
	dtd_.not_empty_ = false;
	dtd_.IncTimeZone();

	mpi_data.echo_waiting_ = false;
	dtd_.ClearAccumFlags();
	dtd_.ClearReduceVars();
}

void MP_LAMP::SendBcastFinish(MPI_Data& mpi_data) {
	DBG(D(2) << "SendBcastFinish" << std::endl
	;);
	int message[1];
	message[0] = 1; // dummy

	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data.bcast_targets_[i] < 0)
			break;
		assert(
				mpi_data.bcast_targets_[i] < mpi_data.nTotalProc_
						&& "SendBcastFinish");
		CallBsend(message, 1, MPI_INT, mpi_data.bcast_targets_[i],
				Tag::BCAST_FINISH);
	}
}

void MP_LAMP::RecvBcastFinish(MPI_Data& mpi_data, int src) {
	DBG(D(2) << "RecvBcastFinish: src=" << src << std::endl
	;);
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::BCAST_FINISH, &recv_status);

	SendBcastFinish(mpi_data);
	dtd_.terminated_ = true;
	DBG(D(2) << "terminated" << std::endl
	;);

	waiting_ = false;
}

//==============================================================================

void MP_LAMP::Search() {
	TreeSearchData* treesearch_data_ = new TreeSearchData(node_stack_);
	GetMinSupData* getminsup_data_;
	GetTestableData* gettestable_data_;
	GetSignificantData* getsignificant_data_;

	log_.d_.search_start_time_ = timer_->Elapsed();
	total_expand_num_ = 0ll;

	// --------
	// preprocess
	CheckPoint();

	expand_num_ = 0ll;
	closed_set_num_ = 0ll;

	lambda_ = 1;

	{
		// push root state to stack
		int * root_itemset;
		node_stack_->PushPre();
		root_itemset = node_stack_->Top();
		node_stack_->SetSup(root_itemset, lambda_max_);
		node_stack_->PushPostNoSort();
	}

	PreProcessRootNode();
	{
		if (mpi_data_.mpiRank_ == 0 && FLAGS_show_progress) {
			std::cout << "# " << "preprocess end\n";
			std::cout << "# " << "lambda=" << lambda_ << "\tcs_thr[lambda]="
					<< std::setw(16) << cs_thr_[lambda_];
			std::cout << "\tpmin_thr[lambda-1]=" << std::setw(12)
					<< d_->PMin(lambda_ - 1);
			std::cout << "\tnum_expand=" << std::setw(12) << expand_num_

			<< "\telapsed_time="
					<< (timer_->Elapsed() - log_.d_.search_start_time_) / GIGA
					<< std::endl;
		}
		DBG(D(2) << "---------------" << std::endl
		;);
		DBG(D(1) << "preprocess phase end" << std::endl
		;);
		DBG(D(2) << "---------------" << std::endl
		;);
		DBG(
				D(2) << "# " << "preprocess end\n" << "# " << "lambda="
						<< lambda_ << "\tcs_thr[lambda]=" << std::setw(16)
						<< cs_thr_[lambda_] << "\tpmin_thr[lambda-1]="
						<< std::setw(12) << d_->PMin(lambda_ - 1)
						<< "\tnum_expand=" << std::setw(12) << expand_num_
						<< "\telapsed_time="
						<< (timer_->Elapsed() - log_.d_.search_start_time_)
								/ GIGA << std::endl
				;);
	}

	// --------
	// prepare phase 1
	phase_ = 1;
	log_.StartPeriodicLog();

	{
		if (mpi_data_.mpiRank_ == 0 && FLAGS_show_progress) {
			std::cout << "# " << "1st phase start\n";
			std::cout << "# " << "lambda=" << lambda_
					<< "\tclosed_set_num[n>=lambda]=" << std::setw(12)
					<< accum_array_[lambda_] << "\tcs_thr[lambda]="
					<< std::setw(16) << cs_thr_[lambda_];
			std::cout << "\tpmin_thr[lambda-1]=" << std::setw(12)
					<< d_->PMin(lambda_ - 1);
			std::cout << "\tnum_expand=" << std::setw(12) << expand_num_

			<< "\telapsed_time="
					<< (timer_->Elapsed() - log_.d_.search_start_time_) / GIGA
					<< std::endl;
		}
		DBG(D(2) << "---------------" << std::endl
		;);
		DBG(D(1) << "1st phase start" << std::endl
		;);
		DBG(D(2) << "---------------" << std::endl
		;);
		DBG(
				D(2) << "lambda=" << lambda_ << "\tclosed_set_num[n>=lambda]="
						<< std::setw(12) << accum_array_[lambda_]
						<< "\tcs_thr[lambda]=" << std::setw(16)
						<< cs_thr_[lambda_] << "\tpmin_thr[lambda-1]="
						<< std::setw(12) << d_->PMin(lambda_ - 1)
						<< "\tnum_expand=" << std::setw(12) << expand_num_
						<< "\telapsed_time="
						<< (timer_->Elapsed() - log_.d_.search_start_time_)
								/ GIGA << std::endl
				;);
		getminsup_data_ = new GetMinSupData(lambda_max_, lambda_, cs_thr_);
//		getminsup_data_ = new GetMinSupData(lambda_max_, lambda_, cs_thr_,
//				dtd_accum_array_base_, accum_array_, dtd_accum_recv_base_,
//				accum_recv_); // The name Phase1 is already so nonsense...
		GetMinimalSupport(mpi_data_, treesearch_data_, getminsup_data_);
		lambda_max_ = getminsup_data_->lambda_max_;
		lambda_ = getminsup_data_->lambda_;
		cs_thr_ = getminsup_data_->cs_thr_; // no needed? trickey.

		// todo: reduce expand_num_
		if (mpi_data_.mpiRank_ == 0 && FLAGS_show_progress) {
			std::cout << "# " << "1st phase end\n";
			std::cout << "# " << "lambda=" << lambda_;
			std::cout << "\tnum_expand=" << std::setw(12) << expand_num_

			<< "\telapsed_time="
					<< (timer_->Elapsed() - log_.d_.search_start_time_) / GIGA
					<< std::endl;
		}
		DBG(D(2) << "---------------" << std::endl
		;);
		DBG(D(1) << "1st phase end" << std::endl
		;);
		DBG(D(2) << "---------------" << std::endl
		;);
		DBG(
				D(2) << "lambda=" << lambda_ << "\tnum_expand=" << std::setw(12)
						<< expand_num_ << "\telapsed_time="
						<< (timer_->Elapsed() - log_.d_.search_start_time_)
								/ GIGA << std::endl
				;);
	}

	log_.d_.dtd_accum_phase_per_sec_ = (double) (log_.d_.dtd_accum_phase_num_)
			/ ((timer_->Elapsed() - log_.d_.search_start_time_) / GIGA);

	// --------
	// prepare phase 2
	phase_ = 2;
	CheckPoint();

	lambda_--;
	CallBcast(&lambda_, 1, MPI_INT);
	final_support_ = lambda_;

	if (!FLAGS_second_phase) {
		log_.d_.search_finish_time_ = timer_->Elapsed();
		log_.GatherLog(mpi_data_.nTotalProc_);
		DBG(D(1) << "log" << std::endl
		;);
		DBG(PrintLog(D(1, false))
		;);
		DBG(PrintPLog(D(1, false))
		;);
		return;
	}

	expand_num_ = 0ll;
	closed_set_num_ = 0ll;

	{
		// push root state to stack
		int * root_itemset;
		node_stack_->Clear();
		node_stack_->PushPre();
		root_itemset = node_stack_->Top();
		node_stack_->SetSup(root_itemset, lambda_max_);
		node_stack_->PushPostNoSort();
	}

	double int_sig_lev = 0.0;
	if (mpi_data_.mpiRank_ == 0)
		int_sig_lev = GetInterimSigLevel(lambda_);
	CallBcast(&int_sig_lev, 1, MPI_DOUBLE);

	// todo: reduce expand_num_

	{
		if (mpi_data_.mpiRank_ == 0 && FLAGS_show_progress) {
			std::cout << "# " << "2nd phase start\n";
			std::cout << "# " << "lambda=" << lambda_ << "\tint_sig_lev="
					<< int_sig_lev << "\telapsed_time="
					<< (timer_->Elapsed() - log_.d_.search_start_time_) / GIGA
					<< std::endl;
		}
		DBG(D(2) << "---------------" << std::endl
		;);
		DBG(D(1) << "2nd phase start" << std::endl
		;);
		DBG(D(2) << "---------------" << std::endl
		;);
		DBG(
				D(2) << "lambda=" << lambda_ << "\tint_sig_lev=" << int_sig_lev
						<< "\telapsed_time="
						<< (timer_->Elapsed() - log_.d_.search_start_time_)
								/ GIGA << std::endl
				;);

		sig_level_ = int_sig_lev;
		gettestable_data_ = new GetTestableData(freq_stack_, freq_map_,
				sig_level_);
		GetTestablePatterns(mpi_data_, treesearch_data_, gettestable_data_);
		freq_stack_ = gettestable_data_->freq_stack_;
		freq_map_ = gettestable_data_->freq_map_;
		sig_level_ = gettestable_data_->sig_level_;
	}

	DBG(D(1) << "closed_set_num=" << closed_set_num_ << std::endl
	;);

	long long int closed_set_num_reduced;
	MPI_Reduce(&closed_set_num_, &closed_set_num_reduced, 1, MPI_LONG_LONG_INT,
	MPI_SUM, 0, MPI_COMM_WORLD);

	DBG(if (mpi_data_.mpiRank_ == 0)
		D(1) << "closed_set_num_reduced=" << closed_set_num_reduced << std::endl
		;);
	if (mpi_data_.mpiRank_ == 0)
		final_closed_set_num_ = closed_set_num_reduced;

	log_.d_.dtd_phase_per_sec_ = (double) (log_.d_.dtd_phase_num_)
			/ ((timer_->Elapsed() - log_.d_.search_start_time_) / GIGA);

	MPI_Barrier( MPI_COMM_WORLD);
	log_.FinishPeriodicLog();

	{
		if (mpi_data_.mpiRank_ == 0 && FLAGS_show_progress) {
			std::cout << "# " << "2nd phase end\n";
			std::cout << "# " << "closed_set_num=" << std::setw(12)
					<< final_closed_set_num_ << "\tsig_lev="
					<< (FLAGS_a / final_closed_set_num_) << "\tnum_expand="
					<< std::setw(12) << expand_num_ << "\telapsed_time="
					<< (timer_->Elapsed() - log_.d_.search_start_time_) / GIGA
					<< std::endl;
		}
	}

	if (!FLAGS_third_phase) {
		log_.d_.search_finish_time_ = timer_->Elapsed();
		log_.GatherLog(mpi_data_.nTotalProc_);
		DBG(D(1) << "log" << std::endl
		;);
		DBG(PrintLog(D(1, false))
		;);
		DBG(PrintPLog(D(1, false))
		;);
		return;
	}

	// prepare 3rd phase
	phase_ = 3;
	//ClearTasks();
	CheckPoint(); // needed for reseting dtd_.terminated_

	if (node_stack_)
		delete node_stack_;
	node_stack_ = NULL;
	significant_stack_ = new VariableLengthItemsetStack(FLAGS_sig_max);
	// significant_stack_ = new VariableLengthItemsetStack(FLAGS_sig_max, lambda_max_);

	final_sig_level_ = FLAGS_a / final_closed_set_num_;
	CallBcast(&final_sig_level_, 1, MPI_DOUBLE);

	{
		getsignificant_data_ = new GetSignificantData();
		GetSignificantPatterns(mpi_data_, getsignificant_data_);
		// TODO: put back to global variables.
	}

	// copy only significant itemset to buffer
	// collect itemset
	//   can reuse the other stack (needs to compute pval again)
	//   or prepare simpler data structure
	if (mpi_data_.mpiRank_ == 0)
		SortSignificantSets();
	log_.d_.search_finish_time_ = timer_->Elapsed();
	MPI_Barrier( MPI_COMM_WORLD);

	{
		if (mpi_data_.mpiRank_ == 0 && FLAGS_show_progress) {
			std::cout << "# " << "3rd phase end\n";
			std::cout << "# " << "sig_lev=" << final_sig_level_
					<< "\telapsed_time="
					<< (log_.d_.search_finish_time_ - log_.d_.search_start_time_)
							/ GIGA << std::endl;
		}
	}

	log_.GatherLog(mpi_data_.nTotalProc_);
	DBG(D(1) << "log" << std::endl
	;);
	DBG(PrintLog(D(1, false))
	;);
	DBG(PrintPLog(D(1, false))
	;);

	//ClearTasks();
}

//void MP_LAMP::SearchStraw1() {
//	log_.d_.search_start_time_ = timer_->Elapsed();
//	total_expand_num_ = 0ll;
//
//	// --------
//	// preprocess
//	CheckPoint();
//
//	expand_num_ = 0ll;
//	closed_set_num_ = 0ll;
//
//	lambda_ = 1;
//
//	{
//		// push root state to stack
//		int * root_itemset;
//		node_stack_->PushPre();
//		root_itemset = node_stack_->Top();
//		node_stack_->SetSup(root_itemset, lambda_max_);
//		node_stack_->PushPostNoSort();
//	}
//
//	PreProcessRootNode();
//	{
//		if (mpi_data_.mpiRank_ == 0 && FLAGS_show_progress) {
//			std::cout << "# " << "preprocess end\n";
//			std::cout << "# " << "lambda=" << lambda_ << "\tcs_thr[lambda]="
//					<< std::setw(16) << cs_thr_[lambda_];
//			std::cout << "\tpmin_thr[lambda-1]=" << std::setw(12)
//					<< d_->PMin(lambda_ - 1);
//			std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
//
//			<< "\telapsed_time="
//					<< (timer_->Elapsed() - log_.d_.search_start_time_) / GIGA
//					<< std::endl;
//		}
//		DBG(D(2) << "---------------" << std::endl
//		;);
//		DBG(D(1) << "preprocess phase end" << std::endl
//		;);
//		DBG(D(2) << "---------------" << std::endl
//		;);
//		DBG(
//				D(2) << "# " << "preprocess end\n" << "# " << "lambda="
//						<< lambda_ << "\tcs_thr[lambda]=" << std::setw(16)
//						<< cs_thr_[lambda_] << "\tpmin_thr[lambda-1]="
//						<< std::setw(12) << d_->PMin(lambda_ - 1)
//						<< "\tnum_expand=" << std::setw(12) << expand_num_
//						<< "\telapsed_time="
//						<< (timer_->Elapsed() - log_.d_.search_start_time_)
//								/ GIGA << std::endl
//				;);
//	}
//
//	// --------
//	// prepare phase 1
//	phase_ = 1;
//	log_.StartPeriodicLog();
//
//	// for strawman
//	stealer_.Finish();
//	waiting_ = false;
//
//	{
//		if (mpi_data_.mpiRank_ == 0 && FLAGS_show_progress) {
//			std::cout << "# " << "1st phase start\n";
//			std::cout << "# " << "lambda=" << lambda_
//					<< "\tclosed_set_num[n>=lambda]=" << std::setw(12)
//					<< accum_array_[lambda_] << "\tcs_thr[lambda]="
//					<< std::setw(16) << cs_thr_[lambda_];
//			std::cout << "\tpmin_thr[lambda-1]=" << std::setw(12)
//					<< d_->PMin(lambda_ - 1);
//			std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
//
//			<< "\telapsed_time="
//					<< (timer_->Elapsed() - log_.d_.search_start_time_) / GIGA
//					<< std::endl;
//		}
//		DBG(D(2) << "---------------" << std::endl
//		;);
//		DBG(D(1) << "1st phase start" << std::endl
//		;);
//		DBG(D(2) << "---------------" << std::endl
//		;);
//		DBG(
//				D(2) << "lambda=" << lambda_ << "\tclosed_set_num[n>=lambda]="
//						<< std::setw(12) << accum_array_[lambda_]
//						<< "\tcs_thr[lambda]=" << std::setw(16)
//						<< cs_thr_[lambda_] << "\tpmin_thr[lambda-1]="
//						<< std::setw(12) << d_->PMin(lambda_ - 1)
//						<< "\tnum_expand=" << std::setw(12) << expand_num_
//						<< "\telapsed_time="
//						<< (timer_->Elapsed() - log_.d_.search_start_time_)
//								/ GIGA << std::endl
//				;);
//
//		MainLoopStraw1();
//
//		// todo: reduce expand_num_
//		if (mpi_data_.mpiRank_ == 0 && FLAGS_show_progress) {
//			std::cout << "# " << "1st phase end\n";
//			std::cout << "# " << "lambda=" << lambda_;
//			std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
//
//			<< "\telapsed_time="
//					<< (timer_->Elapsed() - log_.d_.search_start_time_) / GIGA
//					<< std::endl;
//		}
//		DBG(D(2) << "---------------" << std::endl
//		;);
//		DBG(D(1) << "1st phase end" << std::endl
//		;);
//		DBG(D(2) << "---------------" << std::endl
//		;);
//		DBG(
//				D(2) << "lambda=" << lambda_ << "\tnum_expand=" << std::setw(12)
//						<< expand_num_ << "\telapsed_time="
//						<< (timer_->Elapsed() - log_.d_.search_start_time_)
//								/ GIGA << std::endl
//				;);
//	}
//
//	log_.d_.dtd_accum_phase_per_sec_ = (double) (log_.d_.dtd_accum_phase_num_)
//			/ ((timer_->Elapsed() - log_.d_.search_start_time_) / GIGA);
//
//	// log_.d_.accum_phase_per_sec_ =
//	//     log_.d_.accum_phase_num_ / ((timer_->Elapsed() - log_.d_.search_start_time_) / GIGA);
//
//	// --------
//	// prepare phase 2
//	phase_ = 2;
//	//ClearTasks();
//	CheckPoint();
//
//	// for strawman
//	stealer_.Finish();
//	waiting_ = false;
//
//	lambda_--;
//	CallBcast(&lambda_, 1, MPI_INT);
//	final_support_ = lambda_;
//
//	if (!FLAGS_second_phase) {
//		log_.d_.search_finish_time_ = timer_->Elapsed();
//		log_.GatherLog(mpi_data_.nTotalProc_);
//		DBG(D(1) << "log" << std::endl
//		;);
//		DBG(PrintLog(D(1, false))
//		;);
//		DBG(PrintPLog(D(1, false))
//		;);
//		return;
//	}
//
//	expand_num_ = 0ll;
//	closed_set_num_ = 0ll;
//
//	//if (h_ == 0) // pushing to all rank
//	{
//		// push root state to stack
//		int * root_itemset;
//		node_stack_->Clear();
//		node_stack_->PushPre();
//		root_itemset = node_stack_->Top();
//		node_stack_->SetSup(root_itemset, lambda_max_);
//		node_stack_->PushPostNoSort();
//	}
//
//	double int_sig_lev = 0.0;
//	if (mpi_data_.mpiRank_ == 0)
//		int_sig_lev = GetInterimSigLevel(lambda_);
//	// bcast int_sig_lev
//	CallBcast(&int_sig_lev, 1, MPI_DOUBLE);
//
//	// todo: reduce expand_num_
//
//	{
//		if (mpi_data_.mpiRank_ == 0 && FLAGS_show_progress) {
//			std::cout << "# " << "2nd phase start\n";
//			std::cout << "# " << "lambda=" << lambda_ << "\tint_sig_lev="
//					<< int_sig_lev << "\telapsed_time="
//					<< (timer_->Elapsed() - log_.d_.search_start_time_) / GIGA
//					<< std::endl;
//		}
//		DBG(D(2) << "---------------" << std::endl
//		;);
//		DBG(D(1) << "2nd phase start" << std::endl
//		;);
//		DBG(D(2) << "---------------" << std::endl
//		;);
//		DBG(
//				D(2) << "lambda=" << lambda_ << "\tint_sig_lev=" << int_sig_lev
//						<< "\telapsed_time="
//						<< (timer_->Elapsed() - log_.d_.search_start_time_)
//								/ GIGA << std::endl
//				;);
//
//		sig_level_ = int_sig_lev;
//		MainLoopStraw1();
//	}
//
//	DBG(D(1) << "closed_set_num=" << closed_set_num_ << std::endl
//	;);
//
//	long long int closed_set_num_reduced;
//	MPI_Reduce(&closed_set_num_, &closed_set_num_reduced, 1, MPI_LONG_LONG_INT,
//	MPI_SUM, 0, MPI_COMM_WORLD);
//
//	DBG(if (mpi_data_.mpiRank_ == 0)
//		D(1) << "closed_set_num_reduced=" << closed_set_num_reduced << std::endl
//		;);
//	if (mpi_data_.mpiRank_ == 0)
//		final_closed_set_num_ = closed_set_num_reduced;
//
//	log_.d_.dtd_phase_per_sec_ = (double) (log_.d_.dtd_phase_num_)
//			/ ((timer_->Elapsed() - log_.d_.search_start_time_) / GIGA);
//
//	MPI_Barrier( MPI_COMM_WORLD);
//	log_.FinishPeriodicLog();
//
//	// todo: gather logs
//
//	{
//		if (mpi_data_.mpiRank_ == 0 && FLAGS_show_progress) {
//			std::cout << "# " << "2nd phase end\n";
//			std::cout << "# " << "closed_set_num=" << std::setw(12)
//					<< final_closed_set_num_ << "\tsig_lev="
//					<< (FLAGS_a / final_closed_set_num_) << "\tnum_expand="
//					<< std::setw(12) << expand_num_ << "\telapsed_time="
//					<< (timer_->Elapsed() - log_.d_.search_start_time_) / GIGA
//					<< std::endl;
//		}
//	}
//
//	if (!FLAGS_third_phase) {
//		log_.d_.search_finish_time_ = timer_->Elapsed();
//		log_.GatherLog(mpi_data_.nTotalProc_);
//		DBG(D(1) << "log" << std::endl
//		;);
//		DBG(PrintLog(D(1, false))
//		;);
//		DBG(PrintPLog(D(1, false))
//		;);
//		return;
//	}
//
//	// prepare 3rd phase
//	phase_ = 3;
//	//ClearTasks();
//	CheckPoint(); // needed for reseting dtd_.terminated_
//
//	if (node_stack_)
//		delete node_stack_;
//	node_stack_ = NULL;
//	significant_stack_ = new VariableLengthItemsetStack(FLAGS_sig_max);
//	// significant_stack_ = new VariableLengthItemsetStack(FLAGS_sig_max, lambda_max_);
//
//	final_sig_level_ = FLAGS_a / final_closed_set_num_;
//	CallBcast(&final_sig_level_, 1, MPI_DOUBLE);
//
//	{
//		MainLoopStraw1();
//	}
//
//	// implement this in main loop
//	// copy only significant itemset to buffer
//	// collect itemset
//	//   can reuse the other stack (needs to compute pval again)
//	//   or prepare simpler data structure
//	if (mpi_data_.mpiRank_ == 0)
//		SortSignificantSets();
//	log_.d_.search_finish_time_ = timer_->Elapsed();
//	MPI_Barrier( MPI_COMM_WORLD);
//
//	{
//		if (mpi_data_.mpiRank_ == 0 && FLAGS_show_progress) {
//			std::cout << "# " << "3rd phase end\n";
//			std::cout << "# " << "sig_lev=" << final_sig_level_
//					<< "\telapsed_time="
//					<< (log_.d_.search_finish_time_ - log_.d_.search_start_time_)
//							/ GIGA << std::endl;
//		}
//	}
//
//	log_.GatherLog(mpi_data_.nTotalProc_);
//	DBG(D(1) << "log" << std::endl
//	;);
//	DBG(PrintLog(D(1, false))
//	;);
//	DBG(PrintPLog(D(1, false))
//	;);
//
//	//ClearTasks();
//}

bool MP_LAMP::CheckProcessNodeEnd(int n, bool n_is_ms, int processed,
		long long int start_time) {
	long long int elapsed_time;
	if (n_is_ms) {
		if (processed > 0) {
			elapsed_time = timer_->Elapsed() - start_time;
			if (elapsed_time >= n * 1000000) {  // ms to ns
				return true;
			}
		}
	} else if (processed >= n)
		return true;

	return false;
}

void MP_LAMP::PreProcessRootNode() {
	long long int start_time;
	start_time = timer_->Elapsed();

	expand_num_++;

	// what needed is setting itemset_buf_ to empty itemset
	// can be done faster
	// assuming root itemset is pushed before
	node_stack_->CopyItem(node_stack_->Top(), itemset_buf_);
	node_stack_->Pop();

	// dbg
	DBG(D(2) << "preprocess root node "
	;);
	DBG(node_stack_->Print(D(2), itemset_buf_)
	;);

	// calculate support from itemset_buf_
	bsh_->Set(sup_buf_);
	// skipping rest for root node

	int core_i = g_->CoreIndex(*node_stack_, itemset_buf_);

	int * ppc_ext_buf;
	// todo: use database reduction

	bool is_root_node = true;

	// reverse order
	for (int new_item = d_->NextItemInReverseLoop(is_root_node,
			mpi_data_.mpiRank_, mpi_data_.nTotalProc_, d_->NuItems());
			new_item >= core_i + 1;
			new_item = d_->NextItemInReverseLoop(is_root_node,
					mpi_data_.mpiRank_, mpi_data_.nTotalProc_, new_item)) {
		// skipping not needed because itemset_buf_ if root itemset

		bsh_->Copy(sup_buf_, child_sup_buf_);
		int sup_num = bsh_->AndCountUpdate(d_->NthData(new_item),
				child_sup_buf_);

		if (sup_num < lambda_)
			continue;

		node_stack_->PushPre();
		ppc_ext_buf = node_stack_->Top();

		bool res = g_->PPCExtension(node_stack_, itemset_buf_, child_sup_buf_,
				core_i, new_item, ppc_ext_buf);

		node_stack_->SetSup(ppc_ext_buf, sup_num);
		node_stack_->PushPostNoSort();

		node_stack_->Pop(); // always pop

		if (res) {
			IncCsAccum(sup_num); // increment closed_set_num_array
			assert(sup_num >= lambda_);
			if (ExceedCsThr())
				lambda_ = NextLambdaThr();
		}
	}

	MPI_Reduce(accum_array_, accum_recv_, lambda_max_ + 1, MPI_LONG_LONG_INT,
	MPI_SUM, 0, MPI_COMM_WORLD);

	if (mpi_data_.mpiRank_ == 0) {
		for (int l = 0; l <= lambda_max_; l++) {
			accum_array_[l] = accum_recv_[l]; // overwrite here
			accum_recv_[l] = 0;
		}
	} else { // h_!=0
		for (int l = 0; l <= lambda_max_; l++) {
			accum_array_[l] = 0;
			accum_recv_[l] = 0;
		}
	}

	if (mpi_data_.mpiRank_ == 0)
		if (ExceedCsThr())
			lambda_ = NextLambdaThr();
	CallBcast(&lambda_, 1, MPI_INT);

	// reverse order
	for (int new_item = d_->NextItemInReverseLoop(is_root_node,
			mpi_data_.mpiRank_, mpi_data_.nTotalProc_, d_->NuItems());
			new_item >= core_i + 1;
			new_item = d_->NextItemInReverseLoop(is_root_node,
					mpi_data_.mpiRank_, mpi_data_.nTotalProc_, new_item)) {
		// skipping not needed because itemset_buf_ if root itemset

		bsh_->Copy(sup_buf_, child_sup_buf_);
		int sup_num = bsh_->AndCountUpdate(d_->NthData(new_item),
				child_sup_buf_);

		if (sup_num < lambda_)
			continue;

		node_stack_->PushPre();
		ppc_ext_buf = node_stack_->Top();

		bool res = g_->PPCExtension(node_stack_, itemset_buf_, child_sup_buf_,
				core_i, new_item, ppc_ext_buf);

		node_stack_->SetSup(ppc_ext_buf, sup_num);
		node_stack_->PushPostNoSort();

		if (!res) { // todo: remove this redundancy
			node_stack_->Pop();
		} else {
			node_stack_->SortTop();
			// note: IncCsAccum already done above

			assert(sup_num >= lambda_);
			if (sup_num <= lambda_)
				node_stack_->Pop();
		}
	}

	long long int elapsed_time = timer_->Elapsed() - start_time;
	log_.d_.preprocess_time_ += elapsed_time;

	DBG(
			D(2) << "preprocess root node finished" << "\tlambda=" << lambda_
					<< "\ttime=" << elapsed_time << std::endl
			;);
}

bool MP_LAMP::ProcessNode(MPI_Data& mpi_data, TreeSearchData* treesearch_data,
		GetMinSupData* getminsup_data, GetTestableData* gettestable_data) {
	if (node_stack_->Empty())
		return false;
	long long int start_time, lap_time;
	start_time = timer_->Elapsed();
	lap_time = start_time;

	int processed = 0;
	mpi_data.processing_node_ = true;
	while (!treesearch_data->node_stack_->Empty()) {
		processed++;
		expand_num_++;

		treesearch_data->node_stack_->CopyItem(
				treesearch_data->node_stack_->Top(), itemset_buf_);
		treesearch_data->node_stack_->Pop();

		// dbg
		DBG(D(3) << "expanded "
		;);
		DBG(treesearch_data->node_stack_->Print(D(3), itemset_buf_)
		;);

		// TODO: THIS SHIT IS DEPENDENT ON bsh_
		// calculate support from itemset_buf_
		bsh_->Set(sup_buf_);
		{
			int n = treesearch_data->node_stack_->GetItemNum(itemset_buf_);
			for (int i = 0; i < n; i++) {
				int item = treesearch_data->node_stack_->GetNthItem(
						itemset_buf_, i);
				bsh_->And(d_->NthData(item), sup_buf_);
			}
		}

		int core_i = g_->CoreIndex(*node_stack_, itemset_buf_);

		int * ppc_ext_buf;
		// todo: use database reduction

		assert(
				phase_ != 1
						|| treesearch_data->node_stack_->GetItemNum(
								itemset_buf_) != 0);

		bool is_root_node = (treesearch_data->node_stack_->GetItemNum(
				itemset_buf_) == 0);

		int accum_period_counter_ = 0;
		// reverse order
		// for ( int new_item = d_->NuItems()-1 ; new_item >= core_i+1 ; new_item-- )
		for (int new_item = d_->NextItemInReverseLoop(is_root_node,
				mpi_data.mpiRank_, mpi_data.nTotalProc_, d_->NuItems());
				new_item >= core_i + 1;
				new_item = d_->NextItemInReverseLoop(is_root_node,
						mpi_data.mpiRank_, mpi_data.nTotalProc_, new_item)) {
			// skip existing item
			// todo: improve speed here
			if (treesearch_data->node_stack_->Exist(itemset_buf_, new_item))
				continue;

			{      // Periodic probe. (do in both phases)
				accum_period_counter_++;
				if (FLAGS_probe_period_is_ms) {      // using milli second
					if (accum_period_counter_ >= 64) {
						// to avoid calling timer_ frequently, time is checked once in 64 loops
						// clock_gettime takes 0.3--0.5 micro sec
						accum_period_counter_ = 0;
						long long int elt = timer_->Elapsed();
						if (elt - lap_time >= FLAGS_probe_period * 1000000) {
							log_.d_.process_node_time_ += elt - lap_time;

							Probe(mpi_data);
							Distribute(mpi_data);
							Reject(mpi_data);

							lap_time = timer_->Elapsed();
						}
					}
				} else {            // not using milli second
					if (accum_period_counter_ >= FLAGS_probe_period) {
						accum_period_counter_ = 0;
						log_.d_.process_node_time_ += timer_->Elapsed()
								- lap_time;

						Probe(mpi_data);
						Distribute(mpi_data);
						Reject(mpi_data);

						lap_time = timer_->Elapsed();
					}
				}
				// note: do this before PushPre is called [2015-10-05 21:56]

				// todo: if database reduction is implemented,
				//       do something here for changed lambda_ (skipping new_item value ?)
			}

			bsh_->Copy(sup_buf_, child_sup_buf_);
			int sup_num = bsh_->AndCountUpdate(d_->NthData(new_item),
					child_sup_buf_);

			if (sup_num < lambda_)
				continue;

			treesearch_data->node_stack_->PushPre();
			ppc_ext_buf = treesearch_data->node_stack_->Top();

			bool res = g_->PPCExtension(treesearch_data->node_stack_,
					itemset_buf_, child_sup_buf_, core_i, new_item,
					ppc_ext_buf);

			treesearch_data->node_stack_->SetSup(ppc_ext_buf, sup_num);
			treesearch_data->node_stack_->PushPostNoSort();

			if (!res) {        // todo: remove this redundancy
				treesearch_data->node_stack_->Pop();
			} else {
				treesearch_data->node_stack_->SortTop();

				DBG(if (phase_ == 2) {
					D(3) << "found cs "
					;
					treesearch_data->node_stack_->Print(D(3), ppc_ext_buf)
					;
				});

				if (phase_ == 1)
					IncCsAccum(sup_num); // increment closed_set_num_array
				if (phase_ == 2) {
					closed_set_num_++;
					if (FLAGS_third_phase) {
						int pos_sup_num = bsh_->AndCount(d_->PosNeg(),
								child_sup_buf_);
						double pval = d_->PVal(sup_num, pos_sup_num);
						assert(pval >= 0.0);
						if (pval <= gettestable_data->sig_level_) { // permits == case?
							gettestable_data->freq_stack_->PushPre();
							int * item = gettestable_data->freq_stack_->Top();
							gettestable_data->freq_stack_->CopyItem(ppc_ext_buf,
									item);
							gettestable_data->freq_stack_->PushPostNoSort();

							gettestable_data->freq_map_.insert(
									std::pair<double, int*>(pval, item));
						}
					}
				}

				assert(sup_num >= lambda_);

				// try skipping if supnum_ == sup_threshold,
				// because if sup_num of a node equals to sup_threshold, children will have smaller sup_num
				// therefore no need to check it's children
				// note: skipping node_stack_ full check. allocate enough memory!
				if (sup_num <= lambda_)
					node_stack_->Pop();
			}
		}

		if (CheckProcessNodeEnd(mpi_data_.granularity_,
				mpi_data_.isGranularitySec_, processed, start_time))
			break;
	}

	long long int elapsed_time = timer_->Elapsed() - lap_time;
	log_.d_.process_node_time_ += elapsed_time;
	log_.d_.process_node_num_ += processed;

	DBG(
			D(2) << "processed node num=" << processed << "\ttime="
					<< elapsed_time << std::endl
			;);

	mpi_data_.processing_node_ = false;
	return true;
}

//bool MP_LAMP::ProcessNodeStraw1(int n) {
//	if (node_stack_->Empty())
//		return false;
//	long long int start_time, lap_time;
//	start_time = timer_->Elapsed();
//	lap_time = start_time;
//
//	int processed = 0;
//	mpi_data_.processing_node_ = true;
//	while (!node_stack_->Empty()) {
//		processed++;
//		expand_num_++;
//
//		node_stack_->CopyItem(node_stack_->Top(), itemset_buf_);
//		node_stack_->Pop();
//
//		// dbg
//		DBG(D(3) << "expanded "
//		;);
//		DBG(node_stack_->Print(D(3), itemset_buf_)
//		;);
//
//		// calculate support from itemset_buf_
//		bsh_->Set(sup_buf_);
//		{
//			int n = node_stack_->GetItemNum(itemset_buf_);
//			for (int i = 0; i < n; i++) {
//				int item = node_stack_->GetNthItem(itemset_buf_, i);
//				bsh_->And(d_->NthData(item), sup_buf_);
//			}
//		}
//
//		int core_i = g_->CoreIndex(*node_stack_, itemset_buf_);
//
//		int * ppc_ext_buf;
//		// todo: use database reduction
//
//		assert(phase_ != 1 || node_stack_->GetItemNum(itemset_buf_) != 0);
//
//		bool is_root_node = (node_stack_->GetItemNum(itemset_buf_) == 0);
//
//		int accum_period_counter_ = 0;
//		// reverse order
//		// for ( int new_item = d_->NuItems()-1 ; new_item >= core_i+1 ; new_item-- )
//		for (int new_item = d_->NextItemInReverseLoop(is_root_node,
//				mpi_data_.mpiRank_, mpi_data_.nTotalProc_, d_->NuItems());
//				new_item >= core_i + 1;
//				new_item = d_->NextItemInReverseLoop(is_root_node,
//						mpi_data_.mpiRank_, mpi_data_.nTotalProc_, new_item)) {
//			// skip existing item
//			// todo: improve speed here
//			if (node_stack_->Exist(itemset_buf_, new_item))
//				continue;
//
//			{      // Periodic probe. (do in both phases)
//				accum_period_counter_++;
//				if (FLAGS_probe_period_is_ms) {      // using milli second
//					if (accum_period_counter_ >= 64) {
//						// to avoid calling timer_ frequently, time is checked once in 64 loops
//						// clock_gettime takes 0.3--0.5 micro sec
//						accum_period_counter_ = 0;
//						long long int elt = timer_->Elapsed();
//						if (elt - lap_time >= FLAGS_probe_period * 1000000) {
//							log_.d_.process_node_time_ += elt - lap_time;
//
//							Probe();
//
//							lap_time = timer_->Elapsed();
//						}
//					}
//				} else {            // not using milli second
//					if (accum_period_counter_ >= FLAGS_probe_period) {
//						accum_period_counter_ = 0;
//						log_.d_.process_node_time_ += timer_->Elapsed()
//								- lap_time;
//
//						Probe();
//
//						lap_time = timer_->Elapsed();
//					}
//				}
//				// note: do this before PushPre is called [2015-10-05 21:56]
//
//				// todo: if database reduction is implemented,
//				//       do something here for changed lambda_ (skipping new_item value ?)
//			}
//
//			bsh_->Copy(sup_buf_, child_sup_buf_);
//			int sup_num = bsh_->AndCountUpdate(d_->NthData(new_item),
//					child_sup_buf_);
//
//			if (sup_num < lambda_)
//				continue;
//
//			node_stack_->PushPre();
//			ppc_ext_buf = node_stack_->Top();
//
//			bool res = g_->PPCExtension(node_stack_, itemset_buf_,
//					child_sup_buf_, core_i, new_item, ppc_ext_buf);
//
//			node_stack_->SetSup(ppc_ext_buf, sup_num);
//			node_stack_->PushPostNoSort();
//
//			if (!res) {        // todo: remove this redundancy
//				node_stack_->Pop();
//			} else {
//				node_stack_->SortTop();
//
//				DBG(if (phase_ == 2) {
//					D(3) << "found cs "
//					;
//					node_stack_->Print(D(3), ppc_ext_buf)
//					;
//				});
//
//				if (phase_ == 1)
//					IncCsAccum(sup_num); // increment closed_set_num_array
//				if (phase_ == 2)
//					closed_set_num_++;
//
//				if (phase_ == 2 && FLAGS_third_phase) {
//					int pos_sup_num = bsh_->AndCount(d_->PosNeg(),
//							child_sup_buf_);
//					double pval = d_->PVal(sup_num, pos_sup_num);
//					assert(pval >= 0.0);
//					if (pval <= sig_level_) { // permits == case?
//						freq_stack_->PushPre();
//						int * item = freq_stack_->Top();
//						freq_stack_->CopyItem(ppc_ext_buf, item);
//						freq_stack_->PushPostNoSort();
//
//						freq_map_.insert(std::pair<double, int*>(pval, item));
//					}
//				}
//
//				assert(sup_num >= lambda_);
//
//				// try skipping if supnum_ == sup_threshold,
//				// because if sup_num of a node equals to sup_threshold, children will have smaller sup_num
//				// therefore no need to check it's children
//				// note: skipping node_stack_ full check. allocate enough memory!
//				if (sup_num <= lambda_)
//					node_stack_->Pop();
//			}
//		}
//
//		if (CheckProcessNodeEnd(n, mpi_data_.isGranularitySec_, processed,
//				start_time))
//			break;
//	}
//
//	long long int elapsed_time = timer_->Elapsed() - lap_time;
//	log_.d_.process_node_time_ += elapsed_time;
//	log_.d_.process_node_num_ += processed;
//
//	DBG(
//			D(2) << "processed node num=" << processed << "\ttime="
//					<< elapsed_time << std::endl
//			;);
//
//	mpi_data_.processing_node_ = false;
//	return true;
//}

// lifeline == -1 for random thieves
void MP_LAMP::SendRequest(MPI_Data& mpi_data, int dst, int is_lifeline) {
	assert(dst >= 0);
	int message[2];
	message[0] = dtd_.time_zone_;
	message[1] = is_lifeline; // -1 for random thieves, >=0 for lifeline thieves
	assert(dst < mpi_data.nTotalProc_ && "SendRequest");

	CallBsend(message, 2, MPI_INT, dst, Tag::REQUEST);
	dtd_.OnSend();

	DBG(
			D(2) << "SendRequest: dst=" << dst << "\tis_lifeline="
					<< is_lifeline << "\tdtd_count=" << dtd_.count_ << std::endl
			;);
}

void MP_LAMP::RecvRequest(MPI_Data& mpi_data, int src) {
	DBG(D(2) << "RecvRequest: src=" << src
	;);
	MPI_Status recv_status;
	int message[2];

	CallRecv(&message, 2, MPI_INT, src, Tag::REQUEST, &recv_status);
	dtd_.OnRecv();
	assert(src == recv_status.MPI_SOURCE);

	int timezone = message[0];
	dtd_.UpdateTimeZone(timezone);
	int is_lifeline = message[1]; // -1 for random thieves, >=0 for lifeline thieves
	int thief = src;
	DBG(D(2, false) << "\tis_lifeline=" << is_lifeline
	;);
	DBG(D(2, false) << "\tthief=" << src
	;);
	DBG(D(2, false) << "\tdtd_count=" << dtd_.count_
	;);

	if (node_stack_->Empty()) {
		DBG(D(2, false) << "\tempty and reject" << std::endl
		;);
		if (is_lifeline >= 0) {
			mpi_data.lifeline_thieves_->Push(thief);
			SendReject(mpi_data, thief); // notify
		} else {
			SendReject(mpi_data, thief); // notify
		}
	} else {
		DBG(D(2, false) << "\tpush" << std::endl
		;);
		if (is_lifeline >= 0)
			mpi_data.thieves_->Push(thief);
		else
			mpi_data.thieves_->Push(-thief - 1);
	}
	// todo: take log
}

void MP_LAMP::SendReject(MPI_Data& mpi_data, int dst) {
	assert(dst >= 0);
	int message[1];

	message[0] = dtd_.time_zone_;
	assert(dst < mpi_data.nTotalProc_ && "SendReject");
	CallBsend(message, 1, MPI_INT, dst, Tag::REJECT);
	dtd_.OnSend();

	DBG(
			D(2) << "SendReject: dst=" << dst << "\tdtd_count=" << dtd_.count_
					<< std::endl
			;);
}

void MP_LAMP::RecvReject(MPI_Data& mpi_data, int src) {
	MPI_Status recv_status;
	int message[1];

	CallRecv(&message, 1, MPI_INT, src, Tag::REJECT, &recv_status);
	dtd_.OnRecv();
	assert(src == recv_status.MPI_SOURCE);

	int timezone = message[0];
	dtd_.UpdateTimeZone(timezone);
	DBG(
			D(2) << "RecvReject: src=" << src << "\tdtd_count=" << dtd_.count_
					<< std::endl
			;);

	stealer_.ResetRequesting();
	waiting_ = false;
}

void MP_LAMP::SendGive(MPI_Data& mpi_data, VariableLengthItemsetStack * st,
		int dst, int is_lifeline) {
	assert(dst >= 0);
	st->SetTimestamp(dtd_.time_zone_);
	st->SetFlag(is_lifeline);
	int * message = st->Stack();
	int size = st->UsedCapacity();

	log_.d_.give_stack_max_itm_ = std::max(log_.d_.give_stack_max_itm_,
			(long long int) (st->NuItemset()));
	log_.d_.give_stack_max_cap_ = std::max(log_.d_.give_stack_max_cap_,
			(long long int) (st->UsedCapacity()));

	assert(dst < mpi_data.nTotalProc_ && "SendGive");
	CallBsend(message, size, MPI_INT, dst, Tag::GIVE);
	dtd_.OnSend();

	DBG(
			D(2) << "SendGive: " << "\ttimezone=" << dtd_.time_zone_ << "\tdst="
					<< dst << "\tlfl=" << is_lifeline << "\tsize=" << size
					<< "\tnode=" << st->NuItemset() << "\tdtd_count="
					<< dtd_.count_ << std::endl
			;);
	//st->PrintAll(D(false));
}

void MP_LAMP::RecvGive(MPI_Data& mpi_data, int src, MPI_Status probe_status) {
	int count;
	int error = MPI_Get_count(&probe_status, MPI_INT, &count);
	if (error != MPI_SUCCESS) {
		DBG(D(1) << "error in MPI_Get_count in RecvGive: " << error << std::endl
		;);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	MPI_Status recv_status;

	CallRecv(give_stack_->Stack(), count, MPI_INT, src, Tag::GIVE,
			&recv_status);
	dtd_.OnRecv();
	assert(src == recv_status.MPI_SOURCE);

	int timezone = give_stack_->Timestamp();
	dtd_.UpdateTimeZone(timezone);

	int flag = give_stack_->Flag();
	int orig_nu_itemset = node_stack_->NuItemset();

	node_stack_->MergeStack(
			give_stack_->Stack() + VariableLengthItemsetStack::SENTINEL + 1,
			count - VariableLengthItemsetStack::SENTINEL - 1);
	int new_nu_itemset = node_stack_->NuItemset();

	if (flag >= 0) {
		mpi_data.lifelines_activated_[src] = false;
		log_.d_.lifeline_steal_num_++;
		log_.d_.lifeline_nodes_received_ += (new_nu_itemset - orig_nu_itemset);
	} else {
		log_.d_.steal_num_++;
		log_.d_.nodes_received_ += (new_nu_itemset - orig_nu_itemset);
	}

	DBG(
			D(2) << "RecvGive: src=" << src << "\ttimezone=" << dtd_.time_zone_
					<< "\tlfl=" << flag << "\tsize=" << count << "\tnode="
					<< (new_nu_itemset - orig_nu_itemset) << "\tdtd_count="
					<< dtd_.count_ << std::endl
			;);

	give_stack_->Clear();

	log_.d_.node_stack_max_itm_ = std::max(log_.d_.node_stack_max_itm_,
			(long long int) (node_stack_->NuItemset()));
	log_.d_.node_stack_max_cap_ = std::max(log_.d_.node_stack_max_cap_,
			(long long int) (node_stack_->UsedCapacity()));

	DBG(node_stack_->PrintAll(D(3, false))
	;);

	stealer_.ResetRequesting();
	stealer_.ResetCounters();
	stealer_.SetState(StealState::RANDOM);
	stealer_.SetStealStart();
	waiting_ = false;
}

void MP_LAMP::SendLambda(MPI_Data& mpi_data, int lambda) {
	// send lambda to bcast_targets_
	int message[2];
	message[0] = dtd_.time_zone_;
	message[1] = lambda;

	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data.bcast_targets_[i] < 0)
			break;
		assert(
				mpi_data.bcast_targets_[i] < mpi_data.nTotalProc_
						&& "SendLambda");
		CallBsend(message, 2, MPI_INT, mpi_data.bcast_targets_[i], Tag::LAMBDA);
		dtd_.OnSend();

		DBG(
				D(2) << "SendLambda: dst=" << mpi_data.bcast_targets_[i]
						<< "\tlambda=" << lambda << "\tdtd_count="
						<< dtd_.count_ << std::endl
				;);
	}
}

void MP_LAMP::RecvLambda(MPI_Data& mpi_data, int src) {
	MPI_Status recv_status;
	int message[2];

	CallRecv(&message, 2, MPI_INT, src, Tag::LAMBDA, &recv_status);
	dtd_.OnRecv();
	assert(src == recv_status.MPI_SOURCE);
	int timezone = message[0];
	dtd_.UpdateTimeZone(timezone);

	DBG(
			D(2) << "RecvLambda: src=" << src << "\tlambda=" << message[1]
					<< "\tdtd_count=" << dtd_.count_ << std::endl
			;);

	int new_lambda = message[1];
	if (new_lambda > lambda_) {
		SendLambda(mpi_data, new_lambda);
		lambda_ = new_lambda;
		// todo: do database reduction
	}
}

bool MP_LAMP::AccumCountReady(MPI_Data& mpi_data) const {
	for (int i = 0; i < k_echo_tree_branch; i++)
		if (mpi_data.bcast_targets_[i] >= 0 && !(mpi_data.accum_flag_[i]))
			return false;
	return true;
	// check only valid bcast_targets_
	// always true if leaf
}

void MP_LAMP::CheckCSThreshold(MPI_Data& mpi_data) {
//	assert(mpi_data.mpiRank_ == 0);
	if (ExceedCsThr()) {
		int new_lambda = NextLambdaThr();
		SendLambda(mpi_data, new_lambda);
		lambda_ = new_lambda;
	}
}

//void MP_LAMP::CheckCSThreshold(GetMinSupData* getminsup_data) {
////	assert(mpi_data.mpiRank_ == 0);
//	if (ExceedCsThr()) {
//		int new_lambda = NextLambdaThr();
//		SendLambda(new_lambda);
//		getminsup_data->lambda_ = new_lambda;
//	}
//}

bool MP_LAMP::ExceedCsThr() const {
	// note: > is correct. permit ==
	return (accum_array_[lambda_] > cs_thr_[lambda_]);
}

int MP_LAMP::NextLambdaThr() const {
	int si;
	for (si = lambda_max_; si >= lambda_; si--)
		if (accum_array_[si] > cs_thr_[si])
			break;
	return si + 1;
	// it is safe because lambda_ higher than max results in immediate search finish
}

//int MP_LAMP::NextLambdaThr(GetMinSupData* getminsup_data) const {
//	int si;
//	for (si = getminsup_data->lambda_max_; si >= getminsup_data->lambda_; si--)
//		if (getminsup_data->accum_array_[si] > getminsup_data->cs_thr_[si])
//			break;
//	return si + 1;
//	// it is safe because lambda_ higher than max results in immediate search finish
//}

void MP_LAMP::IncCsAccum(int sup_num) {
	for (int i = sup_num; i >= lambda_ - 1; i--)
		accum_array_[i]++;
}

double MP_LAMP::GetInterimSigLevel(int lambda) const {
	long long int csnum = accum_array_[lambda];
	double lv;
	if (csnum > 0)
		lv = FLAGS_a / (double) csnum;
	else
		lv = FLAGS_a;

	return lv;
}

//==============================================================================

void MP_LAMP::SendResultRequest(MPI_Data& mpi_data) {
	int message[1];
	message[0] = 1; // dummy

	mpi_data.echo_waiting_ = true;

	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data.bcast_targets_[i] < 0)
			break;

		assert(
				mpi_data.bcast_targets_[i] < mpi_data.nTotalProc_
						&& "SendResultRequest");
		CallBsend(message, 1, MPI_INT, mpi_data.bcast_targets_[i],
				Tag::RESULT_REQUEST);
		DBG(
				D(2) << "SendResultRequest: dst=" << mpi_data.bcast_targets_[i]
						<< std::endl
				;);
	}
}

void MP_LAMP::RecvResultRequest(MPI_Data& mpi_data, int src) {
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::RESULT_REQUEST, &recv_status);
	assert(src == recv_status.MPI_SOURCE);

	DBG(D(2) << "RecvResultRequest: src=" << src << std::endl
	;);

	if (IsLeaf(mpi_data))
		SendResultReply(mpi_data);
	else
		SendResultRequest(mpi_data);
}

void MP_LAMP::SendResultReply(MPI_Data& mpi_data) {
	int * message = significant_stack_->Stack();
	int size = significant_stack_->UsedCapacity();
	assert(mpi_data.bcast_source_ < mpi_data.nTotalProc_ && "SendResultReply");
	CallBsend(message, size, MPI_INT, mpi_data.bcast_source_,
			Tag::RESULT_REPLY);

	DBG(D(2) << "SendResultReply: dst=" << mpi_data.bcast_source_ << std::endl
	;);
	DBG(significant_stack_->PrintAll(D(3, false))
	;);

	mpi_data.echo_waiting_ = false;
	dtd_.terminated_ = true;
}

void MP_LAMP::RecvResultReply(MPI_Data& mpi_data, int src,
		MPI_Status probe_status) {
	int count;
	int error = MPI_Get_count(&probe_status, MPI_INT, &count);
	if (error != MPI_SUCCESS) {
		DBG(
				D(1) << "error in MPI_Get_count in RecvResultReply: " << error
						<< std::endl
				;);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	MPI_Status recv_status;
	CallRecv(give_stack_->Stack(), count, MPI_INT, src, Tag::RESULT_REPLY,
			&recv_status);
	assert(src == recv_status.MPI_SOURCE);

	significant_stack_->MergeStack(
			give_stack_->Stack() + VariableLengthItemsetStack::SENTINEL + 1,
			count - VariableLengthItemsetStack::SENTINEL - 1);
	give_stack_->Clear();

	DBG(D(2) << "RecvResultReply: src=" << src << std::endl
	;);
	DBG(significant_stack_->PrintAll(D(3, false))
	;);

	// using the same flags as accum count, should be fixed
	bool flag = false;
	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data.bcast_targets_[i] == src) {
			flag = true;
			mpi_data.accum_flag_[i] = true;
			break;
		}
	}
	assert(flag);

	if (AccumCountReady(mpi_data)) {
		if (mpi_data.mpiRank_ != 0) {
			SendResultReply(mpi_data);
		} else { // root
			mpi_data.echo_waiting_ = false;
			dtd_.terminated_ = true;
		}
	}
}

void MP_LAMP::ExtractSignificantSet() {
	std::multimap<double, int *>::iterator it;
	for (it = freq_map_.begin(); it != freq_map_.end(); ++it) {
		// permits == case
		if ((*it).first <= final_sig_level_) {
			significant_stack_->PushPre();
			int * item = significant_stack_->Top();
			significant_stack_->CopyItem((*it).second, item);
			significant_stack_->PushPostNoSort();
		} else
			break;
	}
}

void MP_LAMP::SortSignificantSets() {
	int * set = significant_stack_->FirstItemset();

	while (set != NULL) {
		// calculate support from set
		bsh_->Set(sup_buf_);
		{
			int n = significant_stack_->GetItemNum(set);
			for (int i = 0; i < n; i++) {
				int item = significant_stack_->GetNthItem(set, i);
				bsh_->And(d_->NthData(item), sup_buf_);
			}
		}

		int sup_num = bsh_->Count(sup_buf_);
		int pos_sup_num = bsh_->AndCount(d_->PosNeg(), sup_buf_);
		double pval = d_->PVal(sup_num, pos_sup_num);

		significant_set_.insert(
				SignificantSetResult(pval, set, sup_num, pos_sup_num,
						significant_stack_));
		set = significant_stack_->NextItemset(set);
	}
}

//==============================================================================

void MP_LAMP::Distribute(MPI_Data& mpi_data) {
	if (mpi_data.nTotalProc_ == 1)
		return;
	DBG(D(3) << "Distribute" << std::endl
	;);
	if (mpi_data.thieves_->Size() > 0
			|| mpi_data.lifeline_thieves_->Size() > 0) {
		int steal_num = node_stack_->Split(give_stack_);
		if (steal_num > 0) {
			DBG(D(3) << "giving" << std::endl
			;);
			DBG(give_stack_->PrintAll(D(3, false))
			;);
			Give(mpi_data, give_stack_, steal_num);
			give_stack_->Clear();
		}
	}
}

void MP_LAMP::Give(MPI_Data& mpi_data, VariableLengthItemsetStack * st,
		int steal_num) {
	DBG(D(3) << "Give: "
	;);
	if (mpi_data.thieves_->Size() > 0) { // random thieves
		int thief = mpi_data.thieves_->Pop();
		if (thief >= 0) { // lifeline thief
			DBG(
					D(3, false) << "thief=" << thief << "\trandom lifeline"
							<< std::endl
					;);
			log_.d_.lifeline_given_num_++;
			log_.d_.lifeline_nodes_given_ += steal_num;
			SendGive(mpi_data, st, thief, mpi_data.mpiRank_);
		} else { // random thief
			DBG(
					D(3, false) << "thief=" << (-thief - 1) << "\trandom"
							<< std::endl
					;);
			log_.d_.given_num_++;
			log_.d_.nodes_given_ += steal_num;
			SendGive(mpi_data, st, (-thief - 1), -1);
		}
	} else {
		assert(mpi_data.lifeline_thieves_->Size() > 0);
		int thief = mpi_data.lifeline_thieves_->Pop();
		assert(thief >= 0);
		DBG(D(3, false) << "thief=" << thief << "\tlifeline" << std::endl
		;);
		log_.d_.lifeline_given_num_++;
		log_.d_.lifeline_nodes_given_ += steal_num;
		SendGive(mpi_data, st, thief, mpi_data.mpiRank_);
	}
}

void MP_LAMP::Reject(MPI_Data& mpi_data) {
	if (mpi_data.nTotalProc_ == 1)
		return;
	DBG(D(3) << "Reject" << std::endl
	;);
	// discard random thieves, move lifeline thieves to lifeline thieves stack
	while (mpi_data.thieves_->Size() > 0) {
		int thief = mpi_data.thieves_->Pop();
		if (thief >= 0) { // lifeline thief
			mpi_data.lifeline_thieves_->Push(thief);
			SendReject(mpi_data, thief);
		} else { // random thief
			SendReject(mpi_data, -thief - 1);
		}
	}
}

void MP_LAMP::Steal(MPI_Data& mpi_data) {
	DBG(D(3) << "Steal" << std::endl
	;);
	if (mpi_data.nTotalProc_ == 1)
		return;
	if (!stealer_.StealStarted())
		return;
	if (stealer_.Requesting())
		return;
	if (!node_stack_->Empty())
		return;

	// reset state, counter

	switch (stealer_.State()) {
	case StealState::RANDOM: {
		int victim = mpi_data.victims_[mpi_data.rand_m_()];
		assert(victim <= mpi_data.nTotalProc_ && "stealrandom");
		SendRequest(mpi_data, victim, -1);
		stealer_.SetRequesting();
		DBG(
				D(2) << "Steal Random:" << "\trequesting="
						<< stealer_.Requesting() << "\trandom counter="
						<< stealer_.RandomCount() << std::endl
				;);
		stealer_.IncRandomCount();
		if (stealer_.RandomCount() == 0)
			stealer_.SetState(StealState::LIFELINE);
	}
		break;
	case StealState::LIFELINE: {
		int lifeline = mpi_data.lifelines_[stealer_.LifelineVictim()];
		assert(lifeline >= 0); // at least lifelines_[0] must be 0 or higher
		if (!mpi_data.lifelines_activated_[lifeline]) {
			mpi_data.lifelines_activated_[lifeline] = true;
			// becomes false only at RecvGive
			assert(lifeline <= mpi_data.nTotalProc_ && "lifeline");
			SendRequest(mpi_data, lifeline, 1);
			DBG(
					D(2) << "Steal Lifeline:" << "\trequesting="
							<< stealer_.Requesting() << "\tlifeline counter="
							<< stealer_.LifelineCounter()
							<< "\tlifeline victim=" << stealer_.LifelineVictim()
							<< "\tz_=" << mpi_data.hypercubeDimension_
							<< std::endl
					;);
		}

		stealer_.IncLifelineCount();
		stealer_.IncLifelineVictim();
		if (stealer_.LifelineCounter() >= mpi_data.hypercubeDimension_
				|| mpi_data.lifelines_[stealer_.LifelineVictim()] < 0) { // hack fix
// note: this resets next_lifeline_victim_, is this OK?
			stealer_.Finish();
			DBG(D(2) << "Steal finish:" << std::endl
			;);
			// note: one steal phase finished, don't restart until receiving give?

			// in x10 glb, if steal() finishes and empty==true,
			// processStack() will finish, but will be restarted by deal()
			// for the same behavior, reseting the states to initial state should be correct
		}
	}
		break;
	default:
		assert(0);
		break;
	}
}

//void MP_LAMP::Steal2() {
//	DBG(D(3) << "Steal" << std::endl
//	;);
//	if (mpi_data.nTotalProc_ == 1)
//		return;
//
//	// random steal
//	for (int i = 0; i < mpi_data.nRandStealTrials_ && node_stack_->Empty();
//			i++) {
//		// todo: log steal trial
//		int victim = mpi_data.victims_[mpi_data.rand_m_()];
//		SendRequest(victim, -1);
//		waiting_ = true;
//		DBG(D(2) << "Steal Random:" << "\trandom counter=" << i << std::endl
//		;);
//		// todo: measure idle time
//		while (waiting_)
//			Probe();
//	}
//	if (dtd_.terminated_)
//		return;
//
//	// lifeline steal
//	for (int i = 0; i < mpi_data.hypercubeDimension_ && node_stack_->Empty();
//			i++) {
//		int lifeline = mpi_data.lifelines_[i];
//		if (lifeline < 0)
//			break; // break if -1
//		if (!mpi_data.lifelines_activated_[lifeline]) {
//			mpi_data.lifelines_activated_[lifeline] = true;
//			// becomes false only at RecvGive
//			SendRequest(lifeline, 1);
//			waiting_ = true;
//			DBG(
//					D(2) << "Steal Lifeline:" << "\ti=" << i << "\tz_="
//							<< mpi_data.hypercubeDimension_ << std::endl
//					;);
//			// todo: measure idle time
//			while (waiting_)
//				Probe();
//		}
//	}
//}

// mainloop
void MP_LAMP::GetMinimalSupport(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, GetMinSupData* getminsup_data) {
	phase_ = 1;
	DBG(D(1) << "MainLoop" << std::endl
	;);
	while (!dtd_.terminated_) {
		while (!dtd_.terminated_) {
			if (ProcessNode(mpi_data, treesearch_data, getminsup_data,
					(GetTestableData*) NULL)) {
				log_.d_.node_stack_max_itm_ =
						std::max(log_.d_.node_stack_max_itm_,
								(long long int) (treesearch_data->node_stack_->NuItemset()));
				log_.d_.node_stack_max_cap_ =
						std::max(log_.d_.node_stack_max_cap_,
								(long long int) (treesearch_data->node_stack_->UsedCapacity()));
				Probe(mpi_data);
				if (dtd_.terminated_)
					break;
				Distribute(mpi_data);
				Reject(mpi_data); // distribute finished, reject remaining requests
				if (mpi_data.mpiRank_ == 0)
					CheckCSThreshold(mpi_data);
			} else
				break;
		}
		if (dtd_.terminated_)
			break;

		log_.idle_start_ = timer_->Elapsed();
		Reject(mpi_data); // node_stack_ empty. reject requests
		Steal(mpi_data); // request steal
		if (dtd_.terminated_) {
			log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
			break;
		}

		Probe(mpi_data);
		if (dtd_.terminated_) {
			log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
			break;
		}
		if (mpi_data.mpiRank_ == 0)
			CheckCSThreshold(mpi_data);
		log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
	}

}

void MP_LAMP::GetTestablePatterns(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, GetTestableData* gettestable_data) {
	phase_ = 2;
	DBG(D(1) << "MainLoop" << std::endl
	;);
	while (!dtd_.terminated_) {
		while (!dtd_.terminated_) {
			if (ProcessNode(mpi_data, treesearch_data, (GetMinSupData*) NULL,
					gettestable_data)) {
				log_.d_.node_stack_max_itm_ = std::max(
						log_.d_.node_stack_max_itm_,
						(long long int) (node_stack_->NuItemset()));
				log_.d_.node_stack_max_cap_ = std::max(
						log_.d_.node_stack_max_cap_,
						(long long int) (node_stack_->UsedCapacity()));
				Probe(mpi_data);
				if (dtd_.terminated_)
					break;
				Distribute(mpi_data);
				Reject(mpi_data); // distribute finished, reject remaining requests
			} else
				break;
		}
		if (dtd_.terminated_)
			break;

		log_.idle_start_ = timer_->Elapsed();
		Reject(mpi_data); // node_stack_ empty. reject requests
		Steal(mpi_data); // request steal
		if (dtd_.terminated_) {
			log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
			break;
		}

		Probe(mpi_data);
		if (dtd_.terminated_) {
			log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
			break;
		}

		log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
	}

}

void MP_LAMP::GetSignificantPatterns(MPI_Data& mpi_data,
		GetSignificantData* getsignificant_data) {
	DBG(D(1) << "MainLoop" << std::endl
	;);

	ExtractSignificantSet();
	if (mpi_data.mpiRank_ == 0)
		SendResultRequest(mpi_data);

	while (!dtd_.terminated_)
		Probe(mpi_data);

}

//// mainloop
//void MP_LAMP::MainLoop() {
//	DBG(D(1) << "MainLoop" << std::endl
//	;);
//	if (phase_ == 1 || phase_ == 2) {
//		while (!dtd_.terminated_) {
//			while (!dtd_.terminated_) {
//				if (ProcessNode(mpi_data.granularity_)) {
//					log_.d_.node_stack_max_itm_ = std::max(
//							log_.d_.node_stack_max_itm_,
//							(long long int) (node_stack_->NuItemset()));
//					log_.d_.node_stack_max_cap_ = std::max(
//							log_.d_.node_stack_max_cap_,
//							(long long int) (node_stack_->UsedCapacity()));
//					Probe();
//					if (dtd_.terminated_)
//						break;
//					Distribute();
//					Reject(); // distribute finished, reject remaining requests
//					if (mpi_data.mpiRank_ == 0 && phase_ == 1)
//						CheckCSThreshold();
//				} else
//					break;
//			}
//			if (dtd_.terminated_)
//				break;
//
//			log_.idle_start_ = timer_->Elapsed();
//			Reject(); // node_stack_ empty. reject requests
//			Steal(); // request steal
//			if (dtd_.terminated_) {
//				log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
//				break;
//			}
//
//			Probe();
//			if (dtd_.terminated_) {
//				log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
//				break;
//			}
//			if (mpi_data.mpiRank_ == 0 && phase_ == 1)
//				CheckCSThreshold();
//			log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
//		}
//	} else if (phase_ == 3) {
//		ExtractSignificantSet();
//		if (mpi_data.mpiRank_ == 0)
//			SendResultRequest();
//
//		while (!dtd_.terminated_)
//			Probe();
//	}
//
//}

// mainloop for straw man 1
//void MP_LAMP::MainLoopStraw1() {
//	DBG(D(1) << "MainLoop" << std::endl
//	;);
//	if (phase_ == 1 || phase_ == 2) {
//		while (!dtd_.terminated_) {
//			while (!dtd_.terminated_) {
//				if (ProcessNodeStraw1(mpi_data.granularity_)) {
//					log_.d_.node_stack_max_itm_ = std::max(
//							log_.d_.node_stack_max_itm_,
//							(long long int) (node_stack_->NuItemset()));
//					log_.d_.node_stack_max_cap_ = std::max(
//							log_.d_.node_stack_max_cap_,
//							(long long int) (node_stack_->UsedCapacity()));
//					Probe();
//					if (dtd_.terminated_)
//						break;
//					if (mpi_data.mpiRank_ == 0 && phase_ == 1)
//						CheckCSThreshold();
//				} else
//					break;
//			}
//			if (dtd_.terminated_)
//				break;
//
//			log_.idle_start_ = timer_->Elapsed();
//
//			Probe();
//			if (dtd_.terminated_) {
//				log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
//				break;
//			}
//			if (mpi_data.mpiRank_ == 0 && phase_ == 1)
//				CheckCSThreshold();
//			log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
//		}
//	} else if (phase_ == 3) {
//		ExtractSignificantSet();
//		if (mpi_data.mpiRank_ == 0)
//			SendResultRequest();
//
//		while (!dtd_.terminated_)
//			Probe();
//	}
//
//}

//==============================================================================

bool MP_LAMP::IsLeaf(MPI_Data& mpi_data) const {
	for (int i = 0; i < k_echo_tree_branch; i++)
		if (mpi_data.bcast_targets_[i] >= 0)
			return false;
	return true;
}

int MP_LAMP::CallIprobe(MPI_Status * status, int * src, int * tag) {
	long long int start_time;
	long long int end_time;
	log_.d_.iprobe_num_++;

	// todo: prepare non-log mode to remove measurement
	// clock_gettime takes 0.3--0.5 micro sec
	LOG(start_time = timer_->Elapsed(););

	int flag;
	int error = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag,
			status);
	if (error != MPI_SUCCESS) {
		DBG(D(1) << "error in MPI_Iprobe: " << error << std::endl
		;);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	LOG(
			end_time = timer_->Elapsed(); log_.d_.iprobe_time_ += end_time - start_time; log_.d_.iprobe_time_max_ = std::max(end_time - start_time, log_.d_.iprobe_time_max_););

	if (flag) {
		log_.d_.iprobe_succ_num_++;
		LOG(
				log_.d_.iprobe_succ_time_ += end_time - start_time; log_.d_.iprobe_succ_time_max_ = std::max(end_time - start_time, log_.d_.iprobe_succ_time_max_););

		*tag = status->MPI_TAG;
		*src = status->MPI_SOURCE;
// #ifdef NDEBUG
//     MPI_Get_count(status, type, count);
// #else
//     {
//       DBG( D(4) << "calling MPI_Get_count in CallIprobe" << std::endl; );
//       int ret = MPI_Get_count(status, type, count);
//       DBG( D(4) << "returnd from MPI_Get_count in CallIprobe" << std::endl; );
//       assert(ret == MPI_SUCCESS);
//     }
// #endif
	} else {
		log_.d_.iprobe_fail_num_++;
		LOG(
				log_.d_.iprobe_fail_time_ += end_time - start_time; log_.d_.iprobe_fail_time_max_ = std::max(end_time - start_time, log_.d_.iprobe_fail_time_max_););
	}

	return flag;
}

int MP_LAMP::CallRecv(void * buffer, int data_count, MPI_Datatype type, int src,
		int tag, MPI_Status * status) {
	long long int start_time;
	long long int end_time;

	// todo: prepare non-log mode to remove measurement
	// clock_gettime takes 0.3--0.5 micro sec
	log_.d_.recv_num_++;
	LOG(start_time = timer_->Elapsed(););

	int error = MPI_Recv(buffer, data_count, type, src, tag, MPI_COMM_WORLD,
			status);

	LOG(
			end_time = timer_->Elapsed(); log_.d_.recv_time_ += end_time - start_time; log_.d_.recv_time_max_ = std::max(end_time - start_time, log_.d_.recv_time_max_););
	return error;
}

int MP_LAMP::CallBsend(void * buffer, int data_count, MPI_Datatype type,
		int dest, int tag) {
	assert(0 <= dest && dest < mpi_data_.nTotalProc_);
	long long int start_time;
	long long int end_time;
	log_.d_.bsend_num_++;
	start_time = timer_->Elapsed();

	int error = MPI_Bsend(buffer, data_count, type, dest, tag, MPI_COMM_WORLD);

	end_time = timer_->Elapsed();
	log_.d_.bsend_time_ += end_time - start_time;
	log_.d_.bsend_time_max_ = std::max(end_time - start_time,
			log_.d_.bsend_time_max_);
	return error;
}

int MP_LAMP::CallBcast(void * buffer, int data_count, MPI_Datatype type) {
	long long int start_time;
	long long int end_time;
	log_.d_.bcast_num_++;
	start_time = timer_->Elapsed();

	int error = MPI_Bcast(buffer, data_count, type, 0, MPI_COMM_WORLD);

	end_time = timer_->Elapsed();
	log_.d_.bcast_time_ += end_time - start_time;
	log_.d_.bcast_time_max_ = std::max(end_time - start_time,
			log_.d_.bcast_time_max_);
	return error;
}

//==============================================================================

std::ostream & MP_LAMP::PrintDBInfo(std::ostream & out) const {
	std::stringstream s;

	s << "# ";
	d_->ShowInfo(s);

	out << s.str() << std::flush;
	return out;
}

std::ostream & MP_LAMP::PrintResults(std::ostream & out) const {
	std::stringstream s;

	s << "# min. sup=" << final_support_;
	if (FLAGS_second_phase)
		s << "\tcorrection factor=" << final_closed_set_num_;
	s << std::endl;

	if (FLAGS_third_phase)
		PrintSignificantSet(s);

	out << s.str() << std::flush;
	return out;
}

std::ostream & MP_LAMP::PrintSignificantSet(std::ostream & out) const {
	std::stringstream s;

	s << "# number of significant patterns=" << significant_set_.size()
			<< std::endl;
	s
			<< "# pval (raw)    pval (corr)         freq     pos        # items items\n";
	for (std::set<SignificantSetResult, sigset_compare>::const_iterator it =
			significant_set_.begin(); it != significant_set_.end(); ++it) {

		s << "" << std::setw(16) << std::left << (*it).pval_ << std::right << ""

		<< std::setw(16) << std::left << (*it).pval_ * final_closed_set_num_
				<< std::right << "" << std::setw(8) << (*it).sup_num_ << ""
				<< std::setw(8) << (*it).pos_sup_num_ << "";
		// s << "pval (raw)="   << std::setw(16) << std::left << (*it).pval_ << std::right
		//   << "pval (corr)="  << std::setw(16) << std::left << (*it).pval_ * final_closed_set_num_ << std::right
		//   << "\tfreq=" << std::setw(8)  << (*it).sup_num_
		//   << "\tpos="  << std::setw(8)  << (*it).pos_sup_num_
		//   << "\titems";

		const int * item = (*it).set_;
		significant_stack_->Print(s, d_->ItemNames(), item);
	}

	out << s.str() << std::flush;
	return out;
}

//==============================================================================

std::ostream & MP_LAMP::PrintAggrLog(std::ostream & out) {
	std::stringstream s;

	log_.Aggregate(mpi_data_.nTotalProc_);

	long long int total_time = log_.d_.search_finish_time_
			- log_.d_.search_start_time_;

	s << std::setprecision(6) << std::setiosflags(std::ios::fixed)
			<< std::right;

	s << "# variable name     :" << "          rank 0" << "       total/max"
			<< "             avg" << std::endl;

	s << "# time              =" << std::setw(16) << total_time / GIGA // s
	<< std::setprecision(3) << std::setw(16) << total_time / MEGA // ms
	<< "(s), (ms)" << std::endl;

	s << std::setprecision(3);

	s << "# process_node_num  =" << std::setw(16) << log_.d_.process_node_num_
			<< std::setw(16) << log_.a_.process_node_num_
			// sum
			<< std::setw(16)
			<< log_.a_.process_node_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# process_node_time =" << std::setw(16)

	<< log_.d_.process_node_time_ / MEGA << std::setw(16)
			<< log_.a_.process_node_time_ / MEGA // sum
			<< std::setw(16)
			<< log_.a_.process_node_time_ / MEGA / mpi_data_.nTotalProc_ // avg
			<< "(ms)" << std::endl;
	s << "# node / second     =" << std::setw(16)

	<< (log_.d_.process_node_num_) / (log_.d_.process_node_time_ / GIGA)
			<< std::setw(16) << "---" << std::setw(16)
			<< (log_.a_.process_node_num_) / (log_.a_.process_node_time_ / GIGA) // avg
			<< std::endl;

	s << "# idle_time         =" << std::setw(16) << log_.d_.idle_time_ / MEGA
			<< std::setw(16) << log_.a_.idle_time_ / MEGA // sum
			<< std::setw(16)
			<< log_.a_.idle_time_ / MEGA / mpi_data_.nTotalProc_ // avg
			<< "(ms)" << std::endl;

	s << "# pval_table_time   =" << std::setw(16)
			<< log_.d_.pval_table_time_ / MEGA << std::setw(16)
			<< log_.a_.pval_table_time_ / MEGA // sum
			<< std::setw(16)
			<< log_.a_.pval_table_time_ / MEGA / mpi_data_.nTotalProc_ // avg
			<< "(ms)" << std::endl;

	s << "# probe_num         =" << std::setw(16) << log_.d_.probe_num_
			<< std::setw(16) << log_.a_.probe_num_ // sum
			<< std::setw(16) << log_.a_.probe_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# probe_time        =" << std::setw(16) << log_.d_.probe_time_ / MEGA
			<< std::setw(16) << log_.a_.probe_time_ / MEGA // sum
			<< std::setw(16)
			<< log_.a_.probe_time_ / MEGA / mpi_data_.nTotalProc_ // avg
			<< "(ms)" << std::endl;
	s << "# probe_time_max    =" << std::setw(16)
			<< log_.d_.probe_time_max_ / MEGA << std::setw(16)
			<< log_.a_.probe_time_max_ / MEGA // max
			<< "(ms)" << std::endl;

	s << "# preprocess_time_  =" << std::setw(16)
			<< log_.d_.preprocess_time_ / MEGA << std::setw(16)
			<< log_.a_.preprocess_time_ / MEGA // sum
			<< std::setw(16)
			<< log_.a_.preprocess_time_ / MEGA / mpi_data_.nTotalProc_ // avg
			<< "(ms)" << std::endl;

	s << "# iprobe_num        =" << std::setw(16) << log_.d_.iprobe_num_ // rank 0
			<< std::setw(16) << log_.a_.iprobe_num_ // sum
			<< std::setw(16) << log_.a_.iprobe_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	LOG(
			s << "# iprobe_time       =" << std::setw(16) << log_.d_.iprobe_time_ / KILO // rank 0
			<< std::setw(16) << log_.a_.iprobe_time_ / KILO// sum
			<< std::setw(16) << log_.a_.iprobe_time_ / KILO / mpi_data_.nTotalProc_// avg
			<< "(us)" << std::endl; s << "# iprobe_time_max   =" << std::setw(16) << log_.d_.iprobe_time_max_ / KILO << std::setw(16) << log_.a_.iprobe_time_max_ / KILO// max
			<< "(us)" << std::endl; s << "# us / (one iprobe) =" << std::setw(16) << log_.d_.iprobe_time_ / log_.d_.iprobe_num_ / KILO// rank 0
			<< std::setw(16) << log_.a_.iprobe_time_ / log_.a_.iprobe_num_ / KILO// global
			<< "(us)" << std::endl;);

	s << "# iprobe_succ_num_  =" << std::setw(16) << log_.d_.iprobe_succ_num_ // rank 0
			<< std::setw(16) << log_.a_.iprobe_succ_num_ // sum
			<< std::setw(16) << log_.a_.iprobe_succ_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	LOG(
			s << "# iprobe_succ_time  =" << std::setw(16) << log_.d_.iprobe_succ_time_ / KILO // rank 0
			<< std::setw(16) << log_.a_.iprobe_succ_time_ / KILO// sum
			<< std::setw(16) << log_.a_.iprobe_succ_time_ / KILO / mpi_data_.nTotalProc_// avg
			<< "(us)" << std::endl; s << "# iprobe_succ_time_m=" << std::setw(16) << log_.d_.iprobe_succ_time_max_ / KILO << std::setw(16) << log_.a_.iprobe_succ_time_max_ / KILO// max
			<< "(us)" << std::endl; s << "# us / (one iprobe) =" << std::setw(16) << log_.d_.iprobe_succ_time_ / log_.d_.iprobe_succ_num_ / KILO// rank 0
			<< std::setw(16) << log_.a_.iprobe_succ_time_ / log_.a_.iprobe_succ_num_ / KILO// global
			<< "(us)" << std::endl;);

	s << "# iprobe_fail_num   =" << std::setw(16) << log_.d_.iprobe_fail_num_ // rank 0
			<< std::setw(16) << log_.a_.iprobe_fail_num_ // sum
			<< std::setw(16) << log_.a_.iprobe_fail_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	LOG(
			s << "# iprobe_fail_time  =" << std::setw(16) << log_.d_.iprobe_fail_time_ / KILO // rank 0
			<< std::setw(16) << log_.a_.iprobe_fail_time_ / KILO// sum
			<< std::setw(16) << log_.a_.iprobe_fail_time_ / KILO / mpi_data_.nTotalProc_// avg
			<< "(us)" << std::endl; s << "# iprobe_fail_time_m=" << std::setw(16) << log_.d_.iprobe_fail_time_max_ / KILO << std::setw(16) << log_.a_.iprobe_fail_time_max_ / KILO// max
			<< "(us)" << std::endl; s << "# us / (one iprobe) =" << std::setw(16) << log_.d_.iprobe_fail_time_ / log_.d_.iprobe_fail_num_ / KILO// rank 0
			<< std::setw(16) << log_.a_.iprobe_fail_time_ / log_.a_.iprobe_fail_num_ / KILO// global
			<< "(us)" << std::endl;);

	s << "# recv_num          =" << std::setw(16) << log_.d_.recv_num_
			<< std::setw(16) << log_.a_.recv_num_ // sum
			<< std::setw(16) << log_.a_.recv_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	LOG(
			s << "# recv_time         =" << std::setw(16) << log_.d_.recv_time_ / MEGA << std::setw(16) << log_.a_.recv_time_ / MEGA // sum
			<< std::setw(16) << log_.a_.recv_time_ / MEGA / mpi_data_.nTotalProc_// avg
			<< "(ms)" << std::endl; s << "# recv_time_max     =" << std::setw(16) << log_.d_.recv_time_max_ / MEGA << std::setw(16) << log_.a_.recv_time_max_ / MEGA// max
			<< "(ms)" << std::endl;);

	s << "# bsend_num         =" << std::setw(16) << log_.d_.bsend_num_
			<< std::setw(16) << log_.a_.bsend_num_ // sum
			<< std::setw(16) << log_.a_.bsend_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# bsend_time        =" << std::setw(16) << log_.d_.bsend_time_ / MEGA
			<< std::setw(16) << log_.a_.bsend_time_ / MEGA // sum
			<< std::setw(16)
			<< log_.a_.bsend_time_ / MEGA / mpi_data_.nTotalProc_ // avg
			<< "(ms)" << std::endl;
	s << "# bsend_time_max    =" << std::setw(16)
			<< log_.d_.bsend_time_max_ / MEGA << std::setw(16)
			<< log_.a_.bsend_time_max_ / MEGA // max
			<< "(ms)" << std::endl;

	s << "# bcast_num         =" << std::setw(16) << log_.d_.bcast_num_
			<< std::setw(16) << log_.a_.bcast_num_ // sum
			<< std::setw(16) << log_.a_.bcast_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# bcast_time        =" << std::setw(16) << log_.d_.bcast_time_ / MEGA
			<< std::setw(16) << log_.a_.bcast_time_ / MEGA // sum
			<< std::setw(16)
			<< log_.a_.bcast_time_ / MEGA / mpi_data_.nTotalProc_ // max
			<< "(ms)" << std::endl;
	s << "# bcast_time_max    =" << std::setw(16)
			<< log_.d_.bcast_time_max_ / MEGA << std::setw(16)
			<< log_.a_.bcast_time_max_ / MEGA // max
			<< "(ms)" << std::endl;

	// s << "# accum_task_time   ="
	//   << std::setw(16) << log_.d_.accum_task_time_ / MEGA
	//   << std::setw(16) << log_.a_.accum_task_time_ / MEGA // sum
	//   << std::setw(16) << log_.a_.accum_task_time_ / MEGA / p_ // avg
	//   << "(ms)" << std::endl;
	// s << "# accum_task_num    ="
	//   << std::setw(16) << log_.d_.accum_task_num_
	//   << std::setw(16) << log_.a_.accum_task_num_ // sum
	//   << std::setw(16) << log_.a_.accum_task_num_ / p_ // avg
	//   << std::endl;
	// s << "# basic_task_time   ="
	//   << std::setw(16) << log_.d_.basic_task_time_ / MEGA
	//   << std::setw(16) << log_.a_.basic_task_time_ / MEGA // sum
	//   << std::setw(16) << log_.a_.basic_task_time_ / MEGA / p_ // avg
	//   << "(ms)" << std::endl;
	// s << "# basic_task_num    ="
	//   << std::setw(16) << log_.d_.basic_task_num_
	//   << std::setw(16) << log_.a_.basic_task_num_ // sum
	//   << std::setw(16) << log_.a_.basic_task_num_ / p_ // avg
	//   << std::endl;
	// s << "# control_task_time ="
	//   << std::setw(16) << log_.d_.control_task_time_ / MEGA
	//   << std::setw(16) << log_.a_.control_task_time_ / MEGA // sum
	//   << std::setw(16) << log_.a_.control_task_time_ / MEGA / p_ // avg
	//   << "(ms)" << std::endl;
	// s << "# control_task_num  ="
	//   << std::setw(16) << log_.d_.control_task_num_
	//   << std::setw(16) << log_.a_.control_task_num_ // sum
	//   << std::setw(16) << log_.a_.control_task_num_ / p_ // avg
	//   << std::endl;

	s << "# dtd_phase_num_    =" << std::setw(16) << log_.d_.dtd_phase_num_
			<< std::endl;
	s << "# dtd_phase_per_sec_=" << std::setw(16) << log_.d_.dtd_phase_per_sec_
			<< std::endl;
	s << "# dtd_accum_ph_num_ =" << std::setw(16)
			<< log_.d_.dtd_accum_phase_num_ << std::endl;
	s << "# dtd_ac_ph_per_sec_=" << std::setw(16)
			<< log_.d_.dtd_accum_phase_per_sec_ << std::endl;

	// s << "# dtd_request_num   ="
	//   << std::setw(16) << log_.d_.dtd_request_num_
	//   << std::endl;
	// s << "# dtd_reply_num     ="
	//   << std::setw(16) << log_.d_.dtd_reply_num_
	//   << std::endl;

	// s << "# accum_phase_num_  ="
	//   << std::setw(16) << log_.d_.accum_phase_num_
	//   << std::endl;
	// s << "# accum_ph_per_sec_ ="
	//   << std::setw(16) << log_.d_.accum_phase_per_sec_
	//   << std::endl;

	s << "# lfl_steal_num     =" << std::setw(16) << log_.d_.lifeline_steal_num_
			<< std::setw(16) << log_.a_.lifeline_steal_num_
			// sum
			<< std::setw(16)
			<< log_.a_.lifeline_steal_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# lfl_nodes_recv    =" << std::setw(16)

	<< log_.d_.lifeline_nodes_received_ << std::setw(16)
			<< log_.a_.lifeline_nodes_received_
			// sum
			<< std::setw(16)
			<< log_.a_.lifeline_nodes_received_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# steal_num         =" << std::setw(16) << log_.d_.steal_num_
			<< std::setw(16) << log_.a_.steal_num_ // sum
			<< std::setw(16) << log_.a_.steal_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# nodes_received    =" << std::setw(16) << log_.d_.nodes_received_
			<< std::setw(16) << log_.a_.nodes_received_ // sum
			<< std::setw(16) << log_.a_.nodes_received_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;

	s << "# lfl_given_num     =" << std::setw(16) << log_.d_.lifeline_given_num_
			<< std::setw(16) << log_.a_.lifeline_given_num_
			// sum
			<< std::setw(16)
			<< log_.a_.lifeline_given_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# lfl_nodes_given   =" << std::setw(16)
			<< log_.d_.lifeline_nodes_given_ << std::setw(16)
			<< log_.a_.lifeline_nodes_given_
			// sum
			<< std::setw(16)
			<< log_.a_.lifeline_nodes_given_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# given_num         =" << std::setw(16) << log_.d_.given_num_
			<< std::setw(16) << log_.a_.given_num_ // sum
			<< std::setw(16) << log_.a_.given_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# nodes_given       =" << std::setw(16) << log_.d_.nodes_given_
			<< std::setw(16) << log_.a_.nodes_given_ // sum
			<< std::setw(16) << log_.a_.nodes_given_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;

	s << "# node_stack_max_itm=" << std::setw(16) << log_.d_.node_stack_max_itm_ // rank 0
			<< std::setw(16) << log_.a_.node_stack_max_itm_ // global
			<< std::endl;
	s << "# give_stack_max_itm=" << std::setw(16) << log_.d_.give_stack_max_itm_ // rank 0
			<< std::setw(16) << log_.a_.give_stack_max_itm_ // global
			<< std::endl;

	s << "# node_stack_max_cap=" << std::setw(16) << log_.d_.node_stack_max_cap_ // rank 0
			<< std::setw(16) << log_.a_.node_stack_max_cap_ // global
			<< std::endl;
	s << "# give_stack_max_cap=" << std::setw(16) << log_.d_.give_stack_max_cap_ // rank 0
			<< std::setw(16) << log_.a_.give_stack_max_cap_ // global
			<< std::endl;

	s << "# cleared_tasks_    =" << std::setw(16) << log_.d_.cleared_tasks_
			<< std::setw(16) << log_.a_.cleared_tasks_ // sum
			<< std::setw(16) << log_.a_.cleared_tasks_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;

	for (int l = 1; l <= final_support_ + 1; l++) {
		s << "# " << "lambda=" << l << "\tclosed_set_num[" << std::setw(3) << l
				<< "]=" << std::setw(12) << accum_array_[l] << "\tcs_thr["
				<< std::setw(3) << l << "]=" << std::setw(16) << cs_thr_[l]
				<< "\tpmin_thr[" << std::setw(3) << l - 1 << "]="
				<< std::setw(12) << d_->PMin(lambda_ - 1) << std::endl;
	}

	out << s.str() << std::flush;
	return out;
}

std::ostream & MP_LAMP::PrintPLog(std::ostream & out) {
	std::stringstream s;

	s << "# periodic log of node stack capacity" << std::endl;
	s << "# phase nano_sec seconds lambda capacity" << std::endl;
	for (std::size_t i = 0; i < log_.plog_.size(); i++) {
		s << "# " << std::setw(1) << log_.plog_[i].phase_ << std::setw(12)

		<< log_.plog_[i].seconds_ << std::setw(6)
				<< (int) (log_.plog_[i].seconds_ / 1000000000) << std::setw(5)
				<< log_.plog_[i].lambda_ << " " << std::setw(13)
				<< log_.plog_[i].capacity_ << std::endl;
	}

	out << s.str() << std::flush;
	return out;
}

std::ostream & MP_LAMP::PrintAggrPLog(std::ostream & out) {
	std::stringstream s;

	s << "# periodic log of node stack capacity" << std::endl;
	s << "# phase nano_sec seconds lambda min max mean sd" << std::endl;
	for (int si = 0; si < log_.sec_max_; si++) {
		long long int sum = 0ll;
		long long int max = -1;
		long long int min = std::numeric_limits<long long int>::max();

		for (int p = 0; p < mpi_data_.nTotalProc_; p++) {
			long long int cap =
					log_.plog_gather_buf_[p * log_.sec_max_ + si].capacity_;
			sum += cap;
			max = std::max(max, cap);
			min = std::min(min, cap);
		}
		double mean = sum / (double) (mpi_data_.nTotalProc_);
		double sq_diff_sum = 0.0;
		for (int p = 0; p < mpi_data_.nTotalProc_; p++) {
			long long int cap =
					log_.plog_gather_buf_[p * log_.sec_max_ + si].capacity_;
			double diff = cap - mean;
			sq_diff_sum += diff * diff;
		}
		double variance = sq_diff_sum / (double) mpi_data_.nTotalProc_;
		double sd = sqrt(variance);

		s << "# " << std::setw(1) << log_.plog_buf_[si].phase_ << std::setw(16)

		<< log_.plog_buf_[si].seconds_ << std::setw(6)
				<< (int) (log_.plog_buf_[si].seconds_ / 1000000000)
				<< std::setw(5) << log_.plog_buf_[si].lambda_ << " "
				<< std::setw(13) << min << " " << std::setw(13) << max
				<< std::setprecision(3) << " " << std::setw(17) << mean << " "
				<< std::setw(13) << sd << std::endl;
	}

	out << s.str() << std::flush;
	return out;
}

std::ostream & MP_LAMP::PrintLog(std::ostream & out) const {
	std::stringstream s;

	long long int total_time = log_.d_.search_finish_time_
			- log_.d_.search_start_time_;

	s << std::setprecision(6) << std::setiosflags(std::ios::fixed)
			<< std::right;

	s << "# time           =" << std::setw(16) << total_time / GIGA
			<< std::setprecision(3) << std::setw(16) << total_time / MEGA // ms
			<< "(s), (ms)" << std::endl;

	s << std::setprecision(3);

	s << "# process_node_num  =" << std::setw(16) << log_.d_.process_node_num_
			<< std::endl;
	s << "# process_node_time =" << std::setw(16)
			<< log_.d_.process_node_time_ / MEGA << "(ms)" << std::endl;
	s << "# node / second     =" << std::setw(16)

	<< (log_.d_.process_node_num_) / (log_.d_.process_node_time_ / GIGA)
			<< std::endl;

	s << "# idle_time         =" << std::setw(16) << log_.d_.idle_time_ / MEGA
			<< "(ms)" << std::endl;

	s << "# pval_table_time   =" << std::setw(16)
			<< log_.d_.pval_table_time_ / MEGA << "(ms)" << std::endl;

	s << "# probe_num         =" << std::setw(16) << log_.d_.probe_num_
			<< std::endl;
	s << "# probe_time        =" << std::setw(16) << log_.d_.probe_time_ / MEGA
			<< "(ms)" << std::endl;
	s << "# probe_time_max    =" << std::setw(16)
			<< log_.d_.probe_time_max_ / MEGA << "(ms)" << std::endl;

	s << "# preprocess_time_  =" << std::setw(16)
			<< log_.d_.preprocess_time_ / MEGA << "(ms)" << std::endl;

	s << "# iprobe_num        =" << std::setw(16) << log_.d_.iprobe_num_ // rank 0
			<< std::endl;
	LOG(
			s << "# iprobe_time       =" << std::setw(16) << log_.d_.iprobe_time_ / MEGA // rank 0
			<< "(ms)" << std::endl; s << "# iprobe_time_max   =" << std::setw(16) << log_.d_.iprobe_time_max_ / MEGA << "(ms)" << std::endl; s << "# ms / (one iprobe) =" << std::setw(16) << log_.d_.iprobe_time_ / log_.d_.iprobe_num_ / MEGA// rank 0
			<< "(ms)" << std::endl;);

	s << "# iprobe_succ_num   =" << std::setw(16) << log_.d_.iprobe_succ_num_ // rank 0
			<< std::endl;
	LOG(
			s << "# iprobe_succ_time  =" << std::setw(16) << log_.d_.iprobe_succ_time_ / MEGA // rank 0
			<< "(ms)" << std::endl; s << "# iprobe_succ_time_m=" << std::setw(16) << log_.d_.iprobe_succ_time_max_ / MEGA << "(ms)" << std::endl; s << "# ms / (one iprobe) =" << std::setw(16) << log_.d_.iprobe_succ_time_ / log_.d_.iprobe_succ_num_ / MEGA// rank 0
			<< "(ms)" << std::endl;);

	s << "# iprobe_fail_num   =" << std::setw(16) << log_.d_.iprobe_fail_num_ // rank 0
			<< std::endl;
	LOG(
			s << "# iprobe_fail_time  =" << std::setw(16) << log_.d_.iprobe_fail_time_ / MEGA // rank 0
			<< "(ms)" << std::endl; s << "# iprobe_fail_time_m=" << std::setw(16) << log_.d_.iprobe_fail_time_max_ / MEGA << "(ms)" << std::endl; s << "# ms / (one iprobe) =" << std::setw(16) << log_.d_.iprobe_fail_time_ / log_.d_.iprobe_fail_num_ / MEGA// rank 0
			<< "(ms)" << std::endl;);

	s << "# recv_num          =" << std::setw(16) << log_.d_.recv_num_
			<< std::endl;
	LOG(
			s << "# recv_time         =" << std::setw(16) << log_.d_.recv_time_ / MEGA << "(ms)" << std::endl; s << "# recv_time_max     =" << std::setw(16) << log_.d_.recv_time_max_ / MEGA << "(ms)" << std::endl;);

	s << "# bsend_num         =" << std::setw(16) << log_.d_.bsend_num_
			<< std::endl;
	s << "# bsend_time        =" << std::setw(16) << log_.d_.bsend_time_ / MEGA
			<< "(ms)" << std::endl;
	s << "# bsend_time_max    =" << std::setw(16)
			<< log_.d_.bsend_time_max_ / MEGA << "(ms)" << std::endl;

	s << "# bcast_num         =" << std::setw(16) << log_.d_.bcast_num_
			<< std::endl;
	s << "# bcast_time        =" << std::setw(16) << log_.d_.bcast_time_ / MEGA
			<< "(ms)" << std::endl;
	s << "# bcast_time_max    =" << std::setw(16)
			<< log_.d_.bcast_time_max_ / MEGA << "(ms)" << std::endl;

	// s << "# accum_task_time   ="
	//   << std::setw(16) << log_.d_.accum_task_time_ / MEGA
	//   << "(ms)" << std::endl;
	// s << "# accum_task_num    ="
	//   << std::setw(16) << log_.d_.accum_task_num_
	//   << std::endl;
	// s << "# basic_task_time   ="
	//   << std::setw(16) << log_.d_.basic_task_time_ / MEGA
	//   << "(ms)" << std::endl;
	// s << "# basic_task_num    ="
	//   << std::setw(16) << log_.d_.basic_task_num_
	//   << std::endl;
	// s << "# control_task_time ="
	//   << std::setw(16) << log_.d_.control_task_time_ / MEGA
	//   << "(ms)" << std::endl;
	// s << "# control_task_num  ="
	//   << std::setw(16) << log_.d_.control_task_num_
	//   << std::endl;

	s << "# dtd_phase_num_    =" << std::setw(16) << log_.d_.dtd_phase_num_
			<< std::endl;
	s << "# dtd_phase_per_sec_=" << std::setw(16) << log_.d_.dtd_phase_per_sec_
			<< std::endl;
	s << "# dtd_accum_ph_num_ =" << std::setw(16)
			<< log_.d_.dtd_accum_phase_num_ << std::endl;
	s << "# dtd_ac_ph_per_sec_=" << std::setw(16)
			<< log_.d_.dtd_accum_phase_per_sec_ << std::endl;

	// s << "# dtd_request_num   ="
	//   << std::setw(16) << log_.d_.dtd_request_num_
	//   << std::endl;
	// s << "# dtd_reply_num     ="
	//   << std::setw(16) << log_.d_.dtd_reply_num_
	//   << std::endl;

	// s << "# accum_phase_num_  ="
	//   << std::setw(16) << log_.d_.accum_phase_num_
	//   << std::endl;
	// s << "# accum_phase_per_sec_  ="
	//   << std::setw(16) << log_.d_.accum_phase_per_sec_
	//   << std::endl;

	s << "# lfl_steal_num     =" << std::setw(16) << log_.d_.lifeline_steal_num_
			<< std::endl;
	s << "# lfl_nodes_recv    =" << std::setw(16)
			<< log_.d_.lifeline_nodes_received_ << std::endl;
	s << "# steal_num         =" << std::setw(16) << log_.d_.steal_num_
			<< std::endl;
	s << "# nodes_received    =" << std::setw(16) << log_.d_.nodes_received_
			<< std::endl;

	s << "# lfl_given_num     =" << std::setw(16) << log_.d_.lifeline_given_num_
			<< std::endl;
	s << "# lfl_nodes_given   =" << std::setw(16)
			<< log_.d_.lifeline_nodes_given_ << std::endl;
	s << "# given_num         =" << std::setw(16) << log_.d_.given_num_
			<< std::endl;
	s << "# nodes_given       =" << std::setw(16) << log_.d_.nodes_given_
			<< std::endl;

	s << "# node_stack_max_itm=" << std::setw(16) << log_.d_.node_stack_max_itm_
			<< std::endl;
	s << "# give_stack_max_itm=" << std::setw(16) << log_.d_.give_stack_max_itm_
			<< std::endl;

	s << "# node_stack_max_cap=" << std::setw(16) << log_.d_.node_stack_max_cap_
			<< std::endl;
	s << "# give_stack_max_cap=" << std::setw(16) << log_.d_.give_stack_max_cap_
			<< std::endl;

	s << "# cleared_tasks_    =" << std::setw(16)
			<< log_.d_.cleared_tasks_ / MEGA << std::endl;

	out << s.str() << std::flush;
	return out;
}

//==============================================================================

void MP_LAMP::SetLogFileName() {
	std::stringstream ss;
	ss << FLAGS_debuglogfile << "_log";
	ss << std::setw(4) << std::setfill('0') << mpi_data_.mpiRank_ << ".txt";
	log_file_name_ = ss.str();
}

std::ofstream MP_LAMP::null_stream_;

std::ostream& MP_LAMP::D(int level, bool show_phase) {
	if (FLAGS_d == 0)
		return null_stream_;

	if (level <= FLAGS_d) {
		if (show_phase)
			lfs_ << std::setw(4) << phase_ << ": ";
		return lfs_;
	} else
		return null_stream_;
}

std::ostream& operator<<(std::ostream & out, const FixedSizeStack & st) {
	std::stringstream s;
	s << "entry:";
	for (int i = 0; i < st.Size(); i++)
		s << " " << st.data(i);
	s << std::endl;
	out << s.str() << std::flush;
	return out;
}

} // namespace lamp_search

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
