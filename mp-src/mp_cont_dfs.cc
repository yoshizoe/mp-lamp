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
#include "mp_cont_dfs.h"
#include "ParallelContinuousPM.h"

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
DEFINE_int32(give_size_max, 1024 * 1024 * 4,
		"maximum size of one give");
DECLARE_int32(freq_max); // 1024*1024*64, "stack size for holding freq sets", lamp.cc
DEFINE_int32(sig_max, 1024 * 1024 * 64,
		"stack size for holding significant sets");

DEFINE_int32(bsend_buffer_size, 1024 * 1024 * 64,
		"size of bsend buffer");

DEFINE_int32(d, 10,
		"debug level. 0: none, higher level produce more log");
DEFINE_string(debuglogfile, "d", "base filename for debug log");
DEFINE_bool(log, false, "show log");

DEFINE_int32(probe_period, 128, "probe period during process node");
DEFINE_bool(probe_period_is_ms, false,
		"true: probe period is milli sec, false: num loops");

namespace lamp_search {

const int MP_CONT_LAMP::k_echo_tree_branch = 3;

const int MP_CONT_LAMP::k_int_max = std::numeric_limits<int>::max();
// assuming (digits in long long int) > (bits of double mantissa)
const long long int MP_CONT_LAMP::k_cs_max = 1ll
		<< (std::numeric_limits<double>::digits - 1);

const int MP_CONT_LAMP::k_probe_period = 128;

MP_CONT_LAMP::MP_CONT_LAMP(ContDatabase* d, int rank, int nu_proc,
		int n, bool n_is_ms, int w, int l, int m) :
		d_(d), dtd_(k_echo_tree_branch), mpi_data_(
				FLAGS_bsend_buffer_size, rank, nu_proc, n, n_is_ms, w,
				l, m, k_echo_tree_branch, &dtd_), timer_(
				Timer::GetInstance()), give_stack_(NULL), stealer_(
				mpi_data_.nRandStealTrials_,
				mpi_data_.hypercubeDimension_), phase_(0), freq_stack_(
		NULL), significant_stack_(
		NULL), total_expand_num_(0ll), expand_num_(0ll), closed_set_num_(
				0ll), num_final_testable_patterns(0ll), final_sig_level_(
				0.0), last_bcast_was_dtd_(false) {
	printf("initializing MP_LAMP\n");
	if (FLAGS_d > 0) {
		SetLogFileName();
		DBG(lfs_.open(log_file_name_.c_str(), std::ios::out)
		;);
	}

	dtd_.Init();

	printf("dtd\n");

	InitTreeRequest(); // TODO: This should be in MPI_Data.h

	Init();

	DBG(D(1) << "lifelines:"
	;);
	DBG(for (int i = 0; i < mpi_data_.hypercubeDimension_; i++)
		D(1, false) << "\t" << mpi_data_.lifelines_[i]
		;);
	DBG(D(1, false) << std::endl
	;);

	DBG(
			D(1) << "thieves_ size=" << mpi_data_.thieves_->Size()
					<< std::endl
			;);
	DBG(
			D(1) << "mpi_data_.lifeline_thieves_ size="
					<< mpi_data_.lifeline_thieves_->Size()
					<< std::endl
			;);

	DBG(
			D(1) << "bcast_targets:" << std::endl; for (std::size_t i = 0; i < k_echo_tree_branch; i++) { D(1) << "\t[" << i << "]=" << mpi_data_.bcast_targets_[i] << std::endl; }
			// << "\t[0]=" << bcast_targets_[0]
			// << "\t[1]=" << bcast_targets_[1]
			// << "\t[2]=" << bcast_targets_[2]
			// D(1) << std::endl;
			);

	node_stack_ = new VariableLengthItemsetStack(FLAGS_stack_size);
	give_stack_ = new VariableLengthItemsetStack(FLAGS_give_size_max);
	freq_stack_ = new VariableLengthItemsetStack(FLAGS_freq_max);

	printf("initialized MP_LAMP\n");
}

void MP_CONT_LAMP::InitTreeRequest() {
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

	if (mpi_data_.mpiRank_ > 0) {
		mpi_data_.lifelines_activated_[(mpi_data_.mpiRank_ - 1)
				/ k_echo_tree_branch] = true;
		mpi_data_.bcast_source_ = (mpi_data_.mpiRank_ - 1)
				/ k_echo_tree_branch;
	}
}

void MP_CONT_LAMP::SetTreeRequest() {
	for (std::size_t i = 1; i <= k_echo_tree_branch; i++) {
		if (k_echo_tree_branch * mpi_data_.mpiRank_ + i
				< mpi_data_.nTotalProc_)
			mpi_data_.lifeline_thieves_->Push(
					k_echo_tree_branch * mpi_data_.mpiRank_ + i);
	}

	if (mpi_data_.mpiRank_ > 0)
		mpi_data_.lifelines_activated_[(mpi_data_.mpiRank_ - 1)
				/ k_echo_tree_branch] = true;
}

MP_CONT_LAMP::~MP_CONT_LAMP() {
	DBG(D(2) << "MP_LAMP destructor begin" << std::endl
	;);

	if (node_stack_)
		delete node_stack_;
	if (give_stack_)
		delete give_stack_;
	if (freq_stack_)
		delete freq_stack_;
	if (significant_stack_)
		delete significant_stack_;

//	if (d_)
//		delete d_;
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
//
//void MP_CONT_LAMP::InitDatabaseRoot(std::istream & is1,
//		std::istream & is2) {
//	uint64 * data = NULL;
//	uint64 * positive = NULL;
//	boost::array<int, 5> counters; // nu_bits, nu_trans, nu_items, max_item_in_transaction
//	counters.assign(-1);
//
//	int nu_trans = 0;
//	int nu_items = 0;
//	int nu_pos_total = 0;
//	int max_item_in_transaction = 0;
//
//	std::vector<std::string> * item_names = NULL;
//	std::vector<std::string> * transaction_names = NULL;
//
//	DatabaseReader<uint64> reader;
//
//	assert(mpi_data_.mpiRank_ == 0);
//	{
//		item_names = new std::vector<std::string>;
//		transaction_names = new std::vector<std::string>;
//
//		if (FLAGS_lcm) {
//			reader.ReadFilesLCM(&bsh_, is1, &data, &nu_trans,
//					&nu_items, is2, &positive, &nu_pos_total,
//					item_names, &max_item_in_transaction);
//		} else {
//			reader.ReadFiles(&bsh_, is1, &data, &nu_trans, &nu_items,
//					is2, &positive, &nu_pos_total, item_names,
//					transaction_names, &max_item_in_transaction);
//		}
//
//		counters[0] = (int) (bsh_->nu_bits);
//		counters[1] = nu_trans;
//		counters[2] = nu_items;
//		counters[3] = nu_pos_total;
//		counters[4] = max_item_in_transaction;
//
//		CallBcast(&counters, 5, MPI_INT);
//	}
//
//	CallBcast(data, bsh_->NewArraySize(nu_items),
//	MPI_UNSIGNED_LONG_LONG);
//	CallBcast(positive, bsh_->NuBlocks(), MPI_UNSIGNED_LONG_LONG);
//
//	long long int start_time = timer_->Elapsed();
//	d_ = new Database<uint64>(bsh_, data, nu_trans, nu_items,
//			positive, nu_pos_total, max_item_in_transaction,
//			item_names, transaction_names);
//	log_.d_.pval_table_time_ = timer_->Elapsed() - start_time;
//
////	g_ = new LampGraph<uint64>(*d_);
//
//	// Lambda max is initially the number of total transactions in the database.
//	lambda_max_ = d_->MaxX(); // used for getMinSup
//	// D() << "max_x=lambda_max=" << lambda_max_ << std::endl;
//	// D() << "pos_total=" << d_->PosTotal() << std::endl;
//	// d_->DumpItems(D(false));
//	// d_->DumpPMinTable(D(false));
//
//	double *pmin_thr_ = new double[lambda_max_ + 1]; // used for getMinSup
//	cs_thr_ = new long long int[lambda_max_ + 1]; // used for getMinSup
//
//	dtd_accum_array_base_ = new long long int[lambda_max_ + 4];
//	dtd_accum_recv_base_ = new long long int[lambda_max_ + 4];
//	accum_array_ = &(dtd_accum_array_base_[3]); // TODO: ???
//	accum_recv_ = &(dtd_accum_recv_base_[3]);
//
//	node_stack_ = new VariableLengthItemsetStack(FLAGS_stack_size);
//	give_stack_ = new VariableLengthItemsetStack(FLAGS_give_size_max);
//	freq_stack_ = new VariableLengthItemsetStack(FLAGS_freq_max);
//	// node_stack_ =  new VariableLengthItemsetStack(FLAGS_stack_size, lambda_max_);
//	// give_stack_ =  new VariableLengthItemsetStack(FLAGS_give_size_max, lambda_max_);
//	// freq_stack_ = new VariableLengthItemsetStack(FLAGS_freq_max, lambda_max_);
//
//	for (int i = 0; i <= lambda_max_; i++)
//		pmin_thr_[i] = d_->PMin(i);
//
//	cs_thr_[0] = 0ll; // should not be used ???
//	for (int i = 1; i <= lambda_max_; i++) {
//		cs_thr_[i] = (long long int) (std::min(
//				std::floor(FLAGS_a / pmin_thr_[i - 1]),
//				(double) (k_cs_max)));
//		DBG(
//				D(2) << "cs_thr[" << i << "]=\t" << cs_thr_[i]
//						<< "\tpmin_thr=" << pmin_thr_[i - 1]
//						<< std::endl
//				;);
//	}
//
//	for (int i = 0; i <= lambda_max_; i++)
//		accum_array_[i] = 0ll;
//
//	sup_buf_ = bsh_->New();
//	child_sup_buf_ = bsh_->New();
//
//	delete pmin_thr_;
//}

void MP_CONT_LAMP::Init() {
	// todo: move accum cs count variable to inner class
	mpi_data_.echo_waiting_ = false;
	for (int i = 0; i < k_echo_tree_branch; i++)
		mpi_data_.accum_flag_[i] = false;

	dtd_.Init();

	log_.Init();

	stealer_.Init();
	mpi_data_.waiting_ = false;

	last_bcast_was_dtd_ = false;
}

void MP_CONT_LAMP::CheckPoint() {
	DBG(D(4) << "CheckPoint Started" << std::endl
	;);
	assert(mpi_data_.echo_waiting_ == false);

	// does this hold?
	for (int i = 0; i < k_echo_tree_branch; i++)
		assert(mpi_data_.accum_flag_[i] == false);

	dtd_.CheckPoint();
	stealer_.Init();
	mpi_data_.waiting_ = false;

	mpi_data_.thieves_->Clear();
	mpi_data_.lifeline_thieves_->Clear();
	for (int pi = 0; pi < mpi_data_.nTotalProc_; pi++)
		mpi_data_.lifelines_activated_[pi] = false;

	SetTreeRequest();
	DBG(D(4) << "CheckPoint Finished" << std::endl
	;);
}

void MP_CONT_LAMP::ClearTasks() {
	MPI_Status probe_status, recv_status;
	int data_count, src, tag;
	int flag;
	int error;

	while (true) {
		error = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG,
		MPI_COMM_WORLD, &flag, &probe_status);
		if (error != MPI_SUCCESS) {
			DBG(
					D(1) << "error in MPI_Iprobe in ClearTasks: "
							<< error << std::endl
					;);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		if (!flag)
			break;

		log_.d_.cleared_tasks_++;

		error = MPI_Get_count(&probe_status, MPI_INT, &data_count);
		if (error != MPI_SUCCESS) {
			DBG(
					D(1) << "error in MPI_Get_count in ClearTasks: "
							<< error << std::endl
					;);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		src = probe_status.MPI_SOURCE;
		tag = probe_status.MPI_TAG;

		error = MPI_Recv(give_stack_->Stack(), data_count, MPI_INT,
				src, tag,
				MPI_COMM_WORLD, &recv_status);
		if (error != MPI_SUCCESS) {
			DBG(
					D(1) << "error in MPI_Recv in ClearTasks: "
							<< error << std::endl
					;);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		give_stack_->Clear(); // drop messages
	}
}

//==============================================================================

void MP_CONT_LAMP::Search() {
	printf("MP_CONT_LAMP::Search\n");

	// TODO: itemset_buf_ is not allocated.
	TreeSearchData* treesearch_data_ = new TreeSearchData(node_stack_,
			give_stack_, &stealer_, itemset_buf_, &log_, timer_);
//	BinaryPatternMiningData* bpm_data_ = new BinaryPatternMiningData(
//			d_, bsh_, sup_buf_, child_sup_buf_);
	ContinuousPatternMiningData* cpm_data_ =
			new ContinuousPatternMiningData(d_);\
	ParallelContinuousPM* psearch = new ParallelContinuousPM(
			cpm_data_, mpi_data_, treesearch_data_, FLAGS_a, &log_,
			timer_, lfs_);
	GetTestableData* gettestable_data_;
	GetContSignificantData* getsignificant_data_;

	log_.d_.search_start_time_ = timer_->Elapsed();
	total_expand_num_ = 0ll;

	expand_num_ = 0ll;
	closed_set_num_ = 0ll;

	double lambda_ = 0.0;

	// --------
	// prepare phase 2
	phase_ = 2;
	CheckPoint();

//	CallBcast(&lambda_, 1, MPI_INT); // Rank-0 process broadcasts its lambda to all the other processes
//	final_support_ = lambda_;

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
//		node_stack_->SetSup(root_itemset, lambda_max_);
		node_stack_->PushPostNoSort();
	}

//	double int_sig_lev = 0.0;
//	if (mpi_data_.mpiRank_ == 0) {
//		int_sig_lev = GetInterimSigLevel(lambda_);
//	}
	// todo: reduce expand_num_

	{
//		if (mpi_data_.mpiRank_ == 0 && FLAGS_show_progress) {
//			std::cout << "# " << "2nd phase start\n";
//			std::cout << "# " << "lambda=" << lambda_
//					<< "\tint_sig_lev=" << int_sig_lev
//					<< "\telapsed_time="
//					<< (timer_->Elapsed() - log_.d_.search_start_time_)
//							/ GIGA << std::endl;
//		}
		DBG(D(2) << "---------------" << std::endl
		;);
		DBG(D(1) << "2nd phase start" << std::endl
		;);
		DBG(D(2) << "---------------" << std::endl
		;);
//		DBG(
//				D(2) << "lambda=" << lambda_ << "\tint_sig_lev="
//						<< int_sig_lev << "\telapsed_time="
//						<< (timer_->Elapsed()
//								- log_.d_.search_start_time_) / GIGA
//						<< std::endl
//				;);

//		CallBcast(&int_sig_lev, 1, MPI_DOUBLE);
		double sig_level_ = FLAGS_a;
		gettestable_data_ = new GetTestableData(lambda_, freq_stack_,
				&freq_map_, sig_level_);

		psearch->GetTestablePatterns(gettestable_data_);
//		GetTestablePatterns(mpi_data_, treesearch_data_, gettestable_data_);
//		freq_stack_ = gettestable_data_->freq_stack_; // No need: two are the same.
		//		freq_map_ = gettestable_data_->freq_map_; // No need?
		sig_level_ = gettestable_data_->sig_level_;
		assert(freq_stack_ == gettestable_data_->freq_stack_);
		assert(
				freq_map_.size()
						== gettestable_data_->freq_map_->size());
		closed_set_num_ = freq_map_.size();
	}

	DBG(D(1) << "closed_set_num=" << closed_set_num_ << std::endl
	;);

	long long int closed_set_num_reduced;
	MPI_Reduce(&closed_set_num_, &closed_set_num_reduced, 1,
	MPI_LONG_LONG_INT,
	MPI_SUM, 0, MPI_COMM_WORLD); // error?

	DBG(
			if (mpi_data_.mpiRank_ == 0)
				D(1) << "closed_set_num_reduced="
						<< closed_set_num_reduced << std::endl
				;);
	if (mpi_data_.mpiRank_ == 0)
		num_final_testable_patterns = closed_set_num_reduced;

	log_.d_.dtd_phase_per_sec_ =
			(double) (log_.d_.dtd_phase_num_)
					/ ((timer_->Elapsed() - log_.d_.search_start_time_)
							/ GIGA);

	MPI_Barrier( MPI_COMM_WORLD);
	log_.FinishPeriodicLog();

	{
		if (mpi_data_.mpiRank_ == 0 && FLAGS_show_progress) {
			std::cout << "# " << "2nd phase end\n";
			std::cout << "# " << "closed_set_num=" << std::setw(12)
					<< num_final_testable_patterns << "\tsig_lev="
					<< (FLAGS_a / num_final_testable_patterns)
					<< "\tnum_expand=" << std::setw(12) << expand_num_
					<< "\telapsed_time="
					<< (timer_->Elapsed() - log_.d_.search_start_time_)
							/ GIGA << std::endl;
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
	significant_stack_ = new VariableLengthItemsetStack(
			FLAGS_sig_max);
	// significant_stack_ = new VariableLengthItemsetStack(FLAGS_sig_max, lambda_max_);

//	final_sig_level_ = FLAGS_a / num_final_testable_patterns;
	final_sig_level_ = gettestable_data_->sig_level_;
	CallBcast(&final_sig_level_, 1, MPI_DOUBLE);

	// TODO: ?
//	printf("freq_stack_->PrintAll()\n");
//	freq_stack_->PrintAll(std::cout);

	{
		// TODO: FLAGS_a should be final_sig_level. For testing purpose we put alpha.
		getsignificant_data_ = new GetContSignificantData(freq_stack_,
				&freq_map_, final_sig_level_, significant_stack_,
				&significant_set_);
//		getsignificant_data_ = new GetSignificantData(freq_stack_,
//				&freq_map_, final_sig_level_, significant_stack_,
//				&significant_set_);
//		GetSignificantPatterns(mpi_data_, getsignificant_data_);
		psearch->GetSignificantPatterns(mpi_data_,
				getsignificant_data_);
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
					<< (log_.d_.search_finish_time_
							- log_.d_.search_start_time_) / GIGA
					<< std::endl;
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
	delete psearch;
	delete treesearch_data_;
	delete cpm_data_;
	delete gettestable_data_;
	delete getsignificant_data_;

}

// TODO: For continuous features we cannot get significance.
double MP_CONT_LAMP::GetInterimSigLevel(int lambda) const {
	return FLAGS_a;
//	long long int csnum = accum_array_[lambda];
//	double lv;
//	if (csnum > 0)
//		lv = FLAGS_a / (double) csnum;
//	else
//		lv = FLAGS_a;
//
//	return lv;
}

// TODO: Put inside of ParallelContinuousPM
void MP_CONT_LAMP::SortSignificantSets() {
	printf("SortSignificantSets: %d items\n",
			significant_stack_->NuItemset());

	int * set = significant_stack_->FirstItemset();
	while (set != NULL) {
		// calculate support from set
		std::vector<int> itemset = significant_stack_->getItems(set);
		std::vector<double> freqs = d_->GetFreqArray(itemset);

		double freq = d_->GetFreq(freqs);
		double pos_freq = d_->GetPositiveFreq(freqs);
		double pval = d_->CalculatePValue(freq, pos_freq);

		significant_set_.insert(
				ContSignificantSetResult(pval, set, freq, pos_freq,
						significant_stack_));
		set = significant_stack_->NextItemset(set);
	}
	printf("%d significant sets\n", significant_set_.size());
}

// TODO: Ideally, this should also be hidden in other class.
int MP_CONT_LAMP::CallBcast(void * buffer, int data_count,
		MPI_Datatype type) {
	long long int start_time;
	long long int end_time;
	log_.d_.bcast_num_++;
	start_time = timer_->Elapsed();

	int error = MPI_Bcast(buffer, data_count, type, 0,
	MPI_COMM_WORLD);

	end_time = timer_->Elapsed();
	log_.d_.bcast_time_ += end_time - start_time;
	log_.d_.bcast_time_max_ = std::max(end_time - start_time,
			log_.d_.bcast_time_max_);
	return error;
}

//==============================================================================

std::ostream & MP_CONT_LAMP::PrintDBInfo(std::ostream & out) const {
	std::stringstream s;

	s << "# ";
	d_->ShowInfo();

	out << s.str() << std::flush;
	return out;
}

std::ostream & MP_CONT_LAMP::PrintResults(std::ostream & out) const {
	std::stringstream s;

	if (FLAGS_second_phase)
		s << "\tcorrection factor=" << num_final_testable_patterns;
	s << std::endl;

	if (FLAGS_third_phase)
		PrintSignificantSet(s);

	out << s.str() << std::flush;
	return out;
}

std::ostream & MP_CONT_LAMP::PrintSignificantSet(
		std::ostream & out) const {
	std::stringstream s;
//	printf("PrintSignificantSet not implemented\n");

	s << "# number of significant patterns="
			<< significant_set_.size() << std::endl;
	s
			<< "# pval (raw)    pval (corr)         freq     pos        # items items\n";
	for (std::set<ContSignificantSetResult, cont_sigset_compare>::const_iterator it =
			significant_set_.begin(); it != significant_set_.end();
			++it) {

		s << "" << std::setw(16) << std::left << (*it).pval_
				<< std::right << ""

				<< std::setw(16) << std::left
				<< (*it).pval_ * num_final_testable_patterns << std::right
				<< "" << std::setw(12) << (*it).sup_num_ << ""
				<< std::setw(12) << (*it).pos_sup_num_ << " ";
//				<< std::endl;
//		 s << "pval (raw)="   << std::setw(16) << std::left << (*it).pval_ << std::right
//		   << "pval (corr)="  << std::setw(16) << std::left << (*it).pval_ * final_closed_set_num_ << std::right
//		   << "\tfreq=" << std::setw(8)  << (*it).sup_num_
//		   << "\tpos="  << std::setw(8)  << (*it).pos_sup_num_
//		   << "\titems";

		const int * item = (*it).set_;
		significant_stack_->Print(s, item);
	}

	out << s.str() << std::flush;
	return out;
}

//==============================================================================

std::ostream & MP_CONT_LAMP::PrintAggrLog(std::ostream & out) {
	std::stringstream s;

	log_.Aggregate(mpi_data_.nTotalProc_);

	long long int total_time = log_.d_.search_finish_time_
			- log_.d_.search_start_time_;

	s << std::setprecision(6) << std::setiosflags(std::ios::fixed)
			<< std::right;

	s << "# variable name     :" << "          rank 0"
			<< "       total/max" << "             avg" << std::endl;

	s << "# time              =" << std::setw(16) << total_time / GIGA // s
	<< std::setprecision(3) << std::setw(16) << total_time / MEGA // ms
	<< "(s), (ms)" << std::endl;

	s << std::setprecision(3);

	s << "# process_node_num  =" << std::setw(16)
			<< log_.d_.process_node_num_ << std::setw(16)
			<< log_.a_.process_node_num_
			// sum
			<< std::setw(16)
			<< log_.a_.process_node_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# process_node_time =" << std::setw(16)

	<< log_.d_.process_node_time_ / MEGA << std::setw(16)
			<< log_.a_.process_node_time_ / MEGA // sum
			<< std::setw(16)
			<< log_.a_.process_node_time_ / MEGA
					/ mpi_data_.nTotalProc_ // avg
			<< "(ms)" << std::endl;
	s << "# node / second     =" << std::setw(16)

			<< (log_.d_.process_node_num_)
					/ (log_.d_.process_node_time_ / GIGA)
			<< std::setw(16) << "---" << std::setw(16)
			<< (log_.a_.process_node_num_)
					/ (log_.a_.process_node_time_ / GIGA) // avg
			<< std::endl;

	s << "# idle_time         =" << std::setw(16)
			<< log_.d_.idle_time_ / MEGA << std::setw(16)
			<< log_.a_.idle_time_ / MEGA // sum
			<< std::setw(16)
			<< log_.a_.idle_time_ / MEGA / mpi_data_.nTotalProc_ // avg
			<< "(ms)" << std::endl;

	s << "# pval_table_time   =" << std::setw(16)
			<< log_.d_.pval_table_time_ / MEGA << std::setw(16)
			<< log_.a_.pval_table_time_ / MEGA // sum
			<< std::setw(16)
			<< log_.a_.pval_table_time_ / MEGA / mpi_data_.nTotalProc_ // avg
			<< "(ms)" << std::endl;

	s << "# probe_num         =" << std::setw(16)
			<< log_.d_.probe_num_ << std::setw(16)
			<< log_.a_.probe_num_
			// sum
			<< std::setw(16)
			<< log_.a_.probe_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# probe_time        =" << std::setw(16)
			<< log_.d_.probe_time_ / MEGA << std::setw(16)
			<< log_.a_.probe_time_ / MEGA // sum
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

	s << "# iprobe_num        =" << std::setw(16)
			<< log_.d_.iprobe_num_
			// rank 0
			<< std::setw(16) << log_.a_.iprobe_num_
			// sum
			<< std::setw(16)
			<< log_.a_.iprobe_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	LOG(
			s << "# iprobe_time       =" << std::setw(16) << log_.d_.iprobe_time_ / KILO // rank 0
			<< std::setw(16) << log_.a_.iprobe_time_ / KILO// sum
			<< std::setw(16) << log_.a_.iprobe_time_ / KILO / mpi_data_.nTotalProc_// avg
			<< "(us)" << std::endl; s << "# iprobe_time_max   =" << std::setw(16) << log_.d_.iprobe_time_max_ / KILO << std::setw(16) << log_.a_.iprobe_time_max_ / KILO// max
			<< "(us)" << std::endl; s << "# us / (one iprobe) =" << std::setw(16) << log_.d_.iprobe_time_ / log_.d_.iprobe_num_ / KILO// rank 0
			<< std::setw(16) << log_.a_.iprobe_time_ / log_.a_.iprobe_num_ / KILO// global
			<< "(us)" << std::endl;);

	s << "# iprobe_succ_num_  =" << std::setw(16)
			<< log_.d_.iprobe_succ_num_
			// rank 0
			<< std::setw(16) << log_.a_.iprobe_succ_num_
			// sum
			<< std::setw(16)
			<< log_.a_.iprobe_succ_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	LOG(
			s << "# iprobe_succ_time  =" << std::setw(16) << log_.d_.iprobe_succ_time_ / KILO // rank 0
			<< std::setw(16) << log_.a_.iprobe_succ_time_ / KILO// sum
			<< std::setw(16) << log_.a_.iprobe_succ_time_ / KILO / mpi_data_.nTotalProc_// avg
			<< "(us)" << std::endl; s << "# iprobe_succ_time_m=" << std::setw(16) << log_.d_.iprobe_succ_time_max_ / KILO << std::setw(16) << log_.a_.iprobe_succ_time_max_ / KILO// max
			<< "(us)" << std::endl; s << "# us / (one iprobe) =" << std::setw(16) << log_.d_.iprobe_succ_time_ / log_.d_.iprobe_succ_num_ / KILO// rank 0
			<< std::setw(16) << log_.a_.iprobe_succ_time_ / log_.a_.iprobe_succ_num_ / KILO// global
			<< "(us)" << std::endl;);

	s << "# iprobe_fail_num   =" << std::setw(16)
			<< log_.d_.iprobe_fail_num_
			// rank 0
			<< std::setw(16) << log_.a_.iprobe_fail_num_
			// sum
			<< std::setw(16)
			<< log_.a_.iprobe_fail_num_ / mpi_data_.nTotalProc_ // avg
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
			<< std::setw(16) << log_.a_.recv_num_
			// sum
			<< std::setw(16)
			<< log_.a_.recv_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	LOG(
			s << "# recv_time         =" << std::setw(16) << log_.d_.recv_time_ / MEGA << std::setw(16) << log_.a_.recv_time_ / MEGA // sum
			<< std::setw(16) << log_.a_.recv_time_ / MEGA / mpi_data_.nTotalProc_// avg
			<< "(ms)" << std::endl; s << "# recv_time_max     =" << std::setw(16) << log_.d_.recv_time_max_ / MEGA << std::setw(16) << log_.a_.recv_time_max_ / MEGA// max
			<< "(ms)" << std::endl;);

	s << "# bsend_num         =" << std::setw(16)
			<< log_.d_.bsend_num_ << std::setw(16)
			<< log_.a_.bsend_num_
			// sum
			<< std::setw(16)
			<< log_.a_.bsend_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# bsend_time        =" << std::setw(16)
			<< log_.d_.bsend_time_ / MEGA << std::setw(16)
			<< log_.a_.bsend_time_ / MEGA // sum
			<< std::setw(16)
			<< log_.a_.bsend_time_ / MEGA / mpi_data_.nTotalProc_ // avg
			<< "(ms)" << std::endl;
	s << "# bsend_time_max    =" << std::setw(16)
			<< log_.d_.bsend_time_max_ / MEGA << std::setw(16)
			<< log_.a_.bsend_time_max_ / MEGA // max
			<< "(ms)" << std::endl;

	s << "# bcast_num         =" << std::setw(16)
			<< log_.d_.bcast_num_ << std::setw(16)
			<< log_.a_.bcast_num_
			// sum
			<< std::setw(16)
			<< log_.a_.bcast_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# bcast_time        =" << std::setw(16)
			<< log_.d_.bcast_time_ / MEGA << std::setw(16)
			<< log_.a_.bcast_time_ / MEGA // sum
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

	s << "# dtd_phase_num_    =" << std::setw(16)
			<< log_.d_.dtd_phase_num_ << std::endl;
	s << "# dtd_phase_per_sec_=" << std::setw(16)
			<< log_.d_.dtd_phase_per_sec_ << std::endl;
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

	s << "# lfl_steal_num     =" << std::setw(16)
			<< log_.d_.lifeline_steal_num_ << std::setw(16)
			<< log_.a_.lifeline_steal_num_
			// sum
			<< std::setw(16)
			<< log_.a_.lifeline_steal_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# lfl_nodes_recv    =" << std::setw(16)

	<< log_.d_.lifeline_nodes_received_ << std::setw(16)
			<< log_.a_.lifeline_nodes_received_
			// sum
			<< std::setw(16)
			<< log_.a_.lifeline_nodes_received_
					/ mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# steal_num         =" << std::setw(16)
			<< log_.d_.steal_num_ << std::setw(16)
			<< log_.a_.steal_num_
			// sum
			<< std::setw(16)
			<< log_.a_.steal_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# nodes_received    =" << std::setw(16)
			<< log_.d_.nodes_received_ << std::setw(16)
			<< log_.a_.nodes_received_
			// sum
			<< std::setw(16)
			<< log_.a_.nodes_received_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;

	s << "# lfl_given_num     =" << std::setw(16)
			<< log_.d_.lifeline_given_num_ << std::setw(16)
			<< log_.a_.lifeline_given_num_
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
	s << "# given_num         =" << std::setw(16)
			<< log_.d_.given_num_ << std::setw(16)
			<< log_.a_.given_num_
			// sum
			<< std::setw(16)
			<< log_.a_.given_num_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;
	s << "# nodes_given       =" << std::setw(16)
			<< log_.d_.nodes_given_ << std::setw(16)
			<< log_.a_.nodes_given_
			// sum
			<< std::setw(16)
			<< log_.a_.nodes_given_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;

	s << "# node_stack_max_itm=" << std::setw(16)
			<< log_.d_.node_stack_max_itm_ // rank 0
			<< std::setw(16) << log_.a_.node_stack_max_itm_ // global
			<< std::endl;
	s << "# give_stack_max_itm=" << std::setw(16)
			<< log_.d_.give_stack_max_itm_ // rank 0
			<< std::setw(16) << log_.a_.give_stack_max_itm_ // global
			<< std::endl;

	s << "# node_stack_max_cap=" << std::setw(16)
			<< log_.d_.node_stack_max_cap_ // rank 0
			<< std::setw(16) << log_.a_.node_stack_max_cap_ // global
			<< std::endl;
	s << "# give_stack_max_cap=" << std::setw(16)
			<< log_.d_.give_stack_max_cap_ // rank 0
			<< std::setw(16) << log_.a_.give_stack_max_cap_ // global
			<< std::endl;

	s << "# cleared_tasks_    =" << std::setw(16)
			<< log_.d_.cleared_tasks_ << std::setw(16)
			<< log_.a_.cleared_tasks_
			// sum
			<< std::setw(16)
			<< log_.a_.cleared_tasks_ / mpi_data_.nTotalProc_ // avg
			<< std::endl;

//	for (int l = 1; l <= final_support_ + 1; l++) {
//		s << "# " << "lambda=" << l << "\tclosed_set_num["
//				<< std::setw(3) << l << "]=" << std::setw(12)
//				<< accum_array_[l] << "\tcs_thr[" << std::setw(3) << l
//				<< "]=" << std::setw(16) << cs_thr_[l]
//				<< "\tpmin_thr[" << std::setw(3) << l - 1 << "]="
//				<< std::setw(12) << d_->PMin(lambda_ - 1)
//				<< std::endl;
//	}

	out << s.str() << std::flush;
	return out;
}

std::ostream & MP_CONT_LAMP::PrintPLog(std::ostream & out) {
	std::stringstream s;

	s << "# periodic log of node stack capacity" << std::endl;
	s << "# phase nano_sec seconds lambda capacity" << std::endl;
	for (std::size_t i = 0; i < log_.plog_.size(); i++) {
		s << "# " << std::setw(1) << log_.plog_[i].phase_
				<< std::setw(12)

				<< log_.plog_[i].seconds_ << std::setw(6)
				<< (int) (log_.plog_[i].seconds_ / 1000000000)
				<< std::setw(5) << log_.plog_[i].lambda_ << " "
				<< std::setw(13) << log_.plog_[i].capacity_
				<< std::endl;
	}

	out << s.str() << std::flush;
	return out;
}

std::ostream & MP_CONT_LAMP::PrintAggrPLog(std::ostream & out) {
	std::stringstream s;

	s << "# periodic log of node stack capacity" << std::endl;
	s << "# phase nano_sec seconds lambda min max mean sd"
			<< std::endl;
	for (int si = 0; si < log_.sec_max_; si++) {
		long long int sum = 0ll;
		long long int max = -1;
		long long int min = std::numeric_limits<long long int>::max();

		for (int p = 0; p < mpi_data_.nTotalProc_; p++) {
			long long int cap = log_.plog_gather_buf_[p
					* log_.sec_max_ + si].capacity_;
			sum += cap;
			max = std::max(max, cap);
			min = std::min(min, cap);
		}
		double mean = sum / (double) (mpi_data_.nTotalProc_);
		double sq_diff_sum = 0.0;
		for (int p = 0; p < mpi_data_.nTotalProc_; p++) {
			long long int cap = log_.plog_gather_buf_[p
					* log_.sec_max_ + si].capacity_;
			double diff = cap - mean;
			sq_diff_sum += diff * diff;
		}
		double variance = sq_diff_sum
				/ (double) mpi_data_.nTotalProc_;
		double sd = sqrt(variance);

		s << "# " << std::setw(1) << log_.plog_buf_[si].phase_
				<< std::setw(16)

				<< log_.plog_buf_[si].seconds_ << std::setw(6)
				<< (int) (log_.plog_buf_[si].seconds_ / 1000000000)
				<< std::setw(5) << log_.plog_buf_[si].lambda_ << " "
				<< std::setw(13) << min << " " << std::setw(13) << max
				<< std::setprecision(3) << " " << std::setw(17)
				<< mean << " " << std::setw(13) << sd << std::endl;
	}

	out << s.str() << std::flush;
	return out;
}

std::ostream & MP_CONT_LAMP::PrintLog(std::ostream & out) const {
	std::stringstream s;

	long long int total_time = log_.d_.search_finish_time_
			- log_.d_.search_start_time_;

	s << std::setprecision(6) << std::setiosflags(std::ios::fixed)
			<< std::right;

	s << "# time           =" << std::setw(16) << total_time / GIGA
			<< std::setprecision(3) << std::setw(16)
			<< total_time / MEGA // ms
			<< "(s), (ms)" << std::endl;

	s << std::setprecision(3);

	s << "# process_node_num  =" << std::setw(16)
			<< log_.d_.process_node_num_ << std::endl;
	s << "# process_node_time =" << std::setw(16)
			<< log_.d_.process_node_time_ / MEGA << "(ms)"
			<< std::endl;
	s << "# node / second     =" << std::setw(16)

			<< (log_.d_.process_node_num_)
					/ (log_.d_.process_node_time_ / GIGA)
			<< std::endl;

	s << "# idle_time         =" << std::setw(16)
			<< log_.d_.idle_time_ / MEGA << "(ms)" << std::endl;

	s << "# pval_table_time   =" << std::setw(16)
			<< log_.d_.pval_table_time_ / MEGA << "(ms)" << std::endl;

	s << "# probe_num         =" << std::setw(16)
			<< log_.d_.probe_num_ << std::endl;
	s << "# probe_time        =" << std::setw(16)
			<< log_.d_.probe_time_ / MEGA << "(ms)" << std::endl;
	s << "# probe_time_max    =" << std::setw(16)
			<< log_.d_.probe_time_max_ / MEGA << "(ms)" << std::endl;

	s << "# preprocess_time_  =" << std::setw(16)
			<< log_.d_.preprocess_time_ / MEGA << "(ms)" << std::endl;

	s << "# iprobe_num        =" << std::setw(16)
			<< log_.d_.iprobe_num_ // rank 0
			<< std::endl;
	LOG(
			s << "# iprobe_time       =" << std::setw(16) << log_.d_.iprobe_time_ / MEGA // rank 0
			<< "(ms)" << std::endl; s << "# iprobe_time_max   =" << std::setw(16) << log_.d_.iprobe_time_max_ / MEGA << "(ms)" << std::endl; s << "# ms / (one iprobe) =" << std::setw(16) << log_.d_.iprobe_time_ / log_.d_.iprobe_num_ / MEGA// rank 0
			<< "(ms)" << std::endl;);

	s << "# iprobe_succ_num   =" << std::setw(16)
			<< log_.d_.iprobe_succ_num_ // rank 0
			<< std::endl;
	LOG(
			s << "# iprobe_succ_time  =" << std::setw(16) << log_.d_.iprobe_succ_time_ / MEGA // rank 0
			<< "(ms)" << std::endl; s << "# iprobe_succ_time_m=" << std::setw(16) << log_.d_.iprobe_succ_time_max_ / MEGA << "(ms)" << std::endl; s << "# ms / (one iprobe) =" << std::setw(16) << log_.d_.iprobe_succ_time_ / log_.d_.iprobe_succ_num_ / MEGA// rank 0
			<< "(ms)" << std::endl;);

	s << "# iprobe_fail_num   =" << std::setw(16)
			<< log_.d_.iprobe_fail_num_ // rank 0
			<< std::endl;
	LOG(
			s << "# iprobe_fail_time  =" << std::setw(16) << log_.d_.iprobe_fail_time_ / MEGA // rank 0
			<< "(ms)" << std::endl; s << "# iprobe_fail_time_m=" << std::setw(16) << log_.d_.iprobe_fail_time_max_ / MEGA << "(ms)" << std::endl; s << "# ms / (one iprobe) =" << std::setw(16) << log_.d_.iprobe_fail_time_ / log_.d_.iprobe_fail_num_ / MEGA// rank 0
			<< "(ms)" << std::endl;);

	s << "# recv_num          =" << std::setw(16) << log_.d_.recv_num_
			<< std::endl;
	LOG(
			s << "# recv_time         =" << std::setw(16) << log_.d_.recv_time_ / MEGA << "(ms)" << std::endl; s << "# recv_time_max     =" << std::setw(16) << log_.d_.recv_time_max_ / MEGA << "(ms)" << std::endl;);

	s << "# bsend_num         =" << std::setw(16)
			<< log_.d_.bsend_num_ << std::endl;
	s << "# bsend_time        =" << std::setw(16)
			<< log_.d_.bsend_time_ / MEGA << "(ms)" << std::endl;
	s << "# bsend_time_max    =" << std::setw(16)
			<< log_.d_.bsend_time_max_ / MEGA << "(ms)" << std::endl;

	s << "# bcast_num         =" << std::setw(16)
			<< log_.d_.bcast_num_ << std::endl;
	s << "# bcast_time        =" << std::setw(16)
			<< log_.d_.bcast_time_ / MEGA << "(ms)" << std::endl;
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

	s << "# dtd_phase_num_    =" << std::setw(16)
			<< log_.d_.dtd_phase_num_ << std::endl;
	s << "# dtd_phase_per_sec_=" << std::setw(16)
			<< log_.d_.dtd_phase_per_sec_ << std::endl;
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

	s << "# lfl_steal_num     =" << std::setw(16)
			<< log_.d_.lifeline_steal_num_ << std::endl;
	s << "# lfl_nodes_recv    =" << std::setw(16)
			<< log_.d_.lifeline_nodes_received_ << std::endl;
	s << "# steal_num         =" << std::setw(16)
			<< log_.d_.steal_num_ << std::endl;
	s << "# nodes_received    =" << std::setw(16)
			<< log_.d_.nodes_received_ << std::endl;

	s << "# lfl_given_num     =" << std::setw(16)
			<< log_.d_.lifeline_given_num_ << std::endl;
	s << "# lfl_nodes_given   =" << std::setw(16)
			<< log_.d_.lifeline_nodes_given_ << std::endl;
	s << "# given_num         =" << std::setw(16)
			<< log_.d_.given_num_ << std::endl;
	s << "# nodes_given       =" << std::setw(16)
			<< log_.d_.nodes_given_ << std::endl;

	s << "# node_stack_max_itm=" << std::setw(16)
			<< log_.d_.node_stack_max_itm_ << std::endl;
	s << "# give_stack_max_itm=" << std::setw(16)
			<< log_.d_.give_stack_max_itm_ << std::endl;

	s << "# node_stack_max_cap=" << std::setw(16)
			<< log_.d_.node_stack_max_cap_ << std::endl;
	s << "# give_stack_max_cap=" << std::setw(16)
			<< log_.d_.give_stack_max_cap_ << std::endl;

	s << "# cleared_tasks_    =" << std::setw(16)
			<< log_.d_.cleared_tasks_ / MEGA << std::endl;

	out << s.str() << std::flush;
	return out;
}

//==============================================================================

void MP_CONT_LAMP::SetLogFileName() {
	std::stringstream ss;
	ss << FLAGS_debuglogfile << "_log";
	ss << std::setw(4) << std::setfill('0') << mpi_data_.mpiRank_
			<< ".txt";
	log_file_name_ = ss.str();
}

std::ofstream MP_CONT_LAMP::null_stream_;

std::ostream& MP_CONT_LAMP::D(int level, bool show_phase) {
	if (FLAGS_d == 0)
		return null_stream_;

	if (level <= FLAGS_d) {
		if (show_phase)
			lfs_ << std::setw(4) << phase_ << ": ";
		return lfs_;
	} else
		return null_stream_;
}

std::ostream& operator<<(std::ostream & out,
		const FixedSizeStack & st) {
	std::stringstream s;
	s << "entry:";
	for (int i = 0; i < st.Size(); i++)
		s << " " << st.data(i);
	s << std::endl;
	out << s.str() << std::flush;
	return out;
}

}	// namespace lamp_search

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
