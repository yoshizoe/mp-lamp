/*
 * ParallelSearch.cpp
 *
 *  Created on: Apr 16, 2017
 *      Author: yuu
 */
#include "ParallelSearch.h"

#include "gflags/gflags.h"
#include "mpi.h"

#include "Log.h"
#include "DTD.h"
#include <assert.h>

#ifdef __CDT_PARSER__
#undef DBG
#define DBG(a)  ;
#endif

#ifdef __CDT_PARSER__
#undef LOG
#define LOG(a)  ;
#endif

// Algorithms
DECLARE_bool(second_phase);// true, "do second phase"
DECLARE_bool(third_phase);// true, "do third phase"
DECLARE_double(alpha);// significance level alpha

// MPI
DEFINE_int32(bsend_buffer_size, 1024 * 1024 * 64, "size of bsend buffer");
DEFINE_int32(probe_period, 128, "probe period during process node");
DEFINE_bool(probe_period_is_ms, false,
		"true: probe period is milli sec, false: num loops");

// Utils
DEFINE_int32(d, 0, "debug level. 0: none, higher level produce more log");

namespace lamp_search {

ParallelSearch::ParallelSearch(MPI_Data& mpi_data) :
		mpi_data(mpi_data), rng_(mpi_data.mpiRank_), dst_p_(0,
				mpi_data.nTotalProc - 1), dst_m_(0,
				mpi_data.nRandStealCands - 1), rand_p_(rng_, dst_p_), rand_m_(
				rng_, dst_m_), d_(NULL), g_(NULL), bsh_(NULL), timer_(
				Timer::GetInstance()), pmin_thr_(NULL), cs_thr_(NULL), dtd_accum_array_base_(
		NULL), accum_array_(NULL), dtd_accum_recv_base_(NULL), accum_recv_(
		NULL), node_stack_(NULL), give_stack_(NULL), processing_node_(false), waiting_(
				false), stealer_(mpi_data.nRandStealTrials,
				mpi_data.hypercubeDimension), phase_(0), sup_buf_(
		NULL), child_sup_buf_(NULL), freq_stack_(NULL), significant_stack_(
		NULL), total_expand_num_(0ll), expand_num_(0ll), closed_set_num_(0ll), final_sig_level_(
				0.0) {

	SetupTopology();
	InitTreeRequest();
	Init();
}

ParallelSearch::~ParallelSearch() {
	DBG(D(2) << "MP_LAMP destructor begin" << std::endl;);

	delete[] lifelines_activated_;
	delete[] accum_flag_;
	delete[] bcast_targets_;

	if (pmin_thr_)
		delete[] pmin_thr_;
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

	DBG(D(2) << "MP_LAMP destructor end" << std::endl;);
	if (FLAGS_d > 0) {
		DBG(lfs_.close(););
	}
}

int ParallelSearch::ComputeZ(int p, int l) {
	int z0 = 1;
	int zz = l;
	while (zz < p) {
		z0++;
		zz *= l;
	}
	return z0;
}

// TODO: not sure what it is doing
void ParallelSearch::InitTreeRequest() {
	// pushing tree to lifeline for the 1st wave
	int num = 0;

	for (std::size_t i = 1; i <= k_echo_tree_branch; i++) {
		if (k_echo_tree_branch * mpiRank_ + i < nTotalProc) {
			lifeline_thieves_->Push(k_echo_tree_branch * mpiRank_ + i);
			bcast_targets_[num++] = k_echo_tree_branch * mpiRank_ + i;
		}
	}
	if (mpiRank_ > 0) {
		lifelines_activated_[(mpiRank_ - 1) / k_echo_tree_branch] = true;
		bcast_source_ = (mpiRank_ - 1) / k_echo_tree_branch;
	}
}

void ParallelSearch::Init() {
	// todo: move accum cs count variable to inner class
	echo_waiting_ = false;
	for (int i = 0; i < k_echo_tree_branch; i++)
		accum_flag_[i] = false;

	dtd_.Init();

	log_.Init();

	stealer_.Init();
	waiting_ = false;

	last_bcast_was_dtd_ = false;
}

void ParallelSearch::SetTreeRequest() {
	for (std::size_t i = 1; i <= k_echo_tree_branch; i++) {
		if (k_echo_tree_branch * mpiRank_ + i < nTotalProc)
			lifeline_thieves_->Push(k_echo_tree_branch * mpiRank_ + i);
	}
	// if (3*h_+1 < p_) lifeline_thieves_->Push(3*h_+1);
	// if (3*h_+2 < p_) lifeline_thieves_->Push(3*h_+2);
	// if (3*h_+3 < p_) lifeline_thieves_->Push(3*h_+3);
	if (mpiRank_ > 0)
		lifelines_activated_[(mpiRank_ - 1) / k_echo_tree_branch] = true;
}

void ParallelSearch::CheckAndClear() {
	DBG(D(4) << "CheckPoint Started"<< std::endl;);
	assert(echo_waiting_ == false);

	// does this hold?
	for (int i = 0; i < k_echo_tree_branch; i++)
		assert(accum_flag_[i] == false);

	dtd_.CheckPoint();
	stealer_.Init();
	waiting_ = false;

	thieves_->Clear();
	lifeline_thieves_->Clear();
	for (int pi = 0; pi < nTotalProc; pi++)
		lifelines_activated_[pi] = false;

	SetTreeRequest();
	DBG(D(4) << "CheckPoint Finished"<< std::endl;);
}

void ParallelSearch::ClearTasks() {
	MPI_Status probe_status, recv_status;
	int data_count, src, tag;
	int flag;
	int error;

	while (true) {
		error = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag,
				&probe_status);
		if (error != MPI_SUCCESS) {
			DBG(
					D(1) << "error in MPI_Iprobe in ClearTasks: " << error << std::endl;);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		if (!flag)
			break;

		log_.d_.cleared_tasks_++;

		error = MPI_Get_count(&probe_status, MPI_INT, &data_count);
		if (error != MPI_SUCCESS) {
			DBG(
					D(1) << "error in MPI_Get_count in ClearTasks: " << error << std::endl;);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		src = probe_status.MPI_SOURCE;
		tag = probe_status.MPI_TAG;

		error = MPI_Recv(give_stack_->Stack(), data_count, MPI_INT, src, tag,
		MPI_COMM_WORLD, &recv_status);
		if (error != MPI_SUCCESS) {
			DBG(
					D(1) << "error in MPI_Recv in ClearTasks: " << error << std::endl;);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		give_stack_->Clear(); // drop messages
	}
}

// TODO: way too awful.
void ParallelSearch::setPhase(int phase) {
	phase_ = phase;
}

void ParallelSearch::Search(int phase) {
	setPhase(phase);
	Search();
}

void ParallelSearch::Search() {
	while (!dtd_.terminated_) {
		while (!dtd_.terminated_) {
			if (ProcessNode(granularity)) {
				log_.d_.node_stack_max_itm_ = std::max(
						log_.d_.node_stack_max_itm_,
						(long long int) (node_stack_->NuItemset()));
				log_.d_.node_stack_max_cap_ = std::max(
						log_.d_.node_stack_max_cap_,
						(long long int) (node_stack_->UsedCapacity()));
				Probe();
				if (dtd_.terminated_)
					break;
				Distribute();
				Reject(); // distribute finished, reject remaining requests
				if (mpiRank_ == 0 && phase_ == 1)
					CheckCSThreshold();
			} else
				break;
		}
		if (dtd_.terminated_)
			break;

		log_.idle_start_ = timer_->Elapsed();
		Reject(); // node_stack_ empty. reject requests
		Steal(); // request steal
		if (dtd_.terminated_) {
			log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
			break;
		}

		Probe();
		if (dtd_.terminated_) {
			log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
			break;
		}
		if (mpiRank_ == 0 && phase_ == 1)
			CheckCSThreshold();
		log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
	}
}

void ParallelSearch::Probe() {
	DBG(D(3) << "Probe" << std::endl;);
	MPI_Status probe_status;
	int probe_src, probe_tag;

	long long int start_time;
	start_time = timer_->Elapsed();

	log_.d_.probe_num_++;

	while (CallIprobe(&probe_status, &probe_src, &probe_tag)) {
		DBG(
				D(4) << "CallIprobe returned src=" << probe_src << "\ttag=" << probe_tag << std::endl;);
		switch (probe_tag) {
		// control tasks
		case Tag::DTD_REQUEST:
			RecvDTDRequest(probe_src);
			break;
		case Tag::DTD_REPLY:
			RecvDTDReply(probe_src);
			break;

		case Tag::DTD_ACCUM_REQUEST:
			assert(phase_ == 1);
			RecvDTDAccumRequest(probe_src);
			break;
		case Tag::DTD_ACCUM_REPLY:
			assert(phase_ == 1);
			RecvDTDAccumReply(probe_src);
			break;

		case Tag::BCAST_FINISH:
			RecvBcastFinish(probe_src);
			break;

// basic tasks
		case Tag::LAMBDA:
			assert(phase_ == 1);
			RecvLambda(probe_src);
			break;

		case Tag::REQUEST:
			RecvRequest(probe_src);
			break;
		case Tag::REJECT:
			RecvReject(probe_src);
			break;
		case Tag::GIVE:
			RecvGive(probe_src, probe_status);
			break;

// third phase tasks
		case Tag::RESULT_REQUEST:
			RecvResultRequest(probe_src);
			break;
		case Tag::RESULT_REPLY:
			RecvResultReply(probe_src, probe_status);
			break;
		default:
			DBG(
					D(1) << "unknown Tag=" << probe_tag << " received in Probe: " << std::endl;)
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

	if (mpiRank_ == 0) {
		// initiate termination detection

		// note: for phase_ 1, accum request and dtd request are unified
		if (!echo_waiting_ && !dtd_.terminated_) {
			if (phase_ == 1) {
				SendDTDAccumRequest();
				// log_.d_.dtd_accum_request_num_++;
			} else if (phase_ == 2) {
				if (node_stack_->Empty()) {
					SendDTDRequest();
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

void ParallelSearch::Distribute() {
	if (nTotalProc == 1)
		return;DBG(D(3) << "Distribute" << std::endl;);
	if (thieves_->Size() > 0 || lifeline_thieves_->Size() > 0) {
		int steal_num = node_stack_->Split(give_stack_);
		if (steal_num > 0) {
			DBG(D(3) << "giving" << std::endl;);DBG(
					give_stack_->PrintAll(D(3,false)););
			Give(give_stack_, steal_num);
			give_stack_->Clear();
		}
	}
}

void ParallelSearch::Reject() {
}

void ParallelSearch::CheckCSThreshold() {
}

bool ParallelSearch::ProcessNode(int n) {
	if (node_stack_->Empty())
		return false;
	long long int start_time, lap_time;
	start_time = timer_->Elapsed();
	lap_time = start_time;

	int processed = 0;
	processing_node_ = true;
	while (!node_stack_->Empty()) {
		processed++;
		expand_num_++;

		node_stack_->CopyItem(node_stack_->Top(), itemset_buf_);
		node_stack_->Pop();

		// dbg
		DBG(D(3) << "expanded ";);DBG(node_stack_->Print(D(3), itemset_buf_););

		// calculate support from itemset_buf_
		bsh_->Set(sup_buf_);
		{
			int n = node_stack_->GetItemNum(itemset_buf_);
			for (int i = 0; i < n; i++) {
				int item = node_stack_->GetNthItem(itemset_buf_, i);
				bsh_->And(d_->NthData(item), sup_buf_);
			}
		}

		int core_i = g_->CoreIndex(*node_stack_, itemset_buf_);

		int * ppc_ext_buf;
		// todo: use database reduction

		assert(phase_ != 1 || node_stack_->GetItemNum(itemset_buf_) != 0);

		bool is_root_node = (node_stack_->GetItemNum(itemset_buf_) == 0);

		int accum_period_counter_ = 0;
		// reverse order
		// for ( int new_item = d_->NuItems()-1 ; new_item >= core_i+1 ; new_item-- )
		for (int new_item = d_->NextItemInReverseLoop(is_root_node, mpiRank_,
				nTotalProc, d_->NuItems()); new_item >= core_i + 1;
				new_item = d_->NextItemInReverseLoop(is_root_node, mpiRank_,
						nTotalProc, new_item)) {
			// skip existing item
			// todo: improve speed here
			if (node_stack_->Exist(itemset_buf_, new_item))
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

							Probe();
							Distribute();
							Reject();

							lap_time = timer_->Elapsed();
						}
					}
				} else {            // not using milli second
					if (accum_period_counter_ >= FLAGS_probe_period) {
						accum_period_counter_ = 0;
						log_.d_.process_node_time_ += timer_->Elapsed()
								- lap_time;

						Probe();
						Distribute();
						Reject();

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

			node_stack_->PushPre();
			ppc_ext_buf = node_stack_->Top();

			bool res = g_->PPCExtension(node_stack_, itemset_buf_,
					child_sup_buf_, core_i, new_item, ppc_ext_buf);

			node_stack_->SetSup(ppc_ext_buf, sup_num);
			node_stack_->PushPostNoSort();

			if (!res) {        // todo: remove this redundancy
				node_stack_->Pop();
			} else {
				node_stack_->SortTop();

				DBG(
						if (phase_ == 2) {D(3) << "found cs "; node_stack_->Print(D(3), ppc_ext_buf);});

				if (phase_ == 1)
					IncCsAccum(sup_num); // increment closed_set_num_array
				if (phase_ == 2)
					closed_set_num_++;

				if (phase_ == 2 && FLAGS_third_phase) {
					int pos_sup_num = bsh_->AndCount(d_->PosNeg(),
							child_sup_buf_);
					double pval = d_->PVal(sup_num, pos_sup_num);
					assert(pval >= 0.0);
					if (pval <= sig_level_) { // permits == case?
						freq_stack_->PushPre();
						int * item = freq_stack_->Top();
						freq_stack_->CopyItem(ppc_ext_buf, item);
						freq_stack_->PushPostNoSort();

						freq_map_.insert(std::pair<double, int*>(pval, item));
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

		if (CheckProcessNodeEnd(n, isGranularitySec, processed, start_time))
			break;
	}

	long long int elapsed_time = timer_->Elapsed() - lap_time;
	log_.d_.process_node_time_ += elapsed_time;
	log_.d_.process_node_num_ += processed;

	DBG(
			D(2) << "processed node num=" << processed << "\ttime=" << elapsed_time << std::endl;);

	processing_node_ = false;
	return true;
}

bool ParallelSearch::CheckProcessNodeEnd(int n, bool n_is_ms, int processed,
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

void ParallelSearch::Steal() {
}

void ParallelSearch::Give(VariableLengthItemsetStack* st, int steal_num) {
	DBG(D(3) << "Give: ";);
	if (thieves_->Size() > 0) { // random thieves
		int thief = thieves_->Pop();
		if (thief >= 0) { // lifeline thief
			DBG(
					D(3,false) << "thief=" << thief << "\trandom lifeline" << std::endl;);
			log_.d_.lifeline_given_num_++;
			log_.d_.lifeline_nodes_given_ += steal_num;
			SendGive(st, thief, mpiRank_);
		} else { // random thief
			DBG(
					D(3,false) << "thief=" << (-thief-1) << "\trandom" << std::endl;);
			log_.d_.given_num_++;
			log_.d_.nodes_given_ += steal_num;
			SendGive(st, (-thief - 1), -1);
		}
	} else {
		assert(lifeline_thieves_->Size() > 0);
		int thief = lifeline_thieves_->Pop();
		assert(thief >= 0);DBG(
				D(3,false) << "thief=" << thief << "\tlifeline" << std::endl;);
		log_.d_.lifeline_given_num_++;
		log_.d_.lifeline_nodes_given_ += steal_num;
		SendGive(st, thief, mpiRank_);
	}
}

//==============================================================================

bool ParallelSearch::IsLeaf() const {
	for (int i = 0; i < k_echo_tree_branch; i++)
		if (bcast_targets_[i] >= 0)
			return false;
	return true;
}

int ParallelSearch::CallIprobe(MPI_Status * status, int * src, int * tag) {
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
		DBG(D(1) << "error in MPI_Iprobe: " << error << std::endl;);
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

int ParallelSearch::CallRecv(void * buffer, int data_count, MPI_Datatype type,
		int src, int tag, MPI_Status * status) {
	long long int start_time;
	long long int end_time;

	// todo: prepare non-log mode to remove measurement
	// clock_gettime takes 0.3--0.5 micro sec
	log_.d_.recv_num_++;
	LOG(start_time = timer_->Elapsed(););

	int error = MPI_Recv(buffer, data_count, type, src, tag,
	MPI_COMM_WORLD, status);

	LOG(
			end_time = timer_->Elapsed(); log_.d_.recv_time_ += end_time - start_time; log_.d_.recv_time_max_ = std::max( end_time - start_time, log_.d_.recv_time_max_););
	return error;
}

int ParallelSearch::CallBsend(void * buffer, int data_count, MPI_Datatype type,
		int dest, int tag) {
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

int ParallelSearch::CallBcast(void * buffer, int data_count,
		MPI_Datatype type) {
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

// TODO: Why are these methods not handled by DTD???
// call this from root rank to start DTD
void ParallelSearch::SendDTDRequest() {
	int message[1];
	message[0] = 1; // dummy

	echo_waiting_ = true;

	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (bcast_targets_[i] < 0)
			break;
		CallBsend(message, 1, MPI_INT, bcast_targets_[i], Tag::DTD_REQUEST);
		DBG(
				D(3) << "SendDTDRequest: dst=" << bcast_targets_[i] << "\ttimezone=" << dtd_.time_zone_ << std::endl;);
	}
}

void ParallelSearch::RecvDTDRequest(int src) {
	DBG(
			D(3) << "RecvDTDRequest: src=" << src << "\ttimezone=" << dtd_.time_zone_ << std::endl;);
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::DTD_REQUEST, &recv_status);

	if (IsLeaf())
		SendDTDReply();
	else
		SendDTDRequest();
}

bool ParallelSearch::DTDReplyReady() const {
	for (int i = 0; i < k_echo_tree_branch; i++)
		if (bcast_targets_[i] >= 0 && !(dtd_.accum_flag_[i]))
			return false;
	return true;
	// check only valid bcast_targets_
	// always true if leaf
}

void ParallelSearch::SendDTDReply() {
	int message[3];
	// using reduced vars
	message[0] = dtd_.count_ + dtd_.reduce_count_;
	bool tw_flag = dtd_.time_warp_ || dtd_.reduce_time_warp_;
	message[1] = (tw_flag ? 1 : 0);
	// for Steal
	dtd_.not_empty_ = !(node_stack_->Empty()) || (thieves_->Size() > 0)
			|| stealer_.StealStarted() || processing_node_; // thieves_ and stealer state check
	// for Steal2
	// dtd_.not_empty_ =
	//     !(node_stack_->Empty()) || (thieves_->Size() > 0) ||
	//     waiting_ || processing_node_;
	bool em_flag = dtd_.not_empty_ || dtd_.reduce_not_empty_;
	message[2] = (em_flag ? 1 : 0);
	DBG(
			D(3) << "SendDTDReply: dst = " << bcast_source_ << "\tcount=" << message[0] << "\ttw=" << tw_flag << "\tem=" << em_flag << std::endl;);

	CallBsend(message, 3, MPI_INT, bcast_source_, Tag::DTD_REPLY);

	dtd_.time_warp_ = false;
	dtd_.not_empty_ = false;
	dtd_.IncTimeZone();

	echo_waiting_ = false;
	dtd_.ClearAccumFlags();
	dtd_.ClearReduceVars();
}

void ParallelSearch::RecvDTDReply(int src) {
	MPI_Status recv_status;
	int message[3];
	CallRecv(&message, 3, MPI_INT, src, Tag::DTD_REPLY, &recv_status);
	assert(src == recv_status.MPI_SOURCE);

	// reduce reply (count, time_warp, not_empty)
	dtd_.Reduce(message[0], (message[1] != 0), (message[2] != 0));

	DBG(
			D(3) << "RecvDTDReply: src=" << src << "\tcount=" << message[0] << "\ttw=" << message[1] << "\tem=" << message[2] << "\treduced_count=" << dtd_.reduce_count_ << "\treduced_tw=" << dtd_.reduce_time_warp_ << "\treduced_em=" << dtd_.reduce_not_empty_ << std::endl;);

	bool flag = false;
	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (bcast_targets_[i] == src) {
			flag = true;
			dtd_.accum_flag_[i] = true;
			break;
		}
	}
	assert(flag);

	if (DTDReplyReady()) {
		if (mpiRank_ == 0) {
			DTDCheck(); // at root
			log_.d_.dtd_phase_num_++;
		} else
			SendDTDReply();
	}
}

void ParallelSearch::SendDTDAccumRequest() {
	int message[1];
	message[0] = 1; // dummy

	echo_waiting_ = true;

	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (bcast_targets_[i] < 0)
			break;
		CallBsend(message, 1, MPI_INT, bcast_targets_[i],
				Tag::DTD_ACCUM_REQUEST);

		DBG(
				D(3) << "SendDTDAccumRequest: dst=" << bcast_targets_[i] << "\ttimezone=" << dtd_.time_zone_ << std::endl;);
	}
}

void ParallelSearch::RecvDTDAccumRequest(int src) {
	DBG(
			D(3) << "RecvDTDAccumRequest: src=" << src << "\ttimezone=" << dtd_.time_zone_ << std::endl;);
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::DTD_ACCUM_REQUEST, &recv_status);

	if (IsLeaf())
		SendDTDAccumReply();
	else
		SendDTDAccumRequest();
}

void ParallelSearch::SendDTDAccumReply() {
	dtd_accum_array_base_[0] = dtd_.count_ + dtd_.reduce_count_;
	bool tw_flag = dtd_.time_warp_ || dtd_.reduce_time_warp_;
	dtd_accum_array_base_[1] = (tw_flag ? 1 : 0);
	// for Steal
	dtd_.not_empty_ = !(node_stack_->Empty()) || (thieves_->Size() > 0)
			|| stealer_.StealStarted() || processing_node_; // thieves_ and stealer state check
	// for Steal2
	// dtd_.not_empty_ =
	//     !(node_stack_->Empty()) || (thieves_->Size() > 0) ||
	//     waiting_ || processing_node_;
	bool em_flag = dtd_.not_empty_ || dtd_.reduce_not_empty_;
	dtd_accum_array_base_[2] = (em_flag ? 1 : 0);

	DBG(
			D(3) << "SendDTDAccumReply: dst = " << bcast_source_ << "\tcount=" << dtd_accum_array_base_[0] << "\ttw=" << tw_flag << "\tem=" << em_flag << std::endl;);

	CallBsend(dtd_accum_array_base_, lambda_max_ + 4,
	MPI_LONG_LONG_INT, bcast_source_, Tag::DTD_ACCUM_REPLY);

	dtd_.time_warp_ = false;
	dtd_.not_empty_ = false;
	dtd_.IncTimeZone();

	echo_waiting_ = false;
	dtd_.ClearAccumFlags();
	dtd_.ClearReduceVars();

	for (int l = 0; l <= lambda_max_; l++)
		accum_array_[l] = 0ll;
}

void ParallelSearch::RecvDTDAccumReply(int src) {
	MPI_Status recv_status;

	CallRecv(dtd_accum_recv_base_, lambda_max_ + 4,
	MPI_LONG_LONG_INT, src, Tag::DTD_ACCUM_REPLY, &recv_status);
	assert(src == recv_status.MPI_SOURCE);

	int count = (int) (dtd_accum_recv_base_[0]);
	bool time_warp = (dtd_accum_recv_base_[1] != 0);
	bool not_empty = (dtd_accum_recv_base_[2] != 0);

	dtd_.Reduce(count, time_warp, not_empty);

	DBG(
			D(3) << "RecvDTDAccumReply: src=" << src << "\tcount=" << count << "\ttw=" << time_warp << "\tem=" << not_empty << "\treduced_count=" << dtd_.reduce_count_ << "\treduced_tw=" << dtd_.reduce_time_warp_ << "\treduced_em=" << dtd_.reduce_not_empty_ << std::endl;);

	for (int l = lambda_ - 1; l <= lambda_max_; l++)
		accum_array_[l] += accum_recv_[l];

	bool flag = false;
	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (bcast_targets_[i] == src) {
			flag = true;
			dtd_.accum_flag_[i] = true;
			break;
		}
	}
	assert(flag);

	if (mpiRank_ == 0) {
		if (ExceedCsThr()) {
			int new_lambda = NextLambdaThr();
			SendLambda(new_lambda);
			lambda_ = new_lambda;
		}
		// if SendLambda is called, dtd_.count_ is incremented and DTDCheck will always fail
		if (DTDReplyReady()) {
			DTDCheck();
			log_.d_.dtd_accum_phase_num_++;
		}
	} else {  // not root
		if (DTDReplyReady())
			SendDTDAccumReply();
	}
}

void ParallelSearch::DTDCheck() {
	assert(mpiRank_ == 0);
	// (count, time_warp, not_empty)
	dtd_.Reduce(dtd_.count_, dtd_.time_warp_, dtd_.not_empty_);

	if (dtd_.reduce_count_ == 0 && dtd_.reduce_time_warp_ == false
			&& dtd_.reduce_not_empty_ == false) {
		// termination
		SendBcastFinish();
		dtd_.terminated_ = true;
		DBG(D(1) << "terminated" << std::endl;);
	}
	// doing same thing as SendDTDReply
	dtd_.time_warp_ = false;
	dtd_.not_empty_ = false;
	dtd_.IncTimeZone();

	echo_waiting_ = false;
	dtd_.ClearAccumFlags();
	dtd_.ClearReduceVars();
}

void ParallelSearch::SendBcastFinish() {
	DBG(D(2) << "SendBcastFinish" << std::endl;);
	int message[1];
	message[0] = 1; // dummy

	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (bcast_targets_[i] < 0)
			break;
		CallBsend(message, 1, MPI_INT, bcast_targets_[i], Tag::BCAST_FINISH);
	}
}

void ParallelSearch::RecvBcastFinish(int src) {
	DBG(D(2) << "RecvBcastFinish: src=" << src << std::endl;);
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::BCAST_FINISH, &recv_status);

	SendBcastFinish();
	dtd_.terminated_ = true;
	DBG(D(2) << "terminated" << std::endl;);

	waiting_ = false;
}

// lifeline == -1 for random thieves
void ParallelSearch::SendRequest(int dst, int is_lifeline) {
	assert(dst >= 0);
	int message[2];
	message[0] = dtd_.time_zone_;
	message[1] = is_lifeline; // -1 for random thieves, >=0 for lifeline thieves

	CallBsend(message, 2, MPI_INT, dst, Tag::REQUEST);
	dtd_.OnSend();

	DBG(
			D(2) << "SendRequest: dst=" << dst << "\tis_lifeline=" << is_lifeline << "\tdtd_count=" << dtd_.count_ << std::endl;);
}

void ParallelSearch::RecvRequest(int src) {
	DBG(D(2) << "RecvRequest: src=" << src;);
	MPI_Status recv_status;
	int message[2];

	CallRecv(&message, 2, MPI_INT, src, Tag::REQUEST, &recv_status);
	dtd_.OnRecv();
	assert(src == recv_status.MPI_SOURCE);

	int timezone = message[0];
	dtd_.UpdateTimeZone(timezone);
	int is_lifeline = message[1]; // -1 for random thieves, >=0 for lifeline thieves
	int thief = src;
	DBG(D(2,false) << "\tis_lifeline=" << is_lifeline;);DBG(
			D(2,false) << "\tthief=" << src;);DBG(
			D(2,false) << "\tdtd_count=" << dtd_.count_;);

	if (node_stack_->Empty()) {
		DBG(D(2,false) << "\tempty and reject" << std::endl;);
		if (is_lifeline >= 0) {
			lifeline_thieves_->Push(thief);
			SendReject(thief); // notify
		} else {
			SendReject(thief); // notify
		}
	} else {
		DBG(D(2,false) << "\tpush" << std::endl;);
		if (is_lifeline >= 0)
			thieves_->Push(thief);
		else
			thieves_->Push(-thief - 1);
	}
	// todo: take log
}

void ParallelSearch::SendReject(int dst) {
	assert(dst >= 0);
	int message[1];

	message[0] = dtd_.time_zone_;
	CallBsend(message, 1, MPI_INT, dst, Tag::REJECT);
	dtd_.OnSend();

	DBG(
			D(2) << "SendReject: dst=" << dst << "\tdtd_count=" << dtd_.count_ << std::endl;);
}

void ParallelSearch::RecvReject(int src) {
	MPI_Status recv_status;
	int message[1];

	CallRecv(&message, 1, MPI_INT, src, Tag::REJECT, &recv_status);
	dtd_.OnRecv();
	assert(src == recv_status.MPI_SOURCE);

	int timezone = message[0];
	dtd_.UpdateTimeZone(timezone);
	DBG(
			D(2) << "RecvReject: src=" << src << "\tdtd_count=" << dtd_.count_ << std::endl;);

	stealer_.ResetRequesting();
	waiting_ = false;
}

void ParallelSearch::SendGive(VariableLengthItemsetStack * st, int dst,
		int is_lifeline) {
	assert(dst >= 0);
	st->SetTimestamp(dtd_.time_zone_);
	st->SetFlag(is_lifeline);
	int * message = st->Stack();
	int size = st->UsedCapacity();

	log_.d_.give_stack_max_itm_ = std::max(log_.d_.give_stack_max_itm_,
			(long long int) (st->NuItemset()));
	log_.d_.give_stack_max_cap_ = std::max(log_.d_.give_stack_max_cap_,
			(long long int) (st->UsedCapacity()));

	CallBsend(message, size, MPI_INT, dst, Tag::GIVE);
	dtd_.OnSend();

	DBG(
			D(2) << "SendGive: " << "\ttimezone=" << dtd_.time_zone_ << "\tdst=" << dst << "\tlfl=" << is_lifeline << "\tsize=" << size << "\tnode=" << st->NuItemset() << "\tdtd_count=" << dtd_.count_ << std::endl;);
	//st->PrintAll(D(false));
}

void ParallelSearch::RecvGive(int src, MPI_Status probe_status) {
	int count;
	int error = MPI_Get_count(&probe_status, MPI_INT, &count);
	if (error != MPI_SUCCESS) {
		DBG(
				D(1) << "error in MPI_Get_count in RecvGive: " << error << std::endl;);
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
		lifelines_activated_[src] = false;
		log_.d_.lifeline_steal_num_++;
		log_.d_.lifeline_nodes_received_ += (new_nu_itemset - orig_nu_itemset);
	} else {
		log_.d_.steal_num_++;
		log_.d_.nodes_received_ += (new_nu_itemset - orig_nu_itemset);
	}

	DBG(
			D(2) << "RecvGive: src=" << src << "\ttimezone=" << dtd_.time_zone_ << "\tlfl=" << flag << "\tsize=" << count << "\tnode=" << (new_nu_itemset - orig_nu_itemset) << "\tdtd_count=" << dtd_.count_ << std::endl;);

	give_stack_->Clear();

	log_.d_.node_stack_max_itm_ = std::max(log_.d_.node_stack_max_itm_,
			(long long int) (node_stack_->NuItemset()));
	log_.d_.node_stack_max_cap_ = std::max(log_.d_.node_stack_max_cap_,
			(long long int) (node_stack_->UsedCapacity()));

	DBG(node_stack_->PrintAll(D(3,false)););

	stealer_.ResetRequesting();
	stealer_.ResetCounters();
	stealer_.SetState(StealState::RANDOM);
	stealer_.SetStealStart();
	waiting_ = false;
}

void ParallelSearch::SendLambda(int lambda) {
	// send lambda to bcast_targets_
	int message[2];
	message[0] = dtd_.time_zone_;
	message[1] = lambda;

	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (bcast_targets_[i] < 0)
			break;
		CallBsend(message, 2, MPI_INT, bcast_targets_[i], Tag::LAMBDA);
		dtd_.OnSend();

		DBG(
				D(2) << "SendLambda: dst=" << bcast_targets_[i] << "\tlambda=" << lambda << "\tdtd_count=" << dtd_.count_ << std::endl;);
	}
}

void ParallelSearch::RecvLambda(int src) {
	MPI_Status recv_status;
	int message[2];

	CallRecv(&message, 2, MPI_INT, src, Tag::LAMBDA, &recv_status);
	dtd_.OnRecv();
	assert(src == recv_status.MPI_SOURCE);
	int timezone = message[0];
	dtd_.UpdateTimeZone(timezone);

	DBG(
			D(2) << "RecvLambda: src=" << src << "\tlambda=" << message[1] << "\tdtd_count=" << dtd_.count_ << std::endl;);

	int new_lambda = message[1];
	if (new_lambda > lambda_) {
		SendLambda(new_lambda);
		lambda_ = new_lambda;
		// todo: do database reduction
	}
}

bool ParallelSearch::AccumCountReady() const {
	for (int i = 0; i < k_echo_tree_branch; i++)
		if (bcast_targets_[i] >= 0 && !(accum_flag_[i]))
			return false;
	return true;
	// check only valid bcast_targets_
	// always true if leaf
}

bool ParallelSearch::ExceedCsThr() const {
	// note: > is correct. permit ==
	return (accum_array_[lambda_] > cs_thr_[lambda_]);
}

int ParallelSearch::NextLambdaThr() const {
	int si;
	for (si = lambda_max_; si >= lambda_; si--)
		if (accum_array_[si] > cs_thr_[si])
			break;
	return si + 1;
	// it is safe because lambda_ higher than max results in immediate search finish
}

void ParallelSearch::IncCsAccum(int sup_num) {
	for (int i = sup_num; i >= lambda_ - 1; i--)
		accum_array_[i]++;
}

double ParallelSearch::GetInterimSigLevel(int lambda) const {
	long long int csnum = accum_array_[lambda];
	double lv;
	if (csnum > 0)
		lv = FLAGS_alpha / (double) csnum;
	else
		lv = FLAGS_alpha;

	return lv;
}

//==============================================================================

void ParallelSearch::SendResultRequest() {
	int message[1];
	message[0] = 1; // dummy

	echo_waiting_ = true;

	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (bcast_targets_[i] < 0)
			break;
		CallBsend(message, 1, MPI_INT, bcast_targets_[i], Tag::RESULT_REQUEST);
		DBG(
				D(2) << "SendResultRequest: dst=" << bcast_targets_[i] << std::endl;);
	}
}

void ParallelSearch::RecvResultRequest(int src) {
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::RESULT_REQUEST, &recv_status);
	assert(src == recv_status.MPI_SOURCE);

	DBG(D(2) << "RecvResultRequest: src=" << src << std::endl;);

	if (IsLeaf())
		SendResultReply();
	else
		SendResultRequest();
}

void ParallelSearch::SendResultReply() {
	int * message = significant_stack_->Stack();
	int size = significant_stack_->UsedCapacity();

	CallBsend(message, size, MPI_INT, bcast_source_, Tag::RESULT_REPLY);

	DBG(D(2) << "SendResultReply: dst=" << bcast_source_ << std::endl;);DBG(
			significant_stack_->PrintAll(D(3,false)););

	echo_waiting_ = false;
	dtd_.terminated_ = true;
}

void ParallelSearch::RecvResultReply(int src, MPI_Status probe_status) {
	int count;
	int error = MPI_Get_count(&probe_status, MPI_INT, &count);
	if (error != MPI_SUCCESS) {
		DBG(
				D(1) << "error in MPI_Get_count in RecvResultReply: " << error << std::endl;);
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

	DBG(D(2) << "RecvResultReply: src=" << src << std::endl;);DBG(
			significant_stack_->PrintAll(D(3,false)););

// using the same flags as accum count, should be fixed
	bool flag = false;
	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (bcast_targets_[i] == src) {
			flag = true;
			accum_flag_[i] = true;
			break;
		}
	}
	assert(flag);

	if (AccumCountReady()) {
		if (mpiRank_ != 0) {
			SendResultReply();
		} else { // root
			echo_waiting_ = false;
			dtd_.terminated_ = true;
		}
	}
}

void ParallelSearch::ExtractSignificantSet() {
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

void ParallelSearch::SortSignificantSets() {
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

} /* namespace lamp_search */
