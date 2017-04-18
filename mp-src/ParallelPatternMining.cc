/*
 * ParallelSearch.cpp
 *
 *  Created on: Apr 17, 2017
 *      Author: yuu
 */

#include "ParallelPatternMining.h"

#include "google/gflags.h"

#include "mpi_tag.h"
#include "Log.h"
#include "FixedSizeStack.h"
#include "DTD.h"
#include "StealState.h"
//#include "SignificantSetResults.h"

#include "../src/timer.h"

DEFINE_int32(probe_period_, 128, "probe period during process node");
DEFINE_bool(probe_period_is_ms_, false,
		"true: probe period is milli sec, false: num loops");
DECLARE_bool(third_phase_); // true, "do third phase"

#ifdef __CDT_PARSER__
#undef DBG
#define DBG(a)  a
#endif

#ifdef __CDT_PARSER__
#undef LOG
#define LOG(a)  ;
#endif

namespace lamp_search {

ParallelPatternMining::ParallelPatternMining(Database<uint64> * d_,
		LampGraph<uint64> * g_, VariableBitsetHelper<uint64> * bsh_,
		MPI_Data& mpi_data, TreeSearchData* treesearch_data, Log* log,
		Timer* timer) :
		ParallelSearch(mpi_data, treesearch_data, log, timer), d_(d_), g_(g_), bsh_(
				bsh_), expand_num_(0), closed_set_num_(0), phase_(0), getminsup_data(
		NULL), gettestable_data(NULL) {

}

ParallelPatternMining::~ParallelPatternMining() {
	// TODO: lots of things to delete
}

void ParallelPatternMining::GetMinimalSupport(GetMinSupData* getminsup_data) {
	this->getminsup_data = getminsup_data;
//	CheckInit();
	phase_ = 1;
	Search();

}

void ParallelPatternMining::GetTestablePatterns(
		GetTestableData* gettestable_data) {
	this->gettestable_data = gettestable_data; // TODO: not sure if this is a good idea.
	this->getminsup_data->lambda_ = gettestable_data->freqThreshold_; // TODO: not the best way.
	CheckInitTestable();
	phase_ = 2;
	DBG(D(1) << "MainLoop" << std::endl
	;);
	printf("GetTestablePatterns\n");
	Search();

}

//==============================================================================

void ParallelPatternMining::GetSignificantPatterns(MPI_Data& mpi_data,
		GetSignificantData* getsignificant_data) {
	this->getsignificant_data = getsignificant_data;
	DBG(D(1) << "MainLoop" << std::endl
	;);
	printf("extract\n");
	ExtractSignificantSet();
	if (mpi_data.mpiRank_ == 0) {
		printf("sendresults\n");
		SendResultRequest(mpi_data);
	}
	printf("probe\n");
	while (!mpi_data.dtd_->terminated_) {
		Probe(mpi_data, treesearch_data);
	}
}

void ParallelPatternMining::Probe(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data) {
	DBG(D(3) << "Probe" << std::endl
	;);
	MPI_Status probe_status;
	int probe_src, probe_tag;

	long long int start_time;
	start_time = timer_->Elapsed();

	log_->d_.probe_num_++;

	while (CallIprobe(&probe_status, &probe_src, &probe_tag)) {
		DBG(
				D(4) << "CallIprobe returned src=" << probe_src << "\ttag="
						<< probe_tag << std::endl
				;);
		switch (probe_tag) {
		// control tasks
		case Tag::DTD_REQUEST:
			RecvDTDRequest(mpi_data, treesearch_data, probe_src);
			break;
		case Tag::DTD_REPLY:
			RecvDTDReply(mpi_data, treesearch_data, probe_src);
			break;

		case Tag::BCAST_FINISH:
			RecvBcastFinish(mpi_data, probe_src);
			break;

// basic tasks
		case Tag::REQUEST:
			RecvRequest(mpi_data, treesearch_data, probe_src);
			break;
		case Tag::REJECT:
			RecvReject(mpi_data, probe_src);
			break;
		case Tag::GIVE:
			RecvGive(mpi_data, treesearch_data, probe_src, probe_status);
			break;

			/**
			 * Messages used for SIGNIFICANCE TEST
			 */
		case Tag::RESULT_REQUEST:
			RecvResultRequest(mpi_data, probe_src);
			break;
		case Tag::RESULT_REPLY:
			RecvResultReply(mpi_data, probe_src, probe_status);
			break;
			/**
			 * Messages used ONLY for Get Minimal Support
			 */
		case Tag::DTD_ACCUM_REQUEST:
			assert(phase_ == 1);
			RecvDTDAccumRequest(mpi_data, probe_src);
			break;
		case Tag::DTD_ACCUM_REPLY:
			assert(phase_ == 1);
			RecvDTDAccumReply(mpi_data, probe_src);
			break;
		case Tag::LAMBDA:
			assert(phase_ == 1);
			RecvLambda(mpi_data, probe_src);
			break;
		default:
			printf("UNKNONW TAG\n");
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
	assert(treesearch_data->node_stack_);
	if (phase_ == 1) {
		log_->TakePeriodicLog(treesearch_data->node_stack_->NuItemset(),
				getminsup_data->lambda_, phase_);
	} else {
		log_->TakePeriodicLog(treesearch_data->node_stack_->NuItemset(),
				gettestable_data->freqThreshold_, phase_);
	}

	if (mpi_data.mpiRank_ == 0) {
		// initiate termination detection

		// note: for phase_ 1, accum request and dtd request are unified
		if (!mpi_data.echo_waiting_ && !mpi_data.dtd_->terminated_) {
			if (phase_ == 1) {
				SendDTDAccumRequest(mpi_data);
				// log_->d_.dtd_accum_request_num_++;
			} else if (phase_ == 2) {
				if (treesearch_data->node_stack_->Empty()) {
					SendDTDRequest(mpi_data);
					// log_->d_.dtd_request_num_++;
				}
			} else {
				// unknown phase
				assert(0);
			}
		}
	}

	long long int elapsed_time = timer_->Elapsed() - start_time;
	log_->d_.probe_time_ += elapsed_time;
	log_->d_.probe_time_max_ = std::max(elapsed_time, log_->d_.probe_time_max_);
}

//==============================================================================

void ParallelPatternMining::Distribute(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data) {
	if (mpi_data.nTotalProc_ == 1)
		return;
	DBG(D(3) << "Distribute" << std::endl
	;);
	if (mpi_data.thieves_->Size() > 0
			|| mpi_data.lifeline_thieves_->Size() > 0) {
		int steal_num = treesearch_data->node_stack_->Split(
				treesearch_data->give_stack_);
		if (steal_num > 0) {
			DBG(D(3) << "giving" << std::endl
			;);
			DBG(treesearch_data->give_stack_->PrintAll(D(3, false))
			;);
			Give(mpi_data, treesearch_data->give_stack_, steal_num);
			treesearch_data->give_stack_->Clear();
		}
	}
}

void ParallelPatternMining::Give(MPI_Data& mpi_data,
		VariableLengthItemsetStack * st, int steal_num) {
	DBG(D(3) << "Give: "
	;);
	if (mpi_data.thieves_->Size() > 0) { // random thieves
		int thief = mpi_data.thieves_->Pop();
		if (thief >= 0) { // lifeline thief
			DBG(
					D(3, false) << "thief=" << thief << "\trandom lifeline"
							<< std::endl
					;);
			log_->d_.lifeline_given_num_++;
			log_->d_.lifeline_nodes_given_ += steal_num;
			SendGive(mpi_data, st, thief, mpi_data.mpiRank_);
		} else { // random thief
			DBG(
					D(3, false) << "thief=" << (-thief - 1) << "\trandom"
							<< std::endl
					;);
			log_->d_.given_num_++;
			log_->d_.nodes_given_ += steal_num;
			SendGive(mpi_data, st, (-thief - 1), -1);
		}
	} else {
		assert(mpi_data.lifeline_thieves_->Size() > 0);
		int thief = mpi_data.lifeline_thieves_->Pop();
		assert(thief >= 0);
		DBG(D(3, false) << "thief=" << thief << "\tlifeline" << std::endl
		;);
		log_->d_.lifeline_given_num_++;
		log_->d_.lifeline_nodes_given_ += steal_num;
		SendGive(mpi_data, st, thief, mpi_data.mpiRank_);
	}
}

void ParallelPatternMining::Reject(MPI_Data& mpi_data) {
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

void ParallelPatternMining::Steal(MPI_Data& mpi_data) {
	DBG(D(3) << "Steal" << std::endl
	;);
	if (mpi_data.nTotalProc_ == 1)
		return;
	if (!treesearch_data->stealer_->StealStarted())
		return;
	if (treesearch_data->stealer_->Requesting())
		return;
	if (!treesearch_data->node_stack_->Empty())
		return;

	// reset state, counter

	switch (treesearch_data->stealer_->State()) {
	case StealState::RANDOM: {
		int victim = mpi_data.victims_[mpi_data.rand_m_()];
		assert(victim <= mpi_data.nTotalProc_ && "stealrandom");
		SendRequest(mpi_data, victim, -1);
		treesearch_data->stealer_->SetRequesting();
		DBG(
				D(2) << "Steal Random:" << "\trequesting="
						<< treesearch_data->stealer_->Requesting()
						<< "\trandom counter="
						<< treesearch_data->stealer_->RandomCount() << std::endl
				;);
		treesearch_data->stealer_->IncRandomCount();
		if (treesearch_data->stealer_->RandomCount() == 0)
			treesearch_data->stealer_->SetState(StealState::LIFELINE);
	}
		break;
	case StealState::LIFELINE: {
		int lifeline =
				mpi_data.lifelines_[treesearch_data->stealer_->LifelineVictim()];
		assert(lifeline >= 0); // at least lifelines_[0] must be 0 or higher
		if (!mpi_data.lifelines_activated_[lifeline]) {
			mpi_data.lifelines_activated_[lifeline] = true;
			// becomes false only at RecvGive
			assert(lifeline <= mpi_data.nTotalProc_ && "lifeline");
			SendRequest(mpi_data, lifeline, 1);
			DBG(
					D(2) << "Steal Lifeline:" << "\trequesting="
							<< treesearch_data->stealer_->Requesting()
							<< "\tlifeline counter="
							<< treesearch_data->stealer_->LifelineCounter()
							<< "\tlifeline victim="
							<< treesearch_data->stealer_->LifelineVictim()
							<< "\tz_=" << mpi_data.hypercubeDimension_
							<< std::endl
					;);
		}

		treesearch_data->stealer_->IncLifelineCount();
		treesearch_data->stealer_->IncLifelineVictim();
		if (treesearch_data->stealer_->LifelineCounter()
				>= mpi_data.hypercubeDimension_
				|| mpi_data.lifelines_[treesearch_data->stealer_->LifelineVictim()]
						< 0) { // hack fix
// note: this resets next_lifeline_victim_, is this OK?
			treesearch_data->stealer_->Finish();
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

void ParallelPatternMining::Check(MPI_Data& mpi_data) {
	if (mpi_data.mpiRank_ == 0 && phase_ == 1) {
		CheckCSThreshold(mpi_data);
	}
	return;
}

bool ParallelPatternMining::ProcessNode(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data) {
	if (treesearch_data->node_stack_->Empty())
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
				treesearch_data->node_stack_->Top(),
				treesearch_data->itemset_buf_);
		treesearch_data->node_stack_->Pop();

		// dbg
		DBG(D(3) << "expanded "
		;);
		DBG(
				treesearch_data->node_stack_->Print(D(3),
						treesearch_data->itemset_buf_)
				;);

		// TODO: THIS SHIT IS DEPENDENT ON bsh_
		// calculate support from itemset_buf_
		bsh_->Set(treesearch_data->sup_buf_);
		{
			int n = treesearch_data->node_stack_->GetItemNum(
					treesearch_data->itemset_buf_);
			for (int i = 0; i < n; i++) {
				int item = treesearch_data->node_stack_->GetNthItem(
						treesearch_data->itemset_buf_, i);
				bsh_->And(d_->NthData(item), treesearch_data->sup_buf_);
			}
		}

		int core_i = g_->CoreIndex(*treesearch_data->node_stack_,
				treesearch_data->itemset_buf_);

		int * ppc_ext_buf;
		// todo: use database reduction

		assert(
				phase_ != 1
						|| treesearch_data->node_stack_->GetItemNum(
								treesearch_data->itemset_buf_) != 0);

		bool is_root_node = (treesearch_data->node_stack_->GetItemNum(
				treesearch_data->itemset_buf_) == 0);

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
			if (treesearch_data->node_stack_->Exist(
					treesearch_data->itemset_buf_, new_item))
				continue;

			{      // Periodic probe. (do in both phases)
				accum_period_counter_++;
				if (FLAGS_probe_period_is_ms_) {      // using milli second
					if (accum_period_counter_ >= 64) {
						// to avoid calling timer_ frequently, time is checked once in 64 loops
						// clock_gettime takes 0.3--0.5 micro sec
						accum_period_counter_ = 0;
						long long int elt = timer_->Elapsed();
						if (elt - lap_time >= FLAGS_probe_period_ * 1000000) {
							log_->d_.process_node_time_ += elt - lap_time;

							Probe(mpi_data, treesearch_data);
							Distribute(mpi_data, treesearch_data);
							Reject(mpi_data);

							lap_time = timer_->Elapsed();
						}
					}
				} else {            // not using milli second
					if (accum_period_counter_ >= FLAGS_probe_period_) {
						accum_period_counter_ = 0;
						log_->d_.process_node_time_ += timer_->Elapsed()
								- lap_time;

						Probe(mpi_data, treesearch_data);
						Distribute(mpi_data, treesearch_data);
						Reject(mpi_data);

						lap_time = timer_->Elapsed();
					}
				}
				// note: do this before PushPre is called [2015-10-05 21:56]

				// todo: if database reduction is implemented,
				//       do something here for changed lambda_ (skipping new_item value ?)
			}

			bsh_->Copy(treesearch_data->sup_buf_,
					treesearch_data->child_sup_buf_);
			int sup_num = bsh_->AndCountUpdate(d_->NthData(new_item),
					treesearch_data->child_sup_buf_);

			if (sup_num < getminsup_data->lambda_)
				continue;

			treesearch_data->node_stack_->PushPre();
			ppc_ext_buf = treesearch_data->node_stack_->Top();

			bool res = g_->PPCExtension(treesearch_data->node_stack_,
					treesearch_data->itemset_buf_,
					treesearch_data->child_sup_buf_, core_i, new_item,
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
					if (true) { // XXX: FLAGS_third_phase_
						int pos_sup_num = bsh_->AndCount(d_->PosNeg(),
								treesearch_data->child_sup_buf_);
						double pval = d_->PVal(sup_num, pos_sup_num);
						assert(pval >= 0.0);
						if (pval <= gettestable_data->sig_level_) { // permits == case?
							gettestable_data->freq_stack_->PushPre();
							int * item = gettestable_data->freq_stack_->Top();
							gettestable_data->freq_stack_->CopyItem(ppc_ext_buf,
									item);
							gettestable_data->freq_stack_->PushPostNoSort();

							gettestable_data->freq_map_->insert(
									std::pair<double, int*>(pval, item));
						}
					}
				}

				assert(sup_num >= getminsup_data->lambda_);

				// try skipping if supnum_ == sup_threshold,
				// because if sup_num of a node equals to sup_threshold, children will have smaller sup_num
				// therefore no need to check it's children
				// note: skipping node_stack_ full check. allocate enough memory!
				if (sup_num <= getminsup_data->lambda_)
					treesearch_data->node_stack_->Pop();
			}
		}

		if (CheckProcessNodeEnd(mpi_data.granularity_,
				mpi_data.isGranularitySec_, processed, start_time))
			break;
	}

	long long int elapsed_time = timer_->Elapsed() - lap_time;
	log_->d_.process_node_time_ += elapsed_time;
	log_->d_.process_node_num_ += processed;

	DBG(
			D(2) << "processed node num=" << processed << "\ttime="
					<< elapsed_time << std::endl
			;);

	mpi_data.processing_node_ = false;
	return true;
}

bool ParallelPatternMining::CheckProcessNodeEnd(int n, bool n_is_ms,
		int processed, long long int start_time) {
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

// call this from root rank to start DTD
void ParallelPatternMining::SendDTDRequest(MPI_Data& mpi_data) {
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
						<< "\ttimezone=" << mpi_data.dtd_->time_zone_
						<< std::endl
				;);
	}
}

void ParallelPatternMining::RecvDTDRequest(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, int src) {
	DBG(
			D(3) << "RecvDTDRequest: src=" << src << "\ttimezone="
					<< mpi_data.dtd_->time_zone_ << std::endl
			;);
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::DTD_REQUEST, &recv_status);

	if (IsLeafInTopology(mpi_data))
		SendDTDReply(mpi_data, treesearch_data);
	else
		SendDTDRequest(mpi_data);
}

bool ParallelPatternMining::DTDReplyReady(MPI_Data& mpi_data) const {
	for (int i = 0; i < k_echo_tree_branch; i++)
		if (mpi_data.bcast_targets_[i] >= 0 && !(mpi_data.dtd_->accum_flag_[i]))
			return false;
	return true;
// check only valid bcast_targets_
// always true if leaf
}

void ParallelPatternMining::SendDTDReply(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data) {
	int message[3];
// using reduced vars
	message[0] = mpi_data.dtd_->count_ + mpi_data.dtd_->reduce_count_;
	bool tw_flag = mpi_data.dtd_->time_warp_
			|| mpi_data.dtd_->reduce_time_warp_;
	message[1] = (tw_flag ? 1 : 0);
// for Steal
	mpi_data.dtd_->not_empty_ = !(treesearch_data->node_stack_->Empty())
			|| (mpi_data.thieves_->Size() > 0)
			|| treesearch_data->stealer_->StealStarted()
			|| mpi_data.processing_node_; // thieves_ and stealer state check
// for Steal2
// dtd_.not_empty_ =
//     !(node_stack_->Empty()) || (thieves_->Size() > 0) ||
//     waiting_ || mpi_data.processing_node_;
	bool em_flag = mpi_data.dtd_->not_empty_
			|| mpi_data.dtd_->reduce_not_empty_;
	message[2] = (em_flag ? 1 : 0);
	DBG(
			D(3) << "SendDTDReply: dst = " << mpi_data.bcast_source_
					<< "\tcount=" << message[0] << "\ttw=" << tw_flag << "\tem="
					<< em_flag << std::endl
			;);

	assert(mpi_data.bcast_source_ < mpi_data.nTotalProc_ && "SendDTDReply");
	CallBsend(message, 3, MPI_INT, mpi_data.bcast_source_, Tag::DTD_REPLY);

	mpi_data.dtd_->time_warp_ = false;
	mpi_data.dtd_->not_empty_ = false;
	mpi_data.dtd_->IncTimeZone();

	mpi_data.echo_waiting_ = false;
	mpi_data.dtd_->ClearAccumFlags();
	mpi_data.dtd_->ClearReduceVars();
}

void ParallelPatternMining::RecvDTDReply(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, int src) {
	MPI_Status recv_status;
	int message[3];
	CallRecv(&message, 3, MPI_INT, src, Tag::DTD_REPLY, &recv_status);
	assert(src == recv_status.MPI_SOURCE);

// reduce reply (count, time_warp, not_empty)
	mpi_data.dtd_->Reduce(message[0], (message[1] != 0), (message[2] != 0));

	DBG(
			D(3) << "RecvDTDReply: src=" << src << "\tcount=" << message[0]
					<< "\ttw=" << message[1] << "\tem=" << message[2]
					<< "\treduced_count=" << mpi_data.dtd_->reduce_count_
					<< "\treduced_tw=" << mpi_data.dtd_->reduce_time_warp_
					<< "\treduced_em=" << mpi_data.dtd_->reduce_not_empty_
					<< std::endl
			;);

	bool flag = false;
	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data.bcast_targets_[i] == src) {
			flag = true;
			mpi_data.dtd_->accum_flag_[i] = true;
			break;
		}
	}
	assert(flag);

	if (DTDReplyReady(mpi_data)) {
		if (mpi_data.mpiRank_ == 0) {
			DTDCheck(mpi_data); // at root
			log_->d_.dtd_phase_num_++;
		} else
			SendDTDReply(mpi_data, treesearch_data);
	}
}

void ParallelPatternMining::SendDTDAccumRequest(MPI_Data& mpi_data) {
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
						<< mpi_data.dtd_->time_zone_ << std::endl
				;);
	}
}

void ParallelPatternMining::RecvDTDAccumRequest(MPI_Data& mpi_data, int src) {
	DBG(
			D(3) << "RecvDTDAccumRequest: src=" << src << "\ttimezone="
					<< mpi_data.dtd_->time_zone_ << std::endl
			;);
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::DTD_ACCUM_REQUEST, &recv_status);

	if (IsLeafInTopology(mpi_data))
		SendDTDAccumReply(mpi_data);
	else
		SendDTDAccumRequest(mpi_data);
}

void ParallelPatternMining::SendDTDAccumReply(MPI_Data& mpi_data) {
	getminsup_data->dtd_accum_array_base_[0] = mpi_data.dtd_->count_
			+ mpi_data.dtd_->reduce_count_;
	bool tw_flag = mpi_data.dtd_->time_warp_
			|| mpi_data.dtd_->reduce_time_warp_;
	getminsup_data->dtd_accum_array_base_[1] = (tw_flag ? 1 : 0);
	// for Steal
	mpi_data.dtd_->not_empty_ = !(treesearch_data->node_stack_->Empty())
			|| (mpi_data.thieves_->Size() > 0)
			|| treesearch_data->stealer_->StealStarted()
			|| mpi_data.processing_node_; // thieves_ and stealer state check
	// for Steal2
	// mpi_data.dtd_->not_empty_ =
	//     !(node_stack_->Empty()) || (thieves_->Size() > 0) ||
	//     waiting_ || mpi_data.processing_node_;
	bool em_flag = mpi_data.dtd_->not_empty_
			|| mpi_data.dtd_->reduce_not_empty_;
	getminsup_data->dtd_accum_array_base_[2] = (em_flag ? 1 : 0);

	DBG(
			D(3) << "SendDTDAccumReply: dst = " << mpi_data.bcast_source_
					<< "\tcount=" << getminsup_data->dtd_accum_array_base_[0]
					<< "\ttw=" << tw_flag << "\tem=" << em_flag << std::endl
			;);

	assert(
			mpi_data.bcast_source_ < mpi_data.nTotalProc_
					&& "SendDTDAccumReply");
	CallBsend(getminsup_data->dtd_accum_array_base_,
			getminsup_data->lambda_max_ + 4, MPI_LONG_LONG_INT,
			mpi_data.bcast_source_, Tag::DTD_ACCUM_REPLY);

	mpi_data.dtd_->time_warp_ = false;
	mpi_data.dtd_->not_empty_ = false;
	mpi_data.dtd_->IncTimeZone();

	mpi_data.echo_waiting_ = false;
	mpi_data.dtd_->ClearAccumFlags();
	mpi_data.dtd_->ClearReduceVars();

	for (int l = 0; l <= getminsup_data->lambda_max_; l++)
		getminsup_data->accum_array_[l] = 0ll;
}

// getminsup_data
void ParallelPatternMining::RecvDTDAccumReply(MPI_Data& mpi_data, int src) {
	MPI_Status recv_status;

	CallRecv(getminsup_data->dtd_accum_recv_base_,
			getminsup_data->lambda_max_ + 4, MPI_LONG_LONG_INT, src,
			Tag::DTD_ACCUM_REPLY, &recv_status);
	assert(src == recv_status.MPI_SOURCE);

	int count = (int) (getminsup_data->dtd_accum_recv_base_[0]);
	bool time_warp = (getminsup_data->dtd_accum_recv_base_[1] != 0);
	bool not_empty = (getminsup_data->dtd_accum_recv_base_[2] != 0);

	mpi_data.dtd_->Reduce(count, time_warp, not_empty);

	DBG(
			D(3) << "RecvDTDAccumReply: src=" << src << "\tcount=" << count
					<< "\ttw=" << time_warp << "\tem=" << not_empty
					<< "\treduced_count=" << mpi_data.dtd_->reduce_count_
					<< "\treduced_tw=" << mpi_data.dtd_->reduce_time_warp_
					<< "\treduced_em=" << mpi_data.dtd_->reduce_not_empty_
					<< std::endl
			;);

	for (int l = getminsup_data->lambda_ - 1; l <= getminsup_data->lambda_max_;
			l++)
		getminsup_data->accum_array_[l] += getminsup_data->accum_recv_[l];

	bool flag = false;
	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data.bcast_targets_[i] == src) {
			flag = true;
			mpi_data.dtd_->accum_flag_[i] = true;
			break;
		}
	}
	assert(flag);

	if (mpi_data.mpiRank_ == 0) {
		if (ExceedCsThr()) {
			int new_lambda = NextLambdaThr();
			SendLambda(mpi_data, new_lambda);
			getminsup_data->lambda_ = new_lambda;
		}
// if SendLambda is called, dtd_.count_ is incremented and DTDCheck will always fail
		if (DTDReplyReady(mpi_data)) {
			DTDCheck(mpi_data);
			log_->d_.dtd_accum_phase_num_++;
		}
	} else {  // not root
		if (DTDReplyReady(mpi_data))
			SendDTDAccumReply(mpi_data);
	}
}

void ParallelPatternMining::DTDCheck(MPI_Data& mpi_data) {
	assert(mpi_data.mpiRank_ == 0);
	// (count, time_warp, not_empty)
	mpi_data.dtd_->Reduce(mpi_data.dtd_->count_, mpi_data.dtd_->time_warp_,
			mpi_data.dtd_->not_empty_);

	if (mpi_data.dtd_->reduce_count_ == 0
			&& mpi_data.dtd_->reduce_time_warp_ == false
			&& mpi_data.dtd_->reduce_not_empty_ == false) {
// termination
		SendBcastFinish(mpi_data);
		mpi_data.dtd_->terminated_ = true;
		DBG(D(1) << "terminated" << std::endl
		;);
	}
	// doing same thing as SendDTDReply
	mpi_data.dtd_->time_warp_ = false;
	mpi_data.dtd_->not_empty_ = false;
	mpi_data.dtd_->IncTimeZone();

	mpi_data.echo_waiting_ = false;
	mpi_data.dtd_->ClearAccumFlags();
	mpi_data.dtd_->ClearReduceVars();
}

void ParallelPatternMining::SendBcastFinish(MPI_Data& mpi_data) {
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

void ParallelPatternMining::RecvBcastFinish(MPI_Data& mpi_data, int src) {
	DBG(D(2) << "RecvBcastFinish: src=" << src << std::endl
	;);
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::BCAST_FINISH, &recv_status);

	SendBcastFinish(mpi_data);
	mpi_data.dtd_->terminated_ = true;
	DBG(D(2) << "terminated" << std::endl
	;);

	mpi_data.waiting_ = false;
}

// lifeline == -1 for random thieves
void ParallelPatternMining::SendRequest(MPI_Data& mpi_data, int dst,
		int is_lifeline) {
	assert(dst >= 0);
	int message[2];
	message[0] = mpi_data.dtd_->time_zone_;
	message[1] = is_lifeline; // -1 for random thieves, >=0 for lifeline thieves
	assert(dst < mpi_data.nTotalProc_ && "SendRequest");

	CallBsend(message, 2, MPI_INT, dst, Tag::REQUEST);
	mpi_data.dtd_->OnSend();

	DBG(
			D(2) << "SendRequest: dst=" << dst << "\tis_lifeline="
					<< is_lifeline << "\tdtd_count=" << mpi_data.dtd_->count_
					<< std::endl
			;);
}

void ParallelPatternMining::RecvRequest(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, int src) {
	DBG(D(2) << "RecvRequest: src=" << src
	;);
	MPI_Status recv_status;
	int message[2];

	CallRecv(&message, 2, MPI_INT, src, Tag::REQUEST, &recv_status);
	mpi_data.dtd_->OnRecv();
	assert(src == recv_status.MPI_SOURCE);

	int timezone = message[0];
	mpi_data.dtd_->UpdateTimeZone(timezone);
	int is_lifeline = message[1]; // -1 for random thieves, >=0 for lifeline thieves
	int thief = src;
	DBG(D(2, false) << "\tis_lifeline=" << is_lifeline
	;);
	DBG(D(2, false) << "\tthief=" << src
	;);
	DBG(D(2, false) << "\tdtd_count=" << mpi_data.dtd_->count_
	;);

	if (treesearch_data->node_stack_->Empty()) {
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

void ParallelPatternMining::SendReject(MPI_Data& mpi_data, int dst) {
	assert(dst >= 0);
	int message[1];

	message[0] = mpi_data.dtd_->time_zone_;
	assert(dst < mpi_data.nTotalProc_ && "SendReject");
	CallBsend(message, 1, MPI_INT, dst, Tag::REJECT);
	mpi_data.dtd_->OnSend();

	DBG(
			D(2) << "SendReject: dst=" << dst << "\tdtd_count="
					<< mpi_data.dtd_->count_ << std::endl
			;);
}

void ParallelPatternMining::RecvReject(MPI_Data& mpi_data, int src) {
	MPI_Status recv_status;
	int message[1];

	CallRecv(&message, 1, MPI_INT, src, Tag::REJECT, &recv_status);
	mpi_data.dtd_->OnRecv();
	assert(src == recv_status.MPI_SOURCE);

	int timezone = message[0];
	mpi_data.dtd_->UpdateTimeZone(timezone);
	DBG(
			D(2) << "RecvReject: src=" << src << "\tdtd_count="
					<< mpi_data.dtd_->count_ << std::endl
			;);

	treesearch_data->stealer_->ResetRequesting();
	mpi_data.waiting_ = false;
}

void ParallelPatternMining::SendGive(MPI_Data& mpi_data,
		VariableLengthItemsetStack * st, int dst, int is_lifeline) {
	assert(dst >= 0);
	st->SetTimestamp(mpi_data.dtd_->time_zone_);
	st->SetFlag(is_lifeline);
	int * message = st->Stack();
	int size = st->UsedCapacity();

	log_->d_.give_stack_max_itm_ = std::max(log_->d_.give_stack_max_itm_,
			(long long int) (st->NuItemset()));
	log_->d_.give_stack_max_cap_ = std::max(log_->d_.give_stack_max_cap_,
			(long long int) (st->UsedCapacity()));

	assert(dst < mpi_data.nTotalProc_ && "SendGive");
	CallBsend(message, size, MPI_INT, dst, Tag::GIVE);
	mpi_data.dtd_->OnSend();

	DBG(
			D(2) << "SendGive: " << "\ttimezone=" << mpi_data.dtd_->time_zone_
					<< "\tdst=" << dst << "\tlfl=" << is_lifeline << "\tsize="
					<< size << "\tnode=" << st->NuItemset() << "\tdtd_count="
					<< mpi_data.dtd_->count_ << std::endl
			;);
	//st->PrintAll(D(false));
}

void ParallelPatternMining::RecvGive(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, int src, MPI_Status probe_status) {
	int count;
	int error = MPI_Get_count(&probe_status, MPI_INT, &count);
	if (error != MPI_SUCCESS) {
		DBG(D(1) << "error in MPI_Get_count in RecvGive: " << error << std::endl
		;);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	MPI_Status recv_status;

	CallRecv(treesearch_data->give_stack_->Stack(), count, MPI_INT, src,
			Tag::GIVE, &recv_status);
	mpi_data.dtd_->OnRecv();
	assert(src == recv_status.MPI_SOURCE);

	int timezone = treesearch_data->give_stack_->Timestamp();
	mpi_data.dtd_->UpdateTimeZone(timezone);

	int flag = treesearch_data->give_stack_->Flag();
	int orig_nu_itemset = treesearch_data->node_stack_->NuItemset();

	treesearch_data->node_stack_->MergeStack(
			treesearch_data->give_stack_->Stack()
					+ VariableLengthItemsetStack::SENTINEL + 1,
			count - VariableLengthItemsetStack::SENTINEL - 1);
	int new_nu_itemset = treesearch_data->node_stack_->NuItemset();

	if (flag >= 0) {
		mpi_data.lifelines_activated_[src] = false;
		log_->d_.lifeline_steal_num_++;
		log_->d_.lifeline_nodes_received_ += (new_nu_itemset - orig_nu_itemset);
	} else {
		log_->d_.steal_num_++;
		log_->d_.nodes_received_ += (new_nu_itemset - orig_nu_itemset);
	}

	DBG(
			D(2) << "RecvGive: src=" << src << "\ttimezone="
					<< mpi_data.dtd_->time_zone_ << "\tlfl=" << flag
					<< "\tsize=" << count << "\tnode="
					<< (new_nu_itemset - orig_nu_itemset) << "\tdtd_count="
					<< mpi_data.dtd_->count_ << std::endl
			;);

	treesearch_data->give_stack_->Clear();

	log_->d_.node_stack_max_itm_ = std::max(log_->d_.node_stack_max_itm_,
			(long long int) (treesearch_data->node_stack_->NuItemset()));
	log_->d_.node_stack_max_cap_ = std::max(log_->d_.node_stack_max_cap_,
			(long long int) (treesearch_data->node_stack_->UsedCapacity()));

	DBG(treesearch_data->node_stack_->PrintAll(D(3, false))
	;);

	treesearch_data->stealer_->ResetRequesting();
	treesearch_data->stealer_->ResetCounters();
	treesearch_data->stealer_->SetState(StealState::RANDOM);
	treesearch_data->stealer_->SetStealStart();
	mpi_data.waiting_ = false;
}

//==============================================================================
//
//void ParallelSearch::SendResultRequest(MPI_Data& mpi_data) {
//	int message[1];
//	message[0] = 1; // dummy
//
//	mpi_data.echo_waiting_ = true;
//
//	for (int i = 0; i < k_echo_tree_branch; i++) {
//		if (mpi_data.bcast_targets_[i] < 0)
//			break;
//
//		assert(
//				mpi_data.bcast_targets_[i] < mpi_data.nTotalProc_
//						&& "SendResultRequest");
//		CallBsend(message, 1, MPI_INT, mpi_data.bcast_targets_[i],
//				Tag::RESULT_REQUEST);
//		DBG(
//				D(2) << "SendResultRequest: dst=" << mpi_data.bcast_targets_[i]
//						<< std::endl
//				;);
//	}
//}

//void ParallelSearch::RecvResultRequest(MPI_Data& mpi_data, int src) {
//	MPI_Status recv_status;
//	int message[1];
//	CallRecv(&message, 1, MPI_INT, src, Tag::RESULT_REQUEST, &recv_status);
//	assert(src == recv_status.MPI_SOURCE);
//
//	DBG(D(2) << "RecvResultRequest: src=" << src << std::endl
//	;);
//
//	if (IsLeafInTopology(mpi_data))
//		SendResultReply(mpi_data);
//	else
//		SendResultRequest(mpi_data);
//}
//
//void ParallelSearch::SendResultReply(MPI_Data& mpi_data) {
//	int * message = significant_stack_->Stack();
//	int size = significant_stack_->UsedCapacity();
//	assert(mpi_data.bcast_source_ < mpi_data.nTotalProc_ && "SendResultReply");
//	CallBsend(message, size, MPI_INT, mpi_data.bcast_source_,
//			Tag::RESULT_REPLY);
//
//	DBG(D(2) << "SendResultReply: dst=" << mpi_data.bcast_source_ << std::endl
//	;);
//	DBG(significant_stack_->PrintAll(D(3, false))
//	;);
//
//	mpi_data.echo_waiting_ = false;
//	mpi_data.dtd_->terminated_ = true;
//}
//
//void ParallelSearch::RecvResultReply(MPI_Data& mpi_data, int src,
//		MPI_Status probe_status) {
//	int count;
//	int error = MPI_Get_count(&probe_status, MPI_INT, &count);
//	if (error != MPI_SUCCESS) {
//		DBG(
//				D(1) << "error in MPI_Get_count in RecvResultReply: " << error
//						<< std::endl
//				;);
//		MPI_Abort(MPI_COMM_WORLD, 1);
//	}
//
//	MPI_Status recv_status;
//	CallRecv(give_stack_->Stack(), count, MPI_INT, src, Tag::RESULT_REPLY,
//			&recv_status);
//	assert(src == recv_status.MPI_SOURCE);
//
//	significant_stack_->MergeStack(
//			give_stack_->Stack() + VariableLengthItemsetStack::SENTINEL + 1,
//			count - VariableLengthItemsetStack::SENTINEL - 1);
//	give_stack_->Clear();
//
//	DBG(D(2) << "RecvResultReply: src=" << src << std::endl
//	;);
//	DBG(significant_stack_->PrintAll(D(3, false))
//	;);
//
//	// using the same flags as accum count, should be fixed
//	bool flag = false;
//	for (int i = 0; i < k_echo_tree_branch; i++) {
//		if (mpi_data.bcast_targets_[i] == src) {
//			flag = true;
//			mpi_data.accum_flag_[i] = true;
//			break;
//		}
//	}
//	assert(flag);
//
//	if (AccumCountReady(mpi_data)) {
//		if (mpi_data.mpiRank_ != 0) {
//			SendResultReply(mpi_data);
//		} else { // root
//			mpi_data.echo_waiting_ = false;
//			mpi_data.dtd_->terminated_ = true;
//		}
//	}
//}

//==============================================================================

bool ParallelPatternMining::IsLeafInTopology(MPI_Data& mpi_data) const {
	for (int i = 0; i < k_echo_tree_branch; i++)
		if (mpi_data.bcast_targets_[i] >= 0)
			return false;
	return true;
}

int ParallelPatternMining::CallIprobe(MPI_Status * status, int * src,
		int * tag) {
	long long int start_time;
	long long int end_time;
	log_->d_.iprobe_num_++;

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
			end_time = timer_->Elapsed(); log_->d_.iprobe_time_ += end_time - start_time; log_->d_.iprobe_time_max_ = std::max(end_time - start_time, log_->d_.iprobe_time_max_););

	if (flag) {
		log_->d_.iprobe_succ_num_++;
		LOG(
				log_->d_.iprobe_succ_time_ += end_time - start_time; log_->d_.iprobe_succ_time_max_ = std::max(end_time - start_time, log_->d_.iprobe_succ_time_max_););

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
		log_->d_.iprobe_fail_num_++;
		LOG(
				log_->d_.iprobe_fail_time_ += end_time - start_time; log_->d_.iprobe_fail_time_max_ = std::max(end_time - start_time, log_->d_.iprobe_fail_time_max_););
	}

	return flag;
}

int ParallelPatternMining::CallRecv(void * buffer, int data_count,
		MPI_Datatype type, int src, int tag, MPI_Status * status) {
	long long int start_time;
	long long int end_time;

	// todo: prepare non-log mode to remove measurement
	// clock_gettime takes 0.3--0.5 micro sec
	log_->d_.recv_num_++;
	LOG(start_time = timer_->Elapsed(););

	int error = MPI_Recv(buffer, data_count, type, src, tag, MPI_COMM_WORLD,
			status);

	LOG(
			end_time = timer_->Elapsed(); log_->d_.recv_time_ += end_time - start_time; log_->d_.recv_time_max_ = std::max(end_time - start_time, log_->d_.recv_time_max_););
	return error;
}

int ParallelPatternMining::CallBsend(void * buffer, int data_count,
		MPI_Datatype type, int dest, int tag) {
//	assert(0 <= dest && dest < mpi_data.nTotalProc_);
	long long int start_time;
	long long int end_time;
	log_->d_.bsend_num_++;
	start_time = timer_->Elapsed();

	int error = MPI_Bsend(buffer, data_count, type, dest, tag, MPI_COMM_WORLD);

	end_time = timer_->Elapsed();
	log_->d_.bsend_time_ += end_time - start_time;
	log_->d_.bsend_time_max_ = std::max(end_time - start_time,
			log_->d_.bsend_time_max_);
	return error;
}

int ParallelPatternMining::CallBcast(void * buffer, int data_count,
		MPI_Datatype type) {
	long long int start_time;
	long long int end_time;
	log_->d_.bcast_num_++;
	start_time = timer_->Elapsed();

	int error = MPI_Bcast(buffer, data_count, type, 0, MPI_COMM_WORLD);

	end_time = timer_->Elapsed();
	log_->d_.bcast_time_ += end_time - start_time;
	log_->d_.bcast_time_max_ = std::max(end_time - start_time,
			log_->d_.bcast_time_max_);
	return error;
}

/**
 * GetMinSup Functions
 *
 *
 */
void ParallelPatternMining::SendLambda(MPI_Data& mpi_data, int lambda) {
	// send lambda to bcast_targets_
	int message[2];
	message[0] = mpi_data.dtd_->time_zone_;
	message[1] = lambda;

	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data.bcast_targets_[i] < 0)
			break;
		assert(
				mpi_data.bcast_targets_[i] < mpi_data.nTotalProc_
						&& "SendLambda");
		CallBsend(message, 2, MPI_INT, mpi_data.bcast_targets_[i], Tag::LAMBDA);
		mpi_data.dtd_->OnSend();

		DBG(
				D(2) << "SendLambda: dst=" << mpi_data.bcast_targets_[i]
						<< "\tlambda=" << lambda << "\tdtd_count="
						<< mpi_data.dtd_->count_ << std::endl
				;);
	}
}

void ParallelPatternMining::RecvLambda(MPI_Data& mpi_data, int src) {
	MPI_Status recv_status;
	int message[2];

	CallRecv(&message, 2, MPI_INT, src, Tag::LAMBDA, &recv_status);
	mpi_data.dtd_->OnRecv();
	assert(src == recv_status.MPI_SOURCE);
	int timezone = message[0];
	mpi_data.dtd_->UpdateTimeZone(timezone);

	DBG(
			D(2) << "RecvLambda: src=" << src << "\tlambda=" << message[1]
					<< "\tdtd_count=" << mpi_data.dtd_->count_ << std::endl
			;);

	int new_lambda = message[1];
	if (new_lambda > getminsup_data->lambda_) {
		SendLambda(mpi_data, new_lambda);
		getminsup_data->lambda_ = new_lambda;
		// todo: do database reduction
	}
}
bool ParallelPatternMining::AccumCountReady(MPI_Data& mpi_data) const {
	for (int i = 0; i < k_echo_tree_branch; i++)
		if (mpi_data.bcast_targets_[i] >= 0 && !(mpi_data.accum_flag_[i]))
			return false;
	return true;
	// check only valid bcast_targets_
	// always true if leaf
}

void ParallelPatternMining::CheckCSThreshold(MPI_Data& mpi_data) {
//	assert(mpi_data.mpiRank_ == 0);
	if (ExceedCsThr()) {
		int new_lambda = NextLambdaThr();
		SendLambda(mpi_data, new_lambda);
		getminsup_data->lambda_ = new_lambda;
	}
}

void ParallelPatternMining::IncCsAccum(int sup_num) {
	for (int i = sup_num; i >= getminsup_data->lambda_ - 1; i--)
		getminsup_data->accum_array_[i]++;
}

bool ParallelPatternMining::ExceedCsThr() const {
	// note: > is correct. permit ==
	return (getminsup_data->accum_array_[getminsup_data->lambda_]
			> getminsup_data->cs_thr_[getminsup_data->lambda_]);
}

int ParallelPatternMining::NextLambdaThr() const {
	int si;
	for (si = getminsup_data->lambda_max_; si >= getminsup_data->lambda_; si--)
		if (getminsup_data->accum_array_[si] > getminsup_data->cs_thr_[si])
			break;
	return si + 1;
	// it is safe because lambda_ higher than max results in immediate search finish
}

/**
 * GETSIGNIFICANT PATTERNS
 *
 *
 *
 *
 *
 *
 *
 *
 */

//==============================================================================
void ParallelPatternMining::SendResultRequest(MPI_Data& mpi_data) {
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

void ParallelPatternMining::RecvResultRequest(MPI_Data& mpi_data, int src) {
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::RESULT_REQUEST, &recv_status);
	assert(src == recv_status.MPI_SOURCE);

	DBG(D(2) << "RecvResultRequest: src=" << src << std::endl
	;);

	if (IsLeafInTopology(mpi_data))
		SendResultReply(mpi_data);
	else
		SendResultRequest(mpi_data);
}

void ParallelPatternMining::SendResultReply(MPI_Data& mpi_data) {
	int * message = getsignificant_data->significant_stack_->Stack();
	int size = getsignificant_data->significant_stack_->UsedCapacity();
	assert(mpi_data.bcast_source_ < mpi_data.nTotalProc_ && "SendResultReply");
	CallBsend(message, size, MPI_INT, mpi_data.bcast_source_,
			Tag::RESULT_REPLY);

	DBG(D(2) << "SendResultReply: dst=" << mpi_data.bcast_source_ << std::endl
	;);
	DBG(getsignificant_data->significant_stack_->PrintAll(D(3, false))
	;);

	mpi_data.echo_waiting_ = false;
	mpi_data.dtd_->terminated_ = true;
}

void ParallelPatternMining::RecvResultReply(MPI_Data& mpi_data, int src,
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
	CallRecv(treesearch_data->give_stack_->Stack(), count, MPI_INT, src,
			Tag::RESULT_REPLY, &recv_status);
	assert(src == recv_status.MPI_SOURCE);

	getsignificant_data->significant_stack_->MergeStack(
			treesearch_data->give_stack_->Stack()
					+ VariableLengthItemsetStack::SENTINEL + 1,
			count - VariableLengthItemsetStack::SENTINEL - 1);
	treesearch_data->give_stack_->Clear();

	DBG(D(2) << "RecvResultReply: src=" << src << std::endl
	;);
	DBG(getsignificant_data->significant_stack_->PrintAll(D(3, false))
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
			mpi_data.dtd_->terminated_ = true;
		}
	}
}

void ParallelPatternMining::ExtractSignificantSet() {
	std::multimap<double, int *>::iterator it;
	for (it = getsignificant_data->freq_map_->begin();
			it != getsignificant_data->freq_map_->end(); ++it) {
		// permits == case
		if ((*it).first <= getsignificant_data->final_sig_level_) {
			getsignificant_data->significant_stack_->PushPre();
			int * item = getsignificant_data->significant_stack_->Top();
			getsignificant_data->significant_stack_->CopyItem((*it).second,
					item);
			getsignificant_data->significant_stack_->PushPostNoSort();
		} else
			break;
	}
}

void ParallelPatternMining::SortSignificantSets() {
	int * set = getsignificant_data->significant_stack_->FirstItemset();

	while (set != NULL) {
		// calculate support from set
		bsh_->Set(treesearch_data->sup_buf_);
		{
			int n = getsignificant_data->significant_stack_->GetItemNum(set);
			for (int i = 0; i < n; i++) {
				int item = getsignificant_data->significant_stack_->GetNthItem(
						set, i);
				bsh_->And(d_->NthData(item), treesearch_data->sup_buf_);
			}
		}

		int sup_num = bsh_->Count(treesearch_data->sup_buf_);
		int pos_sup_num = bsh_->AndCount(d_->PosNeg(),
				treesearch_data->sup_buf_);
		double pval = d_->PVal(sup_num, pos_sup_num);

		getsignificant_data->significant_set_->insert(
				SignificantSetResult(pval, set, sup_num, pos_sup_num,
						getsignificant_data->significant_stack_));
		set = getsignificant_data->significant_stack_->NextItemset(set);
	}
}

void ParallelPatternMining::CheckInit() {
	assert(mpi_data.echo_waiting_ == false);
	assert(mpi_data.waiting_ == false);
	assert(mpi_data.thieves_->Size() == 0);
	assert(mpi_data.lifeline_thieves_->Size() == 0);
	for (int pi = 0; pi < mpi_data.nTotalProc_; pi++)
		assert(mpi_data.lifelines_activated_[pi] == false);
}

// TODO: For debugging: Not sure how to show the caller in assert.
void ParallelPatternMining::CheckInitTestable() {
	mpi_data.lifeline_thieves_->Clear();
	for (int pi = 0; pi < mpi_data.nTotalProc_; pi++) {
		mpi_data.lifelines_activated_[pi] = false;
	}

	assert(mpi_data.echo_waiting_ == false);
	assert(mpi_data.waiting_ == false);
	assert(mpi_data.thieves_->Size() == 0);
	assert(mpi_data.lifeline_thieves_->Size() == 0);
	for (int pi = 0; pi < mpi_data.nTotalProc_; pi++) {
		assert(mpi_data.lifelines_activated_[pi] == false);
	}
}

} /* namespace lamp_search */
