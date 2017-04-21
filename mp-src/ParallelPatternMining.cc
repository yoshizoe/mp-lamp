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

ParallelPatternMining::ParallelPatternMining(BinaryPatternMiningData* bpm_data,
		MPI_Data& mpi_data, TreeSearchData* treesearch_data, Log* log,
		Timer* timer) :
		ParallelDFS(mpi_data, treesearch_data, log, timer), d_(bpm_data->d_), bsh_(
				bpm_data->bsh_), sup_buf_(bpm_data->sup_buf_), child_sup_buf_(
				bpm_data->child_sup_buf_), expand_num_(0), closed_set_num_(0), phase_(
				0), getminsup_data(
		NULL), gettestable_data(NULL) {
	g_ = new LampGraph<uint64>(*d_); // No overhead to generate LampGraph.
}

ParallelPatternMining::~ParallelPatternMining() {
	// TODO: lots of things to delete
	if (g_)
		delete g_;
}

// TODO: This should be under GetMinimumSupport
void ParallelPatternMining::PreProcessRootNode(GetMinSupData* getminsup_data) {
	this->getminsup_data = getminsup_data;
	phase_ = 1;
	long long int start_time;
	start_time = timer_->Elapsed();

	expand_num_++;

	// what needed is setting itemset_buf_ to empty itemset
	// can be done faster
	// assuming root itemset is pushed before
	treesearch_data->node_stack_->CopyItem(treesearch_data->node_stack_->Top(),
			treesearch_data->itemset_buf_);
	treesearch_data->node_stack_->Pop();

	// dbg
	DBG(D(2) << "preprocess root node "
	;);
	DBG(treesearch_data->node_stack_->Print(D(2), treesearch_data->itemset_buf_)
	;);

	// calculate support from itemset_buf_
	bsh_->Set(sup_buf_);
	// skipping rest for root node

	int core_i = g_->CoreIndex(*treesearch_data->node_stack_,
			treesearch_data->itemset_buf_);

	int * ppc_ext_buf;
	// todo: use database reduction

	bool is_root_node = true;

	// reverse order
	for (int new_item = d_->NextItemInReverseLoop(is_root_node,
			mpi_data.mpiRank_, mpi_data.nTotalProc_, d_->NuItems());
			new_item >= core_i + 1;
			new_item = d_->NextItemInReverseLoop(is_root_node,
					mpi_data.mpiRank_, mpi_data.nTotalProc_, new_item)) {
		// skipping not needed because itemset_buf_ if root itemset

		bsh_->Copy(sup_buf_, child_sup_buf_);
		int sup_num = bsh_->AndCountUpdate(d_->NthData(new_item),
				child_sup_buf_);

		if (sup_num < getminsup_data->lambda_)
			continue;

		treesearch_data->node_stack_->PushPre();
		ppc_ext_buf = treesearch_data->node_stack_->Top();
//		treesearch_data->node_stack_->CopyItem(treesearch_data->itemset_buf_,
//				ppc_ext_buf);
		bool res = g_->PPCExtension(treesearch_data->node_stack_,
				treesearch_data->itemset_buf_, child_sup_buf_, core_i, new_item,
				ppc_ext_buf);

		treesearch_data->node_stack_->SetSup(ppc_ext_buf, sup_num);
		treesearch_data->node_stack_->PushPostNoSort();

		treesearch_data->node_stack_->Pop(); // always pop

		if (res) {
			IncCsAccum(sup_num); // increment closed_set_num_array
			assert(sup_num >= getminsup_data->lambda_);
			if (ExceedCsThr())
				getminsup_data->lambda_ = NextLambdaThr();
		}
	}
//	 TODO: lambda_max_ is wrong???
	printf("lambda_max_ = %d\n", getminsup_data->lambda_max_);
	MPI_Reduce(getminsup_data->accum_array_, getminsup_data->accum_recv_,
			getminsup_data->lambda_max_ + 1, MPI_LONG_LONG_INT,
			MPI_SUM, 0, MPI_COMM_WORLD); // error?

	if (mpi_data.mpiRank_ == 0) {
		for (int l = 0; l <= getminsup_data->lambda_max_; l++) {
			getminsup_data->accum_array_[l] = getminsup_data->accum_recv_[l]; // overwrite here
			getminsup_data->accum_recv_[l] = 0;
		}
	} else { // h_!=0
		for (int l = 0; l <= getminsup_data->lambda_max_; l++) {
			getminsup_data->accum_array_[l] = 0;
			getminsup_data->accum_recv_[l] = 0;
		}
	}

	if (mpi_data.mpiRank_ == 0)
		if (ExceedCsThr())
			getminsup_data->lambda_ = NextLambdaThr();
	CallBcast(&getminsup_data->lambda_, 1, MPI_INT);

	// reverse order
	for (int new_item = d_->NextItemInReverseLoop(is_root_node,
			mpi_data.mpiRank_, mpi_data.nTotalProc_, d_->NuItems());
			new_item >= core_i + 1;
			new_item = d_->NextItemInReverseLoop(is_root_node,
					mpi_data.mpiRank_, mpi_data.nTotalProc_, new_item)) {
		// skipping not needed because itemset_buf_ if root itemset

		bsh_->Copy(sup_buf_, child_sup_buf_);
		int sup_num = bsh_->AndCountUpdate(d_->NthData(new_item),
				child_sup_buf_);

		if (sup_num < getminsup_data->lambda_)
			continue;

		treesearch_data->node_stack_->PushPre();
		ppc_ext_buf = treesearch_data->node_stack_->Top();
//		treesearch_data->node_stack_->CopyItem(treesearch_data->itemset_buf_,
//				ppc_ext_buf);

		bool res = g_->PPCExtension(treesearch_data->node_stack_,
				treesearch_data->itemset_buf_, child_sup_buf_, core_i, new_item,
				ppc_ext_buf);

		treesearch_data->node_stack_->SetSup(ppc_ext_buf, sup_num);
		treesearch_data->node_stack_->PushPostNoSort();

		if (!res) { // todo: remove this redundancy
			treesearch_data->node_stack_->Pop();
		} else {
			treesearch_data->node_stack_->SortTop();
			// note: IncCsAccum already done above

			assert(sup_num >= getminsup_data->lambda_);
			if (sup_num <= getminsup_data->lambda_)
				treesearch_data->node_stack_->Pop();
		}
	}

	long long int elapsed_time = timer_->Elapsed() - start_time;
	log_->d_.preprocess_time_ += elapsed_time;

	DBG(
			D(2) << "preprocess root node finished" << "\tlambda="
					<< getminsup_data->lambda_ << "\ttime=" << elapsed_time
					<< std::endl
			;);
}

void ParallelPatternMining::GetMinimalSupport(GetMinSupData* getminsup_data) {
	this->getminsup_data = getminsup_data;
//	CheckInit();
	phase_ = 1; // TODO: remove dependency on this
	Search();

	// return lambda?
}

void ParallelPatternMining::GetTestablePatterns(
		GetTestableData* gettestable_data) {
	this->gettestable_data = gettestable_data; // TODO: not sure if this is a good idea.
	this->getminsup_data->lambda_ = gettestable_data->freqThreshold_; // TODO: not the best way.
	CheckInitTestable();
	phase_ = 2;  // TODO: remove dependency on this
	DBG(D(1) << "MainLoop" << std::endl
	;);
	printf("GetTestablePatterns\n");
	Search();

}

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

//==============================================================================
/**
 * Procedure to run after Probe().
 */
void ParallelPatternMining::ProcAfterProbe() {
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
}

void ParallelPatternMining::ProbeExecute(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, MPI_Status* probe_status,
		int probe_src, int probe_tag) {
	switch (probe_tag) {
	/**
	 * MINIMALSUPPORT
	 */
	case Tag::DTD_ACCUM_REQUEST:
		assert(phase_ == 1);
		RecvDTDAccumRequest(mpi_data, probe_src); // TODO: This can be solved by polymorphism.
		break;
	case Tag::DTD_ACCUM_REPLY:
		assert(phase_ == 1);
		RecvDTDAccumReply(mpi_data, probe_src);
		break;
	case Tag::LAMBDA:
		assert(phase_ == 1);
		RecvLambda(mpi_data, probe_src);
		break;

		/**
		 * SIGNIFICANTSET
		 */
	case Tag::RESULT_REQUEST:
		RecvResultRequest(mpi_data, probe_src);
		break;
	case Tag::RESULT_REPLY:
		RecvResultReply(mpi_data, probe_src, *probe_status);
		break;

	default:
		DBG(
				D(1) << "unknown Tag for ParallelPatternMining=" << probe_tag
						<< " received in Probe: " << std::endl
				;
		);
		ParallelDFS::ProbeExecute(mpi_data, treesearch_data, probe_status,
				probe_src, probe_tag);
		break;
	}
	return;
}

/**
 * Domain specifics
 */

/**
 * GetMinSup specific
 */
void ParallelPatternMining::Check(MPI_Data& mpi_data) {
	if (mpi_data.mpiRank_ == 0 && phase_ == 1) {
		CheckCSThreshold(mpi_data);
	}
	return;
}

// TODO: What does it do after all?
bool ParallelPatternMining::ExpandNode(MPI_Data& mpi_data,
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

		PopNodeFromStack();

		assert(
				phase_ != 1
						|| treesearch_data->node_stack_->GetItemNum(
								treesearch_data->itemset_buf_) != 0);

		int accum_period_counter_ = 0;

		// TODO: core_i should be hidden
		int core_i = g_->CoreIndex(*treesearch_data->node_stack_,
				treesearch_data->itemset_buf_);
		std::vector<int> children = GetChildren(core_i);
		for (std::vector<int>::iterator it = children.begin();
				it != children.end(); ++it) {
			int new_item = *it;

			if (treesearch_data->node_stack_->Exist(
					treesearch_data->itemset_buf_, new_item)) {
				return false;
			}
			// skip existing item
			// todo: improve speed here
			CheckProbe(accum_period_counter_, lap_time);

			if (!TestAndPushNode(new_item, core_i)) {
				continue;
			}

			int* ppc_ext_buf = treesearch_data->node_stack_->Top();
			int sup_num = treesearch_data->node_stack_->GetSup(ppc_ext_buf);

			DBG(if (phase_ == 2) {
				D(3) << "found cs "
				;
				treesearch_data->node_stack_->Print(D(3), ppc_ext_buf)
				;
			});

			ProcessNode(sup_num, ppc_ext_buf);

			assert(sup_num >= getminsup_data->lambda_);

			// try skipping if supnum_ == sup_threshold,
			// because if sup_num of a node equals to sup_threshold, children will have smaller sup_num
			// therefore no need to check it's children
			// note: skipping node_stack_ full check. allocate enough memory!
			if (sup_num == getminsup_data->lambda_) {
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

// TODO: Dependent on d_. This function should be responsible of Domain class.
// TODO: Should make clear distinction of LAMP children and immediate children?
std::vector<int> ParallelPatternMining::GetChildren(int core_i) {
	bool is_root_node = (treesearch_data->node_stack_->GetItemNum(
			treesearch_data->itemset_buf_) == 0);

	// reverse order
	// for ( int new_item = d_->NuItems()-1 ; new_item >= core_i+1 ; new_item-- )
	std::vector<int> children;
	for (int new_item = d_->NextItemInReverseLoop(is_root_node,
			mpi_data.mpiRank_, mpi_data.nTotalProc_, d_->NuItems());
			new_item >= core_i + 1;
			new_item = d_->NextItemInReverseLoop(is_root_node,
					mpi_data.mpiRank_, mpi_data.nTotalProc_, new_item)) {
		children.push_back(new_item);
	}
	return children;
}

/**
 * itemset_buf_ := pointer to the index of itemset
 * sup_buf_     := support* of the itemset (takes AND for each item)
 * support*     := transactions which include all items in the itemset
 */
void ParallelPatternMining::PopNodeFromStack() {
	/**
	 * Pop item index and put into itemset_buf.
	 * This part is not dependent on the feature.
	 */
	treesearch_data->node_stack_->CopyItem(treesearch_data->node_stack_->Top(),
			treesearch_data->itemset_buf_);
	treesearch_data->node_stack_->Pop();

// dbg
	DBG(D(3) << "expanded "
	;);
	DBG(treesearch_data->node_stack_->Print(D(3), treesearch_data->itemset_buf_)
	;);

	/**
	 * Get the support (frequency) of the item into sup_buf_
	 */
// TODO: This part is dependent on feature.
	// Why is it uint64???
	bsh_->Set(sup_buf_);
	{
		// This should definitely be inside the Domain class.
		// n := number of items in the itemset.
		int n = treesearch_data->node_stack_->GetItemNum(
				treesearch_data->itemset_buf_);
		for (int i = 0; i < n; i++) {
			int item = treesearch_data->node_stack_->GetNthItem(
					treesearch_data->itemset_buf_, i);
			bsh_->And(d_->NthData(item), sup_buf_);
		}
	}
}

/**
 * Test two pruning criteria before pushing
 * 1. The support of the itemset should be larger or equal to lambda (minimal support).
 * 2. New node should be a PPC extension of the parent itemset.
 */
bool ParallelPatternMining::TestAndPushNode(int new_item, int core_i) {
//	// todo: if database reduction is implemented,
//	//       do something here for changed lambda_ (skipping new_item value ?)

	bsh_->Copy(sup_buf_, child_sup_buf_);
	int sup_num = bsh_->AndCountUpdate(d_->NthData(new_item), child_sup_buf_);
	// If the support is smaller than the required minimal support for
	// significant pattern (=lambda), then prune it.
	if (sup_num < getminsup_data->lambda_) {
		return false;
	}
	// TODO: Here it is inserting a new item into the node_stack_.
	treesearch_data->node_stack_->PushPre();
	int* ppc_ext_buf = treesearch_data->node_stack_->Top();

	// TODO: return true if child_sup_buf_ is PPC extension???
	bool res = g_->PPCExtension(treesearch_data->node_stack_,
			treesearch_data->itemset_buf_, child_sup_buf_, core_i, new_item,
			ppc_ext_buf);

	// TODO: Those things should be done in other class.
	treesearch_data->node_stack_->SetSup(ppc_ext_buf, sup_num);
	treesearch_data->node_stack_->PushPostNoSort();

	if (!res) {        // todo: remove this redundancy
		treesearch_data->node_stack_->Pop();
		return false;
	}

	treesearch_data->node_stack_->SortTop();

	return true;
}

void ParallelPatternMining::ProcessNode(int sup_num, int* ppc_ext_buf) {
	if (phase_ == 1)
		IncCsAccum(sup_num); // increment closed_set_num_array
	if (phase_ == 2) {
		closed_set_num_++;
		if (true) { // XXX: FLAGS_third_phase_
			int pos_sup_num = bsh_->AndCount(d_->PosNeg(), child_sup_buf_);
			double pval = d_->PVal(sup_num, pos_sup_num);
			assert(pval >= 0.0);
			if (pval <= gettestable_data->sig_level_) { // permits == case?
				gettestable_data->freq_stack_->PushPre();
				int * item = gettestable_data->freq_stack_->Top();
				gettestable_data->freq_stack_->CopyItem(ppc_ext_buf, item);
				gettestable_data->freq_stack_->PushPostNoSort();

				gettestable_data->freq_map_->insert(
						std::pair<double, int*>(pval, item));
			}
		}
	}
}

void ParallelPatternMining::CheckProbe(int& accum_period_counter_,
		long long int lap_time) {
// TODO: whatever this is trying to do, it should be factored into a function.
//       Why is it Probing while in the expansion loop?
	accum_period_counter_++;
	if (FLAGS_probe_period_is_ms_) {      // using milli second
		if (accum_period_counter_ >= 64) { // TODO: what is this magic number?
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
			log_->d_.process_node_time_ += timer_->Elapsed() - lap_time;

			Probe(mpi_data, treesearch_data);
			Distribute(mpi_data, treesearch_data);
			Reject(mpi_data);

			lap_time = timer_->Elapsed();
		}
	}
// note: do this before PushPre is called [2015-10-05 21:56]
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

// TODO: polymorphism to override SendDTDRequest
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
			getminsup_data->lambda_max_ + 4,
			MPI_LONG_LONG_INT, src, Tag::DTD_ACCUM_REPLY, &recv_status);
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
//	mpi_data.lifeline_thieves_->Clear();
//	for (int pi = 0; pi < mpi_data.nTotalProc_; pi++) {
//		mpi_data.lifelines_activated_[pi] = false;
//	}

	assert(mpi_data.echo_waiting_ == false);
	assert(mpi_data.waiting_ == false);
//	assert(mpi_data.thieves_->Size() == 0);
//	assert(mpi_data.lifeline_thieves_->Size() == 0);
//	for (int pi = 0; pi < mpi_data.nTotalProc_; pi++) {
//		assert(mpi_data.lifelines_activated_[pi] == false);
//	}
}

} /* namespace lamp_search */
