/*
 * ParallelSearch.cpp
 *
 *  Created on: Apr 17, 2017
 *      Author: yuu
 */

#include "ParallelContinuousPM.h"

#include "google/gflags.h"

#include "mpi_tag.h"
#include "Log.h"
#include "FixedSizeStack.h"
#include "DTD.h"
#include "StealState.h"
//#include "SignificantSetResults.h"

#include "../src/timer.h"

#include <vector>
#include <algorithm>

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

ParallelContinuousPM::ParallelContinuousPM(
		ContinuousPatternMiningData* bpm_data, MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, double alpha, Log* log,
		Timer* timer, std::ostream& ofs) :
		ParallelDFS(mpi_data, treesearch_data, log, timer, ofs), d_(
				bpm_data->d_), alpha_(alpha), expand_num_(0), closed_set_num_(
				0), phase_(0), gettestable_data(NULL), getsignificant_data(
		NULL), sup_buf_(d_->NumTransactions(), 1.0), thre_freq_(0.0), thre_pmin_(
				alpha_) {
//	g_ = new LampGraph<uint64>(*d_); // No overhead to generate LampGraph.
}

ParallelContinuousPM::~ParallelContinuousPM() {
	// TODO: lots of things to delete
//	if (g_)
//		delete g_;
}

void ParallelContinuousPM::GetTestablePatterns(
		GetTestableData* gettestable_data) {
	this->gettestable_data = gettestable_data; // TODO: not sure if this is a good idea.
//	this->getminsup_data->lambda_ = gettestable_data->freqThreshold_; // TODO: not the best way.
//	CheckInitTestable();
	phase_ = 2;  // TODO: remove dependency on this
	DBG(D(1) << "MainLoop" << std::endl
	;);
	printf("GetTestablePatterns\n");
	Search();

}

void ParallelContinuousPM::GetSignificantPatterns(MPI_Data& mpi_data,
		GetContSignificantData* getsignificant_data) {
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
void ParallelContinuousPM::ProcAfterProbe() {
//	printf("ProcAfterProbe\n");

	if (mpi_data.mpiRank_ == 0) {
		// initiate termination detection

		// TODO: isn't too often to consider calling DTDRequest after every Probe?
		// note: for phase_ 1, accum request and dtd request are unified
		if (!mpi_data.echo_waiting_ && !mpi_data.dtd_->terminated_) {
			if (phase_ == 2) {
				// TODO: This sends too many DTDs?
				// Delay this for testing purpose...?
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

void ParallelContinuousPM::ProbeExecute(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, MPI_Status* probe_status,
		int probe_src, int probe_tag) {
//	printf("ProbeExecute\n");

	switch (probe_tag) {
	/**
	 * MINIMALSUPPORT
	 */
	case Tag::DTD_ACCUM_REQUEST:
//		assert(false);
		RecvMinPValueRequest(mpi_data, probe_src); // TODO: This can be solved by polymorphism.
		break;
	case Tag::DTD_ACCUM_REPLY:
//		assert(false);
		RecvMinPValueReply(mpi_data, probe_src, probe_status);
		break;
	case Tag::LAMBDA:
//		assert(false);
		RecvNewSigLevel(probe_src);
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
				D(1) << "unknown Tag for ParallelPatternMining="
						<< probe_tag << " received in Probe: "
						<< std::endl
				;
		);
		ParallelDFS::ProbeExecute(mpi_data, treesearch_data,
				probe_status, probe_src, probe_tag);
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
void ParallelContinuousPM::Check(MPI_Data& mpi_data) {
//	printf("Check\n");
	if (mpi_data.mpiRank_ == 0 && phase_ == 2) {
//		assert(false);
		if (!mpi_data.echo_waiting_) {
			SendMinPValueRequest(mpi_data);
		} else {
//			printf("MinPValue not sent due to echo_waiting_.\n");
		}
	}

	return;
}

// TODO: What does it do after all?
bool ParallelContinuousPM::ExpandNode(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data) {
//	printf("ExpandNode\n");
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

//		printf("Stack has %d itemsets\n",
//				treesearch_data->node_stack_->NuItemset());
		if (!PopNodeFromStack()) {
			// Discard a node if pmin(node) > pmin_threshould_.
			continue;
		}
//		printf("Stack has %d itemsets\n",
//				treesearch_data->node_stack_->NuItemset());
		int n = treesearch_data->node_stack_->GetItemNum(
				treesearch_data->itemset_buf_);
//		int* init = treesearch_data->node_stack_->GetItemArray(
//				treesearch_data->itemset_buf_);
//		std::vector<int> items(init, init + n);
//		printf("Poped %d\n", n);
//		PrintItemset(treesearch_data->itemset_buf_, sup_buf_);

		assert(
				phase_ != 1
						|| treesearch_data->node_stack_->GetItemNum(
								treesearch_data->itemset_buf_) != 0);
		bool isRoot = (n == 0);

		int accum_period_counter_ = 0;

		// TODO: Implement ContTable to return list of children.

//		int core_i = g_->CoreIndex(*treesearch_data->node_stack_,
//				treesearch_data->itemset_buf_);
		std::vector<int> current_items =
				treesearch_data->node_stack_->getItems(
						treesearch_data->itemset_buf_);
		std::vector<int> children = GetChildren(current_items);

		for (std::vector<int>::iterator it = children.begin();
				it != children.end(); ++it) {
			int new_item = *it;
			assert(0 <= new_item && new_item < d_->NumItems());

			// TODO: no need for duplicate detection
			if (treesearch_data->node_stack_->Exist(
					treesearch_data->itemset_buf_, new_item)) {
				printf("Node Pruned because of duplicated items");
				return false;
			}
			// skip existing item
			// todo: improve speed here
			if (!isRoot) {
				CheckProbe(accum_period_counter_, lap_time);
			}

			if (!TestAndPushNode(new_item)) {
				continue;
			}

			int* ppc_ext_buf = treesearch_data->node_stack_->Top();

//			std::vector<int> child =
//					treesearch_data->node_stack_->getItems(
//							ppc_ext_buf);
//			printf("Pushed %d items:",
//					treesearch_data->node_stack_->GetItemNum(
//							ppc_ext_buf));
//			PrintItemset(ppc_ext_buf, child_sup_buf_);

//			int sup_num = treesearch_data->node_stack_->GetSup(
//					ppc_ext_buf);
//			std::vector<int> items = treesearch_data->node_stack_->getItems(ppc_ext_buf);
			double freq = d_->GetFreq(child_sup_buf_);
//			d_->GetFreq()
//			double pp

			DBG(if (phase_ == 2) {
				D(3) << "found cs "
				;
				treesearch_data->node_stack_->Print(D(3), ppc_ext_buf)
				;
			});

			ProcessNode(freq, ppc_ext_buf);

//			assert(freq >= pmin_thre_);

			// try skipping if supnum_ == sup_threshold,
			// because if sup_num of a node equals to sup_threshold, children will have smaller sup_num
			// therefore no need to check it's children
			// note: skipping node_stack_ full check. allocate enough memory!
//			if (freq == pmin_thre_) {
//				printf("All Children Pruned");
//				treesearch_data->node_stack_->Pop();
//			}
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
std::vector<int> ParallelContinuousPM::GetChildren(
		std::vector<int> items) {
//	printf("GetChildren\n");
	if (items.empty()) {
		// Edge partitioning for root node
		std::vector<int> children = d_->GetChildren(items);
		std::vector<int> responsible;
		for (int i = mpi_data.mpiRank_; i < children.size(); i +=
				mpi_data.nTotalProc_) {
			responsible.push_back(children[i]);
		}
		return responsible;
	} else {
		return d_->GetChildren(items);
	}
}

/**
 * itemset_buf_ := pointer to the index of itemset
 * sup_buf_     := support* of the itemset (takes AND for each item)
 * support*     := transactions which include all items in the itemset
 */
bool ParallelContinuousPM::PopNodeFromStack() {
//	printf("PopNodeFromStack\n");
	/**
	 * Pop item index and put into itemset_buf.
	 * This part is not dependent on the feature.
	 */
	treesearch_data->node_stack_->CopyItem(
			treesearch_data->node_stack_->Top(),
			treesearch_data->itemset_buf_);
//	printf("copied top itemset into itemset_buf_\n");

	treesearch_data->node_stack_->Pop();
//	printf("poped top itemset\n");

// dbg
	DBG(D(3) << "expanded "
	;);
	DBG(
			treesearch_data->node_stack_->Print(D(3),
					treesearch_data->itemset_buf_)
			;);

	/**
	 * Get the support (frequency) of the item into sup_buf_
	 */
// TODO: This part is dependent on feature.
	// Why is it uint64???
	// TODO:
	int n = treesearch_data->node_stack_->GetItemNum(
			treesearch_data->itemset_buf_);
//	printf("%d items in the itemset\n", n);
	int* array = treesearch_data->node_stack_->GetItemArray(
			treesearch_data->itemset_buf_);
	std::vector<int> buf(array, array + n);
	sup_buf_ = d_->GetFreqArray(buf); // TODO: std::move
//	TODO: Put freq into sup_buf

	double freq = d_->GetFreq(sup_buf_);
//	double pmin = d_->CalculatePMin(freq);
	if (freq < thre_freq_) {
		printf(
				"Discarded node as freq = %.4f < %.4f (= thre_freq_)\n",
				freq, thre_freq_);
		return false;
	}
	return true;
}

/**
 * Test two pruning criteria before pushing
 * 1. The support of the itemset should be larger or equal to lambda (minimal support).
 * 2. New node should be a PPC extension of the parent itemset.
 */
bool ParallelContinuousPM::TestAndPushNode(int new_item) {
//	printf("TestAndPushNode\n");
//	bsh_->Copy(sup_buf_, child_sup_buf_);
//	memcpy(sup_buf_, child_sup_buf_, );
// TODO: COPY sup_buf to child_sup_buf
	child_sup_buf_ = d_->GetChildrenFreq(sup_buf_, new_item);
//	printf("GetChildrenFreq done\n");
	Feature child_freq = d_->GetFreq(child_sup_buf_);
//	printf("GetFreq done\n");

//	int sup_num = bsh_->AndCountUpdate(d_->NthData(new_item), child_sup_buf_);
// If the support is smaller than the required minimal support for
// significant pattern (=lambda), then prune it.

// TODO: Can we inverse calculate the frequency from a pmin???
//	double pmin = d_->CalculatePMin(child_freq);
	if (child_freq < thre_freq_) {
//		printf("Node Pruned as pmin %.8f < pmin_thre_ %.8f\n",
//				pmin, pmin_thre_);
		return false;
	}

// TODO: Here it is inserting a new item into the node_stack_.
	treesearch_data->node_stack_->PushPre();
	int* ppc_ext_buf = treesearch_data->node_stack_->Top();
//	printf("treesearch_data->node_stack_->Top() done\n");
// TODO: Not sure what PPCExtension does
	bool res = d_->PPCExtension(treesearch_data->node_stack_,
			treesearch_data->itemset_buf_, new_item, ppc_ext_buf);

// TODO: Those things should be done in other class.
//	treesearch_data->node_stack_->SetSup(ppc_ext_buf, child_freq); // TODO: we can use the stack to hold double?
	treesearch_data->node_stack_->PushPostNoSort();
//	printf("treesearch_data->node_stack_->PushPostNoSort() done\n");

	if (!res) {        // todo: remove this redundancy
		treesearch_data->node_stack_->Pop();
//		printf("treesearch_data->node_stack_->Pop() done\n");

		printf("Node Pruned as node is not PPCExtension\n");
		return false;
	}

	treesearch_data->node_stack_->SortTop();
//	printf("treesearch_data->node_stack_->SortTop() done\n");

	return true;
}

// TODO: what we need here?
void ParallelContinuousPM::ProcessNode(double freq,
		int* ppc_ext_buf) {
//	printf("ProcessNode\n");
	closed_set_num_++;

// TODO: For continuous pattern mining calculating p value should be put later?
	if (true) { // XXX: FLAGS_third_phase_
		// CalculatePValue;
//		printf("Frequency = %.2f -> Pmin = %.2f\n",
//				d_->GetFreq(child_sup_buf_), pmin);
//		int pos_sup_num = bsh_->AndCount(d_->PosNeg(), child_sup_buf_);
//		double pval = d_->PVal(sup_num, pos_sup_num);

		double pmin = d_->CalculatePMin(freq);
		// TODO: Store into freq_stack_
//		assert(freq);
//		assert(pmin < pmin_thre_);
		if (pmin < thre_pmin_) { // permits == case?
			gettestable_data->freq_stack_->PushPre();
			int * item = gettestable_data->freq_stack_->Top();
			gettestable_data->freq_stack_->CopyItem(ppc_ext_buf,
					item);
			gettestable_data->freq_stack_->PushPostNoSort();

			gettestable_data->freq_map_->insert(
					std::pair<double, int*>(freq, item));
//			double pmin = d_->CalculatePMin(freq);
			freq_stack_.push_back(freq);
		}

	}

}

void ParallelContinuousPM::CheckProbe(int& accum_period_counter_,
		long long int lap_time) {
//	printf("CheckProbe\n");
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
			log_->d_.process_node_time_ += timer_->Elapsed()
					- lap_time;

			Probe(mpi_data, treesearch_data);
			Distribute(mpi_data, treesearch_data);
			Reject(mpi_data);

			lap_time = timer_->Elapsed();
		}
	}
// note: do this before PushPre is called [2015-10-05 21:56]
}

bool ParallelContinuousPM::CheckProcessNodeEnd(int n, bool n_is_ms,
		int processed, long long int start_time) {
//	printf("CheckProcessNodeEnd\n");
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

bool ParallelContinuousPM::AccumCountReady(MPI_Data& mpi_data) const {
	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data.bcast_targets_[i] >= 0
				&& !(mpi_data.accum_flag_[i])) {
			printf("AccumCount is Not Ready\n");
			return false;
		}
	}
	printf("AccumCount is Ready\n");
	return true;
// check only valid bcast_targets_
// always true if leaf
}

/**
 * Methods for Maintaining threshold value
 */
void ParallelContinuousPM::SendMinPValueRequest(MPI_Data& mpi_data) {
	printf("SendMinPValueRequest\n");
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

void ParallelContinuousPM::RecvMinPValueRequest(MPI_Data& mpi_data,
		int src) {
	printf("RecvMinPValueRequest\n");

//	DBG(
//			D(3) << "RecvDTDAccumRequest: src=" << src
//					<< "\ttimezone=" << mpi_data.dtd_->time_zone_
//					<< std::endl
//			;);
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::DTD_ACCUM_REQUEST,
			&recv_status);

	if (IsLeafInTopology(mpi_data))
		SendMinPValueReply(mpi_data);
	else
		SendMinPValueRequest(mpi_data);
}

void ParallelContinuousPM::SendMinPValueReply(MPI_Data& mpi_data) {
	printf("SendMinPValueReply\n");
//	std::sort(freq_stack_.begin(), freq_stack_.end());
// TODO: Apply pruning items with pmin higher than alpha/k.

	int size = freq_stack_.size();
	printf("SendMinPValueReply: send %d items\n", size);
	CallBsend(freq_stack_.data(), size, MPI_DOUBLE,
			mpi_data.bcast_source_, Tag::DTD_ACCUM_REPLY);

// TODO: Integrate with terminate detection or not?
//	mpi_data.dtd_->time_warp_ = false;
//	mpi_data.dtd_->not_empty_ = false;
//	mpi_data.dtd_->IncTimeZone();
//
	mpi_data.echo_waiting_ = false;
	mpi_data.dtd_->ClearAccumFlags();
//	mpi_data.dtd_->ClearReduceVars();

	freq_stack_.clear();
}

// getminsup_data
void ParallelContinuousPM::RecvMinPValueReply(MPI_Data& mpi_data,
		int src, MPI_Status* probe_status) {

//	MPI_Status recv_status;
// TODO: get status from caller
	int num = 0;
	MPI_Get_count(probe_status, MPI_DOUBLE, &num);
	printf("RecvMinPValueReply: recved %d items\n", num);

	int prev_size = freq_stack_.size();
	freq_stack_.resize(prev_size + num);

//	std::vector<double> buffer(num);
//	pmins_stack_
// TODO: single data is not enough!
	CallRecv(&(freq_stack_[prev_size]), num,
	MPI_DOUBLE, src, Tag::DTD_ACCUM_REPLY, probe_status);
	assert(src == probe_status->MPI_SOURCE);

// TODO?
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
		if (DTDReplyReady(mpi_data)) {
			CalculateThreshold();
			mpi_data.echo_waiting_ = false;
			mpi_data.dtd_->ClearAccumFlags();
		}
		// TODO: Broadcast Threshold!

// if SendLambda is called, dtd_.count_ is incremented and DTDCheck will always fail
//		if (DTDReplyReady(mpi_data)) {
//			DTDCheck(mpi_data);
//			log_->d_.dtd_accum_phase_num_++;
//		}
	} else {  // not root
		if (DTDReplyReady(mpi_data)) {
			// TODO: Merge and sort pmins_stack here for efficiency.
			SendMinPValueReply(mpi_data);
		}
	}
}

void ParallelContinuousPM::CalculateThreshold() {
	printf("CalculateThreshold\n");
// pmin_stack_ has a stack of pmins from all the processes.
// Our goal here is to find an improved threshold which satisfies pmin(k) = alpha/k.
// TODO: Binary search may be more efficient.
// Sort is already inefficient.
	int numItemsets = freq_stack_.size();

	std::vector<std::pair<double, double>> freq_pmin;
//	std::vector<double> pmins;
	for (int i = 0; i < freq_stack_.size(); ++i) {
		double pmin = d_->CalculatePMin(freq_stack_[i]);
		if (pmin < thre_pmin_) {
			freq_pmin.push_back(
					std::pair<double, double>(freq_stack_[i], pmin));
		} else {
			freq_stack_.erase(freq_stack_.begin() + i);
		}
	}
	printf("%d itemsets out of %d itemsets is pmin(X) < thre_pmin\n",
			freq_pmin.size(), numItemsets);
	auto cmp =
			[](std::pair<double,double> const & a, std::pair<double,double> const & b)
			{
				return a.first < b.first;
			};
	std::sort(freq_pmin.begin(), freq_pmin.end(), cmp);
//	std::sort(pmins.begin(), pmins.end());
	double prev_thre_freq = thre_freq_;
	while (freq_pmin.size() * thre_pmin_ >= alpha_) {
//		printf(
//				"thre_freq_ %.4f is too low!: %d (|T|) * %.4f (pmin(X)) >= %.2f\n",
//				thre_freq_, freq_pmin.size(), thre_pmin_, alpha_);
		double min_freq = freq_pmin.begin()->first;
		thre_freq_ = std::min(d_->NumPosRatio(), min_freq);
		thre_pmin_ = d_->CalculatePMin(thre_freq_);
//		printf("thre_freq updated to %.4f, thre_pmin_ = %.4f\n",
//				thre_freq_, thre_pmin_);
		freq_pmin.erase(freq_pmin.begin());
		// TODO: should also remove from freq_stack_.
	}

	assert(prev_thre_freq <= thre_freq_);
	if (prev_thre_freq < thre_freq_) {
		SendNewSigLevel(thre_freq_);
	}
	printf("thre_pmin_ = %.8f for thre_freq = %.5f\n", thre_pmin_,
			thre_freq_);
}

/**
 * GetMinSup Functions
 *
 *
 */
void ParallelContinuousPM::SendNewSigLevel(double sig_level) {
// send lambda to bcast_targets_
//	printf("SendNewSigLevel is not implemented yet.\n");
//	assert(false);
//	return;

	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data.bcast_targets_[i] < 0)
			break;
		assert(
				mpi_data.bcast_targets_[i] < mpi_data.nTotalProc_
						&& "SendLambda");
		CallBsend(&sig_level, 1, MPI_DOUBLE,
				mpi_data.bcast_targets_[i], Tag::LAMBDA);

		DBG(
				D(2) << "SendLambda: dst="
						<< mpi_data.bcast_targets_[i] << "\tlambda="
						<< sig_level << "\tdtd_count="
						<< mpi_data.dtd_->count_ << std::endl
				;);
	}
}

void ParallelContinuousPM::RecvNewSigLevel(int src) {
	printf("RecvNewSigLevel\n");
	MPI_Status recv_status;
	double msg = -1.0;

	CallRecv(&msg, 1, MPI_DOUBLE, src, Tag::LAMBDA, &recv_status);
	assert(0.0 <= msg);
// TODO: Each process may have its own threshold updated.
	assert(msg > thre_freq_);
	if (msg > thre_freq_) {
		SendNewSigLevel(msg);
		thre_freq_ = msg;
		thre_pmin_ = d_->CalculatePMin(thre_freq_);
// todo: do database reduction
	}
// TODO: remove freqs in freq_stacks_ smaller than thre_freq_.
	std::sort(freq_stack_.begin(), freq_stack_.end());
	std::vector<double>::iterator it = std::lower_bound(
			freq_stack_.begin(), freq_stack_.end(), thre_freq_);
	if (it != freq_stack_.end()) {
		int prev_size = freq_stack_.size();
		freq_stack_.erase(it, freq_stack_.end());
		printf(
				"With new thre_freq_=%.8f, freq_stack_ shrank from %d to %d.\n",
				thre_freq_, prev_size, freq_stack_.size());
	} else {
		printf("All freqs in stack is smaller than thre_freq_\n");
	}
}

/**
 * GETSIGNIFICANT PATTERNS
 *
 */

//==============================================================================
void ParallelContinuousPM::SendResultRequest(MPI_Data& mpi_data) {
	printf("SendResultRequest\n");
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
				D(2) << "SendResultRequest: dst="
						<< mpi_data.bcast_targets_[i] << std::endl
				;);
	}
}

void ParallelContinuousPM::RecvResultRequest(MPI_Data& mpi_data,
		int src) {
	printf("RecvResultRequest\n");
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::RESULT_REQUEST,
			&recv_status);
	assert(src == recv_status.MPI_SOURCE);

	DBG(D(2) << "RecvResultRequest: src=" << src << std::endl
	;);

	if (IsLeafInTopology(mpi_data))
		SendResultReply(mpi_data);
	else
		SendResultRequest(mpi_data);
}

void ParallelContinuousPM::SendResultReply(MPI_Data& mpi_data) {
	printf("SendResultReply\n");
	int * message = getsignificant_data->significant_stack_->Stack();
	int size =
			getsignificant_data->significant_stack_->UsedCapacity();
	assert(
			mpi_data.bcast_source_ < mpi_data.nTotalProc_
					&& "SendResultReply");
	CallBsend(message, size, MPI_INT, mpi_data.bcast_source_,
			Tag::RESULT_REPLY);

	DBG(
			D(2) << "SendResultReply: dst=" << mpi_data.bcast_source_
					<< std::endl
			;);
	DBG(getsignificant_data->significant_stack_->PrintAll(D(3, false))
	;);

	mpi_data.echo_waiting_ = false;
	mpi_data.dtd_->terminated_ = true;
}

void ParallelContinuousPM::RecvResultReply(MPI_Data& mpi_data,
		int src, MPI_Status probe_status) {
	printf("RecvResultReply\n");
	int count;
	int error = MPI_Get_count(&probe_status, MPI_INT, &count);
	if (error != MPI_SUCCESS) {
		DBG(
				D(1) << "error in MPI_Get_count in RecvResultReply: "
						<< error << std::endl
				;);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	MPI_Status recv_status;
	CallRecv(treesearch_data->give_stack_->Stack(), count, MPI_INT,
			src, Tag::RESULT_REPLY, &recv_status);
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

void ParallelContinuousPM::ExtractSignificantSet() {
	printf("ExtractSignificantSet\n");
//	double thre_bonferroni = d_->CalculatePMin(pmin_thre_);
	std::multimap<double, int *>::iterator it;
	printf("bonferroni corrected threshold = %.8f\n", thre_pmin_);
	for (it = getsignificant_data->freq_map_->begin();
			it != getsignificant_data->freq_map_->end(); ++it) {
		// TODO: Here we should implement calculating p-value.
		std::vector<int> itemset =
				treesearch_data->node_stack_->getItems((*it).second);
		double actual_pvalue = d_->CalculatePValue(itemset);

		// TODO: equal??
		if (actual_pvalue <= thre_pmin_) {
//			printf("Significant Itemset = ");
//			for (int i = 0; i < itemset.size(); ++i) {
//				printf("%d ", itemset[i]);
//			}
//			printf("\n Pvalue = %.2f, Pmin = %.2f \n", actual_pvalue,
//					(*it).first);
			getsignificant_data->significant_stack_->PushPre();
			int * item =
					getsignificant_data->significant_stack_->Top();

			getsignificant_data->significant_stack_->CopyItem(
					(*it).second, item);
			getsignificant_data->significant_stack_->PushPostNoSort();
		} else {
			//
		}
	}
}

void ParallelContinuousPM::PrintItemset(int* itembuf,
		std::vector<Feature> freqs) {
	std::vector<int> itemset = treesearch_data->node_stack_->getItems(
			itembuf);
	printf("node = ");
	for (int i = 0; i < itemset.size(); ++i) {
		printf("%d ", itemset[i]);
	}
//	printf("\nFreq = ");
//	for (int i = 0; i < freqs.size(); ++i) {
//		printf("%.2f ", freqs[i]);
//	}
	printf("\n");
}

} /* namespace lamp_search */
