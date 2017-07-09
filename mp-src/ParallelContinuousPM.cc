/*
 * ParallelSearch.cpp
 *
 *  Created on: Apr 17, 2017
 *      Author: yuu
 */

#include "ParallelContinuousPM.h"

#include "gflags/gflags.h"

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
				alpha_), freq_received(true) {
//	g_ = new LampGraph<uint64>(*d_); // No overhead to generate LampGraph.
}

ParallelContinuousPM::~ParallelContinuousPM() {
	// TODO: lots of things to delete
//	if (g_)
//		delete g_;
}

void ParallelContinuousPM::GetMinimalSupport() {
//	this->getminsup_data = getminsup_data;
//	CheckInit();
	phase_ = 1; // TODO: remove dependency on this
	Search();

	// return lambda?
}

void ParallelContinuousPM::GetDiscretizedMinimalSupport(
		double freqRatio) {
	printf("GetDiscretizedMinimalSupport\n");
//	CheckInit();
	phase_ = 4; // TODO: remove dependency on this

	// TODO: Edit getminsup_data here.

	// TODO: How do we select the discretization???
	thresholds = InitializeThresholdTable(freqRatio, 128, alpha_);
//	printf("thresholds.size() = %d\n", thresholds.size());
	long long int* dtd_accum_array_base_ =
			new long long int[thresholds.size() + 4];
	long long int* dtd_accum_recv_base_ =
			new long long int[thresholds.size() + 4];
	long long int* accum_array_ = &(dtd_accum_array_base_[3]); // TODO: ???
	long long int* accum_recv_ = &(dtd_accum_recv_base_[3]);
	long long int* cs_thr_; // TODO: moc
	for (int i = 0; i < thresholds.size(); ++i) {
		accum_array_[i] = 0;
		accum_recv_[i] = 0;
	}
//	count.assign(thresholds.size(), 0);
	int lambda_ = 1;
	int lambda_max_ = thresholds.size();
//	printf("getminsup_data\n");
	this->getminsup_data = new GetMinSupData(lambda_max_, lambda_,
			cs_thr_, dtd_accum_array_base_, accum_array_,
			dtd_accum_recv_base_, accum_recv_);

//	printf("Ready GetDiscretizedMinimalSupport\n");
	Search();
//	printf("Done GetDiscretizedMinimalSupport\n");

	// TODO: thre_pmin_ should be alpha / number of items with frequencies higher
	thre_freq_ = thresholds[getminsup_data->lambda_ - 1].first;
	// TODO: Where should we put a threshold?
	thre_pmin_ = thresholds[getminsup_data->lambda_ - 1].second;

	// return lambda?
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

void ParallelContinuousPM::GetSignificantPatterns(
		GetContSignificantData* getsignificant_data) {
	this->getsignificant_data = getsignificant_data;
	DBG(D(1) << "MainLoop" << std::endl
	;);
	printf("extract\n");
	ExtractSignificantSet();
	if (mpi_data.mpiRank_ == 0) {
		printf("sendresults\n");
		SendResultRequest();
	}
	printf("probe\n");
	while (!mpi_data.dtd_->terminated_) {
		Probe(treesearch_data);
	}
}

/**
 * Find the Patterns with K lowest p-values.
 *
 * Sequential version pseudocode
 * Expand(n)
 * 1 calculate min-pvalue
 * 2 if min-pvalue > threshold
 * 3    continue
 * 4 else
 * 5    calculate pvalue
 * 6    if pvalue < threshold
 * 7        insert N in Top-K
 * 8    expand(child in n.children)
 *
 * Process of Top-K: similar to LAMP:
 * 1. Find the threshold pvalue for the Top-K
 * 2. List all patterns with pvalue smaller than the threshold.
 *
 */
void ParallelContinuousPM::GetTopKPvalue(int k) {
	phase_ = 5;
	topk = k;
	printf("GetTopKPvalue\n");
	Search();
}

void ParallelContinuousPM::SearchSignificantPatterns(double pvalue) {
	phase_ = 6;
	printf("SearchSignificantPatterns\n");
	Search();
}

//==============================================================================
/**
 * Procedure to run after Probe().
 */
void ParallelContinuousPM::ProcAfterProbe() {
//	printf("ProcAfterProbe\n");

	if (mpi_data.mpiRank_ == 0) {
		// initiate termination detection
		// note: for phase_ 1, accum request and dtd request are unified
		if (!mpi_data.echo_waiting_ && !mpi_data.dtd_->terminated_) {
			if (phase_ == 1 || phase_ == 5) {
				if (treesearch_data->node_stack_->Empty()) {
					if (freq_received) {
						printf("SendMinPValueRequest because "
								"stack is empty and freq received\n");
						SendMinPValueRequest();
					} else {
						printf(
								"SendDTDRequest because "
										"stack is empty and freq not received\n");
						SendDTDRequest();
					}
				}
			} else if (phase_ == 2 || phase_ == 6) {
				if (treesearch_data->node_stack_->Empty()) {
					SendDTDRequest();
				}
			} else if (phase_ == 4) {
				if (treesearch_data->node_stack_->Empty()) {
					SendDTDAccumRequest();
				}
			} else {
				// unknown phase
				assert(0);
			}
		}
	}
}

void ParallelContinuousPM::ProbeExecute(
		TreeSearchData* treesearch_data, MPI_Status* probe_status,
		int probe_src, int probe_tag) {
//	printf("ProbeExecute\n");

	switch (probe_tag) {
	/**
	 * MINIMALSUPPORT
	 */
	case Tag::CONT_REQUEST:
		assert(phase_ == 1 || phase_ == 5);
		RecvMinPValueRequest(probe_src); // TODO: This can be solved by polymorphism.
		break;
	case Tag::CONT_REPLY:
		assert(phase_ == 1 || phase_ == 5);
		RecvMinPValueReply(probe_src, probe_status);
		break;
	case Tag::CONT_LAMBDA:
		assert(phase_ == 1 || phase_ == 5);
		RecvNewSigLevel(probe_src);
		break;
		/**
		 * Linear Space CPM Minimal Support
		 */
	case Tag::DTD_ACCUM_REQUEST:
		assert(phase_ == 4);
		RecvDTDAccumRequest(probe_src);
		break;
	case Tag::DTD_ACCUM_REPLY:
		assert(phase_ == 4);
		RecvDTDAccumReply(probe_src);
		break;
	case Tag::LAMBDA:
		assert(phase_ == 4);
		RecvLambda(probe_src);
		break;
		/**
		 * SIGNIFICANTSET
		 */
	case Tag::RESULT_REQUEST:
		RecvResultRequest(probe_src);
		break;
	case Tag::RESULT_REPLY:
		RecvResultReply(probe_src, *probe_status);
		break;

	default:
		DBG(
				D(1) << "unknown Tag for ParallelPatternMining="
						<< probe_tag << " received in Probe: "
						<< std::endl
				;
		);
		ParallelDFS::ProbeExecute(treesearch_data, probe_status,
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
void ParallelContinuousPM::Check() {
//	printf("Check\n");
	if (mpi_data.mpiRank_ == 0) {
		if (phase_ == 1 || phase_ == 5) {
			if (!mpi_data.echo_waiting_) {
				SendMinPValueRequest();
			}
		} else if (phase_ == 4) {
			if (!mpi_data.echo_waiting_) {
				CheckCSThreshold();
			}
		} else {
		}
	}

	return;
}

bool ParallelContinuousPM::ExpandNode(
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

			// TODO TOPK: This should check if the pvalue to the current top-k pvalue.
			if (!TestAndPushNode(new_item)) {
				continue;
			}

			int* child_node_buf = treesearch_data->node_stack_->Top();

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
//			d_->GetFreq()
//			double pp

			DBG(
					if (phase_ == 2) {
						D(3) << "found cs "
						;
						treesearch_data->node_stack_->Print(D(3),
								child_node_buf)
						;
					});

			double freq = d_->GetFreq(child_sup_buf_);

			// TODO TOPK: Store pvalue instead of frequency.
			ProcessNode(freq, child_node_buf);

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
	sup_buf_ = d_->GetFreqArray(std::move(buf)); // TODO: std::move
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
	if (phase_ == 1) {
		frequencies.push_back(freq);
	} else if (phase_ == 4) {
		int disFreq = GetDiscretizedFrequency(freq);
		IncCsAccum(disFreq);
	} else if (phase_ == 5) {
		// For Top-K we store pvalue instead of frequencies.
		int n = treesearch_data->node_stack_->GetItemNum(ppc_ext_buf);
		int* array = treesearch_data->node_stack_->GetItemArray(
				ppc_ext_buf);
		std::vector<int> buf(array, array + n);
		double actual_pvalue = d_->CalculatePValue(buf);
		frequencies.push_back(actual_pvalue);
	} else if (phase_ == 6) {
		assert(false && "phase_ 6 not implemented yet");
	} else {
//		assert(freq);
//		assert(pmin < pmin_thre_);
		if (freq > thre_freq_) { // permits == case?
			// TODO: Item is considered testable if freq > threshold.
//			double pmin = d_->CalculatePMin(freq);
			gettestable_data->freq_stack_->PushPre();
			int * item = gettestable_data->freq_stack_->Top();
			gettestable_data->freq_stack_->CopyItem(ppc_ext_buf,
					item);
			gettestable_data->freq_stack_->PushPostNoSort();
			gettestable_data->freq_map_->insert(
					std::pair<double, int*>(freq, item));
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

				Probe(treesearch_data);
				Distribute(treesearch_data);
				Reject();

				lap_time = timer_->Elapsed();
			}
		}
	} else {            // not using milli second
		if (accum_period_counter_ >= FLAGS_probe_period_) {
			accum_period_counter_ = 0;
			log_->d_.process_node_time_ += timer_->Elapsed()
					- lap_time;

			Probe(treesearch_data);
			Distribute(treesearch_data);
			Reject();

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

bool ParallelContinuousPM::HasJobToDo() {
	return !(treesearch_data->node_stack_->Empty())
			|| (mpi_data.thieves_->Size() > 0)
			|| treesearch_data->stealer_->StealStarted()
			|| mpi_data.processing_node_ || !frequencies.empty();
}

//bool ParallelContinuousPM::AccumCountReady(MPI_Data& mpi_data) const {
//	for (int i = 0; i < k_echo_tree_branch; i++) {
//		if (mpi_data.bcast_targets_[i] >= 0
//				&& !(mpi_data.accum_flag_[i])) {
//			printf("AccumCount is Not Ready\n");
//			return false;
//		}
//	}
//	printf("AccumCount is Ready\n");
//	return true;
//// check only valid bcast_targets_
//// always true if leaf
//}

/**
 * Methods for Maintaining threshold value
 */
void ParallelContinuousPM::SendMinPValueRequest() {
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
				Tag::CONT_REQUEST);

		DBG(
				D(3) << "SendDTDAccumRequest: dst="
						<< mpi_data.bcast_targets_[i] << "\ttimezone="
						<< mpi_data.dtd_->time_zone_ << std::endl
				;);
	}
}

void ParallelContinuousPM::RecvMinPValueRequest(int src) {
	printf("RecvMinPValueRequest\n");

//	DBG(
//			D(3) << "RecvDTDAccumRequest: src=" << src
//					<< "\ttimezone=" << mpi_data.dtd_->time_zone_
//					<< std::endl
//			;);
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::CONT_REQUEST,
			&recv_status);

	if (IsLeafInTopology())
		SendMinPValueReply();
	else
		SendMinPValueRequest();
}

void ParallelContinuousPM::SendMinPValueReply() {
	printf("SendMinPValueReply\n");
//	std::sort(freq_stack_.begin(), freq_stack_.end());
// TODO: Apply pruning items with pmin higher than alpha/k.

	int size = frequencies.size();
	printf("SendMinPValueReply: send %d items\n", size);
	CallBsend(frequencies.data(), size, MPI_DOUBLE,
			mpi_data.bcast_source_, Tag::CONT_REPLY);

// TODO: Integrate with terminate detection or not?
//	mpi_data.dtd_->time_warp_ = false;
//	mpi_data.dtd_->not_empty_ = false;
//	mpi_data.dtd_->IncTimeZone();
//
	mpi_data.echo_waiting_ = false;
	mpi_data.dtd_->ClearAccumFlags();
//	mpi_data.dtd_->ClearReduceVars();

	frequencies.clear();
}

// getminsup_data
void ParallelContinuousPM::RecvMinPValueReply(int src,
		MPI_Status* probe_status) {

//	MPI_Status recv_status;
// TODO: get status from caller
	int num = 0;
	MPI_Get_count(probe_status, MPI_DOUBLE, &num);
	printf("RecvMinPValueReply: recved %d items\n", num);

	int prev_size = frequencies.size();
	frequencies.resize(prev_size + num);

//	std::vector<double> buffer(num);
//	pmins_stack_
// TODO: single data is not enough!
	CallRecv(&(frequencies[prev_size]), num,
	MPI_DOUBLE, src, Tag::CONT_REPLY, probe_status);
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
		if (DTDReplyReady()) {
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
		if (DTDReplyReady()) {
			// TODO: Merge and sort pmins_stack here for efficiency.
			SendMinPValueReply();
		}
	}
}

/**
 * In Top-K (phase_ == 5), frequencies actually contains pvalue.
 * To find the threshold pvalue, we find Top-K pvalue.
 * We use that to prune patterns:
 * TODO: curent implementation of top-K relies solely on pvalue pruning.
 * 1. If    pvalue(n) > thre_pvalue then prune the node.
 * 2. If minpvalue(n) > thre_pvalue then prune the node and all the children.
 */
void ParallelContinuousPM::CalculateThreshold() {
	assert(phase_ == 1 || phase_ == 5);

	// TODO: Threshold can be easily pruned by merge sort.
	printf("CalculateThreshold\n");
	double prev_thre_freq = thre_freq_;
//
// pmin_stack_ has a stack of pmins from all the processes.
// Our goal here is to find an improved threshold which satisfies pmin(k) = alpha/k.
// TODO: Binary search may be more efficient.
// Sort is already inefficient.
	int numItemsets = frequencies.size();
	if (numItemsets == 0) {
		printf("NO NEW FREQ ACQUIRED! READY FOR TERMINATION\n");
		// TODO: Put terminate detection.
		freq_received = false;
		return;
	}

	// TODO: This should be more efficient with merge-sort.
	topKFrequencies.insert(topKFrequencies.end(), frequencies.begin(),
			frequencies.end());
	frequencies.clear();

	std::sort(topKFrequencies.rbegin(), topKFrequencies.rend());

	if (phase_ == 1) {
		for (int i = 0; i < topKFrequencies.size(); ++i) {
			double pmin = d_->CalculatePLowerBound(
					topKFrequencies[i]);
			if (pmin * (i + 1) >= alpha_) {
				if (thre_freq_ < topKFrequencies[i]) {
					thre_freq_ = topKFrequencies[i];
				}
				topKFrequencies.erase(topKFrequencies.begin() + i,
						topKFrequencies.end());
				break;
			}
		}
	} else if (phase_ == 5) {
		if (topKFrequencies.size() > topk) {
			printf("Shrink topKFrequencies: %d topK to %d items\n",
					topKFrequencies.size(), topk);
			topKFrequencies.erase(topKFrequencies.begin() + topk,
					topKFrequencies.end());
		}
	} else {
		assert(false && "unknown phase in CalculateThreshold");
	}

//	std::vector<double> pmins;
//	for (int i = 0; i < frequencies.size(); ++i) {
//		double pmin = d_->CalculatePMin(frequencies[i]);
//		if (pmin < thre_pmin_) {
//			freq_pmin.push_back(
//					std::pair<double, double>(frequencies[i], pmin));
//			// TODO: we can check threshold everytime inserting an item.
//		}
//	}
//	frequencies.clear();
//
//	printf("%d itemsets out of %d itemsets is pmin(X) < thre_pmin\n",
//			freq_pmin.size(), numItemsets);
//	auto cmp =
//			[](std::pair<double,double> const & a, std::pair<double,double> const & b)
//			{
//				return a.first < b.first;
//			};
//	std::sort(freq_pmin.begin(), freq_pmin.end(), cmp);
////	std::sort(pmins.begin(), pmins.end());
//	while (freq_pmin.size() * thre_pmin_ >= alpha_) {
////		printf(
////				"thre_freq_ %.4f is too low!: %d (|T|) * %.4f (pmin(X)) >= %.2f\n",
////				thre_freq_, freq_pmin.size(), thre_pmin_, alpha_);
//
//		// TODO: Can we remove multiple items at the same time?
//		double min_freq = freq_pmin.begin()->first;
//
//		if (d_->NumPosRatio() < min_freq) {
//			thre_freq_ = d_->NumPosRatio();
////			thre_pmin_ = d_->PMinWithNumPosRatio();
//		} else {
//			thre_freq_ = min_freq;
////			thre_pmin_ = freq_pmin.begin()->second;
//		}
//		freq_pmin.erase(freq_pmin.begin());
//	}

	if (phase_ == 1) {
		assert(prev_thre_freq <= thre_freq_);
		if (prev_thre_freq < thre_freq_) {
			SendNewSigLevel(thre_freq_);
		}
		freq_received = true;
		thre_pmin_ = alpha_ / (double) topKFrequencies.size();
	} else if (phase_ == 5) {
		thre_pmin_ = *topKFrequencies.rbegin();
		SendNewSigLevel(thre_pmin_);
		freq_received = true;
	} else {
		assert(false);
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
				mpi_data.bcast_targets_[i], Tag::CONT_LAMBDA);

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

	CallRecv(&msg, 1, MPI_DOUBLE, src, Tag::CONT_LAMBDA,
			&recv_status);

	if (phase_ == 1) {
		assert(0.0 <= msg);
// TODO: Each process may have its own threshold updated.
		assert(msg > thre_freq_);
		if (msg > thre_freq_) {
			SendNewSigLevel(msg);
			thre_freq_ = msg;
			thre_pmin_ = d_->CalculatePLowerBound(thre_freq_);
		}

// TODO: remove freqs in freq_stacks_ smaller than thre_freq_.
		std::sort(frequencies.begin(), frequencies.end());
		std::vector<double>::iterator it = std::lower_bound(
				frequencies.begin(), frequencies.end(), thre_freq_);
		if (it != frequencies.end()) {
			int prev_size = frequencies.size();
			frequencies.erase(it, frequencies.end());
//		printf(
//				"With new thre_freq_=%.8f, freq_stack_ shrank from %d to %d.\n",
//				thre_freq_, prev_size, frequencies.size());
		} else {
//		printf("All freqs in stack is smaller than thre_freq_\n");
		}
	} else {
		thre_pmin_ = msg;
		std::sort(frequencies.begin(), frequencies.end());
		std::vector<double>::iterator it = std::lower_bound(
				frequencies.begin(), frequencies.end(), thre_pmin_);
		if (it != frequencies.end()) {
			int prev_size = frequencies.size();
			frequencies.erase(it, frequencies.end());
		}
	}
}
/**
 * LINEAR SPACE CONTINUOUS PATTERN MINING
 *
 */
void ParallelContinuousPM::SendDTDAccumRequest() {
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

void ParallelContinuousPM::RecvDTDAccumRequest(int src) {
	DBG(
			D(3) << "RecvDTDAccumRequest: src=" << src
					<< "\ttimezone=" << mpi_data.dtd_->time_zone_
					<< std::endl
			;);
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::DTD_ACCUM_REQUEST,
			&recv_status);

	if (IsLeafInTopology())
		SendDTDAccumReply();
	else
		SendDTDAccumRequest();
}

void ParallelContinuousPM::SendDTDAccumReply() {
	getminsup_data->dtd_accum_array_base_[0] = mpi_data.dtd_->count_
			+ mpi_data.dtd_->reduce_count_;
	bool tw_flag = mpi_data.dtd_->time_warp_
			|| mpi_data.dtd_->reduce_time_warp_;
	getminsup_data->dtd_accum_array_base_[1] = (tw_flag ? 1 : 0);
// for Steal
	mpi_data.dtd_->not_empty_ =
			!(treesearch_data->node_stack_->Empty())
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
			D(3) << "SendDTDAccumReply: dst = "
					<< mpi_data.bcast_source_ << "\tcount="
					<< getminsup_data->dtd_accum_array_base_[0]
					<< "\ttw=" << tw_flag << "\tem=" << em_flag
					<< std::endl
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
void ParallelContinuousPM::RecvDTDAccumReply(int src) {
	MPI_Status recv_status;

	CallRecv(getminsup_data->dtd_accum_recv_base_,
			getminsup_data->lambda_max_ + 4,
			MPI_LONG_LONG_INT, src, Tag::DTD_ACCUM_REPLY,
			&recv_status);
	assert(src == recv_status.MPI_SOURCE);

	int count = (int) (getminsup_data->dtd_accum_recv_base_[0]);
	bool time_warp = (getminsup_data->dtd_accum_recv_base_[1] != 0);
	bool not_empty = (getminsup_data->dtd_accum_recv_base_[2] != 0);

	mpi_data.dtd_->Reduce(count, time_warp, not_empty);

	DBG(
			D(3) << "RecvDTDAccumReply: src=" << src << "\tcount="
					<< count << "\ttw=" << time_warp << "\tem="
					<< not_empty << "\treduced_count="
					<< mpi_data.dtd_->reduce_count_ << "\treduced_tw="
					<< mpi_data.dtd_->reduce_time_warp_
					<< "\treduced_em="
					<< mpi_data.dtd_->reduce_not_empty_ << std::endl
			;);

	for (int l = getminsup_data->lambda_ - 1;
			l <= getminsup_data->lambda_max_; l++)
		getminsup_data->accum_array_[l] +=
				getminsup_data->accum_recv_[l];

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
			SendLambda(new_lambda);
			getminsup_data->lambda_ = new_lambda;
			thre_freq_ =
					thresholds[getminsup_data->lambda_ - 1].first;
//			thre_pmin_ = thresholds[getminsup_data->lambda_].second;
		}
// if SendLambda is called, dtd_.count_ is incremented and DTDCheck will always fail
		if (DTDReplyReady()) {
			DTDCheck();
			log_->d_.dtd_accum_phase_num_++;
		}
	} else {  // not root
		if (DTDReplyReady())
			SendDTDAccumReply();
	}
}

int ParallelContinuousPM::GetDiscretizedFrequency(double freq) const {
	int i = 0;
	while (freq > thresholds[i].first && i < thresholds.size()) {
		++i;
	}
//	--i;
	if (i == thresholds.size()) {
		printf("Fr_d(%.2f) = %d\n", freq, i);
		--i;
	} else {
		assert(thresholds[i].first > freq);
	}
	return i;
}

void ParallelContinuousPM::CheckCSThreshold() {
	//	printf("CheckCSThreshold\n");
//	assert(mpi_data.mpiRank_ == 0);
	if (ExceedCsThr()) {
		printf("Exceeded!\n");
		int new_lambda = NextLambdaThr();
		SendLambda(new_lambda);
		printf("Lambda updated to %d\n", new_lambda);
		getminsup_data->lambda_ = new_lambda;
		thre_freq_ = thresholds[getminsup_data->lambda_ - 1].first;
	} else {
//		printf("Did not exceeded\n");

	}
}

bool ParallelContinuousPM::ExceedCsThr() const {
//	printf("ExceedCsThr: fr_d = %d, k f(T_k) = %d * %.6f = %.6f\n",
//			getminsup_data->lambda_,
//			getminsup_data->accum_array_[getminsup_data->lambda_],
//			thresholds[getminsup_data->lambda_].second,
//			getminsup_data->accum_array_[getminsup_data->lambda_]
//					* thresholds[getminsup_data->lambda_].second);

	return (getminsup_data->accum_array_[getminsup_data->lambda_]
			* thresholds[getminsup_data->lambda_].second >= alpha_);
}

int ParallelContinuousPM::NextLambdaThr() const {
	printf("NextLambdaThr\n");

	for (int si = getminsup_data->lambda_max_;
			si >= getminsup_data->lambda_; si--) {
		if (getminsup_data->accum_array_[si] * thresholds[si].second
				>= alpha_) {
			return si + 1;
		}
	}
	assert(false && "NextLambdaThr not returning correct threshold.");
	return 0;
// it is safe because lambda_ higher than max results in immediate search finish
}

void ParallelContinuousPM::IncCsAccum(int sup_num) {
//	printf("IncCsAccum\n");
	for (int i = 0; i <= sup_num; i++)
		getminsup_data->accum_array_[i]++;
}

int ParallelContinuousPM::NumberOfTestablePatterns() const {
	return getminsup_data->accum_array_[getminsup_data->lambda_ - 1];
}

void ParallelContinuousPM::SendLambda(int lambda) {
	printf("SendLambda\n");
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
		CallBsend(message, 2, MPI_INT, mpi_data.bcast_targets_[i],
				Tag::LAMBDA);
		mpi_data.dtd_->OnSend();

		DBG(
				D(2) << "SendLambda: dst="
						<< mpi_data.bcast_targets_[i] << "\tlambda="
						<< lambda << "\tdtd_count="
						<< mpi_data.dtd_->count_ << std::endl
				;);
	}
}

void ParallelContinuousPM::RecvLambda(int src) {
//	printf("RecvLambda\n");
	MPI_Status recv_status;
	int message[2];

	CallRecv(&message, 2, MPI_INT, src, Tag::LAMBDA, &recv_status);
	mpi_data.dtd_->OnRecv();
	assert(src == recv_status.MPI_SOURCE);
	int timezone = message[0];
	mpi_data.dtd_->UpdateTimeZone(timezone);

	DBG(
			D(2) << "RecvLambda: src=" << src << "\tlambda="
					<< message[1] << "\tdtd_count="
					<< mpi_data.dtd_->count_ << std::endl
			;);

	int new_lambda = message[1];
	if (new_lambda > getminsup_data->lambda_) {
		SendLambda(new_lambda);
		getminsup_data->lambda_ = new_lambda;
		thre_freq_ = thresholds[getminsup_data->lambda_ - 1].first;
// todo: do database reduction
	}
}

std::vector<std::pair<double, double> > ParallelContinuousPM::InitializeThresholdTable(
		double ratio, int size, double alpha) {
	printf("InitializeThresholdTable\n");
	// TODO: the table should be more efficient with inversed.
	std::vector<double> thresholds;
	double max_freq = 1.0;
	// TODO: Current discretization is way too rough.
	//       Need to find a way to edit the granularity.
	for (int i = 0; i < size; ++i) {
		max_freq = max_freq * ratio; // TODO
		double pbound = d_->CalculatePLowerBound(max_freq);
		if (pbound >= alpha) {
			break;
		} else {
			thresholds.push_back(max_freq);
		}
	}

	std::sort(thresholds.begin(), thresholds.end());

	std::vector<std::pair<double, double> > table(thresholds.size());

	for (int i = 0; i < thresholds.size(); ++i) {
		table[i].first = thresholds[i];
		table[i].second = d_->CalculatePLowerBound(table[i].first);
	}
	printf("The domain of discrete Fr(X) = {0..%d}\n",
			thresholds.size());
//	for (int i = 0; i < thresholds.size(); ++i) {
//		printf("freq/minp = %.2f/%.6f\n", table[i].first,
//				table[i].second);
//	}

	return table;
}

/**
 * GETSIGNIFICANT PATTERNS
 *
 */

//==============================================================================
void ParallelContinuousPM::SendResultRequest() {
//	printf("SendResultRequest\n");
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

void ParallelContinuousPM::RecvResultRequest(int src) {
//	printf("RecvResultRequest\n");
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::RESULT_REQUEST,
			&recv_status);
	assert(src == recv_status.MPI_SOURCE);

	DBG(D(2) << "RecvResultRequest: src=" << src << std::endl
	;);

	if (IsLeafInTopology())
		SendResultReply();
	else
		SendResultRequest();
}

void ParallelContinuousPM::SendResultReply() {
//	printf("SendResultReply\n");
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

void ParallelContinuousPM::RecvResultReply(int src,
		MPI_Status probe_status) {
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

	if (AccumCountReady()) {
		if (mpi_data.mpiRank_ != 0) {
			SendResultReply();
		} else { // root
			mpi_data.echo_waiting_ = false;
			mpi_data.dtd_->terminated_ = true;
		}
	}
}

void ParallelContinuousPM::ExtractSignificantSet() {
//	printf("ExtractSignificantSet\n");

//	double thre_bonferroni = d_->CalculatePMin(pmin_thre_);
	std::multimap<double, int *>::iterator it;
//	printf("bonferroni corrected threshold = %.8f (sig_level_)\n",
//			getsignificant_data->final_sig_level_);
//	printf("bonferroni corrected threshold = %.8f (thre_pmin_)\n",
//			thre_pmin_);
//	printf("#Testable Pattern = %d (freq_pmin.size())\n",
//			freq_pmin.size());
//	printf("#Testable Pattern = %.4f (alpha/thre_pmin_)\n",
//			alpha_ / thre_pmin_);
//	printf("#Testable Pattern = %.4f (alpha/final_sig_level_)\n",
//			alpha_ / getsignificant_data->final_sig_level_);
	for (it = getsignificant_data->freq_map_->begin();
			it != getsignificant_data->freq_map_->end(); ++it) {
		// TODO: Here we should implement calculating p-value.
		std::vector<int> itemset =
				getsignificant_data->freq_stack_->getItems(
						(*it).second);
		double actual_pvalue = d_->CalculatePValue(itemset);
//		double minimal_pvalue = d_->CalculatePMin((*it).first);
		// TODO: equal??
		if (actual_pvalue <= thre_pmin_) {
//			printf("Significant Itemset = ");
//			for (int i = 0; i < itemset.size(); ++i) {
//				printf("%d ", itemset[i]);
//			}
//			printf(": Pvalue = %.8f, Pmin = %.8f, Freq = %.8f \n",
//					actual_pvalue, minimal_pvalue, (*it).first);
			getsignificant_data->significant_stack_->PushPre();
			int * item =
					getsignificant_data->significant_stack_->Top();

			getsignificant_data->significant_stack_->CopyItem(
					(*it).second, item);
			getsignificant_data->significant_stack_->PushPostNoSort();
		} else {
//			printf("Insignificant but Testable Itemset = ");
//			for (int i = 0; i < itemset.size(); ++i) {
//				printf("%d ", itemset[i]);
//			}
//			printf(": Pvalue = %.8f, Pmin = %.8f, Freq = %.8f \n",
//					actual_pvalue, minimal_pvalue, (*it).first);
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
