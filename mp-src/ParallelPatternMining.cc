/*
 * ParallelSearch.cpp
 *
 *  Created on: Apr 17, 2017
 *      Author: yuu
 */

#include "ParallelPatternMining.h"

#include "gflags/gflags.h"

#include "mpi_tag.h"
#include "Log.h"
#include "FixedSizeStack.h"
#include "DTD.h"
#include "StealState.h"
//#include "SignificantSetResults.h"

#include "../src/timer.h"
#include "DFSParallelPatternMiningStack.h"

#ifdef __CDT_PARSER__
#undef DBG
#define DBG(a)  a
#endif

#ifdef __CDT_PARSER__
#undef LOG
#define LOG(a)  ;
#endif

namespace lamp_search {

ParallelPatternMining::ParallelPatternMining(
        BinaryPatternMiningData* bpm_data, MPI_Data& mpi_data,
        TreeSearchData* treesearch_data, Log* log, Timer* timer,
		std::ostream& ofs) :
        ParallelDFS(new DFSParallelPatternMiningStack(this, treesearch_data, bpm_data->d_, mpi_data, ofs, std::bind(&ParallelPatternMining::CallbackForProbe, this)),
					mpi_data, treesearch_data->give_stack_, log, timer, ofs),
        mpi_data_(mpi_data), d_(bpm_data->d_), g_(new LampGraph<uint64>(*d_)), bsh_(bpm_data->bsh_),
        sup_buf_(bpm_data->sup_buf_), child_sup_buf_(bpm_data->child_sup_buf_),
        getminsup_data_(NULL), gettestable_data_(NULL),
        phase_(0), expand_num_(0), closed_set_num_(0) {}

ParallelPatternMining::~ParallelPatternMining() {
	// TODO: lots of things to delete
	if (g_) {
        delete g_;
    }
}

// TODO: This should be under GetMinimumSupport
void ParallelPatternMining::PreProcessRootNode(GetMinSupData* getminsup_data) {
	this->getminsup_data_ = getminsup_data;
	phase_ = 1;
	long long int start_time;
	start_time = timer_->Elapsed();
    expand_num_++;

    DFSWithTreeSearchDataStack* stack = dynamic_cast<DFSWithTreeSearchDataStack*>(dfs_stack_);
    TreeSearchData* treesearch_data = stack->getTreeSearchData();

    // what needed is setting itemset_buf_ to empty itemset
    // can be done faster
    // assuming root itemset is pushed before
    treesearch_data->node_stack_->CopyItem(treesearch_data->node_stack_->Top(), treesearch_data->itemset_buf_);
    treesearch_data->node_stack_->Pop();

    // dbg
    DBG( D(2) << "preprocess root node "; );
    DBG( treesearch_data->node_stack_->Print(D(2), treesearch_data->itemset_buf_); );

	// calculate support from itemset_buf_
	bsh_->Set(sup_buf_);
	// skipping rest for root node

	int core_i = g_->CoreIndex(*treesearch_data->node_stack_, treesearch_data->itemset_buf_);

	int * ppc_ext_buf;
	// todo: use database reduction

	bool is_root_node = true;

	// reverse order
	for (int new_item = d_->NextItemInReverseLoop(is_root_node, mpi_data_.mpiRank_, mpi_data_.nTotalProc_, d_->NuItems());
		new_item >= core_i + 1;
		new_item = d_->NextItemInReverseLoop(is_root_node, mpi_data_.mpiRank_, mpi_data_.nTotalProc_, new_item)) {
		// skipping not needed because itemset_buf_ if root itemset

		bsh_->Copy(sup_buf_, child_sup_buf_);
		int sup_num = bsh_->AndCountUpdate(d_->NthData(new_item), child_sup_buf_);

		if (sup_num < getminsup_data->lambda_) {
			continue;
		}

		treesearch_data->node_stack_->PushPre();
		ppc_ext_buf = treesearch_data->node_stack_->Top();
//		treesearch_data->node_stack_->CopyItem(treesearch_data->itemset_buf_, ppc_ext_buf);
		bool res = g_->PPCExtension(treesearch_data->node_stack_, treesearch_data->itemset_buf_, child_sup_buf_, core_i, new_item, ppc_ext_buf);

		treesearch_data->node_stack_->SetSup(ppc_ext_buf, sup_num);
		treesearch_data->node_stack_->PushPostNoSort();

		treesearch_data->node_stack_->Pop(); // always pop

		if (res) {
			IncCsAccum(sup_num); // increment closed_set_num_array
			assert(sup_num >= getminsup_data->lambda_);
			if (ExceedCsThr()) {
				getminsup_data->lambda_ = NextLambdaThr();
			}
		}
	}
	//	 TODO: lambda_max_ is wrong???
	printf("lambda_max_ = %d\n", getminsup_data->lambda_max_);
	MPI_Reduce(getminsup_data->accum_array_, getminsup_data->accum_recv_, getminsup_data->lambda_max_ + 1,
			   MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD); // error?

	if (mpi_data_.mpiRank_ == 0) {
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

	if (mpi_data_.mpiRank_ == 0) {
        if (ExceedCsThr()) {
            getminsup_data->lambda_ = NextLambdaThr();
        }
    }
	CallBcast(&getminsup_data->lambda_, 1, MPI_INT);

	// reverse order
	for (int new_item = d_->NextItemInReverseLoop(is_root_node,	mpi_data_.mpiRank_, mpi_data_.nTotalProc_, d_->NuItems());
			new_item >= core_i + 1;
			new_item = d_->NextItemInReverseLoop(is_root_node, mpi_data_.mpiRank_, mpi_data_.nTotalProc_, new_item)) {
		// skipping not needed because itemset_buf_ if root itemset

		bsh_->Copy(sup_buf_, child_sup_buf_);
		int sup_num = bsh_->AndCountUpdate(d_->NthData(new_item), child_sup_buf_);

		if (sup_num < getminsup_data->lambda_) {
            continue;
        }

		treesearch_data->node_stack_->PushPre();
		ppc_ext_buf = treesearch_data->node_stack_->Top();
//		treesearch_data->node_stack_->CopyItem(treesearch_data->itemset_buf_,
//				ppc_ext_buf);

		bool res = g_->PPCExtension(treesearch_data->node_stack_,
				treesearch_data->itemset_buf_, child_sup_buf_, core_i,
				new_item, ppc_ext_buf);

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
                << getminsup_data_->lambda_ << "\ttime="
                << elapsed_time << std::endl;
    );
}

void ParallelPatternMining::GetMinimalSupport(GetMinSupData* getminsup_data) {
	this->getminsup_data_ = getminsup_data;
//	CheckInit();
	phase_ = 1; // TODO: remove dependency on this
	Search();

	// return lambda?
}

void ParallelPatternMining::GetTestablePatterns(GetTestableData* gettestable_data) {
	this->gettestable_data_ = gettestable_data; // TODO: not sure if this is a good idea.
	this->getminsup_data_->lambda_ = gettestable_data->freqThreshold_; // TODO: not the best way.
	CheckInitTestable();
	phase_ = 2;  // TODO: remove dependency on this
	DBG(D(1) << "MainLoop" << std::endl
	;);
	printf("GetTestablePatterns\n");
	Search();

}

void ParallelPatternMining::GetSignificantPatterns(GetSignificantData* getsignificant_data) {
	this->getsignificant_data_ = getsignificant_data;
	DBG( D(1) << "MainLoop" << std::endl; );
	printf("extract\n");
	ExtractSignificantSet();
	if (mpi_data_.mpiRank_ == 0) {
		printf("sendresults\n");
		SendResultRequest();
	}
	printf("probe\n");
	while (!mpi_data_.dtd_->terminated_) {
		Probe();
	}
}

//==============================================================================
/**
 * Procedure to run after Probe().
 */
void ParallelPatternMining::ProcAfterProbe() {
	if (phase_ == 1) {
		log_->TakePeriodicLog(dfs_stack_->Count(), getminsup_data_->lambda_, phase_);
	} else {
		log_->TakePeriodicLog(dfs_stack_->Count(), gettestable_data_->freqThreshold_, phase_);
	}

	if (mpi_data_.mpiRank_ == 0) {
		// initiate termination detection

		// note: for phase_ 1, accum request and dtd request are unified
		if (!mpi_data_.echo_waiting_ && !mpi_data_.dtd_->terminated_) {
			if (phase_ == 1) {
				SendDTDAccumRequest();
				// log_->d_.dtd_accum_request_num_++;
			} else if (phase_ == 2) {
				if (dfs_stack_->IsEmpty()) {
					printf("Empty: root sends DTDRequest.\n");
					SendDTDRequest();
					// log_->d_.dtd_request_num_++;
				}
			} else {
				// unknown phase
				assert(0);
			}
		}
	}
}

void ParallelPatternMining::ProbeExecute(MPI_Status* probe_status, int probe_src, int probe_tag) {
	switch (probe_tag) {
	/**
	 * MINIMALSUPPORT
	 */
	case Tag::DTD_ACCUM_REQUEST:
		assert(phase_ == 1);
		RecvDTDAccumRequest(probe_src); // TODO: This can be solved by polymorphism.
		break;
	case Tag::DTD_ACCUM_REPLY:
		assert(phase_ == 1);
		RecvDTDAccumReply(probe_src);
		break;
	case Tag::LAMBDA:
		assert(phase_ == 1);
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
                    << std::endl;
		);
		ParallelDFS::ProbeExecute(probe_status, probe_src, probe_tag);
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
void ParallelPatternMining::Check() {
	if (mpi_data_.mpiRank_ == 0 && phase_ == 1) {
		CheckCSThreshold();
	}
	return;
}

// TODO: polymorphism to override SendDTDRequest
void ParallelPatternMining::SendDTDAccumRequest() {
	int message[1];
	message[0] = 1; // dummy

	mpi_data_.echo_waiting_ = true;

	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data_.bcast_targets_[i] < 0) {
			break;
		}
		assert(mpi_data_.bcast_targets_[i] < mpi_data_.nTotalProc_ && "SendDTDAccumRequest");
		CallBsend(message, 1, MPI_INT, mpi_data.bcast_targets_[i], Tag::DTD_ACCUM_REQUEST);

		DBG(
			D(3) << "SendDTDAccumRequest: dst="
					<< mpi_data_.bcast_targets_[i] << "\ttimezone="
					<< mpi_data_.dtd_->time_zone_ << std::endl;
		);
	}
}

void ParallelPatternMining::RecvDTDAccumRequest(int src) {
	DBG(
		D(3) << "RecvDTDAccumRequest: src=" << src
				<< "\ttimezone=" << mpi_data_.dtd_->time_zone_
				<< std::endl;
	);
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::DTD_ACCUM_REQUEST, &recv_status);

	if (IsLeafInTopology()) {
		SendDTDAccumReply();
	} else {
		SendDTDAccumRequest();
	}
}

void ParallelPatternMining::SendDTDAccumReply() {
	getminsup_data_->dtd_accum_array_base_[0] = mpi_data_.dtd_->count_ + mpi_data_.dtd_->reduce_count_;
	bool tw_flag = mpi_data_.dtd_->time_warp_ || mpi_data_.dtd_->reduce_time_warp_;
	getminsup_data_->dtd_accum_array_base_[1] = (tw_flag ? 1 : 0);

    // for Steal
    // thieves_ and stealer state check
	mpi_data_.dtd_->not_empty_ =	!(dfs_stack_->IsEmpty()) || (mpi_data_.thieves_->Size() > 0) || dfs_stack_->IsStealStarted()	|| mpi_data_.processing_node_;

    // for Steal2
	bool em_flag = mpi_data_.dtd_->not_empty_ || mpi_data_.dtd_->reduce_not_empty_;
	getminsup_data_->dtd_accum_array_base_[2] = (em_flag ? 1 : 0);

	DBG(
        D(3) << "SendDTDAccumReply: dst = "
                << mpi_data_.bcast_source_ << "\tcount="
                << getminsup_data_->dtd_accum_array_base_[0]
                << "\ttw=" << tw_flag << "\tem=" << em_flag
                << std::endl;
    );

	assert(mpi_data_.bcast_source_ < mpi_data_.nTotalProc_ && "SendDTDAccumReply");
	CallBsend(getminsup_data_->dtd_accum_array_base_, getminsup_data_->lambda_max_ + 4, MPI_LONG_LONG_INT, mpi_data_.bcast_source_, Tag::DTD_ACCUM_REPLY);

	mpi_data_.dtd_->time_warp_ = false;
	mpi_data_.dtd_->not_empty_ = false;
	mpi_data_.dtd_->IncTimeZone();

	mpi_data_.echo_waiting_ = false;
	mpi_data_.dtd_->ClearAccumFlags();
	mpi_data_.dtd_->ClearReduceVars();

	for (int i = 0; i <= getminsup_data_->lambda_max_; i++) {
        getminsup_data_->accum_array_[i] = 0ll;
    }
}

// getminsup_data
void ParallelPatternMining::RecvDTDAccumReply(int src) {
	MPI_Status recv_status;

	CallRecv(getminsup_data_->dtd_accum_recv_base_, getminsup_data_->lambda_max_ + 4, MPI_LONG_LONG_INT, src, Tag::DTD_ACCUM_REPLY, &recv_status);
	assert(src == recv_status.MPI_SOURCE);

	int count = (int) (getminsup_data_->dtd_accum_recv_base_[0]);
	bool time_warp = (getminsup_data_->dtd_accum_recv_base_[1] != 0);
	bool not_empty = (getminsup_data_->dtd_accum_recv_base_[2] != 0);

	mpi_data_.dtd_->Reduce(count, time_warp, not_empty);

	DBG(
        D(3) << "RecvDTDAccumReply: src=" << src << "\tcount="
                << count << "\ttw=" << time_warp << "\tem="
                << not_empty << "\treduced_count="
                << mpi_data_.dtd_->reduce_count_ << "\treduced_tw="
                << mpi_data_.dtd_->reduce_time_warp_
                << "\treduced_em="
                << mpi_data_.dtd_->reduce_not_empty_ << std::endl;
    );

	for (int i = getminsup_data_->lambda_ - 1; i <= getminsup_data_->lambda_max_; i++) {
        getminsup_data_->accum_array_[i] += getminsup_data_->accum_recv_[i];
    }

	bool flag = false;
	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data_.bcast_targets_[i] == src) {
			flag = true;
			mpi_data_.dtd_->accum_flag_[i] = true;
			break;
		}
	}
	assert(flag);

	if (mpi_data_.mpiRank_ == 0) {
		if (ExceedCsThr()) {
			int new_lambda = NextLambdaThr();
			SendLambda(new_lambda);
			getminsup_data_->lambda_ = new_lambda;
		}
// if SendLambda is called, dtd_.count_ is incremented and DTDCheck will always fail
		if (DTDReplyReady()) {
			DTDCheck();
			log_->d_.dtd_accum_phase_num_++;
		}
	} else {  // not root
		if (DTDReplyReady()) {
            SendDTDAccumReply();
        }
	}
}

/**
 * GetMinSup Functions
 *
 *
 */
void ParallelPatternMining::SendLambda(int lambda) {
// send lambda to bcast_targets_
	int message[2];
	message[0] = mpi_data_.dtd_->time_zone_;
	message[1] = lambda;

	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data_.bcast_targets_[i] < 0) {
            break;
        }
		assert(mpi_data_.bcast_targets_[i] < mpi_data_.nTotalProc_ && "SendLambda");
		CallBsend(message, 2, MPI_INT, mpi_data_.bcast_targets_[i], Tag::LAMBDA);
		mpi_data_.dtd_->OnSend();

		DBG(
            D(2) << "SendLambda: dst="
                    << mpi_data_.bcast_targets_[i] << "\tlambda="
                    << lambda << "\tdtd_count="
                    << mpi_data_.dtd_->count_ << std::endl;
        );
	}
}

void ParallelPatternMining::RecvLambda(int src) {
	MPI_Status recv_status;
	int message[2];

	CallRecv(&message, 2, MPI_INT, src, Tag::LAMBDA, &recv_status);
	mpi_data_.dtd_->OnRecv();
	assert(src == recv_status.MPI_SOURCE);
	int timezone = message[0];
	mpi_data_.dtd_->UpdateTimeZone(timezone);

	DBG(
        D(2) << "RecvLambda: src=" << src << "\tlambda="
                << message[1] << "\tdtd_count="
                << mpi_data_.dtd_->count_ << std::endl;
    );

	int new_lambda = message[1];
	if (new_lambda > getminsup_data_->lambda_) {
		SendLambda(new_lambda);
		getminsup_data_->lambda_ = new_lambda;
// todo: do database reduction
	}
}

void ParallelPatternMining::CheckCSThreshold() {
//	assert(mpi_data_.mpiRank_ == 0);
	if (ExceedCsThr()) {
		int new_lambda = NextLambdaThr();
		SendLambda(new_lambda);
		getminsup_data_->lambda_ = new_lambda;
	}
}

bool ParallelPatternMining::ExceedCsThr() const {
// note: > is correct. permit ==
	return (getminsup_data_->accum_array_[getminsup_data_->lambda_]
			> getminsup_data_->cs_thr_[getminsup_data_->lambda_]);
}

int ParallelPatternMining::NextLambdaThr() const {
	int si;
	for (si = getminsup_data_->lambda_max_; si >= getminsup_data_->lambda_; si--) {
        if (getminsup_data_->accum_array_[si] > getminsup_data_->cs_thr_[si]) {
            break;
        }
    }
	return si + 1;
// it is safe because lambda_ higher than max results in immediate search finish
}

/**
 * GETSIGNIFICANT PATTERNS
 *
 */

//==============================================================================
void ParallelPatternMining::SendResultRequest() {
	int message[1];
	message[0] = 1; // dummy

	mpi_data_.echo_waiting_ = true;

	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data_.bcast_targets_[i] < 0) {
            break;
        }

		assert(mpi_data_.bcast_targets_[i] < mpi_data_.nTotalProc_ && "SendResultRequest");
		CallBsend(message, 1, MPI_INT, mpi_data_.bcast_targets_[i], Tag::RESULT_REQUEST);
		DBG( D(2) << "SendResultRequest: dst=" << mpi_data_.bcast_targets_[i] << std::endl; );
	}
}

void ParallelPatternMining::RecvResultRequest(int src) {
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::RESULT_REQUEST, &recv_status);
	assert(src == recv_status.MPI_SOURCE);

	DBG( D(2) << "RecvResultRequest: src=" << src << std::endl; );

	if (IsLeafInTopology()) {
        SendResultReply();
    } else {
        SendResultRequest();
    }
}

void ParallelPatternMining::SendResultReply() {
	int * message = getsignificant_data_->significant_stack_->Stack();
	int size = getsignificant_data_->significant_stack_->UsedCapacity();
	assert(mpi_data_.bcast_source_ < mpi_data_.nTotalProc_ && "SendResultReply");
	CallBsend(message, size, MPI_INT, mpi_data_.bcast_source_, Tag::RESULT_REPLY);

	DBG( D(2) << "SendResultReply: dst=" << mpi_data_.bcast_source_ << std::endl; );
	DBG( getsignificant_data_->significant_stack_->PrintAll(D(3, false)); );

	mpi_data_.echo_waiting_ = false;
	mpi_data_.dtd_->terminated_ = true;
}

void ParallelPatternMining::RecvResultReply(int src, MPI_Status probe_status) {
	int count;
	int error = MPI_Get_count(&probe_status, MPI_INT, &count);
	if (error != MPI_SUCCESS) {
		DBG( D(1) << "error in MPI_Get_count in RecvResultReply: " << error << std::endl; );
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	MPI_Status recv_status;
	CallRecv(give_stack_->Stack(), count, MPI_INT, src, Tag::RESULT_REPLY, &recv_status);
	assert(src == recv_status.MPI_SOURCE);

	getsignificant_data_->significant_stack_->MergeStack(
            give_stack_->Stack() + VariableLengthItemsetStack::SENTINEL + 1,
            count - VariableLengthItemsetStack::SENTINEL - 1);
	give_stack_->Clear();

	DBG( D(2) << "RecvResultReply: src=" << src << std::endl; );
	DBG( getsignificant_data_->significant_stack_->PrintAll(D(3, false) )
	;);

// using the same flags as accum count, should be fixed
	bool flag = false;
	for (int i = 0; i < k_echo_tree_branch; i++) {
		if (mpi_data_.bcast_targets_[i] == src) {
			flag = true;
			mpi_data_.accum_flag_[i] = true;
			break;
		}
	}
	assert(flag);

	if (AccumCountReady()) {
		if (mpi_data_.mpiRank_ != 0) {
			SendResultReply();
		} else { // root
			mpi_data_.echo_waiting_ = false;
			mpi_data_.dtd_->terminated_ = true;
		}
	}
}

void ParallelPatternMining::ExtractSignificantSet() {
	std::multimap<double, int *>::iterator it;
	for (it = getsignificant_data_->freq_map_->begin(); it != getsignificant_data_->freq_map_->end(); ++it) {
// permits == case
		if ((*it).first <= getsignificant_data_->final_sig_level_) {
            getsignificant_data_->significant_stack_->PushPre();
			int * item = getsignificant_data_->significant_stack_->Top();
			getsignificant_data_->significant_stack_->CopyItem((*it).second, item);
			getsignificant_data_->significant_stack_->PushPostNoSort();
		} else {
            break;
        }
	}
}

void ParallelPatternMining::CheckInit() {
	assert(mpi_data_.echo_waiting_ == false);
	assert(mpi_data_.waiting_ == false);
	assert(mpi_data_.thieves_->Size() == 0);
	assert(mpi_data_.lifeline_thieves_->Size() == 0);
	for (int pi = 0; pi < mpi_data_.nTotalProc_; pi++) {
        assert(mpi_data_.lifelines_activated_[pi] == false);
    }
}

// TODO: For debugging: Not sure how to show the caller in assert.
void ParallelPatternMining::CheckInitTestable() {
//	mpi_data_.lifeline_thieves_->Clear();
//	for (int pi = 0; pi < mpi_data_.nTotalProc_; pi++) {
//		mpi_data_.lifelines_activated_[pi] = false;
//	}

	assert(mpi_data_.echo_waiting_ == false);
	assert(mpi_data_.waiting_ == false);
//	assert(mpi_data_.thieves_->Size() == 0);
//	assert(mpi_data_.lifeline_thieves_->Size() == 0);
//	for (int pi = 0; pi < mpi_data_.nTotalProc_; pi++) {
//		assert(mpi_data_.lifelines_activated_[pi] == false);
//	}
}

	void ParallelPatternMining::CallbackForProbe() {

	}

    void ParallelPatternMining::UpdateFreqStack(double pval, int* ppc_ext_buf) {
        if (pval <= gettestable_data_->sig_level_) { // permits == case?
            gettestable_data_->freq_stack_->PushPre();
            int *item = gettestable_data_->freq_stack_->Top();
            gettestable_data_->freq_stack_->CopyItem(ppc_ext_buf, item);
            gettestable_data_->freq_stack_->PushPostNoSort();

            gettestable_data_->freq_map_->insert(std::pair<double, int *>(pval, item));
        }
    }

    int ParallelPatternMining::CalcCoreI(VariableLengthItemsetStack* node_stack, const int* itemset_buf) const {
        g_->CoreIndex(*node_stack, itemset_buf);
    }

    int ParallelPatternMining::GetLamba() const {
        return getminsup_data_->lambda_;
    }

    void ParallelPatternMining::IncCsAccum(int sup_num) {
        for (int i = sup_num; i >= getminsup_data_->lambda_ - 1; i--) {
            getminsup_data_->accum_array_[i]++;
        }
    }

    int ParallelPatternMining::UpdateChildSupBuf(int new_item, int core_i) {
        // todo: if database reduction is implemented,
        //       do something here for changed lambda_ (skipping new_item value ?)

        bsh_->Copy(sup_buf_, child_sup_buf_);
        int sup_num = bsh_->AndCountUpdate(d_->NthData(new_item), child_sup_buf_);
        // If the support is smaller than the required minimal support for
        // significant pattern (=lambda), then prune it.
        if (sup_num < getminsup_data_->lambda_) {
            return -1;
        }

        return sup_num;
    }

    bool ParallelPatternMining::PPCExtension(VariableLengthItemsetStack* node_stack, const int* itemset_buf, int core_i, int new_item, int* ppc_ext_buf) {
        return g_->PPCExtension(node_stack, itemset_buf, child_sup_buf_, core_i, new_item, ppc_ext_buf);
    }

} /* namespace lamp_search */
