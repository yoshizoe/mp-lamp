/*
 * ParallelSearch.cc
 *
 *  Created on: Apr 18, 2017
 *      Author: yuu
 */

#include "ParallelDFS.h"

#ifdef __CDT_PARSER__
#undef DBG
#define DBG(a)  a
#endif

#ifdef __CDT_PARSER__
#undef LOG
#define LOG(a)  ;
#endif

namespace lamp_search {

const int ParallelDFS::k_echo_tree_branch = 3;

ParallelDFS::ParallelDFS(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, Log* log, Timer* timer,
		std::ostream& ofs) :
		mpi_data(mpi_data), treesearch_data(treesearch_data), log_(
				log), timer_(timer), lfs_(ofs) {

}

ParallelDFS::~ParallelDFS() {
	// TODO Auto-generated destructor stub
}

void ParallelDFS::Search() {

	printf("ParallelDFS::Search\n");
	DBG(D(1) << "MainLoop" << std::endl
	; );
	while (!mpi_data.dtd_->terminated_) {
		while (!mpi_data.dtd_->terminated_) {
			if (ExpandNode(mpi_data, treesearch_data)) {
				log_->d_.node_stack_max_itm_ =
						std::max(log_->d_.node_stack_max_itm_,
								(long long int) (treesearch_data->node_stack_->NuItemset()));
				log_->d_.node_stack_max_cap_ =
						std::max(log_->d_.node_stack_max_cap_,
								(long long int) (treesearch_data->node_stack_->UsedCapacity()));
				Probe(mpi_data, treesearch_data);
				if (mpi_data.dtd_->terminated_)
					break;
				Distribute(mpi_data, treesearch_data);
				Reject(mpi_data); // distribute finished, reject remaining requests

				Check(mpi_data); // Implement CheckCSThreshold().

			} else
				break;
		}
		if (mpi_data.dtd_->terminated_)
			break;

		log_->idle_start_ = timer_->Elapsed();
		Reject(mpi_data); // node_stack_ empty. reject requests
		Steal(mpi_data); // request steal
		if (mpi_data.dtd_->terminated_) {
			log_->d_.idle_time_ += timer_->Elapsed()
					- log_->idle_start_;
			break;
		}

		Probe(mpi_data, treesearch_data);
		if (mpi_data.dtd_->terminated_) {
			log_->d_.idle_time_ += timer_->Elapsed()
					- log_->idle_start_;
			break;
		}
		Check(mpi_data); // Implement CheckCSThreshold().

		log_->d_.idle_time_ += timer_->Elapsed() - log_->idle_start_;
	}
}

//==============================================================================

bool ParallelDFS::Probe(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data) {
	DBG(D(3) << "Probe" << std::endl
	; );
	MPI_Status probe_status;
	int probe_src, probe_tag;

	long long int start_time;
	start_time = timer_->Elapsed();

	log_->d_.probe_num_++;

	while (CallIprobe(&probe_status, &probe_src, &probe_tag)) {
		DBG(
				D(4) << "CallIprobe returned src=" << probe_src
						<< "\ttag=" << probe_tag << std::endl
				; );
		ProbeExecute(mpi_data, treesearch_data, &probe_status,
				probe_src, probe_tag);
	}

	// capacity, lambda, phase
	assert(treesearch_data->node_stack_);

	/**
	 * Domain Specifics
	 */
	ProcAfterProbe();

	long long int elapsed_time = timer_->Elapsed() - start_time;
	log_->d_.probe_time_ += elapsed_time;
	log_->d_.probe_time_max_ = std::max(elapsed_time,
			log_->d_.probe_time_max_);
}

void ParallelDFS::Distribute(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data) {
	if (mpi_data.nTotalProc_ == 1)
		return;
	DBG(D(3) << "Distribute" << std::endl
	; );
	if (mpi_data.thieves_->Size() > 0
			|| mpi_data.lifeline_thieves_->Size() > 0) {
		int steal_num = treesearch_data->node_stack_->Split(
				treesearch_data->give_stack_);
		if (steal_num > 0) {
			DBG(D(3) << "giving" << std::endl
			; );
			DBG(treesearch_data->give_stack_->PrintAll(D(3, false))
			; );
			Give(mpi_data, treesearch_data->give_stack_, steal_num);
			treesearch_data->give_stack_->Clear();
		}
	}
}

// TODO ParallelDFS
void ParallelDFS::Give(MPI_Data& mpi_data,
		VariableLengthItemsetStack * st, int steal_num) {
	DBG(D(3) << "Give: "
	; );
	if (mpi_data.thieves_->Size() > 0) { // random thieves
		int thief = mpi_data.thieves_->Pop();
		if (thief >= 0) { // lifeline thief
			DBG(
					D(3, false) << "thief=" << thief
							<< "\trandom lifeline" << std::endl
					; );
			log_->d_.lifeline_given_num_++;
			log_->d_.lifeline_nodes_given_ += steal_num;
			SendGive(mpi_data, st, thief, mpi_data.mpiRank_);
		} else { // random thief
			DBG(
					D(3, false) << "thief=" << (-thief - 1)
							<< "\trandom" << std::endl
					; );
			log_->d_.given_num_++;
			log_->d_.nodes_given_ += steal_num;
			SendGive(mpi_data, st, (-thief - 1), -1);
		}
	} else {
		assert(mpi_data.lifeline_thieves_->Size() > 0);
		int thief = mpi_data.lifeline_thieves_->Pop();
		assert(thief >= 0);
		DBG(
				D(3, false) << "thief=" << thief << "\tlifeline"
						<< std::endl
				; );
		log_->d_.lifeline_given_num_++;
		log_->d_.lifeline_nodes_given_ += steal_num;
		SendGive(mpi_data, st, thief, mpi_data.mpiRank_);
	}
}
// TODO ParallelDFS
void ParallelDFS::Reject(MPI_Data& mpi_data) {
	if (mpi_data.nTotalProc_ == 1)
		return;
	DBG(D(3) << "Reject" << std::endl
	; );
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

// TODO ParallelDFS
void ParallelDFS::Steal(MPI_Data& mpi_data) {
	DBG(D(3) << "Steal" << std::endl
	; );
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
						<< treesearch_data->stealer_->RandomCount()
						<< std::endl
				; );
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
					; );
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
			; );
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

//==============================================================================
/**
 * This function should be called if all the subclass ProbeExecute failed.
 *
 */
void ParallelDFS::ProbeExecute(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, MPI_Status* probe_status,
		int probe_src, int probe_tag) {
//	if (ParallelSearch::ProbeExecute(mpi_data, XX)) {
//		return true;
//	}

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
		RecvGive(mpi_data, treesearch_data, probe_src, *probe_status);
		break;

	default:
		printf("UNKNOWN TAG\n");
		DBG(
				D(1) << "unknown Tag=" << probe_tag
						<< " received in Probe: " << std::endl
				;
		);
		MPI_Abort(MPI_COMM_WORLD, 1);
		break;
	}
	return;
}

//==============================================================================
bool ParallelDFS::AccumCountReady(MPI_Data& mpi_data) const {
	for (int i = 0; i < k_echo_tree_branch; i++)
		if (mpi_data.bcast_targets_[i] >= 0 && !(mpi_data.accum_flag_[i]))
			return false;
	return true;
// check only valid bcast_targets_
// always true if leaf
}

// call this from root rank to start DTD
void ParallelDFS::SendDTDRequest(MPI_Data& mpi_data) {
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
				D(3) << "SendDTDRequest: dst="
						<< mpi_data.bcast_targets_[i] << "\ttimezone="
						<< mpi_data.dtd_->time_zone_ << std::endl
				; );
	}
}

void ParallelDFS::RecvDTDRequest(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, int src) {
	DBG(
			D(3) << "RecvDTDRequest: src=" << src << "\ttimezone="
					<< mpi_data.dtd_->time_zone_ << std::endl
			; );
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::DTD_REQUEST,
			&recv_status);

	if (IsLeafInTopology(mpi_data))
		SendDTDReply(mpi_data, treesearch_data);
	else
		SendDTDRequest(mpi_data);
}

bool ParallelDFS::DTDReplyReady(MPI_Data& mpi_data) const {
	for (int i = 0; i < k_echo_tree_branch; i++)
		if (mpi_data.bcast_targets_[i] >= 0
				&& !(mpi_data.dtd_->accum_flag_[i]))
			return false;
	return true;
// check only valid bcast_targets_
// always true if leaf
}

void ParallelDFS::SendDTDReply(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data) {
	int message[3];
// using reduced vars
	message[0] = mpi_data.dtd_->count_ + mpi_data.dtd_->reduce_count_;
	bool tw_flag = mpi_data.dtd_->time_warp_
			|| mpi_data.dtd_->reduce_time_warp_;
	message[1] = (tw_flag ? 1 : 0);
// for Steal
	mpi_data.dtd_->not_empty_ =
			!(treesearch_data->node_stack_->Empty())
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
					<< "\tcount=" << message[0] << "\ttw=" << tw_flag
					<< "\tem=" << em_flag << std::endl
			; );

	assert(
			mpi_data.bcast_source_ < mpi_data.nTotalProc_
					&& "SendDTDReply");
	CallBsend(message, 3, MPI_INT, mpi_data.bcast_source_,
			Tag::DTD_REPLY);

	mpi_data.dtd_->time_warp_ = false;
	mpi_data.dtd_->not_empty_ = false;
	mpi_data.dtd_->IncTimeZone();

	mpi_data.echo_waiting_ = false;
	mpi_data.dtd_->ClearAccumFlags();
	mpi_data.dtd_->ClearReduceVars();
}

void ParallelDFS::RecvDTDReply(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, int src) {
	MPI_Status recv_status;
	int message[3];
	CallRecv(&message, 3, MPI_INT, src, Tag::DTD_REPLY, &recv_status);
	assert(src == recv_status.MPI_SOURCE);

// reduce reply (count, time_warp, not_empty)
	mpi_data.dtd_->Reduce(message[0], (message[1] != 0),
			(message[2] != 0));

	DBG(
			D(3) << "RecvDTDReply: src=" << src << "\tcount="
					<< message[0] << "\ttw=" << message[1] << "\tem="
					<< message[2] << "\treduced_count="
					<< mpi_data.dtd_->reduce_count_ << "\treduced_tw="
					<< mpi_data.dtd_->reduce_time_warp_
					<< "\treduced_em="
					<< mpi_data.dtd_->reduce_not_empty_ << std::endl
			; );

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

void ParallelDFS::DTDCheck(MPI_Data& mpi_data) {
	assert(mpi_data.mpiRank_ == 0);
	// (count, time_warp, not_empty)
	mpi_data.dtd_->Reduce(mpi_data.dtd_->count_,
			mpi_data.dtd_->time_warp_, mpi_data.dtd_->not_empty_);

	if (mpi_data.dtd_->reduce_count_ == 0
			&& mpi_data.dtd_->reduce_time_warp_ == false
			&& mpi_data.dtd_->reduce_not_empty_ == false) {
// termination
		SendBcastFinish(mpi_data);
		mpi_data.dtd_->terminated_ = true;
		DBG(D(1) << "terminated" << std::endl
		; );
	}
	// doing same thing as SendDTDReply
	mpi_data.dtd_->time_warp_ = false;
	mpi_data.dtd_->not_empty_ = false;
	mpi_data.dtd_->IncTimeZone();

	mpi_data.echo_waiting_ = false;
	mpi_data.dtd_->ClearAccumFlags();
	mpi_data.dtd_->ClearReduceVars();
}

void ParallelDFS::SendBcastFinish(MPI_Data& mpi_data) {
	DBG(D(2) << "SendBcastFinish" << std::endl
	; );
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

void ParallelDFS::RecvBcastFinish(MPI_Data& mpi_data, int src) {
	DBG(D(2) << "RecvBcastFinish: src=" << src << std::endl
	; );
	MPI_Status recv_status;
	int message[1];
	CallRecv(&message, 1, MPI_INT, src, Tag::BCAST_FINISH,
			&recv_status);

	SendBcastFinish(mpi_data);
	mpi_data.dtd_->terminated_ = true;
	DBG(D(2) << "terminated" << std::endl
	; );

	mpi_data.waiting_ = false;
}

// lifeline == -1 for random thieves
void ParallelDFS::SendRequest(MPI_Data& mpi_data, int dst,
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
					<< is_lifeline << "\tdtd_count="
					<< mpi_data.dtd_->count_ << std::endl
			; );
}

void ParallelDFS::RecvRequest(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, int src) {
	DBG(D(2) << "RecvRequest: src=" << src
	; );
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
	; );
	DBG(D(2, false) << "\tthief=" << src
	; );
	DBG(D(2, false) << "\tdtd_count=" << mpi_data.dtd_->count_
	; );

	if (treesearch_data->node_stack_->Empty()) {
		DBG(D(2, false) << "\tempty and reject" << std::endl
		; );
		if (is_lifeline >= 0) {
			mpi_data.lifeline_thieves_->Push(thief);
			SendReject(mpi_data, thief); // notify
		} else {
			SendReject(mpi_data, thief); // notify
		}
	} else {
		DBG(D(2, false) << "\tpush" << std::endl
		; );
		if (is_lifeline >= 0)
			mpi_data.thieves_->Push(thief);
		else
			mpi_data.thieves_->Push(-thief - 1);
	}
	// todo: take log
}

void ParallelDFS::SendReject(MPI_Data& mpi_data, int dst) {
	assert(dst >= 0);
	int message[1];

	message[0] = mpi_data.dtd_->time_zone_;
	assert(dst < mpi_data.nTotalProc_ && "SendReject");
	CallBsend(message, 1, MPI_INT, dst, Tag::REJECT);
	mpi_data.dtd_->OnSend();

	DBG(
			D(2) << "SendReject: dst=" << dst << "\tdtd_count="
					<< mpi_data.dtd_->count_ << std::endl
			; );
}

void ParallelDFS::RecvReject(MPI_Data& mpi_data, int src) {
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
			; );

	treesearch_data->stealer_->ResetRequesting();
	mpi_data.waiting_ = false;
}

void ParallelDFS::SendGive(MPI_Data& mpi_data,
		VariableLengthItemsetStack * st, int dst, int is_lifeline) {
	assert(dst >= 0);
	st->SetTimestamp(mpi_data.dtd_->time_zone_);
	st->SetFlag(is_lifeline);
	int * message = st->Stack();
	int size = st->UsedCapacity();

	log_->d_.give_stack_max_itm_ = std::max(
			log_->d_.give_stack_max_itm_,
			(long long int) (st->NuItemset()));
	log_->d_.give_stack_max_cap_ = std::max(
			log_->d_.give_stack_max_cap_,
			(long long int) (st->UsedCapacity()));

	assert(dst < mpi_data.nTotalProc_ && "SendGive");
	CallBsend(message, size, MPI_INT, dst, Tag::GIVE);
	mpi_data.dtd_->OnSend();

	DBG(
			D(2) << "SendGive: " << "\ttimezone="
					<< mpi_data.dtd_->time_zone_

					<< "\tdst="

					<< dst << "\tlfl=" << is_lifeline << "\tsize="
					<< size << "\tnode=" << st->NuItemset()
					<< "\tdtd_count=" << mpi_data.dtd_->count_
					<< std::endl
			; );
	//st->PrintAll(D(false));
}

void ParallelDFS::RecvGive(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, int src,
		MPI_Status probe_status) {
	int count;
	int error = MPI_Get_count(&probe_status, MPI_INT, &count);
	if (error != MPI_SUCCESS) {
		DBG(
				D(1) << "error in MPI_Get_count in RecvGive: "
						<< error << std::endl
				; );
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	MPI_Status recv_status;

	CallRecv(treesearch_data->give_stack_->Stack(), count, MPI_INT,
			src, Tag::GIVE, &recv_status);
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
		log_->d_.lifeline_nodes_received_ += (new_nu_itemset
				- orig_nu_itemset);
	} else {
		log_->d_.steal_num_++;
		log_->d_.nodes_received_ +=
				(new_nu_itemset - orig_nu_itemset);
	}

	DBG(
			D(2) << "RecvGive: src=" << src << "\ttimezone="

			<< mpi_data.dtd_->time_zone_

			<< "\tlfl=" << flag << "\tsize=" << count << "\tnode="
					<< (new_nu_itemset - orig_nu_itemset)
					<< "\tdtd_count=" << mpi_data.dtd_->count_
					<< std::endl
			; );

	treesearch_data->give_stack_->Clear();

	log_->d_.node_stack_max_itm_ =
			std::max(log_->d_.node_stack_max_itm_,
					(long long int) (treesearch_data->node_stack_->NuItemset()));
	log_->d_.node_stack_max_cap_ =
			std::max(log_->d_.node_stack_max_cap_,
					(long long int) (treesearch_data->node_stack_->UsedCapacity()));

	DBG(treesearch_data->node_stack_->PrintAll(D(3, false))
	; );

	treesearch_data->stealer_->ResetRequesting();
	treesearch_data->stealer_->ResetCounters();
	treesearch_data->stealer_->SetState(StealState::RANDOM);
	treesearch_data->stealer_->SetStealStart();
	mpi_data.waiting_ = false;
}

//==============================================================================

bool ParallelDFS::IsLeafInTopology(MPI_Data& mpi_data) const {
	for (int i = 0; i < k_echo_tree_branch; i++)
		if (mpi_data.bcast_targets_[i] >= 0)
			return false;
	return true;
}

int ParallelDFS::CallIprobe(MPI_Status * status, int * src,
		int * tag) {
	long long int start_time;
	long long int end_time;
	log_->d_.iprobe_num_++;

	// todo: prepare non-log mode to remove measurement
	// clock_gettime takes 0.3--0.5 micro sec
	LOG(start_time = timer_->Elapsed() ;);

	int flag;
	int error = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG,
			MPI_COMM_WORLD, &flag, status);
	if (error != MPI_SUCCESS) {
		DBG(D(1) << "error in MPI_Iprobe: " << error << std::endl
		; );
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

int ParallelDFS::CallRecv(void * buffer, int data_count,
		MPI_Datatype type, int src, int tag, MPI_Status * status) {
	long long int start_time;
	long long int end_time;

	// todo: prepare non-log mode to remove measurement
	// clock_gettime takes 0.3--0.5 micro sec
	log_->d_.recv_num_++;
	LOG(start_time = timer_->Elapsed() ;);

	int error = MPI_Recv(buffer, data_count, type, src, tag,
	MPI_COMM_WORLD, status);

	LOG(
			end_time = timer_->Elapsed(); log_->d_.recv_time_ += end_time - start_time; log_->d_.recv_time_max_ = std::max(end_time - start_time, log_->d_.recv_time_max_););
	return error;
}

int ParallelDFS::CallBsend(void * buffer, int data_count,
		MPI_Datatype type, int dest, int tag) {
//	assert(0 <= dest && dest < mpi_data.nTotalProc_);
	long long int start_time;
	long long int end_time;
	log_->d_.bsend_num_++;
	start_time = timer_->Elapsed();

	int error = MPI_Bsend(buffer, data_count, type, dest, tag,
	MPI_COMM_WORLD);

	end_time = timer_->Elapsed();
	log_->d_.bsend_time_ += end_time - start_time;
	log_->d_.bsend_time_max_ = std::max(end_time - start_time,
			log_->d_.bsend_time_max_);
	return error;
}

int ParallelDFS::CallBcast(void * buffer, int data_count,
		MPI_Datatype type) {
	long long int start_time;
	long long int end_time;
	log_->d_.bcast_num_++;
	start_time = timer_->Elapsed();

	int error = MPI_Bcast(buffer, data_count, type, 0,
			MPI_COMM_WORLD);

	end_time = timer_->Elapsed();
	log_->d_.bcast_time_ += end_time - start_time;
	log_->d_.bcast_time_max_ = std::max(end_time - start_time,
			log_->d_.bcast_time_max_);
	return error;
}

std::ofstream ParallelDFS::null_stream_;

// TODO: not sure how to use FLAGS. Is it transferrable?
//       Is it global variable? Either way not sure the purpose of the package.
std::ostream& ParallelDFS::D(int level, bool show_phase) {
//	return null_stream_;
	bool FLAGS_d = true;
	if (FLAGS_d == 0)
		return null_stream_;

	if (level <= FLAGS_d) {
		if (show_phase)
			lfs_ << std::setw(4) << "DFS" << ": ";
		return lfs_;
	} else
		return null_stream_;
}

} /* namespace lamp_search */
