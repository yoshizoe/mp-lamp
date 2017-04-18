/*
 * ParallelSearch.cc
 *
 *  Created on: Apr 18, 2017
 *      Author: yuu
 */

#include "ParallelSearch.h"

namespace lamp_search {

const int ParallelSearch::k_echo_tree_branch = 3;

ParallelSearch::ParallelSearch(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, Log* log, Timer* timer) :
		mpi_data(mpi_data), treesearch_data(treesearch_data), log_(log), timer_(
				timer) {

}

ParallelSearch::~ParallelSearch() {
	// TODO Auto-generated destructor stub
}

void ParallelSearch::Search() {
	DBG(D(1) << "MainLoop" << std::endl
			;);
	while (!mpi_data.dtd_->terminated_) {
		while (!mpi_data.dtd_->terminated_) {
			if (ProcessNode(mpi_data, treesearch_data)) {
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
			log_->d_.idle_time_ += timer_->Elapsed() - log_->idle_start_;
			break;
		}

		Probe(mpi_data, treesearch_data);
		if (mpi_data.dtd_->terminated_) {
			log_->d_.idle_time_ += timer_->Elapsed() - log_->idle_start_;
			break;
		}
		Check(mpi_data); // Implement CheckCSThreshold().

		log_->d_.idle_time_ += timer_->Elapsed() - log_->idle_start_;
	}
}


//==============================================================================

bool ParallelSearch::IsLeafInTopology(MPI_Data& mpi_data) const {
	for (int i = 0; i < k_echo_tree_branch; i++)
		if (mpi_data.bcast_targets_[i] >= 0)
			return false;
	return true;
}

int ParallelSearch::CallIprobe(MPI_Status * status, int * src,
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

int ParallelSearch::CallRecv(void * buffer, int data_count,
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

int ParallelSearch::CallBsend(void * buffer, int data_count,
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

int ParallelSearch::CallBcast(void * buffer, int data_count,
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


std::ofstream ParallelSearch::null_stream_;

std::ostream& ParallelSearch::D(int level, bool show_phase) {
	return null_stream_;
//	bool FLAGS_d = true;
//	if (FLAGS_d == 0)
//		return null_stream_;
//
//	if (level <= FLAGS_d) {
//		if (show_phase)
//			lfs_ << std::setw(4) << phase_ << ": ";
//		return lfs_;
//	} else
//		return null_stream_;
}

} /* namespace lamp_search */
