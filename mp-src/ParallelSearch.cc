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
