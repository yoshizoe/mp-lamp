/*
 * ParallelSearch.cpp
 *
 *  Created on: Apr 17, 2017
 *      Author: yuu
 */

#include "ParallelSearch.h"

namespace lamp_search {

ParallelSearch::ParallelSearch(MPI_Data& mpi_data__,
		TreeSearchData* treesearch_data__, Log* log__, Timer* timer__) :
		mpi_data(mpi_data__), treesearch_data(treesearch_data__), log_(log__), timer_(
				timer__) {
}

ParallelSearch::~ParallelSearch() {
	// TODO Auto-generated destructor stub
}

void ParallelSearch::GetMinimalSupport(MPI_Data& mpi_data,
		TreeSearchData* treesearch_data, GetMinSupData* getminsup_data) {
	phase_ = 1;
	DBG(D(1) << "MainLoop" << std::endl
	;);
	while (!mpi_data.dtd_->terminated_) {
		while (!mpi_data.dtd_->terminated_) {
			if (ProcessNode(mpi_data, treesearch_data, getminsup_data,
					(GetTestableData*) NULL)) {
				log_.d_.node_stack_max_itm_ =
						std::max(log_.d_.node_stack_max_itm_,
								(long long int) (treesearch_data->node_stack_->NuItemset()));
				log_.d_.node_stack_max_cap_ =
						std::max(log_.d_.node_stack_max_cap_,
								(long long int) (treesearch_data->node_stack_->UsedCapacity()));
				Probe(mpi_data, treesearch_data);
				if (mpi_data.dtd_->terminated_)
					break;
				Distribute(mpi_data, treesearch_data);
				Reject(mpi_data); // distribute finished, reject remaining requests
				if (mpi_data.mpiRank_ == 0)
					CheckCSThreshold(mpi_data);
			} else
				break;
		}
		if (mpi_data.dtd_->terminated_)
			break;

		log_.idle_start_ = timer_->Elapsed();
		Reject(mpi_data); // node_stack_ empty. reject requests
		Steal(mpi_data); // request steal
		if (mpi_data.dtd_->terminated_) {
			log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
			break;
		}

		Probe(mpi_data, treesearch_data);
		if (mpi_data.dtd_->terminated_) {
			log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
			break;
		}
		if (mpi_data.mpiRank_ == 0)
			CheckCSThreshold(mpi_data);
		log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
	}

}

} /* namespace lamp_search */
