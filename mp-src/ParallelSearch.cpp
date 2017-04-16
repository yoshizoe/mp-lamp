/*
 * ParallelSearch.cpp
 *
 *  Created on: Apr 16, 2017
 *      Author: yuu
 */

#include "ParallelSearch.h"


namespace lamp_search {

ParallelSearch::ParallelSearch() {
	// TODO Auto-generated constructor stub

}

ParallelSearch::~ParallelSearch() {
	// TODO Auto-generated destructor stub
}

void ParallelSearch::Search() {
	while (!dtd_.terminated_) {
		while (!dtd_.terminated_) {
			if (ProcessNode (n_)) {
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
				if (h_ == 0 && phase_ == 1)
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
		if (h_ == 0 && phase_ == 1)
			CheckCSThreshold();
		log_.d_.idle_time_ += timer_->Elapsed() - log_.idle_start_;
	}
}

void ParallelSearch::Probe() {
}

void ParallelSearch::Distribute() {
}

void ParallelSearch::Reject() {
}

void ParallelSearch::CheckCSThreshold() {
}

void ParallelSearch::Steal() {
}

} /* namespace lamp_search */
