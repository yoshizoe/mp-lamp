/*
 * MPI_Data.h
 *
 *  Created on: Apr 16, 2017
 *      Author: yuu
 */

#ifndef MP_SRC_MPI_DATA_H_
#define MP_SRC_MPI_DATA_H_

#include "mpi.h"
#include <boost/random.hpp>
#include "FixedSizeStack.h"
// GLB variables
namespace lamp_search {

struct MPI_Data {
	MPI_Data(int buffer_size, int rank, int nu_proc, int n, bool n_is_ms, int w,
			int m, int l, int k_echo_tree_branch, DTD* dtd_) :
			dtd_(dtd_), mpiRank_(rank), nTotalProc_(nu_proc), granularity_(
					n), isGranularitySec_(n_is_ms), nRandStealTrials_(w), nRandStealCands_(
					m), lHypercubeEdge_(l), hypercubeDimension_(
					ComputeZ(nTotalProc_, lHypercubeEdge_)), rng_(mpiRank_), dst_p_(
					0, nTotalProc_ - 1), dst_m_(0, nRandStealCands_ - 1), rand_p_(
					rng_, dst_p_), rand_m_(rng_, dst_m_) {
		// Initializing temporary variables in parenths.
		processing_node_ = false;
		bsend_buffer_ = new int[buffer_size];

		int ret = MPI_Buffer_attach(bsend_buffer_, buffer_size * sizeof(int));
		if (ret != MPI_SUCCESS) {
			throw std::bad_alloc();
		}
		printf("bsend_buffer\n");

		prepareVictims();
		prepareLifelines();

		thieves_ = new FixedSizeStack(nTotalProc_);
		lifeline_thieves_ = new FixedSizeStack(
				hypercubeDimension_ + k_echo_tree_branch);
		lifelines_activated_ = new bool[nTotalProc_];
		accum_flag_ = new bool[k_echo_tree_branch];
		bcast_targets_ = new int[k_echo_tree_branch];

		for (int pi = 0; pi < nTotalProc_; pi++)
			lifelines_activated_[pi] = false;

		bcast_source_ = -1;
		echo_waiting_ = false;

		for (std::size_t i = 0; i < k_echo_tree_branch; i++) {
			bcast_targets_[i] = -1;
			accum_flag_[i] = false;
		}
		printf("mpidata\n");
	}
	~MPI_Data() {
		delete[] victims_;
		delete[] lifelines_;
		delete thieves_;
		delete lifeline_thieves_;
		delete[] lifelines_activated_;
		delete[] accum_flag_;
		delete[] bcast_targets_;
		int size;
		MPI_Buffer_detach(&bsend_buffer_, &size);
//		assert(size == FLAGS_bsend_buffer_size * sizeof(int));
		if (bsend_buffer_)
			delete bsend_buffer_;
	}
	void prepareVictims() {
		victims_ = new int[nRandStealCands_];
		if (nTotalProc_ > 1) {
			for (int pi = 0; pi < nRandStealCands_; pi++) {
				int r;
				while (true) {
					r = rand_p_();
					if (r != mpiRank_)
						break;
				}
				victims_[pi] = r;
			}
		}
		printf("victims\n");
	}

	void prepareLifelines() {
		lifelines_ = new int[hypercubeDimension_];
		for (int zi = 0; zi < hypercubeDimension_; zi++)
			lifelines_[zi] = -1;

		// lifeline initialization
		// cf. Saraswat et al. "Lifeline-Based Global Load Balancing", PPoPP 2008.
		int radix = 1;
		int lifeline_buddy_id = 0;
		for (int j = 0; j < hypercubeDimension_; j++) { // for each dimension
			int next_radix = radix * lHypercubeEdge_;

			for (int k = 1; k < lHypercubeEdge_; k++) { // loop over an edge of the hypercube
				int base = mpiRank_ - mpiRank_ % next_radix; // lowest in the current ring
				int index = (mpiRank_ + next_radix - radix) % next_radix;
				int lifeline_buddy = base + index; // previous in the current ring
				if (lifeline_buddy < nTotalProc_) {
					lifelines_[lifeline_buddy_id++] = lifeline_buddy;
					break;
				}
			}
			radix *= lHypercubeEdge_;
		}

		printf("topology\n");
	}

	static int ComputeZ(int p, int l) {
		int z0 = 1;
		int zz = l;
		while (zz < p) {
			z0++;
			zz *= l;
		}
		return z0;
	}
	DTD* dtd_;

	int * bsend_buffer_;
	int mpiRank_; // MPI Rank
	int nTotalProc_; // total proc number

	int granularity_; // granularity of tasks
	bool isGranularitySec_; // false: n_ is number of task, true: n_ is milli sec
	int nRandStealTrials_; // number of random steal trials
	int nRandStealCands_; // number of random steal candidates

	int lHypercubeEdge_; // power of lifeline graph (length the hypercube edge)
	int hypercubeDimension_; // dimension of lifeline (dimension of the hypercube)

	// TODO: should this be here?
	boost::mt19937 rng_; // use seed as rank
	boost::uniform_smallint<int> dst_p_;
	boost::uniform_smallint<int> dst_m_;
	boost::variate_generator<boost::mt19937&, boost::uniform_smallint<int> > rand_p_;
	boost::variate_generator<boost::mt19937&, boost::uniform_smallint<int> > rand_m_;

	// TODO: Initialized by mp_dfs. Need to put in here.
	int bcast_source_;
	int * bcast_targets_;
	bool echo_waiting_; // new! [2015-10-05 22:23]
	bool * accum_flag_;
	int * victims_; // proc id of random victims
	int * lifelines_; // proc id of lifeline buddies

	FixedSizeStack * thieves_; // max size == nu_proc_
	FixedSizeStack * lifeline_thieves_; // size == lifelines_ size + 3
	bool * lifelines_activated_;

	// Temporary variables
	bool processing_node_; // prevent termination while processing node

};
}
#endif /* MP_SRC_MPI_DATA_H_ */
