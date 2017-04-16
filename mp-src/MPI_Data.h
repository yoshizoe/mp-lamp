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
	MPI_Data(int rank, int nu_proc, int n, bool n_is_ms, int w, int m, int l) :
			bsend_buffer_(NULL), mpiRank_(rank), nTotalProc_(nu_proc), granularity_(
					n), isGranularitySec_(n_is_ms), nRandStealTrials_(w), nRandStealCands_(
					m), lHypercubeEdge_(l), hypercubeDimension_(
					ComputeZ(nTotalProc_, lHypercubeEdge_)), rng_(mpiRank_), dst_p_(
					0, nTotalProc_ - 1), dst_m_(0, nRandStealCands_ - 1), rand_p_(rng_, dst_p_), rand_m_(
					rng_, dst_m_) {
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
	static int ComputeZ(int p, int l) {
		int z0 = 1;
		int zz = l;
		while (zz < p) {
			z0++;
			zz *= l;
		}
		return z0;
	}
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

};
}
#endif /* MP_SRC_MPI_DATA_H_ */
