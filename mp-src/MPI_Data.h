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
#include "StealState.h"
#include "Log.h"
#include "DTD.h"
#include "SignificantSetResults.h"
#include "../src/variable_length_itemset.h"
#include "../src/utils.h"
#include "../src/lamp_graph.h"
#include "../src/database.h"
#include "../src/variable_bitset_array.h"
#include "../src/timer.h"
#include "../src/contdatabase.h"

// GLB variables
namespace lamp_search {

struct TreeSearchData {
	TreeSearchData(VariableLengthItemsetStack * nstack,
			VariableLengthItemsetStack * gstack, StealState* stealer_,
			int* itemset_buf_, Log *log_, Timer * timer_) :
			node_stack_(nstack), give_stack_(gstack), stealer_(
					stealer_), itemset_buf_(itemset_buf_), log_(log_), timer_(
					timer_) {
	}
	// Search fields
	VariableLengthItemsetStack * node_stack_;
	VariableLengthItemsetStack * give_stack_;
	StealState* stealer_;
	int* itemset_buf_; // Itemset = state.

//	struct Node {
//		int* itemset_buf_;
//		int getNumberOfItems() {
//			return ((-1) * itemset_buf_[0] - 1);
//		}
//		int getSupport() {
//			return itemset_buf_[1];
//		}
//		std::vector<int> getItems() {
//			// TODO: move?
//			std::vector<int> items(itemset_buf_ + 2,
//					itemset_buf_ + 2 + getNumberOfItems());
//			return items;
//		}
//	};

	// Utils // TODO: put in mpi_data
	Log *log_;
	Timer * timer_;

	// TODO: these should be in Parallel Pattern Mining...
	// Domain fields. TODO: Way to wacky to be dependent on these...
//	uint64 * sup_buf_, *child_sup_buf_;
//	Database<uint64> * d_;
//	LampGraph<uint64> * g_;
//	VariableBitsetHelper<uint64> * bsh_; // bitset helper
};

struct ContinuousPatternMiningData {
	ContinuousPatternMiningData(ContDatabase* d_) :
			d_(d_) {
	}
	ContDatabase * d_;
};

struct BinaryPatternMiningData {
	BinaryPatternMiningData(Database<uint64>* d_,
			VariableBitsetHelper<uint64>* bsh_, uint64* sup_buf_,
			uint64* child_sup_buf_) :
			d_(d_), bsh_(bsh_), sup_buf_(sup_buf_), child_sup_buf_(
					child_sup_buf_) {
	}
	Database<uint64> * d_;
	VariableBitsetHelper<uint64> * bsh_; // bitset helper
	uint64 * sup_buf_, *child_sup_buf_; // TODO: these two should be local
};

struct GetMinSupData {
	GetMinSupData(int lambda_max, int lambda, long long int* cs_thr,
			long long int * dtd_accum_array_base,
			long long int * accum_array,
			long long int * dtd_accum_recv_base,
			long long int * accum_recv) :
			lambda_max_(lambda_max), lambda_(lambda), cs_thr_(cs_thr), dtd_accum_array_base_(
					dtd_accum_array_base), accum_array_(accum_array), dtd_accum_recv_base_(
					dtd_accum_recv_base), accum_recv_(accum_recv) {
	}
	int lambda_max_;
	int lambda_;
	long long int * cs_thr_;
	// 0: count, 1: time warp, 2: empty flag, 3--: array
	long long int * dtd_accum_array_base_; // int array of [-3..lambda_max_] (size lambda_max_+4)
	// -3: count, -2: time warp, -1: empty flag, 0--: array
	long long int * accum_array_; // int array of [0...lambda_max_] (size lambda_max_+1)

	// 0: count, 1: time warp, 2: empty flag, 3--: array
	long long int * dtd_accum_recv_base_; // int array of [-3..lambda_max_] (size lambda_max_+4)
	// -3: count, -2: time warp, -1: empty flag, 0--: array
	long long int * accum_recv_; // int array of [0...lambda_max_] (size lambda_max_+1)
};

struct GetTestableData {
	GetTestableData(int lambda_max_minus_one,
			VariableLengthItemsetStack * freq_stack,
			std::multimap<double, int *>* freq_map, double sig_level) :
			freqThreshold_(lambda_max_minus_one), freq_stack_(
					freq_stack), freq_map_(freq_map), sig_level_(
					sig_level) {
	}
	int freqThreshold_;
	// Retrun variables. Used for GetSignificantPatterns.
	VariableLengthItemsetStack * freq_stack_; // record freq itemsets
	std::multimap<double, int *>* freq_map_; // record (pval, *itemsets)
	double sig_level_;
//	VariableLengthItemsetStack * significant_stack_; // TODO:
};

struct GetSignificantData {
	GetSignificantData(VariableLengthItemsetStack * freq_stack_,
			std::multimap<double, int *>* freq_map_,
			double final_sig_level_,
			VariableLengthItemsetStack * significant_stack_,
			std::set<SignificantSetResult, sigset_compare>* significant_set_) :
			freq_stack_(freq_stack_), freq_map_(freq_map_), final_sig_level_(
					final_sig_level_), significant_stack_(
					significant_stack_), significant_set_(
					significant_set_) {
	}
	VariableLengthItemsetStack * freq_stack_; // record freq itemsets
	std::multimap<double, int *>* freq_map_; // record (pval, *itemsets)
	double final_sig_level_;
	VariableLengthItemsetStack * significant_stack_;
	std::set<SignificantSetResult, sigset_compare>* significant_set_;
};

// TODO: Duplicate. Refactor
struct GetContSignificantData {
	GetContSignificantData(VariableLengthItemsetStack * freq_stack_,
			std::multimap<double, int *>* freq_map_,
			double final_sig_level_,
			VariableLengthItemsetStack * significant_stack_,
			std::set<ContSignificantSetResult, cont_sigset_compare>* significant_set_) :
			freq_stack_(freq_stack_), freq_map_(freq_map_), final_sig_level_(
					final_sig_level_), significant_stack_(
					significant_stack_), significant_set_(
					significant_set_) {
	}
	VariableLengthItemsetStack * freq_stack_; // record freq itemsets
	std::multimap<double, int *>* freq_map_; // record (pval, *itemsets)
	double final_sig_level_;
	VariableLengthItemsetStack * significant_stack_;
	std::set<ContSignificantSetResult, cont_sigset_compare>* significant_set_;
};

struct MPI_Data {
	MPI_Data(int buffer_size, int rank, int nu_proc, int n,
			bool n_is_ms, int w, int m, int l, int k_echo_tree_branch,
			DTD* dtd_) :
			dtd_(dtd_), mpiRank_(rank), nTotalProc_(nu_proc), granularity_(
					n), isGranularitySec_(n_is_ms), nRandStealTrials_(
					w), nRandStealCands_(m), lHypercubeEdge_(l), hypercubeDimension_(
					ComputeZ(nTotalProc_, lHypercubeEdge_)), rng_(
					mpiRank_), dst_p_(0, nTotalProc_ - 1), dst_m_(0,
					nRandStealCands_ - 1), rand_p_(rng_, dst_p_), rand_m_(
					rng_, dst_m_), echo_waiting_(false), waiting_(
					false) {
		// Initializing temporary variables in parenths.
		processing_node_ = false;
		bsend_buffer_ = new int[buffer_size];

		int ret = MPI_Buffer_attach(bsend_buffer_,
				buffer_size * sizeof(int));
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
				int index = (mpiRank_ + next_radix - radix)
						% next_radix;
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
	boost::variate_generator<boost::mt19937&,
			boost::uniform_smallint<int> > rand_p_;
	boost::variate_generator<boost::mt19937&,
			boost::uniform_smallint<int> > rand_m_;

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
	bool waiting_;
};
}
#endif /* MP_SRC_MPI_DATA_H_ */
