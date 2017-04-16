/*
 * ParallelSearch.h
 *
 *  Created on: Apr 16, 2017
 *      Author: yuu
 */

#ifndef MP_SRC_PARALLELSEARCH_H_
#define MP_SRC_PARALLELSEARCH_H_

#include "DTD.h"
#include "Log.h"
#include "../src/variable_length_itemset.h"
#include "../src/timer.h"
#include "FixedSizeStack.h"
#include <boost/random.hpp>

namespace lamp_search {

class ParallelSearch {
public:
	ParallelSearch(MPI_Data& mpi_data); // Put variables required for ALL parallel searches.
	~ParallelSearch();
	void SetupTopology();
	void InitTreeRequest();
	void Init(); // Very basic value assignments.
	void SetTreeRequest();
	void CheckAndClear();
	void ClearTasks();

	void setPhase(int phase); // TODO
	virtual ~ParallelSearch();
	void LCM(LCM_Data& lcm_data); // TODO
	void LAMP(LAMP_Data& lamp_data); // TODO
	void SignificantPM(Significant_Data& sig_data); // TODO

protected:
	bool ProcessNode(int granurality);
	bool CheckProcessNodeEnd(int n, bool n_is_ms, int processed,
			long long int start_time);
	void Probe();
	void Distribute();
    void Reject(); // distribute finished, reject remaining requests
    void CheckCSThreshold();
    void Steal(); // request steal
    void Give(VariableLengthItemsetStack * st, int steal_num);

	static int ComputeZ(int p, int l); // topology


    int phase_; // TODO: this really doesn't make sense for novices.

    // TODO: Packing fields required for all parallel searches would be easier
    // These should be useful for binary pattern mining
    /**
     * MPI variables
     */
    MPI_Data mpi_data;
//    int mpiRank_;
//    int nTotalProc;
//	int * bsend_buffer_;
	/**
	 * Hypercube topology TODO: This should be put in something like ParallelTopology class.
	 */
//	int * victims_; // proc id of random victims
//	int * lifelines_; // proc id of lifeline buddies
//	FixedSizeStack * thieves_; // max size == nu_proc_
//	FixedSizeStack * lifeline_thieves_; // size == lifelines_ size + 3
//	int nRandStealTrials; // number of random steal trials
//	int nRandStealCands; // number of random steal candidates
//	int lHypercubeEdge; // power of lifeline graph (length the hypercube edge)
//	int hypercubeDimension; // dimension of lifeline (dimension of the hypercube)

    // Required for algortihms
    VariableLengthItemsetStack* node_stack_;
    VariableLengthItemsetStack* give_stack_;
	bool * lifelines_activated_; // not even sure what the hell it is.
	DTD& dtd_;
	int granularity; // granularity of tasks
	bool isGranularitySec; // false: n_ is number of task, true: n_ is milli sec
	static const int k_echo_tree_branch; // Used for all tasks?
	int bcast_source_;
	int * bcast_targets_; // Used for all tasks?

	/**
	 * Stealing
	 */
	// set this in steal, reset this in RecvGive and RecvReject
	// small difference from x10 implementation
	bool waiting_;
	StealState stealer_; // TODO: not sure what it does.

    // Used in all parallel algorithms
    bool echo_waiting_; // TODO
	bool * accum_flag_;
	bool last_bcast_was_dtd_;

    // Used in parallel searches
    bool processing_node_; // prevent terminating while processing.
	int itemset_buf_[VariableLengthItemsetStack::kMaxItemsPerSet];

    // LCM and LAMP (binary)
    int lambda_;
    int lambda_max_;
	uint64 * sup_buf_, *child_sup_buf_;
	// TODO: it is terrible that the main code is dependent on this deep. Refactor.
	Database<uint64> * d_;
	LampGraph<uint64> * g_;
	VariableBitsetHelper<uint64> * bsh_; // bitset helper

	/**
	 *  For LAMP (frequent pattern mining for significant pattern mining)
	 */
	double sig_level_; // initially set to 1. set to 0.05 (FLAGS_a) / cs_thr_[global_sup_thr_-1]
	double * pmin_thr_; // pmin_thr[sup] == tbl.PMin(sup), maybe redundant
	long long int * cs_thr_; // cs_thr[sup] shows closed set num threshold
	VariableLengthItemsetStack * freq_stack_; // For frequent pattern mining
	std::multimap<double, int *> freq_map_; // Save p-value and items For FPM for Significant Pattern mining
	VariableLengthItemsetStack * significant_stack_;
	std::set<SignificantSetResult, sigset_compare> significant_set_;
	double final_sig_level_;

	// LCM LAMP method
	bool ExceedCsThr() const;
	int NextLambdaThr() const;
	void IncCsAccum(int sup_num);
	double GetInterimSigLevel(int lambda) const;

	/**
	 * DTD variables: pretty sure these should be factored into DTD.h.
	 */
	// 0: count, 1: time warp, 2: empty flag, 3--: array
	long long int * dtd_accum_array_base_; // int array of [-3..lambda_max_] (size lambda_max_+4)
	// -3: count, -2: time warp, -1: empty flag, 0--: array
	long long int * accum_array_; // int array of [0...lambda_max_] (size lambda_max_+1)

	// 0: count, 1: time warp, 2: empty flag, 3--: array
	long long int * dtd_accum_recv_base_; // int array of [-3..lambda_max_] (size lambda_max_+4)
	// -3: count, -2: time warp, -1: empty flag, 0--: array
	long long int * accum_recv_; // int array of [0...lambda_max_] (size lambda_max_+1)

	/*
	 * Random seeds
	 */
	boost::mt19937 rng_; // use seed as rank
	boost::uniform_smallint<int> dst_p_;
	boost::uniform_smallint<int> dst_m_;
	boost::variate_generator<boost::mt19937&, boost::uniform_smallint<int> > rand_p_;
	boost::variate_generator<boost::mt19937&, boost::uniform_smallint<int> > rand_m_;

    /**
     * Probe Methods
     *
     */
	// return flag. if (flag), it is ready to receive
	bool IsLeaf() const;
	int CallIprobe(MPI_Status * status, int * count, int * src);
	int CallRecv(void * buffer, int count, MPI_Datatype type, int src, int tag,
			MPI_Status * status);
	int CallBsend(void * buffer, int count_int, MPI_Datatype type, int dest,
			int tag);
	int CallBcast(void * buffer, int data_count, MPI_Datatype type);

	// 0: count, 1: time warp flag, 2: empty flag
	void SendDTDRequest();
	void RecvDTDRequest(int src);

	bool DTDReplyReady() const;
	void DTDCheck();

	// 0: count, 1: time warp flag, 2: empty flag
	void SendDTDReply();
	void RecvDTDReply(int src);

	bool DTDAccumReady() const;

	// 0: count, 1: time warp flag, 2: empty flag, 3--: data
	void SendDTDAccumRequest();
	void RecvDTDAccumRequest(int src);

	// 0: count, 1: time warp flag, 2: empty flag, 3--: data
	void SendDTDAccumReply();
	void RecvDTDAccumReply(int src);

	void SendBcastFinish();
	void RecvBcastFinish(int src);

	//--------
	// basic

	// send recv functions
	void SendRequest(int dst, int is_lifeline); // for random thieves, is_lifeline = -1
	void RecvRequest(int src);
	// 0: time zone, 1: is_lifeline
	void SendReject(int dst);
	void RecvReject(int src);
	// 1: time zone
	void SendGive(VariableLengthItemsetStack * st, int dst, int is_lifeline);
	// sets lifelines_activated_ = false
	// lifelines_activated_ becomes false only in this case (reject does NOT)
	void RecvGive(int src, MPI_Status status);

	void SendLambda(int lambda);
	void RecvLambda(int src);

	//--------
	// third phase
	bool AccumCountReady() const;
	void SendResultRequest();
	void RecvResultRequest(int src);
	void SendResultReply();
	void RecvResultReply(int src, MPI_Status status);
	void ExtractSignificantSet();
	// insert pointer into significant_map_ (do not sort the stack itself)
	void SortSignificantSets();

    Log& log_;
    Timer* timer_;

    // TODO: Statistics.
	long long int total_expand_num_;
	long long int expand_num_;
	long long int closed_set_num_;
};

} /* namespace lamp_search */

#endif /* MP_SRC_PARALLELSEARCH_H_ */
