/*
 * ParallelSearch.h
 *
 *  Created on: Apr 17, 2017
 *      Author: yuu
 */

#ifndef MP_SRC_PARALLELPATTERNMINING_H_
#define MP_SRC_PARALLELPATTERNMINING_H_
#include "MPI_Data.h"

#include "ParallelDFS.h"

namespace lamp_search {

/**
 * TODO: This class is intended to be a parent class for pattern mining classes.
 *       For now it takes responsibility to run
 *       - GetMinimalSupport
 *       - GetTestablePatterns
 *       - GetSignificantPatterns
 *       These functions should be factored out as subclasses.
 */
class ParallelPatternMining: public ParallelDFS {
public:
	ParallelPatternMining(BinaryPatternMiningData* bpm_data,
			MPI_Data& mpi_data, TreeSearchData* treesearch_data,
			Log* log, Timer* timer, std::ostream& ofs);
	virtual ~ParallelPatternMining();

//	virtual void Search();
	void GetMinimalSupport(GetMinSupData* getminsup_data);
	void PreProcessRootNode(GetMinSupData* getminsup_data);
	void GetTestablePatterns(GetTestableData* gettestable_data);
	void GetSignificantPatterns(
			GetSignificantData* getsignificant_data);

protected:
	// TODO: How can we hide the dependency on those low level structures?
	//       Let's try to understand the semantics of how these methods are used, and factor out.
	/*
	 * Domain graph
	 */
	Database<uint64> * d_;
	LampGraph<uint64> * g_;
	VariableBitsetHelper<uint64> * bsh_;
	uint64 * sup_buf_, *child_sup_buf_; // TODO: sup_buf_ is only used in ProcessNode and PreProcessRootNode!

	/*
	 * Data structure
	 */
//	MPI_Data& mpi_data;
//	TreeSearchData* treesearch_data;
	GetMinSupData* getminsup_data;
	GetTestableData* gettestable_data;
	GetSignificantData* getsignificant_data;

//	Log* log_;
//	Timer* timer_;

	/**
	 * Methods used for ALL searches: Maybe they should be overrided by other methods.
	 *
	 */
//	bool Probe(MPI_Data& mpi_data, TreeSearchData* treesearch_data);
	virtual void ProbeExecute(TreeSearchData* treesearch_data,
			MPI_Status* probe_status, int probe_src, int probe_tag);
//	bool ProbeExecuteMINSUP(MPI_Data& mpi_data, TreeSearchData* treesearch_data,
//			MPI_Status& probe_status, int probe_src, int probe_tag);

//	void Distribute(MPI_Data& mpi_data, TreeSearchData* treesearch_data);
//	void Give(MPI_Data& mpi_data, VariableLengthItemsetStack * st,
//			int steal_num);
////	void Deal(MPI_Data& mpi_data);
//	void Reject(MPI_Data& mpi_data);
//	void Steal(MPI_Data& mpi_data);
	void ProcAfterProbe(); // DOMAINDEPENDENT
	void Check(); // DOMAINDEPENDENT
	bool ExpandNode(TreeSearchData*treesearch_data);
	std::vector<int> GetChildren(int core_i); // DOMAINDEPENDENT
	void PopNodeFromStack();
	bool TestAndPushNode(int new_item, int core_i);
	void ProcessNode(int sup_num, int* ppc_ext_buf);
	void CheckProbe(int& accum_period_counter_,
			long long int lap_time);
	bool CheckProcessNodeEnd(int n, bool n_is_ms, int processed,
			long long int start_time);

	//--------

	int phase_; // 1, 2, 3

	/**
	 * GetMinSup Functions
	 *
	 */
	// 0: count, 1: time warp flag, 2: empty flag, 3--: data
	void SendDTDAccumRequest();
	void RecvDTDAccumRequest(int src);
	// 0: count, 1: time warp flag, 2: empty flag, 3--: data
	void SendDTDAccumReply();
	void RecvDTDAccumReply(int src);
	void SendLambda(int lambda);
	void RecvLambda(int src);
	void CheckCSThreshold();
	bool ExceedCsThr() const; // getMinSup
	int NextLambdaThr() const; // getMinSup
//	int NextLambdaThr(GetMinSupData* getminsup_data) const; // getMinSup
	void IncCsAccum(int sup_num); // getMinSup

	// TODO: These functions should be factored in Get
	/**
	 * Methods For GetSignificant
	 */
	void SendResultRequest();
	void RecvResultRequest(int src);
	void SendResultReply();
	void RecvResultReply(int src, MPI_Status status);
	void ExtractSignificantSet();

// insert pointer into significant_map_ (do not sort the stack itself)
//	void SortSignificantSets();

	/**
	 * Utils
	 *
	 */
	void CheckInit();
	void CheckInitTestable();

	/**
	 * Statistics
	 */
	int expand_num_;
	int closed_set_num_;
};

} /* namespace lamp_search */

#endif /* MP_SRC_PARALLELPATTERNMINING_H_ */
