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
	ParallelPatternMining(Database<uint64> * d_, LampGraph<uint64> * g_,
			VariableBitsetHelper<uint64> * bsh_, MPI_Data& mpi_data,
			TreeSearchData* treesearch_data, Log* log, Timer* timer);
	virtual ~ParallelPatternMining();

//	virtual void Search();
	void GetMinimalSupport(GetMinSupData* getminsup_data);
	void PreProcessRootNode(GetMinSupData* getminsup_data);
	void GetTestablePatterns(GetTestableData* gettestable_data);
	void GetSignificantPatterns(MPI_Data& mpi_data,
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
	virtual void ProbeExecute(MPI_Data& mpi_data, TreeSearchData* treesearch_data,
			MPI_Status* probe_status, int probe_src, int probe_tag);
//	bool ProbeExecuteMINSUP(MPI_Data& mpi_data, TreeSearchData* treesearch_data,
//			MPI_Status& probe_status, int probe_src, int probe_tag);

//	void Distribute(MPI_Data& mpi_data, TreeSearchData* treesearch_data);
//	void Give(MPI_Data& mpi_data, VariableLengthItemsetStack * st,
//			int steal_num);
////	void Deal(MPI_Data& mpi_data);
//	void Reject(MPI_Data& mpi_data);
//	void Steal(MPI_Data& mpi_data);
	void ProcAfterProbe();
	void Check(MPI_Data& mpi_data);
	bool ProcessNode(MPI_Data& mpi_data, TreeSearchData*treesearch_data);
	std::vector<int> GetChildren(bool is_root_node, int coreindex);
	void CheckProbe(int accum_period_counter_, long long int lap_time);
	bool CheckProcessNodeEnd(int n, bool n_is_ms, int processed,
			long long int start_time);
	/**
	 * Methods for Probe and Send.
	 *
	 */

//	// 0: count, 1: time warp flag, 2: empty flag
//	void SendDTDRequest(MPI_Data& mpi_data);
//	void RecvDTDRequest(MPI_Data& mpi_data, TreeSearchData* treesearch_data,
//			int src);
//
//	bool DTDReplyReady(MPI_Data& mpi_data) const;
//	void DTDCheck(MPI_Data& mpi_data);
//
//	// 0: count, 1: time warp flag, 2: empty flag
//	void SendDTDReply(MPI_Data& mpi_data, TreeSearchData* treesearch_data);
//	void RecvDTDReply(MPI_Data& mpi_data, TreeSearchData* treesearch_data,
//			int src);

//	bool DTDAccumReady(MPI_Data& mpi_data) const;

//	void SendBcastFinish(MPI_Data& mpi_data);
//	void RecvBcastFinish(MPI_Data& mpi_data, int src);

	//--------

	int phase_; // 1, 2, 3

	//--------
	// basic

	// send recv functions
//	void SendRequest(MPI_Data& mpi_data, int dst, int is_lifeline); // for random thieves, is_lifeline = -1
//	void RecvRequest(MPI_Data& mpi_data, TreeSearchData* treesearch_data,
//			int src);
//
//	// 0: time zone, 1: is_lifeline
//	void SendReject(MPI_Data& mpi_data, int dst);
//	void RecvReject(MPI_Data& mpi_data, int src);
//
//	// 1: time zone
//	void SendGive(MPI_Data& mpi_data, VariableLengthItemsetStack * st, int dst,
//			int is_lifeline);
//
//	// sets lifelines_activated_ = false
//	// lifelines_activated_ becomes false only in this case (reject does NOT)
//	void RecvGive(MPI_Data& mpi_data, TreeSearchData* treesearch_data, int src,
//			MPI_Status status);


	/**
	 * GetMinSup Functions
	 *
	 */
	// 0: count, 1: time warp flag, 2: empty flag, 3--: data
	void SendDTDAccumRequest(MPI_Data& mpi_data);
	void RecvDTDAccumRequest(MPI_Data& mpi_data, int src);
	// 0: count, 1: time warp flag, 2: empty flag, 3--: data
	void SendDTDAccumReply(MPI_Data& mpi_data);
	void RecvDTDAccumReply(MPI_Data& mpi_data, int src);
	void SendLambda(MPI_Data& mpi_data, int lambda);
	void RecvLambda(MPI_Data& mpi_data, int src);
	void CheckCSThreshold(MPI_Data& mpi_data);
	bool ExceedCsThr() const; // getMinSup
	int NextLambdaThr() const; // getMinSup
//	int NextLambdaThr(GetMinSupData* getminsup_data) const; // getMinSup
	void IncCsAccum(int sup_num); // getMinSup


	// TODO: These functions should be factored in Get
	/**
	 * Methods For GetSignificant
	 */
	void SendResultRequest(MPI_Data& mpi_data);
	void RecvResultRequest(MPI_Data& mpi_data, int src);
	void SendResultReply(MPI_Data& mpi_data);
	void RecvResultReply(MPI_Data& mpi_data, int src, MPI_Status status);
	bool AccumCountReady(MPI_Data& mpi_data) const;
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
