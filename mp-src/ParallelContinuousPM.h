/*
 * ParallelSearch.h
 *
 *  Created on: Apr 17, 2017
 *      Author: yuu
 */

#ifndef MP_SRC_PARALLELCONTINUOUSPM_H_
#define MP_SRC_PARALLELCONTINUOUSPM_H_
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
class ParallelContinuousPM: public ParallelDFS {
	typedef double Feature;
	typedef int Class;
public:
	ParallelContinuousPM(ContinuousPatternMiningData* bpm_data,
			MPI_Data& mpi_data, TreeSearchData* treesearch_data, double alpha, Log* log,
			Timer* timer, std::ostream& ofs);
	virtual ~ParallelContinuousPM();

//	virtual void Search();
//	void GetMinimalSupport(GetMinSupData* getminsup_data);
//	void PreProcessRootNode(GetMinSupData* getminsup_data);
	void GetTestablePatterns(GetTestableData* gettestable_data);
	void GetSignificantPatterns(MPI_Data& mpi_data,
			GetContSignificantData* getsignificant_data);

protected:
	// TODO: How can we hide the dependency on those low level structures?
	//       Let's try to understand the semantics of how these methods are used, and factor out.
	/*
	 * Domain graph
	 */
	ContDatabase* d_;
	// TODO: sup_buf_ is only used in ProcessNode and PreProcessRootNode!
	std::vector<Feature> sup_buf_;
	std::vector<Feature> child_sup_buf_;

	/*
	 * Data structure
	 */
//	MPI_Data& mpi_data;
//	TreeSearchData* treesearch_data;
//	GetMinSupData* getminsup_data;
	GetTestableData* gettestable_data;
	GetContSignificantData* getsignificant_data;

//	Log* log_;
//	Timer* timer_;

	/**
	 * Methods used for ALL searches: Maybe they should be overrided by other methods.
	 *
	 */
//	bool Probe(MPI_Data& mpi_data, TreeSearchData* treesearch_data);
	virtual void ProbeExecute(MPI_Data& mpi_data,
			TreeSearchData* treesearch_data, MPI_Status* probe_status,
			int probe_src, int probe_tag);
//	bool ProbeExecuteMINSUP(MPI_Data& mpi_data, TreeSearchData* treesearch_data,
//			MPI_Status& probe_status, int probe_src, int probe_tag);

	void ProcAfterProbe(); // DOMAINDEPENDENT
	void Check(MPI_Data& mpi_data); // DOMAINDEPENDENT
	bool ExpandNode(MPI_Data& mpi_data, TreeSearchData*treesearch_data);
	std::vector<int> GetChildren(std::vector<int> items); // DOMAINDEPENDENT
	void PopNodeFromStack();
	bool TestAndPushNode(int new_item);
	void ProcessNode(double freq, int* ppc_ext_buf);
	void CheckProbe(int& accum_period_counter_, long long int lap_time);
	bool CheckProcessNodeEnd(int n, bool n_is_ms, int processed,
			long long int start_time);

	//--------

	int phase_; // 1, 2, 3

	/**
	 * Methods for Maintaining threshold value
	 */
	// 0: count, 1: time warp flag, 2: empty flag, 3--: data
	void SendMinPValueRequest(MPI_Data& mpi_data);
	void RecvMinPValueRequest(MPI_Data& mpi_data, int src);
	// 0: count, 1: time warp flag, 2: empty flag, 3--: data
	void SendMinPValueReply(MPI_Data& mpi_data);
	void RecvMinPValueReply(MPI_Data& mpi_data, int src, MPI_Status* probe_status);
	void CalculateThreshold();
	void SendNewSigLevel(double sig_level);
	void RecvNewSigLevel(int src);
	std::vector<double> freqs_stack_;
	double alpha_;
	double thre_freq_;
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

	void PrintItemset(int* itembuf, std::vector<Feature> freqs);

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
