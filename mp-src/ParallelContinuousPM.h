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
			MPI_Data& mpi_data, TreeSearchData* treesearch_data,
			double alpha, Log* log, Timer* timer, std::ostream& ofs);
	virtual ~ParallelContinuousPM();

//	virtual void Search();
	void GetMinimalSupport();
	void GetDiscretizedMinimalSupport(double freqRatio);
//	void PreProcessRootNode(GetMinSupData* getminsup_data);
	void GetTestablePatterns(GetTestableData* gettestable_data);
	void GetSignificantPatterns(
			GetContSignificantData* getsignificant_data);

	void GetTopKPatterns(int k);

	double GetThreFreq() const {
		return thre_freq_;
	}
	double GetThrePmin() const {
		return thre_pmin_;
	}
	int NumberOfTestablePatterns() const;

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
	GetMinSupData* getminsup_data;
	GetTestableData* gettestable_data;
	GetContSignificantData* getsignificant_data;

//	Log* log_;
//	Timer* timer_;

	/**
	 * Methods used for ALL searches: Maybe they should be overrided by other methods.
	 *
	 */
//	bool Probe(MPI_Data& mpi_data, TreeSearchData* treesearch_data);
	virtual void ProbeExecute(
			TreeSearchData* treesearch_data, MPI_Status* probe_status,
			int probe_src, int probe_tag);
//	bool ProbeExecuteMINSUP(MPI_Data& mpi_data, TreeSearchData* treesearch_data,
//			MPI_Status& probe_status, int probe_src, int probe_tag);

	void ProcAfterProbe(); // DOMAINDEPENDENT
	void Check(); // DOMAINDEPENDENT
	bool ExpandNode(
			TreeSearchData*treesearch_data);
	std::vector<int> GetChildren(std::vector<int> items); // DOMAINDEPENDENT
	bool PopNodeFromStack();
	bool TestAndPushNode(int new_item);
	void ProcessNode(double freq, int* ppc_ext_buf);
	void CheckProbe(int& accum_period_counter_,
			long long int lap_time);
	bool CheckProcessNodeEnd(int n, bool n_is_ms, int processed,
			long long int start_time);

	//--------

	int phase_; // 1, 2, 3

	bool HasJobToDo();

	/**
	 * Methods for Maintaining threshold value
	 */
	// 0: count, 1: time warp flag, 2: empty flag, 3--: data
	void SendMinPValueRequest();
	void RecvMinPValueRequest(int src);
	// 0: count, 1: time warp flag, 2: empty flag, 3--: data
	void SendMinPValueReply();
	void RecvMinPValueReply(int src,
			MPI_Status* probe_status);
	void CalculateThreshold();
	void SendNewSigLevel(double sig_level);
	void RecvNewSigLevel(int src);
	std::vector<double> frequencies;
	std::vector<double> topKFrequencies;

	double alpha_;
	double thre_freq_; // Threshold for itemset-set C.
	double thre_pmin_; // Threshold for itemset-set T (testable pattern)

//	std::vector<std::pair<double, double>> freq_pmin; // only for rank-0.
	bool freq_received;
//	int prev_freq_pmin_size_; // only for rank-0.
	// TODO: These functions should be factored in Get
	/**
	 * Methods For GetSignificant
	 */
	void SendResultRequest();
	void RecvResultRequest(int src);
	void SendResultReply();
	void RecvResultReply(int src,
			MPI_Status status);
//	bool AccumCountReady(MPI_Data& mpi_data) const;
	void ExtractSignificantSet();

	void PrintItemset(int* itembuf, std::vector<Feature> freqs);

// insert pointer into significant_map_ (do not sort the stack itself)
//	void SortSignificantSets();

	/**
	 * Linear Space Continuous Pattern Mining
	 *
	 */
	std::vector<std::pair<double, double>> InitializeThresholdTable(
			 double ratio, int size, double alpha);
	// 0: count, 1: time warp flag, 2: empty flag, 3--: data
	void SendDTDAccumRequest();
	void RecvDTDAccumRequest(int src);
	// 0: count, 1: time warp flag, 2: empty flag, 3--: data
	void SendDTDAccumReply();
	void RecvDTDAccumReply(int src);
	int GetDiscretizedFrequency(double freq) const;
	void CheckCSThreshold();
	bool ExceedCsThr() const;
	int NextLambdaThr() const;
	void IncCsAccum(int sup_num);
	void SendLambda(int lambda);
	void RecvLambda(int src);
	std::vector<std::pair<double, double>> thresholds;
//	std::vector<int> count;

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
