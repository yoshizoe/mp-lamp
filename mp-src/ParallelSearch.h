/*
 * ParallelSearch.h
 *
 *  Created on: Apr 17, 2017
 *      Author: yuu
 */

#ifndef MP_SRC_PARALLELSEARCH_H_
#define MP_SRC_PARALLELSEARCH_H_
#include "MPI_Data.h"

namespace lamp_search {

class ParallelSearch {
public:
	ParallelSearch(MPI_Data& mpi_data, TreeSearchData* treesearch_data,
			Log* log, Timer* timer);
	virtual ~ParallelSearch();
//	void Search();

	void GetMinimalSupport(MPI_Data& mpi_data, TreeSearchData* treesearch_data,
			GetMinSupData* getminsup_data);
private:
	MPI_Data& mpi_data;
	TreeSearchData* treesearch_data;
	Log* log_;
	Timer* timer_;

	int phase_; // TODO: remove it.


	/**
	 * Methods used for ALL searches: Maybe they should be overrided by other methods.
	 *
	 */
	void Probe(MPI_Data& mpi_data, TreeSearchData* treesearch_data);
	void Distribute(MPI_Data& mpi_data, TreeSearchData* treesearch_data);
	void Give(MPI_Data& mpi_data, VariableLengthItemsetStack * st,
			int steal_num);
	void Deal(MPI_Data& mpi_data);
	void Reject(MPI_Data& mpi_data);
};

} /* namespace lamp_search */

#endif /* MP_SRC_PARALLELSEARCH_H_ */
