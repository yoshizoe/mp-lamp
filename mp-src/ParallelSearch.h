/*
 * ParallelSearch.h
 *
 *  Created on: Apr 18, 2017
 *      Author: yuu
 */

#ifndef MP_SRC_PARALLELSEARCH_H_
#define MP_SRC_PARALLELSEARCH_H_

#include "../src/variable_length_itemset.h"
#include "MPI_Data.h"
#include "mpi_tag.h"

namespace lamp_search {

/*
 * TODO: Should Extract more fields here.
 */
class ParallelSearch {
public:
	ParallelSearch(MPI_Data& mpi_data, TreeSearchData* treesearch_data, Log* log,
			Timer* timer);
	virtual ~ParallelSearch();
	virtual void Search();

protected:
	virtual bool Probe(MPI_Data& mpi_data, TreeSearchData* treesearch_data) = 0;
	virtual void Distribute(MPI_Data& mpi_data, TreeSearchData* treesearch_data) = 0;
	virtual void Give(MPI_Data& mpi_data, VariableLengthItemsetStack * st,
			int steal_num) = 0;
	virtual void Reject(MPI_Data& mpi_data) = 0;
	virtual void Steal(MPI_Data& mpi_data) = 0;
	virtual void Check(MPI_Data& mpi_data) = 0;
	virtual bool ProcessNode(MPI_Data& mpi_data, TreeSearchData*treesearch_data) = 0;

	MPI_Data& mpi_data;
	TreeSearchData* treesearch_data;
	static const int k_echo_tree_branch;

	/**
	 * Utility
	 */
	Log* log_;
	Timer* timer_;

	static std::ofstream null_stream_;
	std::ostream& D(int level, bool show_phase = true);
};

} /* namespace lamp_search */

#endif /* MP_SRC_PARALLELSEARCH_H_ */
