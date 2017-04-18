/*
 * ParallelSearch.h
 *
 *  Created on: Apr 18, 2017
 *      Author: yuu
 */

#ifndef MP_SRC_PARALLELDFS_H_
#define MP_SRC_PARALLELDFS_H_

#include "../src/variable_length_itemset.h"
#include "MPI_Data.h"
#include "mpi_tag.h"

namespace lamp_search {

/*
 * TODO: Should Extract more fields here.
 */
class ParallelDFS {
public:
	ParallelDFS(MPI_Data& mpi_data, TreeSearchData* treesearch_data, Log* log,
			Timer* timer);
	virtual ~ParallelDFS();
	virtual void Search();

protected:
	/**
	 * Parallel Search Basic Methods.
	 */
	virtual bool Probe(MPI_Data& mpi_data, TreeSearchData* treesearch_data) = 0;
	virtual void Distribute(MPI_Data& mpi_data, TreeSearchData* treesearch_data) = 0;
	virtual void Give(MPI_Data& mpi_data, VariableLengthItemsetStack * st,
			int steal_num) = 0;
	virtual void Reject(MPI_Data& mpi_data) = 0;
	virtual void Steal(MPI_Data& mpi_data) = 0;
	virtual void Check(MPI_Data& mpi_data) = 0;
	virtual bool ProcessNode(MPI_Data& mpi_data, TreeSearchData*treesearch_data) = 0;

	/**
	 * Network
	 */
	bool IsLeafInTopology(MPI_Data& mpi_data) const;

	/**
	 * Helper methods to encapsulate MPI operations with timer and logger.
	 * TODO: these methods should be generally in Parallel class.
	 */
	int CallIprobe(MPI_Status * status, int * count, int * src);
	int CallRecv(void * buffer, int count, MPI_Datatype type, int src, int tag,
			MPI_Status * status);
	int CallBsend(void * buffer, int count_int, MPI_Datatype type, int dest,
			int tag);
	int CallBcast(void * buffer, int data_count, MPI_Datatype type);


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

#endif /* MP_SRC_PARALLELDFS_H_ */
