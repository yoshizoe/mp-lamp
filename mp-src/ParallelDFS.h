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
			Timer* timer, std::ostream& ofs);
	virtual ~ParallelDFS();
	virtual void Search();

protected:
	/**
	 * Parallel Search Basic Methods.
	 */
	virtual bool Probe(TreeSearchData* treesearch_data);
	virtual void ProbeExecute(TreeSearchData* treesearch_data,
			MPI_Status* probe_status, int probe_src, int probe_tag);
	virtual void Distribute(TreeSearchData* treesearch_data);
	virtual void Give(VariableLengthItemsetStack * st,
			int steal_num) ;
	virtual void Reject();
	virtual void Steal();

	virtual void ProcAfterProbe() = 0;
	virtual void Check() = 0;
	virtual bool ExpandNode(TreeSearchData*treesearch_data) = 0;

	bool AccumCountReady() const;
	/**
	 * ProbeExecute implementation
	 */
	// 0: count, 1: time warp flag, 2: empty flag
	void SendDTDRequest();
	void RecvDTDRequest(TreeSearchData* treesearch_data,
			int src);

	bool DTDReplyReady() const;
	void DTDCheck();

	// 0: count, 1: time warp flag, 2: empty flag
	void SendDTDReply(TreeSearchData* treesearch_data);
	void RecvDTDReply(TreeSearchData* treesearch_data,
			int src);
	virtual bool HasJobToDo();

	void SendBcastFinish();
	void RecvBcastFinish(int src);

	void SendRequest(int dst, int is_lifeline); // for random thieves, is_lifeline = -1
	void RecvRequest(TreeSearchData* treesearch_data,
			int src);

	// 0: time zone, 1: is_lifeline
	void SendReject(int dst);
	void RecvReject(int src);

	// 1: time zone
	void SendGive(VariableLengthItemsetStack * st, int dst,
			int is_lifeline);

	// sets lifelines_activated_ = false
	// lifelines_activated_ becomes false only in this case (reject does NOT)
	void RecvGive(TreeSearchData* treesearch_data, int src,
			MPI_Status status);

	/**
	 * Network
	 */
	bool IsLeafInTopology() const;

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
	std::ostream& lfs_;

	static std::ofstream null_stream_;
	std::ostream& D(int level, bool show_phase = true);
};

} /* namespace lamp_search */

#endif /* MP_SRC_PARALLELDFS_H_ */
