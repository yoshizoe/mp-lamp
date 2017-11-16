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

	class DFSStack;

/*
 * TODO: Should Extract more fields here.
 */
class ParallelDFS {
private:
    int phase_;

public:
	ParallelDFS(DFSStack* stack, MPI_Data& mpi_data, VariableLengthItemsetStack* give_stack, Log* log,	Timer* timer, std::ostream& ofs);
	virtual ~ParallelDFS();
	virtual void Search();

protected:
	/**
	 * Parallel Search Basic Methods.
	 */
	virtual bool Probe();
    virtual void ProbeExecute(MPI_Status* probe_status, int probe_src, int probe_tag);
//	virtual void ProbeExecute(TreeSearchData* treesearch_data,
//			MPI_Status* probe_status, int probe_src, int probe_tag);
	virtual void Distribute();
	virtual void Give(VariableLengthItemsetStack * st,
			int steal_num) ;
	virtual void Reject();
	virtual void Steal();

	virtual void ProcAfterProbe() = 0;
	virtual void Check() = 0;

    void SetPhase(int phase) { phase_ = phase; }
    int GetPhase() const { return phase_; }

    bool ExpandNode(TreeSearchData*treesearch_data);

	bool AccumCountReady() const;
	/**
	 * ProbeExecute implementation
	 */
	// 0: count, 1: time warp flag, 2: empty flag
	void SendDTDRequest();
	void RecvDTDRequest(int src);

	bool DTDReplyReady() const;
	void DTDCheck();

	// 0: count, 1: time warp flag, 2: empty flag
	void SendDTDReply();
	void RecvDTDReply(int src);

	void SendBcastFinish();
	void RecvBcastFinish(int src);

	void SendRequest(int dst, int is_lifeline); // for random thieves, is_lifeline = -1
	void RecvRequest(int src);

	// 0: time zone, 1: is_lifeline
	void SendReject(int dst);
	void RecvReject(int src);

	// 1: time zone
	void SendGive(VariableLengthItemsetStack * st, int dst,
			int is_lifeline);

	// sets lifelines_activated_ = false
	// lifelines_activated_ becomes false only in this case (reject does NOT)
	void RecvGive(int src, MPI_Status status);

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
	DFSStack* dfs_stack_;
    VariableLengthItemsetStack* give_stack_;

//	TreeSearchData* treesearch_data;
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
