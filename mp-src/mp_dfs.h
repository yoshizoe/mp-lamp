// Copyright (c) 2016, Kazuki Yoshizoe
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
// may be used to endorse or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// AREDISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef _LAMP_SEARCH_MP_DFS_H_
#define _LAMP_SEARCH_MP_DFS_H_

#include <vector>
#include <algorithm>

#include <boost/array.hpp>
#include <map>

//#include "mpi.h"

#include "utils.h"
#include "sorted_itemset.h"
#include "topk.h"
#include "variable_bitset_array.h"
#include "variable_length_itemset.h"
#include "table_vba.h"
#include "timer.h"

#include "random.h"
#include "database.h"
#include "lamp_graph.h"

#include "DTD.h"
#include "Log.h"
#include "SignificantSetResults.h"
#include "StealState.h"
#include "FixedSizeStack.h"

namespace lamp_search {

// use these as mpi tag
struct Tag {
	enum TaskType {
		// control tasks
		CONTROL_TASK_BEGIN = 0,
		DTD_REQUEST,
		DTD_REPLY,

		DTD_ACCUM_REQUEST, // request reporting accum count
		DTD_ACCUM_REPLY, // reduce closed set count

		BCAST_FINISH,
		CONTROL_TASK_END,

		// basic tasks
		BASIC_TASK_BEGIN,

		LAMBDA, // for lamp

		REQUEST, // request tasks
		REJECT, // reject requests
		GIVE, // send tasks

		BASIC_TASK_END,

		// third phase tasks
		THIRD_PHASE_BEGIN,

		RESULT_REQUEST,
		RESULT_REPLY,

		THIRD_PHASE_END,
	};
};

struct MPI_Data {
	MPI_Data(int rank, int nu_proc, int n, bool n_is_ms, int w, int l,
			int m) : mpiRank_(rank), nTotalProc(nu_proc), granularity(n), isGranularitySec(n_is_ms), {
	}
	/**
	 * MPI variables
	 */
	int mpiRank_;
	int nTotalProc;
	int * bsend_buffer_;
	/**
	 * Hypercube topology TODO: This should be put in something like ParallelTopology class.
	 */
	int lHypercubeEdge; // power of lifeline graph (length the hypercube edge)
	int hypercubeDimension; // dimension of lifeline (dimension of the hypercube)
	int * victims_; // proc id of random victims
	int * lifelines_; // proc id of lifeline buddies
	bool * lifelines_activated_;
	FixedSizeStack * thieves_; // max size == nu_proc_
	FixedSizeStack * lifeline_thieves_; // size == lifelines_ size + 3
	int bcast_source_;
	int * bcast_targets_;
	int nRandStealTrials; // number of random steal trials
	int nRandStealCands; // number of random steal candidates
	int granularity;
	bool isGranularitySec;

	// TODO: Specific to LCM and LAMP???????
	VariableLengthItemsetStack * node_stack_;
	VariableLengthItemsetStack * give_stack_;

	// TODO: what is it????
	bool * accum_flag_;
};

struct LCM_Data {

};

struct LAMP_Data {

};

struct Significant_Data {

};

class MP_LAMP {
public:
	/**
	 * n: granularity of tasks
	 * m: number of random victim candidates
	 * w: number of random victim tries
	 * l: power of lifeline graph
	 */
	MP_LAMP(int rank, int nu_proc, int n, bool n_is_ms, int w, int l, int m);

	~MP_LAMP();

	// echo operation uses n-array tree (default: 3-ary)
	static const int k_echo_tree_branch;

//	static int ComputeZ(int p, int l);

	// call this before mainloop ?
	// set dtd counters and misc.
//	void Init();

	void CheckPoint();

	// clear unreceived MPI messages
	void ClearTasks();

//	void InitTreeRequest();
	void SetTreeRequest();

	// read file, prepare database, broadcast to all procs
	// for h_==0
	void InitDatabaseRoot(std::istream & is1, std::istream &is2);
	void InitDatabaseRoot(std::istream & is1, int posnum);
	// other
	void InitDatabaseSub(bool pos);

	void Search();
//	void MainLoop();

	void SearchStraw1(); // no workload distribution
	void MainLoopStraw1();
	void SearchStraw2(); // centralized queue
	void MainLoopStraw2();

	const Database<uint64> & GetDatabase() {
		return *d_;
	}

	std::ostream & PrintDBInfo(std::ostream & out) const;

	std::ostream & PrintResults(std::ostream & out) const;
	// std::ostream & PrintSignificantMap(std::ostream & out) const;
	std::ostream & PrintSignificantSet(std::ostream & out) const;

	std::ostream & PrintLog(std::ostream & out) const;
	std::ostream & PrintAggrLog(std::ostream & out);

	std::ostream & PrintPLog(std::ostream & out);
	std::ostream & PrintAggrPLog(std::ostream & out);

private:
	void InitParallel();

	ParallelSearch* search;

	MPI_Data mpi_data;
	DTD dtd_;

	static const int k_int_max;
	// assuming (digits in long long int) > (bits of double mantissa)
	static const long long int k_cs_max;

	static const int k_probe_period;

//	int * bsend_buffer_;

	// GLB variables

	boost::mt19937 rng_; // use seed as rank
	boost::uniform_smallint<int> dst_p_;
	boost::uniform_smallint<int> dst_m_;
	boost::variate_generator<boost::mt19937&, boost::uniform_smallint<int> > rand_p_;
	boost::variate_generator<boost::mt19937&, boost::uniform_smallint<int> > rand_m_;

//	int bcast_source_;
//	int * bcast_targets_;

	// bool dtd_request_;
	// bool dtd_flag_[3];

//	bool echo_waiting_; // new! [2015-10-05 22:23]
	//  bool bcast_requesting_; // new! [2015-10-01 13:49]

	//bool accum_requesting_;
//	bool * accum_flag_;

//	int * victims_; // proc id of random victims
//	int * lifelines_; // proc id of lifeline buddies
//
//	FixedSizeStack * thieves_; // max size == nu_proc_
//	FixedSizeStack * lifeline_thieves_; // size == lifelines_ size + 3
//	bool * lifelines_activated_;

	// proc id of lifeline thieves
	// bool * thieves_requests_;
	// bool * lifeline_requests_;

	// LAMP variables

	Database<uint64> * d_;
	LampGraph<uint64> * g_;
	VariableBitsetHelper<uint64> * bsh_; // bitset helper

	Log log_;
	Timer * timer_;

	// variables for LAMP
	int lambda_max_; // equals to maximum support of single item
	// initially set to 1. will be incremented to N if cs_thr[N] exceeded
	// in 1st phase, search will be pruned if (sup_num < lambda_)
	int lambda_;
	// todo: initial value can be set after checking all single item itemsets
	double sig_level_; // initially set to 1. set to 0.05 (FLAGS_a) / cs_thr_[global_sup_thr_-1]
	double * pmin_thr_; // pmin_thr[sup] == tbl.PMin(sup), maybe redundant
	long long int * cs_thr_; // cs_thr[sup] shows closed set num threshold

//	void CheckCSThreshold();
//
//	bool ExceedCsThr() const;
//	int NextLambdaThr() const;
//	void IncCsAccum(int sup_num);
//	double GetInterimSigLevel(int lambda) const;

	// cs_accum_array is int array of 0..lambda_max_ (size lambda_max_+1)
	// cs_accum_array_[sup] shows closed set num with support higher than or equals to sup

	// note: allocate cs_accum_array_ as
	//   cs_accum_array_base_ = new int[lambda_max_+2]
	// and then
	//   cs_accum_array_ = cs_accum_array_base_+1
	// so that cs_accum_array_[-1] (== cs_accum_array_base_[0])
	// can be used for sending the timestamp counter
	// long long int * cs_accum_array_base_; // int array of -1..lambda_max_ (size lambda_max_+2)

	// 0: count, 1: time warp, 2: empty flag, 3--: array
//	long long int * dtd_accum_array_base_; // int array of [-3..lambda_max_] (size lambda_max_+4)
//	// -3: count, -2: time warp, -1: empty flag, 0--: array
//	long long int * accum_array_; // int array of [0...lambda_max_] (size lambda_max_+1)
//
//	// 0: count, 1: time warp, 2: empty flag, 3--: array
//	long long int * dtd_accum_recv_base_; // int array of [-3..lambda_max_] (size lambda_max_+4)
//	// -3: count, -2: time warp, -1: empty flag, 0--: array
//	long long int * accum_recv_; // int array of [0...lambda_max_] (size lambda_max_+1)
//
//	VariableLengthItemsetStack * node_stack_;
//	// todo: prepare stack with no sup hist
//	// 0: time zone, 1: is_lifeline
//	VariableLengthItemsetStack * give_stack_;

	// periodic closed set count reduce.
//	void Probe();
//	bool processing_node_; // prevent termination while processing node

	// void ProbeAccumTask();
	// void ProbeBasicTask();
	// void ProbeControlTask();

	// first, send to random thief and then to lifeline theives
	// random thief has higher priority
//	void Distribute();

//	void Give(VariableLengthItemsetStack * st, int steal_num);

//	void Deal();

	// send reject to remaining requests
//	void Reject();

	// set this in steal, reset this in RecvGive and RecvReject
	// small difference from x10 implementation
//	bool waiting_;

	// send steal requests
	// two phase, 1, random, 2, lifeline
	// if succeeds, break
	// note: original x10 GLB implementation always probes from i=0..z, is this OK?
	// how about prepare int steal_id_; and do steal_id_++/ steal_id_ %= z ?
	// note:
	// don't send multiple requests at once
//	void Steal();
	// Steal needs change from x10 because of "bool waiting"
	// Steal sends one request each time it is called
	// there should be steal_state and counters c_r and c_l (random and lifeline)
	// c_r=0, c_l=0
	// 0, if (state == RANDOM && w > 0) goto 1, else goto 4
	// 1, send random request, c_r++
	// 2, wait for reject or give
	// 3, if c_r >= w, c_r = 0, set state LIFELINE and return
	//
	// 4, send lifeline request, c_l++
	// 5, wait for reject or give
	// 6, if c_l >= z, c_l = 0, set state RANDOM and return

	// steal with probe
//	void Steal2();

	// will be false if w random steal and all lifeline steal finished
	// will be true if RecvGive
//	StealState stealer_;

	//--------
	// control

	// 0: count, 1: time warp flag, 2: empty flag
//	void SendDTDRequest();
//	void RecvDTDRequest(int src);
//
//	bool DTDReplyReady() const;
//	void DTDCheck();
//
//	// 0: count, 1: time warp flag, 2: empty flag
//	void SendDTDReply();
//	void RecvDTDReply(int src);
//
//	bool DTDAccumReady() const;
//
//	// 0: count, 1: time warp flag, 2: empty flag, 3--: data
//	void SendDTDAccumRequest();
//	void RecvDTDAccumRequest(int src);
//
//	// 0: count, 1: time warp flag, 2: empty flag, 3--: data
//	void SendDTDAccumReply();
//	void RecvDTDAccumReply(int src);
//
//	void SendBcastFinish();
//	void RecvBcastFinish(int src);

	//--------

	int phase_; // 1, 2, 3

	//--------
	// basic

	// send recv functions
//	void SendRequest(int dst, int is_lifeline); // for random thieves, is_lifeline = -1
//	void RecvRequest(int src);
//
//	// 0: time zone, 1: is_lifeline
//	void SendReject(int dst);
//	void RecvReject(int src);
//
//	// 1: time zone
//	void SendGive(VariableLengthItemsetStack * st, int dst, int is_lifeline);
//
//	// sets lifelines_activated_ = false
//	// lifelines_activated_ becomes false only in this case (reject does NOT)
//	void RecvGive(int src, MPI_Status status);
//
//	void SendLambda(int lambda);
//	void RecvLambda(int src);

	// 0: time zone
	// search for depth 1 and get initial lambda
	void PreProcessRootNode();

	// provide int n and bool n_is_ms_ (if n_is_ms_==false, it shows number of nodes)
	bool ProcessNode(int n);
	bool CheckProcessNodeEnd(int n, bool n_is_ms, int processed,
			long long int start_time);

	bool ProcessNodeStraw1(int n);
	bool ProcessNodeStraw2(int n);

	int itemset_buf_[VariableLengthItemsetStack::kMaxItemsPerSet];
	uint64 * sup_buf_, *child_sup_buf_;

	VariableLengthItemsetStack * freq_stack_; // record freq itemsets
	std::multimap<double, int *> freq_map_; // record (pval, *itemsets)

	VariableLengthItemsetStack * significant_stack_;

	/* sorting itemset
	 1, ascending order of pval
	 2, descending order of item numbers
	 3, dictionary order of items
	 */

	std::set<SignificantSetResult, sigset_compare> significant_set_;

	// add pos_sup_num info in significant map
	// add stable sort
	// std::multimap< double, int * > significant_map_;

	long long int total_expand_num_;
	long long int expand_num_;
	long long int closed_set_num_;

	//--------
	// third phase

	// void ProbeThirdPhaseTask();

	bool AccumCountReady() const;

	void SendResultRequest();
	void RecvResultRequest(int src);

	void SendResultReply();
	void RecvResultReply(int src, MPI_Status status);

	void ExtractSignificantSet();

	// insert pointer into significant_map_ (do not sort the stack itself)
	void SortSignificantSets();

	//--------
	// for printing results

	long long int final_closed_set_num_;
	int final_support_;
	double final_sig_level_;

	// true if bcast_targets_ are all -1
//	bool IsLeaf() const;
//
//	// return flag. if (flag), it is ready to receive
//	int CallIprobe(MPI_Status * status, int * count, int * src);
//
//	int CallRecv(void * buffer, int count, MPI_Datatype type, int src, int tag,
//			MPI_Status * status);
//	int CallBsend(void * buffer, int count_int, MPI_Datatype type, int dest,
//			int tag);
//
//	int CallBcast(void * buffer, int data_count, MPI_Datatype type);
//	// todo: implement call reduce, call gather

	//--------
	// for debug

	void SetLogFileName();

	std::string log_file_name_;
	std::ofstream lfs_;

	/** Null stream.
	 This file stream will never be opened and acts as a null stream
	 immitating Fuego's SgDebug(). */
	static std::ofstream null_stream_;
	std::ostream& D(int level, bool show_phase = true);

	// static std::stringstream alert_ss;
	// std::ostream& Alert();

	//--------
	// testing [2015-10-01 13:54]
//	bool last_bcast_was_dtd_;

};

} // namespace lamp_search

#endif // _LAMP_SEARCH_MP_DFS_H_

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
