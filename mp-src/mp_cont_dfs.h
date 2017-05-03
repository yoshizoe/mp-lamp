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

#ifndef _LAMP_SEARCH_MP_CONT_DFS_H_
#define _LAMP_SEARCH_MP_CONT_DFS_H_

#include <vector>
#include <algorithm>

#include <boost/array.hpp>
#include <map>

#include "mpi.h"

#include "utils.h"
#include "topk.h"
#include "variable_length_itemset.h"
#include "timer.h"
#include "../src/contdatabase.h"

#include "random.h"

#include "MPI_Data.h"
#include "StealState.h"
#include "SignificantSetResults.h"
#include "Log.h"
#include "DTD.h"
#include "FixedSizeStack.h"
#include "mpi_tag.h"

namespace lamp_search {

class MP_CONT_LAMP {
public:
	/**
	 * n: granularity of tasks
	 * m: number of random victim candidates
	 * w: number of random victim tries
	 * l: power of lifeline graph
	 */
	MP_CONT_LAMP(ContDatabase* d, int rank, int nu_proc, int n,
			bool n_is_ms, int w, int l, int m, int disretizeFreq = 0);

	~MP_CONT_LAMP();

	// echo operation uses n-array tree (default: 3-ary)
	static const int k_echo_tree_branch;

//	static int ComputeZ(int p, int l);

// call this before mainloop ?
// set dtd counters and misc.
	void Init();

	void CheckPoint();

	// clear unreceived MPI messages
	void ClearTasks();

	void InitTreeRequest();
	void SetTreeRequest();

	void Search();

	std::ostream & PrintDBInfo(std::ostream & out) const;

	std::ostream & PrintResults(std::ostream & out) const;
	// std::ostream & PrintSignificantMap(std::ostream & out) const;
	std::ostream & PrintSignificantSet(std::ostream & out) const;

	std::ostream & PrintLog(std::ostream & out) const;
	std::ostream & PrintAggrLog(std::ostream & out);

	std::ostream & PrintPLog(std::ostream & out);
	std::ostream & PrintAggrPLog(std::ostream & out);

private:
	static const int k_int_max;
	// assuming (digits in long long int) > (bits of double mantissa)
	static const long long int k_cs_max;

	static const int k_probe_period;

	DTD dtd_;
	MPI_Data mpi_data_;

	ContDatabase* d_;

	int disretizeFreq;

	Log log_;
	Timer * timer_;

	// Domain Data is included here for now.

//	// variables for LAMP
//	int lambda_max_; // equals to maximum support of single item // getMinSup
//	// initially set to 1. will be incremented to N if cs_thr[N] exceeded
//	// in 1st phase, search will be pruned if (sup_num < lambda_)
//	int lambda_; // getMinSup
//	// todo: initial value can be set after checking all single item itemsets
//	long long int * cs_thr_; // cs_thr[sup] shows closed set num threshold // getMinSup
//
//	double sig_level_; // initially set to 1. set to 0.05 (FLAGS_a) / cs_thr_[global_sup_thr_-1] // LAMP

	double GetInterimSigLevel(int lambda) const; // LAMP

	VariableLengthItemsetStack * node_stack_;
// todo: prepare stack with no sup hist
// 0: time zone, 1: is_lifeline
	VariableLengthItemsetStack * give_stack_;

// will be false if w random steal and all lifeline steal finished
// will be true if RecvGive
	StealState stealer_; // TODO

	int phase_; // 1, 2, 3

	int itemset_buf_[VariableLengthItemsetStack::kMaxItemsPerSet];
//	uint64 * sup_buf_, *child_sup_buf_;

	VariableLengthItemsetStack * freq_stack_; // record freq itemsets
	std::multimap<double, int *> freq_map_; // record (pval, *itemsets)

	VariableLengthItemsetStack * significant_stack_;

	/* sorting itemset
	 1, ascending order of pval
	 2, descending order of item numbers
	 3, dictionary order of items
	 */
	std::set<ContSignificantSetResult, cont_sigset_compare> significant_set_;

	long long int total_expand_num_;
	long long int expand_num_;
	long long int closed_set_num_;

//--------
// third phase
// insert pointer into significant_map_ (do not sort the stack itself)
	void SortSignificantSets();

//--------
// for printing results

	long long int num_final_testable_patterns;
//	int final_support_;
	double final_sig_level_;

	int CallBcast(void * buffer, int data_count, MPI_Datatype type);
// todo: implement call reduce, call gather

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
	bool last_bcast_was_dtd_;

};

} // namespace lamp_search

#endif // _LAMP_SEARCH_MP_DFS_H_

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
