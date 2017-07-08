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

// MP LAMP
#include <iostream>
#include <iomanip>
#include <fstream>
#include "gflags/gflags.h"

#include <sys/time.h>
#include <sys/resource.h>

#include <cstdio>
//#include <thread> // std::this_thread::sleep_for
//#include <chrono> // std::chrono::seconds

#include "mpi.h"

//#include "../mp-src/ParallelContinuousPM.h"
#include "timer.h"
#include "../mp-src/mp_cont_dfs.h"
#include "../src/contdatabase.h"

//DEFINE_bool(lcm, false, "item file is lcm style");

DEFINE_string(item, "", "filename of item set");
DEFINE_string(pos, "", "filename of positive / negative file");
DEFINE_int32(posnum, 0,
		"positive total (used if for 1st phase only)");

DECLARE_bool(log); // false, "show log", mp-lamp.cc , true, "show log"

//DECLARE_bool(second_phase);// true, "do second phase"
DECLARE_bool(third_phase);// true, "do third phase"
//
DEFINE_int32(n, 1000, "granularity of one Node process");
DEFINE_bool(n_is_ms, true,
		"true: n is milli sec, false: n is num task");
DEFINE_int32(w, 1, "number of random steal attempts");
DEFINE_int32(m, 2, "number of maximum random steal candidates"); // 1024 is default of x10 GLB
DEFINE_int32(l, 2, "power of lifeline graph");

DEFINE_bool(discretize, true, "true if we run awesome");
DEFINE_double(ratio, 0.95, "diminishing ratio of the discretized thresholds");

//
//DEFINE_bool(straw1, false, "use Strawman1 for comparison");
//
//DEFINE_int32(sleep, 0, "sleep in the beinning (for debugger attach)");

using
namespace lamp_search;

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	google::ParseCommandLineFlags(&argc, &argv, true);

	int rank, nu_proc;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nu_proc);

	if (rank == 0) {
	  printf("CONT-LAMP: %d\n", nu_proc);
	}
	
	long long int search_start_time, search_end_time;

	{
		std::ifstream item_file;
		item_file.open(FLAGS_item.c_str(), std::ios::in);
		std::ifstream class_file;
		class_file.open(FLAGS_pos.c_str(), std::ios::in);
		ContDatabase* d = new ContDatabase(item_file, class_file);
		d->ShowInfo();
//		exit(0);
		int discretize = 0;
		if (FLAGS_discretize) {
			discretize = 1;
		}
		printf("discretize = %d\n", discretize);
		MP_CONT_LAMP* search = new MP_CONT_LAMP(d, rank, nu_proc,
				FLAGS_n, FLAGS_n_is_ms, FLAGS_w, FLAGS_l, FLAGS_m,
				discretize, FLAGS_ratio);

		search_start_time = Timer::GetInstance()->Elapsed();
		search->Search();
		search_end_time = Timer::GetInstance()->Elapsed();

		if (rank == 0)
			std::cout << "# time all=" << std::setw(12)
					<< search_end_time / GIGA << "\ttime search="
					<< std::setw(12)
					<< (search_end_time - search_start_time) / GIGA
					<< std::endl;
		if (rank == 0) {
			search->PrintResults(std::cout);
			//if (FLAGS_log) search->PrintLog(std::cout);
			if (FLAGS_log)
				search->PrintAggrLog(std::cout);
			if (FLAGS_log)
				search->PrintAggrPLog(std::cout);
		}

		delete search;
		delete d;
		MPI_Barrier( MPI_COMM_WORLD);

	}

	MPI_Finalize();
	return 0;
}

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
