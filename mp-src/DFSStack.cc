//
// Created by Takuo Doi on 2017/11/05.
//

#include <cstdio>
#include <cassert>
#include "gflags/gflags.h"

#include "DFSStack.h"
#include "MPI_Data.h"

DEFINE_int32(probe_period_, 128, "probe period during process node");
DEFINE_bool(probe_period_is_ms_, false, "true: probe period is milli sec, false: num loops");
DECLARE_bool(third_phase_); // true, "do third phase"

namespace lamp_search {

    std::ofstream DFSStack::null_stream_;

    void DFSStack::CallForProbe() {
        callback_();
    }

    std::ostream& DFSStack::D(int level, bool show_phase) const {
        bool FLAGS_d = true;
        if (FLAGS_d == 0)
            return null_stream_;
        if (level <= FLAGS_d) {
            if (show_phase) {
                lfs_ << std::setw(4) << "DFS" << ": ";
            }
            return lfs_;
        } else {
            return null_stream_;
        }
    }

    void DFSStack::CheckProbe(int& accum_period_counter_, long long int& lap_time, Timer* timer, Log* log) {
        // TODO: whatever this is trying to do, it should be factored into a function.
        //       Why is it Probing while in the expansion loop?
        accum_period_counter_++;
        if (FLAGS_probe_period_is_ms_) {      // using milli second
            if (accum_period_counter_ >= 64) { // TODO: what is this magic number?
                // to avoid calling timer_ frequently, time is checked once in 64 loops
                // clock_gettime takes 0.3--0.5 micro sec
                accum_period_counter_ = 0;
                long long int elt = timer->Elapsed();
                if (elt - lap_time >= FLAGS_probe_period_ * 1000000) {
                    log->d_.process_node_time_ += elt - lap_time;

                    CallForProbe();

                    lap_time = timer->Elapsed();
                }
            }
        } else {            // not using milli second
            if (accum_period_counter_ >= FLAGS_probe_period_) {
                accum_period_counter_ = 0;
                log->d_.process_node_time_ += timer->Elapsed() - lap_time;

                CallForProbe();

                lap_time = timer->Elapsed();
            }
        }
        // note: do this before PushPre is called [2015-10-05 21:56]
    }


};