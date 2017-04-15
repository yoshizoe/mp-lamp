/*
 * Log.cpp
 *
 *  Created on: Apr 14, 2017
 *      Author: yuu
 */
#include "Log.h"

#include "mp_dfs.h"
#include <vector>

namespace lamp_search {

Log::Log() :
		plog_buf_(NULL), plog_gather_buf_(NULL), gather_buf_(NULL) {
	Init();
}

Log::~Log() {
	if (plog_buf_)
		delete[] plog_buf_;
	if (plog_gather_buf_)
		delete[] plog_gather_buf_;
	if (gather_buf_)
		delete[] gather_buf_;
}

void Log::Init() {
	idle_start_ = 0;

	InitPeriodicLog();

	d_.Init();
}

void Log::GatherLog(int nu_proc) {
	int sec = plog_.size();
	MPI_Allreduce(&sec, &sec_max_, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	assert(sec <= sec_max_);

	if (sec_max_ > 0) {
		plog_buf_ = new PeriodicLog_T[sec_max_];
		for (int s = 0; s < sec; s++)
			plog_buf_[s] = plog_[s];
		for (int s = sec; s < sec_max_; s++)
			plog_buf_[s].Clear();

		plog_gather_buf_ = new PeriodicLog_T[nu_proc * sec_max_];
		MPI_Gather((void*) (plog_buf_), sizeof(PeriodicLog_T) * sec_max_,
				MPI_CHAR, (void*) (plog_gather_buf_),
				sizeof(PeriodicLog_T) * sec_max_, MPI_CHAR, 0, MPI_COMM_WORLD);
	}

	gather_buf_ = new LogData[nu_proc];
	for (int i = 0; i < nu_proc; i++)
		gather_buf_->Init();

	MPI_Gather((void*) (&d_), sizeof(LogData), MPI_CHAR, (void*) (gather_buf_),
			sizeof(LogData), MPI_CHAR, 0, MPI_COMM_WORLD);
}

void Log::InitPeriodicLog() {
	periodic_log_start_ = -1; // -1 means not started
	next_log_time_in_second_ = 0;
}

void Log::StartPeriodicLog() {
	periodic_log_start_ = Timer::GetInstance()->Elapsed();
	next_log_time_in_second_ = 0;
}

void Log::FinishPeriodicLog() {
	periodic_log_start_ = -1; // -1 means not started
}

void Log::TakePeriodicLog(long long int capacity, int lambda,
		int phase) {
	if (periodic_log_start_ >= 0) {
		long long int current_time = Timer::GetInstance()->Elapsed();
		long long int elapsed = current_time - periodic_log_start_;
		if (elapsed >= next_log_time_in_second_ * 1000000000) {

			PeriodicLog_T t;
			t.seconds_ = elapsed;
			t.capacity_ = capacity;
			t.lambda_ = lambda;
			t.phase_ = phase;

			plog_.push_back(t);

			next_log_time_in_second_++;
		}
	}
}

void Log::LogData::Init() {
	search_start_time_ = 0ll;

	iprobe_num_ = 0ll;
	iprobe_time_ = 0ll;
	iprobe_time_max_ = 0ll;

	iprobe_succ_num_ = 0ll;
	iprobe_succ_time_ = 0ll;
	iprobe_succ_time_max_ = 0ll;

	iprobe_fail_num_ = 0ll;
	iprobe_fail_time_ = 0ll;
	iprobe_fail_time_max_ = 0ll;

	probe_num_ = 0ll;
	probe_time_ = 0ll;
	probe_time_max_ = 0ll;

	recv_num_ = 0ll;
	recv_time_ = 0ll;
	recv_time_max_ = 0ll;

	bsend_num_ = 0ll;
	bsend_time_ = 0ll;
	bsend_time_max_ = 0ll;

	bcast_num_ = 0ll;
	bcast_time_ = 0ll;
	bcast_time_max_ = 0ll;

	// accum_task_time_ = 0ll;
	// accum_task_num_ = 0ll;
	// basic_task_time_ = 0ll;
	// basic_task_num_ = 0ll;
	// control_task_time_ = 0ll;
	// control_task_num_ = 0ll;

	dtd_phase_num_ = 0ll;
	dtd_phase_per_sec_ = 0.0;
	dtd_accum_phase_num_ = 0ll;
	dtd_accum_phase_per_sec_ = 0.0;
	// dtd_reply_num_  = 0ll;

	// accum_phase_num_  = 0ll;
	// accum_phase_per_sec_  = 0.0;

	lifeline_steal_num_ = 0ll;
	lifeline_nodes_received_ = 0ll;
	steal_num_ = 0ll;
	nodes_received_ = 0ll;

	lifeline_given_num_ = 0ll;
	lifeline_nodes_given_ = 0ll;
	given_num_ = 0ll;
	nodes_given_ = 0ll;

	process_node_num_ = 0ll;
	process_node_time_ = 0ll;

	preprocess_time_ = 0ll;

	idle_time_ = 0ll;

	pval_table_time_ = 0ll;

	node_stack_max_itm_ = 0ll;
	give_stack_max_itm_ = 0ll;
	node_stack_max_cap_ = 0ll;
	give_stack_max_cap_ = 0ll;

	cleared_tasks_ = 0ll;
}

void Log::Aggregate(int nu_proc) {
	a_.Init();
	for (int i = 0; i < nu_proc; i++) {
		a_.iprobe_num_ += gather_buf_[i].iprobe_num_;
		a_.iprobe_time_ += gather_buf_[i].iprobe_time_;
		a_.iprobe_time_max_ = std::max(a_.iprobe_time_max_,
				gather_buf_[i].iprobe_time_max_);

		a_.iprobe_succ_num_ += gather_buf_[i].iprobe_succ_num_;
		a_.iprobe_succ_time_ += gather_buf_[i].iprobe_succ_time_;
		a_.iprobe_succ_time_max_ = std::max(a_.iprobe_succ_time_max_,
				gather_buf_[i].iprobe_succ_time_max_);

		a_.iprobe_fail_num_ += gather_buf_[i].iprobe_fail_num_;
		a_.iprobe_fail_time_ += gather_buf_[i].iprobe_fail_time_;
		a_.iprobe_fail_time_max_ = std::max(a_.iprobe_fail_time_max_,
				gather_buf_[i].iprobe_fail_time_max_);

		a_.probe_num_ += gather_buf_[i].probe_num_;
		a_.probe_time_ += gather_buf_[i].probe_time_;
		a_.probe_time_max_ = std::max(a_.probe_time_max_,
				gather_buf_[i].probe_time_max_);

		a_.recv_num_ += gather_buf_[i].recv_num_;
		a_.recv_time_ += gather_buf_[i].recv_time_;
		a_.recv_time_max_ = std::max(a_.recv_time_max_,
				gather_buf_[i].recv_time_max_);

		a_.bsend_num_ += gather_buf_[i].bsend_num_;
		a_.bsend_time_ += gather_buf_[i].bsend_time_;
		a_.bsend_time_max_ = std::max(a_.bsend_time_max_,
				gather_buf_[i].bsend_time_max_);

		a_.bcast_num_ += gather_buf_[i].bcast_num_;
		a_.bcast_time_ += gather_buf_[i].bcast_time_;
		a_.bcast_time_max_ = std::max(a_.bcast_time_max_,
				gather_buf_[i].bcast_time_max_);

		// a_.accum_task_time_ += gather_buf_[i].accum_task_time_;
		// a_.accum_task_num_ += gather_buf_[i].accum_task_num_;
		// a_.basic_task_time_ += gather_buf_[i].basic_task_time_;
		// a_.basic_task_num_ += gather_buf_[i].basic_task_num_;
		// a_.control_task_time_ += gather_buf_[i].control_task_time_;
		// a_.control_task_num_ += gather_buf_[i].control_task_num_;

		// note: dtd_request_num_ only available at rank 0
		// note: accum_phase_num_ only available at rank 0

		a_.lifeline_steal_num_ += gather_buf_[i].lifeline_steal_num_;
		a_.lifeline_nodes_received_ += gather_buf_[i].lifeline_nodes_received_;
		a_.steal_num_ += gather_buf_[i].steal_num_;
		a_.nodes_received_ += gather_buf_[i].nodes_received_;

		a_.lifeline_given_num_ += gather_buf_[i].lifeline_given_num_;
		a_.lifeline_nodes_given_ += gather_buf_[i].lifeline_nodes_given_;
		a_.given_num_ += gather_buf_[i].given_num_;
		a_.nodes_given_ += gather_buf_[i].nodes_given_;

		a_.process_node_num_ += gather_buf_[i].process_node_num_;
		a_.process_node_time_ += gather_buf_[i].process_node_time_;

		a_.preprocess_time_ += gather_buf_[i].preprocess_time_;

		a_.idle_time_ += gather_buf_[i].idle_time_;

		a_.pval_table_time_ += gather_buf_[i].pval_table_time_;

		a_.node_stack_max_itm_ = std::max(a_.node_stack_max_itm_,
				gather_buf_[i].node_stack_max_itm_);
		a_.give_stack_max_itm_ = std::max(a_.give_stack_max_itm_,
				gather_buf_[i].give_stack_max_itm_);
		a_.node_stack_max_cap_ = std::max(a_.node_stack_max_cap_,
				gather_buf_[i].node_stack_max_cap_);
		a_.give_stack_max_cap_ = std::max(a_.give_stack_max_cap_,
				gather_buf_[i].give_stack_max_cap_);

		a_.cleared_tasks_ += gather_buf_[i].cleared_tasks_;
	}
}

}
