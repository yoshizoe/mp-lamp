/*
 * Log.h
 *
 *  Created on: Apr 14, 2017
 *      Author: yuu
 */

#ifndef MP_SRC_LOG_H_
#define MP_SRC_LOG_H_

#include <vector>

namespace lamp_search {
class Log {
public:
	Log();
	~Log();
	void Init();

	long long int idle_start_;

	void GatherLog(int nu_proc);

	//--------
	// periodic log

	void InitPeriodicLog();
	void StartPeriodicLog();
	void FinishPeriodicLog();

	long long int periodic_log_start_;
	long long int next_log_time_in_second_;

	// show this separately for phase_ 1 and 2
	void TakePeriodicLog(long long int capacity, int lambda, int phase);
	struct PeriodicLog_T {
		void Clear() {
			seconds_ = 0ll;
			capacity_ = 0ll;
			lambda_ = 0;
			phase_ = 0;
		}
		long long int seconds_;
		long long int capacity_;
		int lambda_;
		int phase_;
	};
	std::vector<PeriodicLog_T> plog_;

	void AggregatePLog(int nu_proc);
	int sec_max_;
	PeriodicLog_T * plog_buf_; // gather results for copying local plog
	PeriodicLog_T * plog_gather_buf_; // gather plog of all ranks here

	//--------
	// other log

	struct LogData {
		void Init();

		long long int search_start_time_;
		long long int search_finish_time_;

		long long int iprobe_num_;
		long long int iprobe_time_;
		long long int iprobe_time_max_;

		long long int iprobe_succ_num_;
		long long int iprobe_succ_time_;
		long long int iprobe_succ_time_max_;

		long long int iprobe_fail_num_;
		long long int iprobe_fail_time_;
		long long int iprobe_fail_time_max_;

		long long int probe_num_;
		long long int probe_time_;
		long long int probe_time_max_;

		long long int recv_num_;
		long long int recv_time_;
		long long int recv_time_max_;

		long long int bsend_num_;
		long long int bsend_time_;
		long long int bsend_time_max_;

		long long int bcast_num_;
		long long int bcast_time_;
		long long int bcast_time_max_;

		// long long int accum_task_time_;
		// long long int accum_task_num_;
		// long long int basic_task_time_;
		// long long int basic_task_num_;
		// long long int control_task_time_;
		// long long int control_task_num_;

		long long int dtd_phase_num_;
		double dtd_phase_per_sec_;
		long long int dtd_accum_phase_num_;
		double dtd_accum_phase_per_sec_;
		//long long int dtd_reply_num_;

		// long long int accum_phase_num_; // completed accum phase num
		// double accum_phase_per_sec_; // completed accum phase per second

		long long int lifeline_steal_num_;
		long long int lifeline_nodes_received_;
		long long int steal_num_;
		long long int nodes_received_;

		long long int lifeline_given_num_;
		long long int lifeline_nodes_given_;
		long long int given_num_;
		long long int nodes_given_;

		long long int process_node_num_;
		long long int process_node_time_;

		long long int preprocess_time_;

		long long int idle_time_;

		long long int pval_table_time_;

		long long int node_stack_max_itm_;
		long long int give_stack_max_itm_;

		long long int node_stack_max_cap_;
		long long int give_stack_max_cap_;

		long long int cleared_tasks_;

	};
	LogData d_;

	LogData * gather_buf_;
	LogData a_; // aggregated
	void Aggregate(int nu_proc);
};
}

#endif /* MP_SRC_LOG_H_ */
