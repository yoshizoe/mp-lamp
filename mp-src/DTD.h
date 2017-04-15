/*
 * TerminateDectection.h
 *
 *  Created on: Apr 14, 2017
 *      Author: yuu
 */

#ifndef MP_SRC_DTD_H_
#define MP_SRC_DTD_H_

namespace lamp_search {

// Mattern 1987
// bounded clock-counter version of TIME algorithm
class DTD {
public:
	DTD(int branch) :
			terminated_(false), count_(0), time_zone_(0), time_warp_(false), not_empty_(
					false), k_echo_tree_branch(branch) {
		accum_flag_ = new bool[k_echo_tree_branch];
	}

	~DTD() {
		delete[] accum_flag_;
	}

	void Init() {
		terminated_ = false;

		count_ = 0;
		time_zone_ = 0;
		time_warp_ = false;
		not_empty_ = false;

		// for reduce
		//requesting_ = false;
		ClearAccumFlags();
		ClearReduceVars();
	}

	void CheckPoint() {
		terminated_ = false;
		//assert(requesting_ == false);
		for (int i = 0; i < k_echo_tree_branch; i++)
			assert(accum_flag_[i] == false);
	}

	bool terminated_;

	int count_;
	int time_zone_;
	bool time_warp_;
	bool not_empty_;
	static const int kNuTimezone = 1024 * 1024;

	int k_echo_tree_branch;

	// on receiving basic message
	void UpdateTimeZone(int ts) {
		time_warp_ = time_warp_ || (ts == (time_zone_ + 1) % kNuTimezone);
	}

	void OnSend() {
		count_++;
	}
	void OnRecv() {
		count_--;
	}

	// on receiving control messages
	void ResetTimeWarp() {
		time_warp_ = false;
	}
	void IncTimeZone() {
		time_zone_++;
		time_zone_ = time_zone_ % kNuTimezone;
	}

	// for reduce
	//bool requesting_;
	bool * accum_flag_;

	int reduce_count_;
	bool reduce_time_warp_;
	bool reduce_not_empty_;

	void Reduce(int count, bool time_warp, bool not_empty) {
		reduce_count_ += count;
		reduce_time_warp_ = reduce_time_warp_ || time_warp;
		reduce_not_empty_ = reduce_not_empty_ || not_empty;
	}

	void ClearAccumFlags() {
		for (int i = 0; i < k_echo_tree_branch; i++)
			accum_flag_[i] = false;
	}

	void ClearReduceVars() {
		reduce_count_ = 0;
		reduce_time_warp_ = false;
		reduce_not_empty_ = false;
	}

};

} /* namespace search */

#endif /* MP_SRC_DTD_H_ */
