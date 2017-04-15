/*
 * StealState.h
 *
 *  Created on: Apr 14, 2017
 *      Author: yuu
 */

#ifndef MP_SRC_STEALSTATE_H_
#define MP_SRC_STEALSTATE_H_

namespace lamp_search {
class StealState {
public:
	enum {
		RANDOM = 0, LIFELINE,
	};

	StealState(int w, int z) :
			w_(w), z_(z) {
		Init();
	}

	int w_;
	int z_;

	void Init() {
		steal_phase_started_ = true;

		random_counter_ = 0;
		lifeline_counter_ = 0;
		next_lifeline_victim_ = 0;

		if (w_)
			state_ = RANDOM;
		else
			state_ = LIFELINE;
		requesting_ = false;
	}

	void Finish() {
		steal_phase_started_ = false;

		random_counter_ = 0;
		lifeline_counter_ = 0;

		next_lifeline_victim_ = 0;
		// this is currently needed to resset lifeline victim
		// when lifeline has -1 (if p_ is not powers of 2)

		if (w_)
			state_ = RANDOM;
		else
			state_ = LIFELINE;
		requesting_ = false;
	}

	void ResetCounters() {
		random_counter_ = 0;
		lifeline_counter_ = 0;
	}

	void IncRandomCount() {
		random_counter_++;
		if (random_counter_ >= w_)
			random_counter_ = 0;
	}
	int RandomCount() const {
		return random_counter_;
	}

	void IncLifelineCount() {
		lifeline_counter_++;
	}
	int LifelineCounter() const {
		return lifeline_counter_;
	}

	void IncLifelineVictim() {
		next_lifeline_victim_++;
		if (next_lifeline_victim_ >= z_)
			next_lifeline_victim_ = 0;
	}
	int LifelineVictim() const {
		return next_lifeline_victim_;
	}

	void SetState(int state) {
		state_ = state;
	}
	int State() const {
		if (w_ == 0)
			return LIFELINE;
		else
			return state_;
	}

	bool Requesting() const {
		return requesting_;
	}
	void SetRequesting(bool flag = true) {
		requesting_ = flag;
	}
	void ResetRequesting() {
		requesting_ = false;
	}

	bool StealStarted() const {
		return steal_phase_started_;
	}
	void SetStealStart() {
		steal_phase_started_ = true;
	}
	void ResetStealStart() {
		steal_phase_started_ = true;
	}

	bool steal_phase_started_;

	int random_counter_; // random victim counter
	int lifeline_counter_; // lifeline victim counter
	int next_lifeline_victim_; // next lifeline victim range [0..z-1]

	int state_; // 0: random steal, 1: lifeline steal
	bool requesting_; // true if requesting. become false if rejected or given
	// after sending request, wait until give or reject
};

}

#endif /* MP_SRC_STEALSTATE_H_ */
