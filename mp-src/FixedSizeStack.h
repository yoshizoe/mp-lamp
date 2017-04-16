/*
 * FixedSizeStack.h
 *
 *  Created on: Apr 16, 2017
 *      Author: yuu
 */

#ifndef MP_SRC_FIXEDSIZESTACK_H_
#define MP_SRC_FIXEDSIZESTACK_H_

class FixedSizeStack {
public:
	FixedSizeStack(int capacity) {
		size_ = 0;
		stack_ = new int[capacity];
	}
	~FixedSizeStack() {
		delete[] stack_;
	}
	void Push(int num) {
		stack_[size_] = num;
		size_++;
	}
	int Pop() {
		size_--;
		return stack_[size_];
	}
	int Size() const {
		return size_;
	}

	void Clear() {
		size_ = 0;
	}

	friend std::ostream& operator<<(std::ostream & out,
			const FixedSizeStack & st);

private:
	int size_;
	int * stack_;
};



#endif /* MP_SRC_FIXEDSIZESTACK_H_ */
