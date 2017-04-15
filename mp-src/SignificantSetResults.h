/*
 * SignificantSetResults.h
 *
 *  Created on: Apr 14, 2017
 *      Author: yuu
 */

#ifndef MP_SRC_SIGNIFICANTSETRESULTS_H_
#define MP_SRC_SIGNIFICANTSETRESULTS_H_

namespace lamp_search {
class SignificantSetResult {
public:
	SignificantSetResult(double p, int * s, int nu_sup, int nu_pos,
			VariableLengthItemsetStack * ss) :
			pval_(p), set_(s), sup_num_(nu_sup), pos_sup_num_(nu_pos), ss_(ss) {
	}

	double pval_;
	int * set_;
	int sup_num_;
	int pos_sup_num_;

	const VariableLengthItemsetStack * ss_;
};

struct sigset_compare {
	bool operator()(const SignificantSetResult & lhs,
			const SignificantSetResult & rhs) {
		if (lhs.pval_ < rhs.pval_)
			return true;
		else if (lhs.pval_ > rhs.pval_)
			return false;
		else {
			int l_item_num = lhs.ss_->GetItemNum(lhs.set_);
			int r_item_num = rhs.ss_->GetItemNum(rhs.set_);

			if (l_item_num > r_item_num)
				return true;
			else if (l_item_num < r_item_num)
				return false;
			else {
				// sort based on dictionary order of item
				int n = l_item_num;
				for (int i = 0; i < n; i++) {
					int l_item = lhs.ss_->GetNthItem(lhs.set_, i);
					int r_item = rhs.ss_->GetNthItem(rhs.set_, i);
					if (l_item < r_item)
						return true;
					else if (l_item > r_item)
						return false;
				}
				throw std::runtime_error("identical duplicate itemsets found");
				return false;
			}
		}
		return false;
	}
};

}
#endif /* MP_SRC_SIGNIFICANTSETRESULTS_H_ */
