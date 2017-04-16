/*
 * ParallelSearch.h
 *
 *  Created on: Apr 16, 2017
 *      Author: yuu
 */

#ifndef MP_SRC_PARALLELSEARCH_H_
#define MP_SRC_PARALLELSEARCH_H_

#include "DTD.h"
#include "Log.h"

namespace lamp_search {

class ParallelSearch {
public:
	ParallelSearch();
	virtual ~ParallelSearch();
	void Search();

private:
	void ProcessNode();
	void Probe();
	void Distribute();
    void Reject(); // distribute finished, reject remaining requests
    void CheckCSThreshold();
    void Steal(); // request steal


    DTD dtd_;

    Log log_;
};

} /* namespace lamp_search */

#endif /* MP_SRC_PARALLELSEARCH_H_ */
