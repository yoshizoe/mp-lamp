/*
 * contdatabase.h
 *
 *  Created on: Apr 15, 2017
 *      Author: yuu
 */

#ifndef SRC_CONTDATABASE_H_
#define SRC_CONTDATABASE_H_

#include <vector>

namespace lamp_search {

// TODO: Current implementation of database<Block> is pigeon-holing for optimization.
// This made me make a new class which makes things increasingly complex...

class ContDatabase {
	// TODO: template<typename Feature, typename Class>
	typedef double Feature;
	typedef int Class;
public:
	ContDatabase();
	virtual ~ContDatabase();
//	void readFromCSV(std::ifstream& ifs);
	void readFromCSV(std::istream& ifs, int dim_limit =
			std::numeric_limits<int>::max(), bool reverse = false);
	void readClassFromCSV(std::istream& ifs);

	int NuItems() const {
		return nu_items_;
	}
	int NuTransaction() const {
		return nu_transactions_;
	}
	int PosTotal() const {
		return nu_pos_total_;
	}
private:
	std::vector<std::vector<Feature>> features; // TODO: align values for each feature continuously.
	std::vector<Class> classes;

	// auxilary
	int nu_items_;
	int nu_transactions_;
	int nu_pos_total_;
};

} /* namespace lamp_search */

#endif /* SRC_CONTDATABASE_H_ */
