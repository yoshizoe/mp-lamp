/*
 * contdatabase.h
 *
 *  Created on: Apr 15, 2017
 *      Author: yuu
 */

#ifndef SRC_CONTDATABASE_H_
#define SRC_CONTDATABASE_H_

#include <vector>
#include <stdio.h>
#include <istream>
#include <limits>
#include "variable_length_itemset.h"

namespace lamp_search {

// TODO: Current implementation of database<Block> is pigeon-holing for optimization.
// This made me make a new class which makes things increasingly complex...

class ContDatabase {
	// TODO: template<typename Feature, typename Class>
	typedef double Ftype;
	typedef int Ctype;
public:
	ContDatabase();
	ContDatabase(std::istream& features, std::istream& classes);
	virtual ~ContDatabase();
	void readFromCSV(std::istream& ifs, int dim_limit =
			std::numeric_limits<int>::max(), bool reverse = false);
	void readClassFromCSV(std::istream& ifs);

	// TODO: implement
	std::vector<int> GetChildren(int itemset_id); // TODO: Should this be in database?
	std::vector<int> GetChildren(std::vector<int> items); // TODO: Should this be in database?
	// TODO: implement
	std::vector<Ftype> GetFreqArray(std::vector<int> itemset_item);
	// TODO: implement
	Ftype GetFreq(std::vector<int> itemset_items);
	// TODO: implement
	Ftype GetFreq(std::vector<Ftype> itemset_freqs);
	// TODO: implement
	std::vector<Ftype> GetChildrenFreq(
			std::vector<Ftype>& parent_freq, int new_item);
	bool PPCExtension(VariableLengthItemsetStack * st, int* parent,
			int new_item, int* child);
	double CalculatePValue(Ftype sup_num,
			std::vector<Ftype>& itemset_freqs);

	int NumItems() const {
		return nu_items_;
	}
	int NumTransactions() const {
		return nu_transactions_;
	}
	int NumPositiveItems() const {
		return nu_pos_total_;
	}

	void ShowInfo();
private:
	// TODO: This should be transposed for better memory access.
	std::vector<std::vector<Ftype> > features; // TODO: align values for each feature continuously.
	std::vector<Ctype> classes;

	// auxilary
	int nu_items_;
	int nu_transactions_;
	int nu_pos_total_;
};

} /* namespace lamp_search */

#endif /* SRC_CONTDATABASE_H_ */
