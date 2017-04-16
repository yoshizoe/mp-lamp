/*
 * contdatabase.cpp
 *
 *  Created on: Apr 15, 2017
 *      Author: yuu
 */

#include "contdatabase.h"
#include <iostream>
#include <assert.h>
using namespace std;

namespace lamp_search {

ContDatabase::ContDatabase() {
	// TODO Auto-generated constructor stub

}

ContDatabase::~ContDatabase() {
	// TODO Auto-generated destructor stub
}

// read a database file
void ContDatabase::readFromCSV(istream& ifs, int dim_limit, bool reverse) {
	assert(features.empty());
	assert(typeid(Feature) == typeid(double) && "Feature is not double: need to refactor.");
	string line;
	while (getline(ifs, line)) {
		stringstream lineStream(line);
		string cell;
		vector<Feature> tmp;
		while (getline(lineStream, cell, ',') && tmp.size() <= dim_limit) {
			tmp.push_back(stod(cell)); // TODO: template stod
		}
		// YJ: Add reversed (negated) values for all features so that the method can detect negative correlation.
		if (reverse) {
			int size = tmp.size();
			for (int i = 0; i < size; ++i) {
				tmp.push_back(-tmp[i]);
			}
		}
		features.push_back(tmp);
	}
	nu_transactions_ = features.size();
	nu_items_ = features[0].size();
}

void ContDatabase::readClassFromCSV(istream& ifs) {
	assert(typeid(Feature) == typeid(int) && "Class is not int: need to refactor.");
	string line;
	while (getline(ifs, line)) {
		int cls = stoi(line);
		classes.push_back(cls); // TODO
		++nu_pos_total_;
	}
}

} /* namespace lamp_search */
