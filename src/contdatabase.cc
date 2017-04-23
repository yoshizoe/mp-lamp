/*
 * contdatabase.cpp
 *
 *  Created on: Apr 15, 2017
 *      Author: yuu
 */

#include "contdatabase.h"
#include <iostream>
#include <assert.h>
#include <string>
#include <algorithm>
#include <sstream>
using namespace std;

namespace lamp_search {

ContDatabase::ContDatabase() {
	// TODO Auto-generated constructor stub

}

ContDatabase::~ContDatabase() {
	// TODO Auto-generated destructor stub
}

ContDatabase::ContDatabase(std::istream& features,
		std::istream& classes) {
	readFromCSV(features);
	readClassFromCSV(classes);
}

// read a database file
void ContDatabase::readFromCSV(istream& ifs, int dim_limit,
		bool reverse) {
	std::vector<std::vector<Ftype>> transposed;
	assert(features.empty());
	assert(
			typeid(Ftype) == typeid(double)
					&& "Feature is not double: need to refactor.");
	string line;
	while (getline(ifs, line)) {
		stringstream lineStream(line);
		string cell;
		vector<Ftype> tmp;
		while (getline(lineStream, cell, ',')
				&& tmp.size() <= dim_limit) {
			tmp.push_back(stod(cell)); // TODO: template stod
		}
		// YJ: Add reversed (negated) values for all features so that the method can detect negative correlation.
		if (reverse) {
			int size = tmp.size();
			for (int i = 0; i < size; ++i) {
				tmp.push_back(-tmp[i]);
			}
		}
		transposed.push_back(tmp);
	}

	// TODO: Transpose the feature table
	nu_transactions_ = transposed.size();
	nu_items_ = transposed[0].size();

	features.resize(nu_items_);
	for (int i = 0; i < nu_items_; ++i) {
		features[i].resize(nu_transactions_);
	}
	for (int i = 0; i < nu_items_; ++i) {
		for (int j = 0; j < nu_transactions_; ++j) {
			features[i][j] = transposed[j][i];
		}
	}

}

void ContDatabase::readClassFromCSV(istream& ifs) {
	assert(
			typeid(Ctype) == typeid(int)
					&& "Class is not int: need to refactor.");
	string line;
	while (getline(ifs, line)) {
		int cls = stoi(line);
		classes.push_back(cls); // TODO
		++nu_pos_total_;
	}
	assert(classes.size() == nu_transactions_);
}

// TODO: inefficient
//std::vector<int> ContDatabase::GetChildren(int itemset_id) {
//	std::vector<int> items(nu_items_);
//	for (int i = 0; i < nu_items_; ++i) {
//		items[i] = i;
//	}
//	return items;
//}

std::vector<int> ContDatabase::GetChildren(std::vector<int> items) {
	if (items.empty()) {
		// Expanding root node.
		// TODO: Generate only a responsible nodes for the root.
	}

	// Assert that items are sorted.
	if (items.size() >= 2) {
		for (int i = 0; i < items.size() - 1; ++i) {
			assert(items[i] < items[i + 1]);
		}
	}
	std::vector<int> children;
	int last = 0;
	if (!items.empty()) {
		last = items.back();
		++last;
	}
	for (int i = last; i < nu_items_; ++i) {
		children.push_back(i);
	}
	return children;
}

std::vector<ContDatabase::Ftype> ContDatabase::GetFreqArray(
		std::vector<int> itemset_items) {
	std::vector<ContDatabase::Ftype> freqs(nu_transactions_, 1.0);
	for (int i = 0; i < itemset_items.size(); i++) {
		int item = itemset_items[i];
		// TODO: how is features aligned? Shouldn't it be transposed?
		for (int j = 0; j < nu_transactions_; ++j) {
			freqs[j] = freqs[j] * features[item][j];
		}
	}
	return freqs;
}

ContDatabase::Ftype ContDatabase::GetFreq(
		std::vector<Ftype> itemset_freqs) {
	Ftype freq = 0;
	for (int j = 0; j < nu_transactions_; ++j) {
		freq += itemset_freqs[j];
	}
	return freq;
}

ContDatabase::Ftype ContDatabase::GetFreq(
		std::vector<int> itemset_items) {
	return GetFreq(GetFreqArray(itemset_items));
}

std::vector<ContDatabase::Ftype> ContDatabase::GetChildrenFreq(
		std::vector<Ftype>& parent_freq, int new_item) {
	assert(parent_freq.size() == nu_transactions_);
	std::vector<Ftype> child(nu_transactions_);
	for (int j = 0; j < nu_transactions_; ++j) {
		child[j] = parent_freq[j] * features[new_item][j];
	}
	return child;
}

// TODO: This function is awfully complicated like a spagetti.
bool ContDatabase::PPCExtension(VariableLengthItemsetStack * st, int* parent,
		int new_item, int* child) {
	// TODO: Prune a node if it is not PPCExtension.
	st->CopyItem(parent, child);
	st->PushOneItem(new_item);
	return true;
}

double ContDatabase::CalculatePValue(Ftype total_freqs,
		std::vector<Ftype>& itemset_freqs) {
	// Calculate PValue...
	Ftype positive_freqs = 0;
//	int positive items
	for (int j = 0; j < nu_transactions_; ++j) {
		if (classes[j] == 1) {
			positive_freqs += itemset_freqs[j];
		}
	}
	// TODO: calculate P value using total_freqs and positive_freqs.
	// TODO: For testing let's assume everything is significant.
	return (double) positive_freqs;
}

void ContDatabase::ShowInfo() {
	printf("features:\n");
	for (int i = 0; i < nu_items_; i++) {
		for (int j = 0; j < nu_transactions_; ++j) {
			printf("%.2f ", features[i][j]);
		}
		printf("\n");
	}
	printf("classes:\n");
	for (int i = 0; i < nu_transactions_; i++) {
		printf("%d ", classes[i]);
	}
	printf("\n");
}

} /* namespace lamp_search */
