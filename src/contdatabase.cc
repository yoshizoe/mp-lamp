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
#include <boost/math/distributions/chi_squared.hpp>
using namespace std;

namespace lamp_search {

ContDatabase::ContDatabase() :
		nu_transactions_(0), nu_pos_total_(0), nu_items_(0) {
	// TODO Auto-generated constructor stub

}

ContDatabase::~ContDatabase() {
	// TODO Auto-generated destructor stub
}

ContDatabase::ContDatabase(std::istream& features, std::istream& classes) :
		nu_transactions_(0), nu_pos_total_(0), nu_items_(0) {
	if (!features.good() || !classes.good()) {
		assert(false && "File not found");
		exit(1);
	}
	printf("ContDatabase\n");
	readFromCSV(features);
	readClassFromCSV(classes);
//	ShowInfo();
}

// read a database file
void ContDatabase::readFromCSV(istream& ifs, int dim_limit, bool reverse) {
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
		transposed.push_back(tmp);
	}

	// TODO: Transpose the feature table
	nu_transactions_ = transposed.size();
	nu_items_ = transposed[0].size();

	vector<vector<Ftype>> tmp(nu_items_, vector<Ftype>(nu_transactions_, 0.0));
	for (int i = 0; i < nu_items_; ++i) {
		for (int j = 0; j < nu_transactions_; ++j) {
			tmp[i][j] = transposed[j][i];
		}
	}
//	features = tmp;
	// Convert into Ranking.
	features = vector<vector<Ftype>>(nu_items_,
			vector<Ftype>(nu_transactions_, 0.0));
	vector<size_t> idx(nu_transactions_);
	for (int i = 0; i < nu_items_; ++i) {
		vector<Ftype> freqs = tmp[i];
		// initialize index vector
		iota(idx.begin(), idx.end(), 0);
		// sort indexes based on comparing values in v
		sort(idx.begin(), idx.end(),
				[&freqs](int i1, int i2) {return freqs[i1] < freqs[i2];});
		for (int j = 0; j < nu_transactions_; j++) {
//			rank[i][j] = idx[j] + 1;
			features[i][j] = (double) (idx[j] + 1) / (double) nu_transactions_;
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
		// TODO: Which class should be the positive?
		if (cls == 1) {
			++nu_pos_total_;
		}
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

std::vector<int> ContDatabase::GetChildren(std::vector<int> items) const {
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
		std::vector<int> itemset_items) const {
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
		std::vector<Ftype> itemset_freqs) const {
	Ftype freq = 0;
	for (int j = 0; j < nu_transactions_; ++j) {
		freq += itemset_freqs[j];
	}
	return freq / (double) nu_transactions_;
}

ContDatabase::Ftype ContDatabase::GetPositiveFreq(
		std::vector<Ftype> itemset_freqs) const {
	Ftype freq = 0;
	for (int j = 0; j < nu_transactions_; ++j) {
		if (classes[j] == 1) { // TODO: which is positive? refactor.
			freq += itemset_freqs[j];
		}
	}
	return freq / (double) nu_transactions_;
}

ContDatabase::Ftype ContDatabase::GetFreq(
		std::vector<int> itemset_items) const {
	return GetFreq(GetFreqArray(itemset_items));
}

std::vector<ContDatabase::Ftype> ContDatabase::GetChildrenFreq(
		std::vector<Ftype>& parent_freq, int new_item) const {
	assert(parent_freq.size() == nu_transactions_);
	std::vector<Ftype> child(nu_transactions_);
	for (int j = 0; j < nu_transactions_; ++j) {
		child[j] = parent_freq[j] * features[new_item][j];
	}
	return child;
}

// TODO: This function is awfully complicated like a spagetti.
bool ContDatabase::PPCExtension(VariableLengthItemsetStack * st, int* parent,
		int new_item, int* child) const {
	// TODO: Prune a node if it is not PPCExtension.
	st->CopyItem(parent, child);
	st->PushOneItem(new_item);
	return true;
}

double ContDatabase::CalculatePValue(Ftype total_freqs, Ftype pos_freq) const {
	assert(0 <= pos_freq);
	assert(pos_freq <= total_freqs);
	assert(total_freqs <= 1.0);
	double k = kl(total_freqs, pos_freq);
	return computePvalue(k, nu_transactions_);
}

double ContDatabase::CalculatePMin(Ftype total_freqs) const {
	assert(0.0 <= total_freqs && total_freqs <= 1.0);
	double kl = kl_max_fast(total_freqs, nu_transactions_ - nu_pos_total_,
			nu_transactions_);
	double pmin = computePvalue(kl, nu_transactions_);
//	printf("frequency = %.2f, kl = %.2f, pvalue = %.2f\n",
//			total_freqs, kl, pmin);
	assert(0.0 <= kl);
	assert(0.0 <= pmin && pmin <= 1.00);
	return pmin;
}

double ContDatabase::CalculatePLowerBound(Ftype total_freqs) const {
	assert(0.0 <= total_freqs && total_freqs <= 1.0);
	double kl = kl_max_fast_bound(total_freqs, nu_transactions_ - nu_pos_total_,
			nu_transactions_);
	double pmin = computePvalue(kl, nu_transactions_);
//	printf("frequency = %.2f, kl = %.2f, pvalue = %.2f\n",
//			total_freqs, kl, pmin);
	assert(0.0 <= kl);
	assert(0.0 <= pmin && pmin <= 1.00);
	return pmin;
}

double ContDatabase::CalculatePValue(std::vector<int>& itemset_items) const {
	std::vector<Ftype> freqs = GetFreqArray(itemset_items);
	Ftype tot_freqs = 0;
	Ftype pos_freqs = 0;

	for (int j = 0; j < nu_transactions_; ++j) {
		tot_freqs += freqs[j];
		if (classes[j] == 1) { // TODO: which is positive? refactor.
			pos_freqs += freqs[j];
		}
	}
	tot_freqs /= (double) nu_transactions_;
	pos_freqs /= (double) nu_transactions_;
	return CalculatePValue(tot_freqs, pos_freqs);
}

/**
 * Things I have to understand
 */
// Sugiyama's code
// compute p-value
double ContDatabase::computePvalue(double kl, int N) const {
	assert(0 <= kl);
	assert(0 < N);
	boost::math::chi_squared chisq_dist(1);
	// else pval = 1 - boost::math::cdf(chisq_dist, 2 * (double)N * kl);
	// if (pval > 1) pval = 1.0;
	// if (VERBOSE) cout << "kl: " << kl << endl;
	double pval = 0.0;
	if (kl <= pow(10, -8))
		pval = 1.0;
	else
		pval = 1.0 - boost::math::cdf(chisq_dist, 2 * (double) N * kl);
	return pval;
}

// TODO: ???
// Sugiyama's code
double ContDatabase::kl_max_fast(double freq, int N0, int N) const {
	double r0 = (double) N0 / (double) N;
	if (freq < r0) {
		double d = freq * log(1.0 / r0)
				+ (r0 - freq) * log((r0 - freq) / (r0 - r0 * freq))
				+ (1.0 - r0) * log(1.0 / (1.0 - freq));
		assert(0 <= d);
		return d;
	} else {
		double d = r0 * log(1.0 / freq)
				+ (freq - r0) * log((freq - r0) / (freq - freq * r0))
				+ (1.0 - freq) * log(1.0 / (1.0 - r0));
		assert(0 <= d);
		return d;
	}
}

double ContDatabase::kl_max_fast_bound(double freq, int N0, int N) const {
	double r0 = (double) N0 / (double) N;
//	printf("N0=%d, N=%d, r0=%.2f\n", r0);
	assert(0 <= r0 && r0 <= 1.0);
	assert(0 <= freq && freq <= 1.0);
	if (freq < r0) {
		double d = freq * log(1.0 / r0)
				+ (r0 - freq) * log((r0 - freq) / (r0 - r0 * freq))
				+ (1.0 - r0) * log(1.0 / (1.0 - freq));
		assert(0.0 <= d);
		return d;
	} else {
		double d = r0 * log(1.0 / r0) + (1.0 - r0) * log(1.0 / (1.0 - r0));
		assert(-0.001 <= d);
		return d;
	}
}

// TODO: Why kl < 0?
double ContDatabase::kl(double total_freq, double pos_freq) const {
	double r1 = (double) nu_pos_total_ / (double) nu_transactions_;
	double r0 = (double) (nu_transactions_ - nu_pos_total_)
			/ (double) nu_transactions_;
	double neg_freq = total_freq - pos_freq;

	assert(0 <= pos_freq);
	assert(pos_freq <= total_freq);
	assert(total_freq <= 1.0);
	assert(0.0 <= r0 && r0 <= 1.0);
	assert(0.0 <= r1 && r1 <= 1.0);
	assert(0 <= neg_freq && neg_freq <= total_freq);

	vector<double> po; // { neg_freq, pos_freq, r0 - neg_freq, r1 - pos_freq };
	po.push_back(neg_freq);
	po.push_back(pos_freq);
	po.push_back(r0 - neg_freq);
	po.push_back(r1 - pos_freq);

	vector<double> pe;
	pe.push_back(r0 * total_freq);
	pe.push_back(r1 * total_freq);
	pe.push_back(r0 - r0 * total_freq);
	pe.push_back(r1 - r1 * total_freq);

	double kl = 0.0;
	for (int i = 0; i < po.size(); ++i) {
		assert(0 <= pe[i]);
		assert(0 <= po[i]);
		kl += po[i] * log(po[i] / pe[i]);
	}
//	if (!(0.0 <= kl)) {
//		for (int i = 0; i < po.size(); ++i) {
//			printf(
//					"po[%d] = %.4f, pe[%d] = %.4f, po * log(po/pe) = %.4f\n",
//					i, po[i], i, pe[i], po[i] * log(po[i] / pe[i]));
//		}
//		printf("kl = %.4f\n", kl);
//		kl = 0.0;
//	}
	assert(
			-0.0001 <= kl || (fprintf(stderr, "kl=%f, total_freq=%f, pos_freq=%f\n", kl, total_freq, pos_freq) && 0));
	return kl;
}

void ContDatabase::ShowInfo() const {
	printf("#Database: #Trans= %d #Items= %d #PosTrans= %d\n", nu_transactions_,
			nu_items_, nu_pos_total_);
	assert(0 < nu_transactions_);
	assert(0 < nu_items_);
	assert(0 < nu_pos_total_);
	assert(nu_pos_total_ < nu_transactions_);

//	printf("features:\n");
//	for (int i = 0; i < nu_items_; i++) {
//		for (int j = 0; j < nu_transactions_; ++j) {
//			printf("%.2f ", features[i][j]);
//		}
//		printf("\n");
//	}
//	printf("classes:\n");
//	for (int i = 0; i < nu_transactions_; i++) {
//		printf("%d ", classes[i]);
//	}
//	printf("\n");
//
//	for (int i = 0; i < nu_items_; i++) {
//		for (int j = 0; j < nu_transactions_; ++j) {
//			assert(0.0 <= features[i][j]);
//			assert(features[i][j] <= 1.001);
//		}
//	}
}

} /* namespace lamp_search */
