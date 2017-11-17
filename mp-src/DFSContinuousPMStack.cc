//
// Created by Takuo Doi on 2017/11/05.
//

#include <cstdio>
#include <iostream>
#include "gflags/gflags.h"

#include "DFSContinuousPMStack.h"
#include "ParallelContinuousPM.h"
#include "MPI_Data.h"

namespace lamp_search {

    int DFSContinuousPMStack::Process(int phase, Timer* timer, Log* log) {

        if (treesearch_data_->node_stack_->Empty()) {
            return 0;
        }

        long long int start_time, lap_time;
        start_time = timer->Elapsed();
        lap_time = start_time;

        mpi_data_.processing_node_ = true;
        int processed = 0;
        while (!treesearch_data_->node_stack_->Empty()) {
            processed++;
            expand_num_++;

            DBG( D(3) << "Stack has " << treesearch_data_->node_stack_->NuItemset() << " itemsets" );
            if (!PopNodeFromStack()) {
                // Discard a node if pmin(node) > pmin_threshould_.
                continue;
            }
//		printf("Stack has %d itemsets\n", treesearch_data->node_stack_->NuItemset());
            int n = treesearch_data_->node_stack_->GetItemNum(treesearch_data_->itemset_buf_);
//		int* init = treesearch_data->node_stack_->GetItemArray(treesearch_data->itemset_buf_);
//		std::vector<int> items(init, init + n);
//		printf("Poped %d\n", n);
//		PrintItemset(treesearch_data->itemset_buf_, sup_buf_);

            bool isRoot = (n == 0);
            int accum_period_counter_ = 0;

            // TODO: Implement ContTable to return list of children.
//	        int core_i = g_->CoreIndex(*treesearch_data->node_stack_, treesearch_data->itemset_buf_);
            std::vector<int> current_items = treesearch_data_->node_stack_->getItems(treesearch_data_->itemset_buf_);
            std::vector<int> children = GetChildren(current_items);

            for (std::vector<int>::iterator it = children.begin();
                 it != children.end(); ++it) {
                int new_item = *it;
                assert(0 <= new_item && new_item < d_->NumItems());

                // TODO: no need for duplicate detection
                if (treesearch_data_->node_stack_->Exist(treesearch_data_->itemset_buf_, new_item)) {
                    printf("Node Pruned because of duplicated items");
                    return 0;
                }
                // skip existing item
                // todo: improve speed here
                if (!isRoot) {
                    CheckProbe(accum_period_counter_, lap_time, timer, log);
                }

                // TODO TOPK: This should check if the pvalue to the current top-k pvalue.
                if (!TestAndPushNode(new_item)) {
                    continue;
                }

                int *child_node_buf = treesearch_data_->node_stack_->Top();

//			std::vector<int> child = treesearch_data->node_stack_->getItems(ppc_ext_buf);
//			printf("Pushed %d items:", treesearch_data->node_stack_->GetItemNum(ppc_ext_buf));
//			PrintItemset(ppc_ext_buf, child_sup_buf_);

//			int sup_num = treesearch_data->node_stack_->GetSup(ppc_ext_buf);
//			std::vector<int> items = treesearch_data->node_stack_->getItems(ppc_ext_buf);
//			d_->GetFreq()
//			double pp

                DBG(if (phase == 2) {
                    D(3) << "found cs ";
                    treesearch_data_->node_stack_->Print(D(3), child_node_buf);
                });

                Feature freq = algorithm_->GetFreqFromDatabase();

                // TODO TOPK: Store pvalue instead of frequency.
                ProcessNode(phase, freq, child_node_buf);

//			assert(freq >= pmin_thre_);

                // try skipping if supnum_ == sup_threshold,
                // because if sup_num of a node equals to sup_threshold, children will have smaller sup_num
                // therefore no need to check it's children
                // note: skipping node_stack_ full check. allocate enough memory!
//			if (freq == pmin_thre_) {
//				printf("All Children Pruned");
//				treesearch_data->node_stack_->Pop();
//			}
            }

            if (CheckProcessNodeEnd(mpi_data_.granularity_, mpi_data_.isGranularitySec_, processed, start_time, timer)) {
                break;
            }
        }

        long long int elapsed_time = timer->Elapsed() - lap_time;
        log->d_.process_node_time_ += elapsed_time;
        log->d_.process_node_num_ += processed;

        mpi_data_.processing_node_ = false;

        log->d_.node_stack_max_itm_ = std::max(log->d_.node_stack_max_itm_, (long long int) (treesearch_data_->node_stack_->NuItemset()));
        log->d_.node_stack_max_cap_ = std::max(log->d_.node_stack_max_cap_, (long long int) (treesearch_data_->node_stack_->UsedCapacity()));

        DBG(
                D(2) << "processed node num=" << processed << "\ttime=" << elapsed_time << std::endl;
        );

        return processed;
    }

    int DFSContinuousPMStack::Split(VariableLengthItemsetStack* give_stack) {
        return treesearch_data_->node_stack_->Split(give_stack);
    }

    void DFSContinuousPMStack::Merge(VariableLengthItemsetStack* give_stack, int count) {
        treesearch_data_->node_stack_->MergeStack(
                give_stack->Stack() + VariableLengthItemsetStack::SENTINEL + 1,
                count - VariableLengthItemsetStack::SENTINEL - 1);
    }

    /**
     * itemset_buf_ := pointer to the index of itemset
     * sup_buf_     := support* of the itemset (takes AND for each item)
     * support*     := transactions which include all items in the itemset
     */
    bool DFSContinuousPMStack::AfterPopNodeFromStack() {
        /**
         * Get the support (frequency) of the item into sup_buf_
         */
        // TODO: This part is dependent on feature.
        // Why is it uint64???
        // TODO:
        int n = treesearch_data_->node_stack_->GetItemNum(treesearch_data_->itemset_buf_);
        DBG(
            D(3) << n << " items in the itemset\n";
        );
        int* array = treesearch_data_->node_stack_->GetItemArray(treesearch_data_->itemset_buf_);
        std::vector<int> buf(array, array + n);

        return algorithm_->UpdateSupBuf(buf);
    }

    std::vector<int> DFSContinuousPMStack::GetChildren(std::vector<int> items) {
        if (items.empty()) {
            // Edge partitioning for root node
            std::vector<int> children = d_->GetChildren(items);
            std::vector<int> responsible;
            for (size_t i = mpi_data_.mpiRank_; i < children.size();
                 i += mpi_data_.nTotalProc_) {
                responsible.push_back(children[i]);
            }
            return responsible;
        } else {
            return d_->GetChildren(items);
        }
    }

    // TODO: what we need here?
    void DFSContinuousPMStack::ProcessNode(int phase, double freq, int* ppc_ext_buf) {
        closed_set_num_++;

        // TODO: For continuous pattern mining calculating p value should be put later?
        if (phase == 1) {
            algorithm_->AddFrequence(freq);
        } else if (phase == 4) {
            int disFreq = algorithm_->GetDiscretizedFrequency(freq);
            IncCsAccum(disFreq);
        } else if (phase == 5) {
            // For Top-K we store pvalue instead of frequencies.
            int n = treesearch_data_->node_stack_->GetItemNum(ppc_ext_buf);
            int* array = treesearch_data_->node_stack_->GetItemArray(ppc_ext_buf);
            std::vector<int> buf(array, array + n);
            double actual_pvalue = d_->CalculatePValue(buf);
            int disPvalue = algorithm_->GetDiscretizedFrequency(actual_pvalue);
            // TODO: larger pvalue should come first.
            disPvalue = algorithm_->GetThresholds().size() - disPvalue;

//		printf("ProcessNode:");
//		for (int i = 0; i < buf.size(); ++i) {
//			printf(" %d", buf[i]);
//		}
//		printf(", pval = %.6f, discrete-pval = %d\n", actual_pvalue, disPvalue);

//		frequencies.push_back(actual_pvalue);
            IncCsAccum(disPvalue);
        } else if (phase == 6) {
            // TODO: Item is considered testable if freq > threshold.
//				double pmin = d_->CalculatePMin(freq);
            int n = treesearch_data_->node_stack_->GetItemNum(ppc_ext_buf);
            int* array = treesearch_data_->node_stack_->GetItemArray(ppc_ext_buf);
            std::vector<int> buf(array, array + n);
            double actual_pvalue = d_->CalculatePValue(buf);
//		printf("pval = %f, pmin = %f\n", actual_pvalue, thre_pmin_);
            if (actual_pvalue < algorithm_->GetThrePmin()) {
                algorithm_->UpdateGetTestableData(ppc_ext_buf, freq);
            }
        } else if (phase == 2) {
            if (freq > algorithm_->GetThreFreq()) { // permits == case?
                // TODO: Item is considered testable if freq > threshold.
                algorithm_->UpdateGetTestableData(ppc_ext_buf, freq);
            }
        } else {
            assert(false && "unknown phase");
        }
    }

    bool DFSContinuousPMStack::CheckProcessNodeEnd(int n, bool n_is_ms, int processed, long long int start_time, Timer* timer) {
//	printf("CheckProcessNodeEnd\n");
        long long int elapsed_time;
        if (n_is_ms) {
            if (processed > 0) {
                elapsed_time = timer->Elapsed() - start_time;
                if (elapsed_time >= n * 1000000) {  // ms to ns
                    return true;
                }
            }
        } else if (processed >= n)
            return true;

        return false;
    }

    void DFSContinuousPMStack::IncCsAccum(int sup_num) {
        algorithm_->IncCsAccum(sup_num);
    }

    bool DFSContinuousPMStack::TestAndPushNode(int new_item) {
        //	printf("TestAndPushNode\n");
        //	bsh_->Copy(sup_buf_, child_sup_buf_);
        //	memcpy(sup_buf_, child_sup_buf_, );
        // TODO: COPY sup_buf to child_sup_buf
        algorithm_->UpdateChildSupBuf(new_item);
        //	printf("GetChildrenFreq done\n");
        Feature child_freq = algorithm_->GetFreqFromDatabase();
        //	printf("GetFreq done\n");

        //	int sup_num = bsh_->AndCountUpdate(d_->NthData(new_item), child_sup_buf_);
        // If the support is smaller than the required minimal support for
        // significant pattern (=lambda), then prune it.

        // TODO: Can we inverse calculate the frequency from a pmin???
        //	double pmin = d_->CalculatePMin(child_freq);
        if (child_freq < algorithm_->GetThreFreq()) {
        //		printf("Node Pruned as pmin %.8f < pmin_thre_ %.8f\n",
        //				pmin, pmin_thre_);
            return false;
        }

        // TODO: Here it is inserting a new item into the node_stack_.
        treesearch_data_->node_stack_->PushPre();
        int* ppc_ext_buf = treesearch_data_->node_stack_->Top();
        //	printf("treesearch_data->node_stack_->Top() done\n");
        // TODO: Not sure what PPCExtension does
        bool res = d_->PPCExtension(treesearch_data_->node_stack_,
                                    treesearch_data_->itemset_buf_, new_item, ppc_ext_buf);

        // TODO: Those things should be done in other class.
        //	treesearch_data->node_stack_->SetSup(ppc_ext_buf, child_freq); // TODO: we can use the stack to hold double?
        treesearch_data_->node_stack_->PushPostNoSort();
        //	printf("treesearch_data->node_stack_->PushPostNoSort() done\n");

        if (!res) {        // todo: remove this redundancy
            treesearch_data_->node_stack_->Pop();
//		printf("treesearch_data->node_stack_->Pop() done\n");

            printf("Node Pruned as node is not PPCExtension\n");
            return false;
        }

        treesearch_data_->node_stack_->SortTop();
//	printf("treesearch_data->node_stack_->SortTop() done\n");

        return true;
    }

};