//
// Created by Takuo Doi on 2017/11/15.
//

#include <gflags/gflags.h>
#include "DFSParallelPatternMiningStack.h"
#include "MPI_Data.h"
#include "ParallelPatternMining.h"

namespace lamp_search {

    int DFSParallelPatternMiningStack::Process(int phase, Timer *timer, Log *log) {
        long long int start_time, lap_time;
        start_time = timer->Elapsed();
        lap_time = start_time;

        int processed = 0;
        mpi_data_.processing_node_ = true;
        while (!treesearch_data_->node_stack_->Empty()) {
            processed++;
            expand_num_++;

            PopNodeFromStack();

            assert(phase != 1 || treesearch_data_->node_stack_->GetItemNum(treesearch_data_->itemset_buf_) != 0);

            int accum_period_counter_ = 0;

            // TODO: core_i should be hidden
            int core_i = algorithm_->CalcCoreI(treesearch_data_->node_stack_, treesearch_data_->itemset_buf_);

            std::vector<int> children = GetChildren(core_i);
            for (std::vector<int>::iterator it = children.begin(); it != children.end(); ++it) {
                int new_item = *it;

                if (treesearch_data_->node_stack_->Exist(treesearch_data_->itemset_buf_, new_item)) {
                    return false;
                }
                // skip existing item
                // todo: improve speed here
                CheckProbe(accum_period_counter_, lap_time, timer, log);

                if (!TestAndPushNode(new_item, core_i)) {
                    continue;
                }

                int *ppc_ext_buf = treesearch_data_->node_stack_->Top();
                int sup_num = treesearch_data_->node_stack_->GetSup(ppc_ext_buf);

                DBG(
                    if (phase == 2) {
                        D(3) << "found cs ";
                        treesearch_data_->node_stack_->Print(D(3), ppc_ext_buf);
                    }
                );

                ProcessNode(phase, sup_num, ppc_ext_buf);

                int lambda = algorithm_->GetLamba();

                assert(sup_num >= lambda);

                // try skipping if supnum_ == sup_threshold,
                // because if sup_num of a node equals to sup_threshold, children will have smaller sup_num
                // therefore no need to check it's children
                // note: skipping node_stack_ full check. allocate enough memory!
                if (sup_num == lambda) {
                    treesearch_data_->node_stack_->Pop();
                }
            }

            if (CheckProcessNodeEnd(mpi_data_.granularity_, mpi_data_.isGranularitySec_, processed, start_time, timer)) {
                break;
            }
        }

        long long int elapsed_time = timer->Elapsed() - lap_time;
        log->d_.process_node_time_ += elapsed_time;
        log->d_.process_node_num_ += processed;

        DBG(D(2) << "processed node num=" << processed << "\ttime=" << elapsed_time << std::endl;);

        mpi_data_.processing_node_ = false;
        return true;
    }

    int DFSParallelPatternMiningStack::Split(VariableLengthItemsetStack* give_stack) {
        return treesearch_data_->node_stack_->Split(give_stack);
    }

    void DFSParallelPatternMiningStack::Merge(VariableLengthItemsetStack* give_stack, int count) {
        treesearch_data_->node_stack_->MergeStack(
                give_stack->Stack() + VariableLengthItemsetStack::SENTINEL + 1,
                count - VariableLengthItemsetStack::SENTINEL - 1);
    }

    bool DFSParallelPatternMiningStack::AfterPopNodeFromStack() {
        /**
         * Get the support (frequency) of the item into sup_buf_
         */
        // TODO: This part is dependent on feature.
        // Why is it uint64???
        bsh_->Set(sup_buf_);
        {
            // This should definitely be inside the Domain class.
            // n := number of items in the itemset.
            int n = treesearch_data_->node_stack_->GetItemNum(treesearch_data_->itemset_buf_);
            for (int i = 0; i < n; i++) {
                int item = treesearch_data_->node_stack_->GetNthItem(treesearch_data_->itemset_buf_, i);
                bsh_->And(d_->NthData(item), sup_buf_);
            }
        }
    }

    // TODO: Dependent on d_. This function should be responsible of Domain class.
    // TODO: Should make clear distinction of LAMP children and immediate children?
    std::vector<int> DFSParallelPatternMiningStack::GetChildren(int core_i) {
        bool is_root_node = (treesearch_data_->node_stack_->GetItemNum(treesearch_data_->itemset_buf_) == 0);

        // reverse order
        // for ( int new_item = d_->NuItems()-1 ; new_item >= core_i+1 ; new_item-- )
        std::vector<int> children;
        for (int new_item = d_->NextItemInReverseLoop(is_root_node, mpi_data_.mpiRank_, mpi_data_.nTotalProc_,
                                                      d_->NuItems());
             new_item >= core_i + 1;
             new_item = d_->NextItemInReverseLoop(is_root_node, mpi_data_.mpiRank_, mpi_data_.nTotalProc_, new_item)) {
            children.push_back(new_item);
        }
        return children;
    }

    /**
     * Test two pruning criteria before pushing
     * 1. The support of the itemset should be larger or equal to lambda (minimal support).
     * 2. New node should be a PPC extension of the parent itemset.
     */
    bool DFSParallelPatternMiningStack::TestAndPushNode(int new_item, int core_i) {

        int sup_num = algorithm_->UpdateChildSupBuf(new_item, core_i);
        if (sup_num < 0) {
            return false;
        }

        // TODO: Here it is inserting a new item into the node_stack_.
        treesearch_data_->node_stack_->PushPre();
        int *ppc_ext_buf = treesearch_data_->node_stack_->Top();

        // TODO: return true if child_sup_buf_ is PPC extension???
        bool res = algorithm_->PPCExtension(treesearch_data_->node_stack_, treesearch_data_->itemset_buf_, core_i, new_item, ppc_ext_buf);

        // TODO: Those things should be done in other class.
        treesearch_data_->node_stack_->SetSup(ppc_ext_buf, sup_num);
        treesearch_data_->node_stack_->PushPostNoSort();

        if (!res) {
            // todo: remove this redundancy
            treesearch_data_->node_stack_->Pop();
            return false;
        }

        treesearch_data_->node_stack_->SortTop();

        return true;
    }

    void DFSParallelPatternMiningStack::ProcessNode(int phase, int sup_num, int *ppc_ext_buf) {
        if (phase == 1) {
            IncCsAccum(sup_num); // increment closed_set_num_array
        } else if (phase == 2) {
            closed_set_num_++;
            if (true) { // XXX: FLAGS_third_phase_
                int pos_sup_num = bsh_->AndCount(d_->PosNeg(), child_sup_buf_);
                double pval = d_->PVal(sup_num, pos_sup_num);
                assert(pval >= 0.0);
                algorithm_->UpdateFreqStack(pval, ppc_ext_buf);
            }
        }
    }

    void DFSParallelPatternMiningStack::IncCsAccum(int sup_num) {
        algorithm_->IncCsAccum(sup_num);
    }

    bool DFSParallelPatternMiningStack::CheckProcessNodeEnd(int n, bool n_is_ms, int processed, long long int start_time, Timer* timer) {
        long long int elapsed_time;
        if (n_is_ms) {
            if (processed > 0) {
                elapsed_time = timer->Elapsed() - start_time;
                if (elapsed_time >= n * 1000000) {  // ms to ns
                    return true;
                }
            }
        } else if (processed >= n) {
            return true;
        }

        return false;
    }

}