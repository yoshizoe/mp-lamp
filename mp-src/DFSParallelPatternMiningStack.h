//
// Created by Takuo Doi on 2017/11/15.
//

#ifndef MP_LAMP_PARALLELPATTERNMININGSTACK_H
#define MP_LAMP_PARALLELPATTERNMININGSTACK_H

#include <lamp_graph.h>
#include "DFSStack.h"
#include "DFSWithTreeSearchDataStack.h"

namespace lamp_search {

    class ParallelPatternMining;
    struct TreeSearchData;
    struct MPI_Data;
    struct GetMinSupData;
    struct GetTestableData;
    struct GetContSignificantData;
    struct GetSignificantData;

    class DFSParallelPatternMiningStack : public DFSWithTreeSearchDataStack {
    private:
        ParallelPatternMining* algorithm_;

        TreeSearchData *treesearch_data_;
        Database <uint64> *d_;
        MPI_Data &mpi_data_;

//        LampGraph<uint64> * g_;
        VariableBitsetHelper<uint64> * bsh_;
        uint64 * sup_buf_, *child_sup_buf_; // TODO: sup_buf_ is only used in ProcessNode and PreProcessRootNode!

        /**
         * Statistics
         */
        int expand_num_;
        int closed_set_num_;

    public:
        DFSParallelPatternMiningStack(ParallelPatternMining* algorithm, TreeSearchData *treesearch_data, Database <uint64> *d, MPI_Data &mpi_data, std::ostream &lfs, ProbeCallback callback)
                : DFSWithTreeSearchDataStack(treesearch_data, mpi_data, lfs, callback),
                  algorithm_(algorithm), treesearch_data_(treesearch_data), d_(d), mpi_data_(mpi_data) {}

        virtual int Process(int phase, Timer* timer, Log* log) override;
        virtual int Split(VariableLengthItemsetStack* give_stack) override;
        virtual void Merge(VariableLengthItemsetStack* give_stack, int count) override;

        virtual bool AfterPopNodeFromStack() override;

    private:
        std::vector<int> GetChildren(int core_i);
        bool TestAndPushNode(int new_item, int core_i);
        void ProcessNode(int phase, int sup_num, int* ppc_ext_buf);
        void IncCsAccum(int sup_num);
        bool CheckProcessNodeEnd(int n, bool n_is_ms, int processed, long long int start_time, Timer* timer);
    };

}

#endif //MP_LAMP_PARALLELPATTERNMININGSTACK_H
