//
// Created by Takuo Doi on 2017/11/15.
//

#ifndef MP_LAMP_DFSCONTINUOUSPMSTACK_H
#define MP_LAMP_DFSCONTINUOUSPMSTACK_H

#include "DFSWithTreeSearchDataStack.h"

namespace lamp_search {

    class ContDatabase;
    class ParallelContinuousPM;

    struct TreeSearchData;
    struct MPI_Data;

    typedef double Feature;

    class DFSContinuousPMStack: public DFSWithTreeSearchDataStack {

    private:
        ParallelContinuousPM* algorithm_;

        ContDatabase* d_;

        // Statistics
        int expand_num_;
        int closed_set_num_;

    public:
        DFSContinuousPMStack(ParallelContinuousPM* algorithm, TreeSearchData* treesearch_data, ContDatabase* d, MPI_Data& mpi_data, std::ostream& lfs, ProbeCallback callback)
                : DFSWithTreeSearchDataStack(treesearch_data, mpi_data, lfs, callback), algorithm_(algorithm), d_(d) {}

        virtual int Process(int phase, Timer* timer, Log* log) override;
        virtual int Split(VariableLengthItemsetStack* give_stack) override;
        virtual void Merge(VariableLengthItemsetStack* give_stack, int count) override;

        virtual bool AfterPopNodeFromStack() override;

        std::vector<int> GetChildren(std::vector<int> items);
        void ProcessNode(int phase, double freq, int* ppc_ext_buf);
        bool CheckProcessNodeEnd(int n, bool n_is_ms, int processed, long long int start_time, Timer* timer);
        void ProcAfterProbe();

    private:
        void IncCsAccum(int sup_num);

        bool TestAndPushNode(int new_item);

        int GetDiscretizedFrequency(double freq) const;
    };
}

#endif //MP_LAMP_DFSCONTINUOUSPMSTACK_H
