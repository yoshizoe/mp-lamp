//
// Created by Takuo Doi on 2017/11/15.
//

#ifndef MP_LAMP_DFSSTACKWITHTREESEARCHDATA_H
#define MP_LAMP_DFSSTACKWITHTREESEARCHDATA_H

#include "DFSStack.h"

namespace lamp_search {

    struct TreeSearchData;
    struct MPI_Data;
    struct GetMinSupData;
    struct GetTestableData;
    struct GetContSignificantData;
    struct GetSignificantData;

    class DFSWithTreeSearchDataStack : public DFSStack {

    protected:
        TreeSearchData *treesearch_data_;
        MPI_Data& mpi_data_;

    public:
        DFSWithTreeSearchDataStack(TreeSearchData* treesearch_data, MPI_Data& mpi_data, std::ostream& lfs, ProbeCallback callback)
        : DFSStack(lfs, callback), treesearch_data_(treesearch_data), mpi_data_(mpi_data) {}

        virtual int Count() const override;
        virtual int UsedCapacity() const override;

        virtual bool IsStealStarted() const override;
        virtual bool IsStealAllowed() const override;
        virtual const std::pair<int, int> GetStealInfo() const override;
        virtual void UpdateStealInfo() override;
        virtual void ResetSteal() override;
        virtual bool HasJobToDo() const override;
        virtual void InitSteal() override;
        virtual void CheckStealFinish() override;

        void PrintItemset(int* itembuf) const override;
        void PrintAll(std::ostream& stream) const override;

        // TODO: This function should be removed for capsulization
        TreeSearchData* getTreeSearchData() { return treesearch_data_; }

    protected:
        bool PopNodeFromStack();
        virtual bool AfterPopNodeFromStack();

    };

}


#endif //MP_LAMP_DFSSTACKWITHTREESEARCHDATA_H
