//
// Created by Takuo Doi on 2017/11/15.
//

#include "DFSWithTreeSearchDataStack.h"

#include "MPI_Data.h"

namespace lamp_search {

    bool DFSWithTreeSearchDataStack::PopNodeFromStack() {
        /**
         * Pop item index and put into itemset_buf.
         * This part is not dependent on the feature.
         */
        treesearch_data_->node_stack_->CopyItem(treesearch_data_->node_stack_->Top(), treesearch_data_->itemset_buf_);
        treesearch_data_->node_stack_->Pop();

        DBG(
            D(3) << "expanded ";
            treesearch_data_->node_stack_->Print(D(3), treesearch_data_->itemset_buf_);
        );

        return AfterPopNodeFromStack();
    }

    bool DFSWithTreeSearchDataStack::AfterPopNodeFromStack() {
        return false;
    }

    int DFSWithTreeSearchDataStack::Count() const {
        return treesearch_data_->node_stack_->NuItemset();
    }

    int DFSWithTreeSearchDataStack::UsedCapacity() const {
        return treesearch_data_->node_stack_->UsedCapacity();
    }

    bool DFSWithTreeSearchDataStack::IsStealStarted() const {
        return treesearch_data_->stealer_->StealStarted();
    }

    bool DFSWithTreeSearchDataStack::IsStealAllowed() const {
        return mpi_data_.nTotalProc_ != 1 &&
               treesearch_data_->stealer_->StealStarted() &&
               !treesearch_data_->stealer_->Requesting() &&
               treesearch_data_->node_stack_->Empty();
    }

    const std::pair<int, int> DFSWithTreeSearchDataStack::GetStealInfo() const {
        switch (treesearch_data_->stealer_->State()) {
            case StealState::RANDOM: {
                int victim = mpi_data_.victims_[mpi_data_.rand_m_()];
                assert(victim <= mpi_data_.nTotalProc_ && "stealrandom");
                DBG(
                        D(2) << "Steal Random:" << "\trequesting="
                             << treesearch_data_->stealer_->Requesting()
                             << "\trandom counter="
                             << treesearch_data_->stealer_->RandomCount()
                             << std::endl;
                );
                return std::make_pair(victim, -1);
            }
            case StealState::LIFELINE: {
                int lifeline = mpi_data_.lifelines_[treesearch_data_->stealer_->LifelineVictim()];
                assert(lifeline >= 0); // at least lifelines_[0] must be 0 or higher
                if (!mpi_data_.lifelines_activated_[lifeline]) {
                    mpi_data_.lifelines_activated_[lifeline] = true;
                    // becomes false only at RecvGive
                    assert(lifeline <= mpi_data_.nTotalProc_ && "lifeline");
                    DBG(
                            D(2) << "Steal Lifeline:" << "\trequesting="
                                 << treesearch_data_->stealer_->Requesting()
                                 << "\tlifeline counter="
                                 << treesearch_data_->stealer_->LifelineCounter()
                                 << "\tlifeline victim="
                                 << treesearch_data_->stealer_->LifelineVictim()
                                 << "\tz_=" << mpi_data_.hypercubeDimension_
                                 << std::endl;
                    );
                    return std::make_pair(lifeline, 1);
                }
                break;
            }
            default:
                assert(0);
                break;
        }
        throw new std::exception();
    }

    void DFSWithTreeSearchDataStack::UpdateStealInfo() {
        switch (treesearch_data_->stealer_->State()) {
            case StealState::RANDOM: {
                treesearch_data_->stealer_->SetRequesting();
                treesearch_data_->stealer_->IncRandomCount();
                if (treesearch_data_->stealer_->RandomCount() == 0) {
                    treesearch_data_->stealer_->SetState(StealState::LIFELINE);
                }
                break;
            }
            case StealState::LIFELINE: {
                treesearch_data_->stealer_->IncLifelineCount();
                treesearch_data_->stealer_->IncLifelineVictim();
                if (treesearch_data_->stealer_->LifelineCounter() >= mpi_data_.hypercubeDimension_ ||
                    mpi_data_.lifelines_[treesearch_data_->stealer_->LifelineVictim()] < 0) { // hack fix
                    // note: this resets next_lifeline_victim_, is this OK?
                    treesearch_data_->stealer_->Finish();
                    DBG(D(2) << "Steal finish:" << std::endl;);
                    // note: one steal phase finished, don't restart until receiving give?

                    // in x10 glb, if steal() finishes and empty==true,
                    // processStack() will finish, but will be restarted by deal()
                    // for the same behavior, reseting the states to initial state should be correct
                }
                break;
            }
            default:
                assert(0);
                break;
        }
    }

    void DFSWithTreeSearchDataStack::ResetSteal() {
        treesearch_data_->stealer_->ResetRequesting();
    };

    bool DFSWithTreeSearchDataStack::HasJobToDo() const {
        return !(treesearch_data_->node_stack_->Empty())
               || (mpi_data_.thieves_->Size() > 0)
               || treesearch_data_->stealer_->StealStarted()
               || mpi_data_.processing_node_; // thieves_ and stealer state check
    }

    void DFSWithTreeSearchDataStack::InitSteal() {
        treesearch_data_->stealer_->ResetRequesting();
        treesearch_data_->stealer_->ResetCounters();
        treesearch_data_->stealer_->SetState(StealState::RANDOM);
        treesearch_data_->stealer_->SetStealStart();
    }

    void DFSWithTreeSearchDataStack::PrintItemset(int* itembuf) const {
        std::vector<int> itemset = treesearch_data_->node_stack_->getItems(itembuf);
        printf("node = ");
        for (size_t i = 0; i < itemset.size(); ++i) {
            printf("%d ", itemset[i]);
        }
    }

    void DFSWithTreeSearchDataStack::PrintAll(std::ostream& stream) const {
        treesearch_data_->node_stack_->PrintAll(stream);
    }

}