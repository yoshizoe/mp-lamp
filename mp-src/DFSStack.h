//
// Created by Takuo Doi on 2017/11/05.
//

#ifndef MP_LAMP_DFSSTACK_H
#define MP_LAMP_DFSSTACK_H

#include <vector>
#include <fstream>

#include "variable_length_itemset.h"
#include "utils.h"
#include "database.h"

namespace lamp_search {

    class Timer;
    class Log;

    class DFSStack {
    public:
        typedef std::function<void()> ProbeCallback;

    private:
        std::ostream& lfs_;
        ProbeCallback callback_;

    public:
        DFSStack(std::ostream& lfs, ProbeCallback callback) : lfs_(lfs), callback_(callback) {}

        virtual int Process(int phase, Timer* timer, Log* log) = 0;
        virtual int Split(VariableLengthItemsetStack* give_stack) = 0;
        virtual void Merge(VariableLengthItemsetStack* give_stack, int count) = 0;
        virtual int Count() const = 0;
        virtual int UsedCapacity() const = 0;

        virtual bool IsStealStarted() const = 0;
        virtual bool IsStealAllowed() const = 0;
        virtual const std::pair<int, int> GetStealInfo() const = 0;
        virtual void UpdateStealInfo() = 0;
        virtual void ResetSteal() = 0;
        virtual bool HasJobToDo() const = 0;
        virtual void InitSteal() = 0;
        virtual void CheckStealFinish() = 0;

        virtual void PrintItemset(int* itembuf) const = 0;
        virtual void PrintAll(std::ostream& stream) const = 0;

        virtual void CheckProbe(int& accum_period_counter_, long long int& lap_time, Timer* timer, Log* log);

        bool IsEmpty() const { return Count() == 0; }

    protected:
        void CallForProbe();

        static std::ofstream null_stream_;
        std::ostream& D(int level, bool show_phase = true) const;

    };

}

#endif //MP_LAMP_DFSSTACK_H
