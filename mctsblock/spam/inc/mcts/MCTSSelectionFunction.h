/*---------------------------------------------------------------------------*/
#ifndef SPAM_MCTSSELECTIONFUNCTION_H
#define SPAM_MCTSSELECTIONFUNCTION_H
/*---------------------------------------------------------------------------*/
#include <cfloat>
#include <mcts/MCTSTree.h>
/*---------------------------------------------------------------------------*/
struct ISelectionFunction {
    virtual ~ISelectionFunction() = default;
    virtual MCTSTree* select(MCTSTree* ANode) const =0;
};
/*---------------------------------------------------------------------------*/
struct UCBSelectionFunction: public ISelectionFunction
{
    UCBSelectionFunction(const double AC = 1.42);
    MCTSTree* select(MCTSTree* ANode) const override ;
private:
    double m_c;
};
/*---------------------------------------------------------------------------*/
#endif //SPAM_MCTSSELECTIONFUNCTIONS_H
/*---------------------------------------------------------------------------*/
