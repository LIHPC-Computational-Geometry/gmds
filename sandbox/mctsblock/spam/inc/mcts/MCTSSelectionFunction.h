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
struct SPUCTSelectionFunction: public ISelectionFunction
{
	 SPUCTSelectionFunction(const double AC = 1.42, const double AD = 100);
	 MCTSTree* select(MCTSTree* ANode) const override ;
  private:
	 double m_c;
	 double m_d;
};
/*---------------------------------------------------------------------------*/
#endif //SPAM_MCTSSELECTIONFUNCTIONS_H
/*---------------------------------------------------------------------------*/
