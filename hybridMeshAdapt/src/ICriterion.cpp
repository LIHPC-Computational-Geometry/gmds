/******************************************************************/
#include <gmds/hybridMeshAdapt/ICriterion.h>
/******************************************************************/
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesCell;
using namespace simplicesNode;
/******************************************************************/
#define EPSILON 10E-20
/******************************************************************/
bool VolumeCriterion::execute(SimplexMesh* m_simplex_mesh,const TSimplexID simplexId, const TInt nodeIndexLocal,const gmds::math::Point& nextPt) const
{
  bool flag = false;
  flag = (SimplicesCell(m_simplex_mesh, simplexId).signedBarycentric(nodeIndexLocal, nextPt) <= EPSILON)? true:false;
  return flag;
}
/******************************************************************/
CriterionRAIS::CriterionRAIS(ICriterion* criterion):
m_criterion(criterion)
{
}
/******************************************************************/
CriterionRAIS::~CriterionRAIS()
{
  if(m_criterion != nullptr)
  {
      delete m_criterion;
  }
  else
  {
    /*TODO* exeption criterion* deja supprimé*/
    std::cout << "EXEPTION LEvé dans le destructeur CriterionRAIS" << std::endl;
  }
}
/******************************************************************/
bool CriterionRAIS::execute(SimplexMesh* m_simplex_mesh, const TSimplexID simplexId, const TInt nodeIndexLocal, const gmds::math::Point& nextPt) const
{

  if(m_criterion != nullptr)
  {
      return m_criterion->execute(m_simplex_mesh, simplexId, nodeIndexLocal, nextPt);
  }
  else
  {
    /*TODO exception ...*/
    std::cout << "EXEPTION LEVé dans la fonction execute du fichier ICRITERION.CPP" << std::endl;
  }
  return false;
}
/******************************************************************/
