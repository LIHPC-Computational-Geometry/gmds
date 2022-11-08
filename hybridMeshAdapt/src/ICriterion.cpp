/******************************************************************/
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/ICriterion.h>
#include <gmds/hybridMeshAdapt/SimplicesCell.h>
/******************************************************************/
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesCell;
using namespace simplicesTriangle;
using namespace simplicesNode;
/******************************************************************/
#define EPSILON  10E-10
#define EPSILON_VOLUME 0.0//10E-8
/******************************************************************/
bool VolumeCriterion::execute(SimplexMesh* m_simplex_mesh,const TSimplexID simplexId, const TInt nodeIndexLocal,const gmds::math::Point& nextPt) const
{
  bool flag = false;
  TInt border = std::numeric_limits<TInt>::min();
  math::Orientation::Sign orientation;

  if(simplexId >= 0)
  {
      orientation = SimplicesCell(m_simplex_mesh, simplexId).orientation(nodeIndexLocal, nextPt);
  }
  else
  {
    std::cout << "SimplexId < 0 in VolumeCriterion" << std::endl;
  }

  flag = (orientation < 1)? true:false;

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
    std::cout << "EXEPTION levée dans le destructeur CriterionRAIS" << std::endl;
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
