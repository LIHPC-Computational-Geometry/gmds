#ifndef CRITERION_H_
#define CRITERION_H_
/******************************************************************************/
#include <gmds/sofiane/SimplexMesh.h>
#include <gmds/sofiane/SimplicesNode.h>
#include <gmds/sofiane/SimplicesCell.h>
/******************************************************************************/

namespace gmds
{
    namespace hybrid
    {
      namespace operators
      {
        class ICriterion
        {
        public:

          ICriterion(){};

          virtual ~ICriterion(){};

          virtual bool execute(SimplexMesh* m_simplex_mesh, const TSimplexID simplexId, const TInt nodeIndexLocal, const gmds::math::Point& nextPt) const = 0;
        };

        class VolumeCriterion : public ICriterion
        {
        public:

          VolumeCriterion  (){};

          ~VolumeCriterion (){};

          virtual bool execute(SimplexMesh* m_simplex_mesh, const TSimplexID simplexId, const TInt nodeIndexLocal, const gmds::math::Point& nextPt) const ;
        };

        class metricLengthCriterion : public ICriterion
        {
        public:

          metricLengthCriterion  (){};

          ~metricLengthCriterion (){};

          virtual bool execute(SimplexMesh* m_simplex_mesh, const TSimplexID simplexId, const TInt nodeIndexLocal, const gmds::math::Point& nextPt) const {return false; /*TODO*/}
        };


        class CriterionRAIS
        {

        public:
          CriterionRAIS(ICriterion* criterion = new VolumeCriterion());

          ~CriterionRAIS();

          bool execute(SimplexMesh* m_simplex_mesh, const TSimplexID simplexId, const TInt nodeIndexLocal, const gmds::math::Point& nextPt) const ;

        private:
          ICriterion * m_criterion = nullptr;
        };

    }
  }
}
#endif //CRITERION_H_
