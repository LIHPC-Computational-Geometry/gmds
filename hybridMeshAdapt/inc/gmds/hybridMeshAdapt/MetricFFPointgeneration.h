#ifndef METRIC_ADAPTATION
#define METRIC_ADAPTATION
/*****************************************************************************/
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
#include <gmds/math/Matrix.h>
/*****************************************************************************/
#include "gmds/hybridMeshAdapt/SimplicesCell.h"
#include "gmds/hybridMeshAdapt/SimplicesNode.h"
#include "gmds/hybridMeshAdapt/SimplicesTriangle.h"
#include "gmds/hybridMeshAdapt/Metric.h"
#include <gmds/hybridMeshAdapt/ICriterion.h>
#include <gmds/hybridMeshAdapt/ISimplexMeshIOService.h>
#include "gmds/hybridMeshAdapt/PointInsertion.h"
#include "gmds/hybridMeshAdapt/Octree.h"
/*****************************************************************************/
namespace gmds{
  namespace hybrid{

    class SimplexMesh;

    class MetricFFPointgeneration
    {
      public:
        MetricFFPointgeneration(SimplexMesh* m_simplexMesh);

        ~MetricFFPointgeneration();

        std::map<unsigned int, std::vector<TInt>> buildSortedEdges() const;

        std::vector<std::vector<double>> buildParamEdgeU(const std::map<unsigned int, std::vector<TInt>>& sortedEdge, std::vector<double> & length_edges) const;

        void subdivideEdgeUsingMetric_Relaxation(std::vector<TInt>& nodesAdded, const std::vector<TInt>& edge, const std::vector<double>& edgeU, const double sizeEdge, const unsigned int edgeId) ;

        void subdivideEdgeUsingMetric_Dichotomie(std::vector<TInt>& nodesAdded, const std::vector<TInt>& edge, const std::vector<double>& edgeU, const double sizeEdge) const;

        void nodesSpreading(std::vector<TInt>& nodesAdded);

        void execute();

        void metricSamplingEdge(const unsigned int n, std::vector<double>& res, const std::vector<TInt>& edge, const std::vector<double>& edgeU) const;

        math::Point computeTheEdgeNodeCoordinate(const double u, const std::vector<TInt>& edge, const std::vector<double>& edgeU) const;

        bool toReject(const math::Point& pt) const;

        void findOptimimalPosition(const TInt node, math::Point &newCoord) ;

        //return true if there is no point near (anisotropic length) pt, using octree for fast result
        bool nodeFiltering(const TInt node, const math::Point& pt);

      private:
        SimplexMesh* m_simplexMesh = nullptr;

        SimplexMesh m_nodesMesh;

        Octree m_oc;
    };
  }
}

#endif //METRIC_ADAPTATION
