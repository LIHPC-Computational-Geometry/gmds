#ifndef METRIC_ADAPTATION
#define METRIC_ADAPTATION
/*****************************************************************************/
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
#include <gmds/math/Matrix.h>
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

        void subdivideEdgeUsingMetric(Eigen::Vector3d & dir, std::vector<TInt>& nodesAdded, const std::vector<TInt>& edge, const std::vector<double>& edgeU, const double sizeEdge) const;

        void execute();

      private:
        SimplexMesh* m_simplexMesh = nullptr;
    };
  }
}

#endif //METRIC_ADAPTATION
