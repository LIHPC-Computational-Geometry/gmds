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
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
/*****************************************************************************/
namespace gmds{
  namespace hybrid{

    class SimplexMesh;

    class MetricFFPointgeneration
    {
      public:

        struct DataEdges {
          std::vector<int> subEdge {};
          std::vector<double> subEdgeU {};
          double sizeEdge = 0;
        };

        //structure will permit us to sort the sampling node in order to
        //saple low metric first
        //node is the node where coord came from
        struct nodeSamplingData {
          TInt node;
          double m;
          math::Point coord;
        };


        MetricFFPointgeneration(SimplexMesh* m_simplexMesh);

        ~MetricFFPointgeneration();

        std::map<unsigned int, std::vector<TInt>> buildSortedEdges() const;

        std::vector<TInt> findBoundedNode(const double t, const std::vector<TInt>& edgeNodes) const;

        Eigen::VectorXd findCubicInterpolation(const double t, const std::vector<TInt>& edgeNodes, std::vector<int>& subEdgeIdx) const;

        double curvature(const double t, const std::vector<TInt>& edgeNodes) const;

        Eigen::Matrix3d computeIntersectionMetric(const Eigen::Matrix3d& m1, const Eigen::Matrix3d& m2) const;

        Eigen::VectorXd interpolateCubicCurve(double global_t, unsigned int edgeId);

        std::vector<std::vector<double>> buildParamEdgeU(const std::map<unsigned int, std::vector<TInt>>& sortedEdge, std::vector<double> & length_edges) const;

        std::vector<DataEdges> subdivideEdge(const std::vector<TInt>& edge, const std::vector<double>& edgeU, const double sizeEdge) const;

        void subdivideEdgeUsingMetric_Relaxation(std::vector<TInt>& nodesAdded, const std::vector<TInt>& edge, const std::vector<double>& edgeU, const double sizeEdge, const unsigned int edgeId) ;

        void nodesSpreading(std::vector<TInt>& nodesAdded, bool surfaceFlag = false);

        void execute();

        bool metricSamplingEdge(const unsigned int n, std::vector<double>& res, const std::vector<TInt>& edge, const std::vector<double>& edgeU) const;

        math::Point computeTheEdgeNodeCoordinate(const double u, const std::vector<TInt>& edge, const std::vector<double>& edgeU) const;

        void findOptimimalPosition(const TInt node, math::Point &newCoord) ;

        void nodeFiltering(const math::Point& pt, std::vector<TInt> & neighboorNode, double k = 1.1 * (sqrt(2.0) * 0.5), bool flag = false);

        void computeQuadFaces(std::set<std::vector<TInt>> & faces) const ;

        void computeHexa(std::set<std::vector<TInt>> & hexas) ;

        void correctNodeLabel() ;

        void processNodesStructure();

      private:
        std::unordered_map<TInt, std::vector<TInt>> m_nodeStructure;

        SimplexMesh* m_simplexMesh = nullptr;

        SimplexMesh m_nodesMesh;

        Octree m_oc;
    };
  }
}

#endif //METRIC_ADAPTATION
