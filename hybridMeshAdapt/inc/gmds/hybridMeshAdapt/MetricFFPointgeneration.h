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


        MetricFFPointgeneration(SimplexMesh* m_simplexMesh, const std::string& name);

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

        Eigen::Matrix3d metricInterpolationWithDistorsion(const Eigen::Matrix3d  & metric, const Eigen::Matrix3d  & frameField) const;

        double computeDistorsionMetric(const Eigen::Matrix3d  & m) const;

        void nodesSpreading(std::vector<TInt>& nodesAdded, bool surfaceFlag = false);

        void execute();

        bool metricSamplingEdge(const unsigned int n, std::vector<double>& res, const std::vector<TInt>& edge, const std::vector<double>& edgeU) const;

        math::Point computeTheEdgeNodeCoordinate(const double u, const std::vector<TInt>& edge, const std::vector<double>& edgeU, TInt& NodeA, TInt& NodeB) const;

        void findOptimimalPosition0(const TInt node, math::Point &newCoord, bool surfaceFlag = false, int cpt = 10, double epsilon = 0.1/*0.01*/) ;

        bool findOptimimalPosition(const TInt node, math::Point &newCoord, bool surfaceFlag = false, int cpt = 10, double epsilon = 0.1/*0.01*/) ;

        bool nodeFiltering(const math::Point& pt, const TInt fromNode, const TSimplexID simplex, std::vector<TInt> & neighboorNode);

        void computeQuadFaces(std::set<std::vector<TInt>> & faces) const ;

        void computeHexa(std::set<std::vector<TInt>> & hexas) ;

        void correctNodeLabel() ;

        void addNodeToLayer(const TInt nodeId, const TInt fromNode = -1, bool surfaceFlag = false);

        void incrementLayer();

        void initializeGridWithEdge();

        void correctionNodeStructure();

        void correctionVertexNode();

        void sortBySurfaceNodeAdded(std::vector<TInt>& nodesAdded);

        bool belongToEdge(const math::Point & nodeCoord);

        //This function will order the first curve nodes
        void correctUnwantedConnectionSURFACE();

        void correctUnwantedConnectionVOLUME();

        void connectionWithNeighbor(const std::vector<TInt>& nodesAdded);

      private:
        std::unordered_map<TInt, TSimplexID> nodeBelongingTO;

        unsigned int m_layerNbr;

        std::unordered_map<TInt, std::list<TInt>> m_layers;

        std::unordered_map<TInt, int> m_nodeLayerNbr;

        //this map will help us to corret the different unwanted connection in m_nodeStructure
        std::unordered_map<TInt, std::vector<TInt>> m_nodeGeneratedBy;

        std::unordered_map<TInt, std::vector<TInt>> m_listEdge;

        std::unordered_map<TInt, std::vector<TInt>> m_nodeStructure;

        SimplexMesh* m_simplexMesh = nullptr;

        SimplexMesh m_nodesMesh;

        Octree m_oc;

        double m_minDistance;

        std::string m_name;
    };
  }
}

#endif //METRIC_ADAPTATION
