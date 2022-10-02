#ifndef METRIC_ADAPTATION
#define METRIC_ADAPTATION

namespace gmds{
  namespace hybrid{

    class SimplexMesh;

    class MetricFFPointgeneration
    {
      public:
        MetricFFPointgeneration(SimplexMesh* m_simplexMesh);

        ~MetricFFPointgeneration();

        std::map<unsigned int, std::vector<TInt>> buildSortedEdges() const;

        std::vector<std::vector<double>> buildParamEdgeU(const std::map<unsigned int, std::vector<TInt>>& sortedEdge) const;

        void execute();

      private:
        SimplexMesh* m_simplexMesh = nullptr;
    };
  }
}

#endif //METRIC_ADAPTATION
