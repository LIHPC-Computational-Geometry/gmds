#ifndef METRIC_ADAPTATION
#define METRIC_ADAPTATION

namespace gmds{
  namespace hybrid{

    class SimplexMesh;
    class MetricAdaptation
    {
      public:
        MetricAdaptation(SimplexMesh* m_simplexMesh);

        ~MetricAdaptation();

        //smooth the current mesh's metric
        void metricCorrection();

        //compute the new mesh
        void execute();

        void buildEdgesMap();

      private:
        SimplexMesh* m_simplexMesh = nullptr;

        std::multimap<TInt, TInt> m_edgesMap;
    };
  }
}

#endif //METRIC_ADAPTATION
