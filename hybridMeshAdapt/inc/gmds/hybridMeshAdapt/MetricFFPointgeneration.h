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

        void execute();

      private:
        SimplexMesh* m_simplexMesh = nullptr;
    };
  }
}

#endif //METRIC_ADAPTATION
