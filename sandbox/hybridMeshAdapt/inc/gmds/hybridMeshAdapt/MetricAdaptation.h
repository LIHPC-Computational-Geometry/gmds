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

        //Compute the mesh Slicing
        unsigned int computeSlicing();

        bool computeSlicing(const TInt nodeA, const TInt nodeB) ;

        //Compute the mesh Slicing
        unsigned int computeEdgeRemove();

        bool computeEdgeRemove(const TInt nodeA, const TInt nodeB) const ;
        //Compute the mesh Slicing
        unsigned int computeFaceSwap();

        //Compute the EdgeSwap only on the surface
        unsigned int computeSurfaceEdgeSwap();

        //Compute the mesh Slicing
        unsigned int computePointSmoothing();

        //compute the new mesh
        void execute();

        //compute the new mesh
        void executeCustomMethod();

        void buildEdgesMap();

      private:
        SimplexMesh* m_simplexMesh = nullptr;

        std::multimap<TInt, TInt> m_edgesMap;
    };
  }
}

#endif //METRIC_ADAPTATION
