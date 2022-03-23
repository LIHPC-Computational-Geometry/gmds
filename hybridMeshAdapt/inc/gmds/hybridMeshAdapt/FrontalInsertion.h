#ifndef FRONTAL_INSERTION.H
#define FRONTAL_INSERTION.H

namespace gmds
{
    namespace hybrid
    {
        class SimplexMesh;
        class FrontalInsertion
        {
          public:
            FrontalInsertion(SimplexMesh* simplexMesh):m_simplexMesh(simplexMesh){}

            ~FrontalInsertion(){}

            void execute();

            double lenghtMetric(const simplicesNode::SimplicesNode& nodeA, const simplicesNode::SimplicesNode& nodeB);

          private:
          SimplexMesh* m_simplexMesh = nullptr;
        };
    }
}
#endif // FRONTAL_INSERTION.H
