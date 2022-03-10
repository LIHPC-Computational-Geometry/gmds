#ifndef CAVITY_OPERATOR_H_
#define CAVITY_OPERATOR_H_
/******************************************************************/
#include "gmds/hybridMeshAdapt/SimplicesNode.h"
#include "gmds/hybridMeshAdapt/SimplicesCell.h"
#include "gmds/hybridMeshAdapt/SimplicesTriangle.h"
/******************************************************************/
#include <set>
/******************************************************************/
namespace gmds
{
    namespace hybrid
    {
      class SimplexMesh;
      namespace operators
      {
        class CriterionRAIS;

        class CavityOperator
        {

        public:

          class CavityIO
          {
          public:
            /*Constructor that initialize the quantity of simplex in the cavity In oand Out*/
            CavityIO(SimplexMesh* simplexMesh, const std::vector<TSimplexID>& cavityIn, const std::vector<TSimplexID>& cavityTriangleConnectedToP, const std::vector<TSimplexID>& cavityTriangleNotConnectedToP);

            CavityIO(SimplexMesh* simplexMesh):m_simplex_mesh(simplexMesh)
            {

            }

            /*je supprime les constructeur par copies et d'affectation non semantic pour eviter qu'on ne les utilise sans aire expres...*/
            CavityIO(const CavityIO& cavity);

            CavityIO& operator=(const CavityIO& cavityIO);

            CavityIO(CavityIO&& cavity);

            CavityIO& operator=(CavityIO&& cavityIO);

            ~CavityIO(){};


            /*return a vector of node that are in the cavity*/
            bool                                 nodeInCavity(const TInt node);

            const std::vector<TInt> &            getNodeInCavity(){return m_nodeInCavity;}

            const std::vector<TInt> &            getSurfaceNodeInCavity(){return m_surfaceNodeInCavity;}

            void addNodeInCavity           (const TInt node){m_nodeInCavity.push_back(node);}

            const std::vector<TSimplexID>  &  cellInCavity() const {return m_cavityCellIn;}

            const std::vector<TSimplexID>  &  getTrianglesConnectedToPInCavity() const {return m_cavityTriangleConnectedToP;}

            const std::vector<TSimplexID>  &  getTrianglesNotConnectedToPInCavity() const {return m_cavityTriangleNotConnectedToP;}

            const std::vector<std::vector<TInt>>  &  getNodesToReconnect() const {return m_nodesToReconnect;}

            const std::vector<TSimplexID>  &  getOppositesSimplex() const {return m_oppositeCell;}

            const std::vector<std::vector<TInt>>  &  getTriangleReconstructionInfo() const {return m_localsNodeForReconnectionWithTriangle;}

            const std::vector<std::vector<TSimplexID>>  &  getoppositeTriangle          () const {return m_oppositeTriangle;}

            const std::vector<std::vector<TInt>>  &  getTriangleIndices           () const {return m_triangleIndices;}

            const std::vector<std::vector<TInt>>  &  getBorderEdges               () const {return m_borderSurfaceNode;}

            /*Tableau of node to connect in order to create tetraedre/triangle */
            std::vector<std::vector<TInt>> pointToConnect(std::vector<TSimplexID>& orderedAdjTet, std::set<TSimplexID>& complementarySimplex, const std::vector<TInt>& nodesInsideCavity, std::vector<TInt>& cellOfPointToConnect, std::vector<TInt>& nodeNotConnected, const TInt nodeToInsert) const ;

            void nodesReconnection();

            void findExtSimplex(std::vector<TSimplexID>& extSimplex);

            bool simplexInborderOfCavity(const TSimplexID simplex, std::vector<TInt>& nodesLocal);

            void setSimplexCavity(const std::vector<TSimplexID>& cavityCell, const std::vector<TSimplexID>& cavityTriangleConnectedToP, const std::vector<TSimplexID>& cavityTriangleNotConnectedToP){
              m_cavityCellIn = cavityCell;
              m_cavityTriangleConnectedToP = cavityTriangleConnectedToP;
              m_cavityTriangleNotConnectedToP = cavityTriangleNotConnectedToP;
            }

            bool isTetragonalizableFrom(const TInt nodeToInsert);

          private:
            /*represente the cavity In*/
            std::vector<TSimplexID> m_cavityCellIn;

            std::vector<TSimplexID> m_cavityTriangleConnectedToP;

            std::vector<TSimplexID> m_cavityTriangleNotConnectedToP;

            std::vector<std::vector<TInt>> m_nodesToReconnect;

            std::vector<TSimplexID> m_oppositeCell;

            std::vector<std::vector<TInt>> m_localsNodeForReconnectionWithTriangle;

            std::vector<std::vector<TSimplexID>> m_oppositeTriangle;

            std::vector<std::vector<TInt>> m_triangleIndices;

            std::vector<std::vector<TInt>> m_borderSurfaceNode;

            std::vector<TInt> m_nodeInCavity;

            std::vector<TInt> m_surfaceNodeInCavity;

            SimplexMesh*      m_simplex_mesh   = nullptr;
          };




          CavityOperator                     () = delete;

          CavityOperator                     (hybrid::SimplexMesh* simplexMesh);

          ~CavityOperator                    (){};

          bool cavityEnlargement(CavityIO& cavityIO, std::vector<TSimplexID>& initCavityCell, std::vector<TSimplexID>& initCavityTriangle, const simplicesNode::SimplicesNode& node,
                                  const CriterionRAIS& criterion, const std::multimap<TInt, TInt>& facesAlreadyBuilt, const std::vector<TSimplexID> markedSimplex);

          //void cavityReduction(CavityIO& cavityIO, std::vector<TSimplexID>& initCavity, const simplicesNode::SimplicesNode& node, const CriterionRAIS& criterion, const CavityReduction& cavityReduction, const std::vector<TSimplexID> v = std::vector<TSimplexID>{});

          void selectConnexTriangle(const TSimplexID& firstTriangle, const gmds::BitVector& triangleInCav, gmds::BitVector& connexTriangle);

          void fillSelectedIds(const std::vector<TSimplexID>& cavity);

          void clearSelectedIds();

        private:

          SimplexMesh*      m_simplex_mesh   = nullptr;

          BitVector         m_tetSelectedIds;

        };
    }
  }
}

//#include "../../../src/CavityOperator.cpp"


#endif //CAVITY_OPERATOR_H_
