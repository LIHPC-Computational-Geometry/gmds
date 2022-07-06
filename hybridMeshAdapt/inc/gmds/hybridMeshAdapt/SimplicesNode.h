#ifndef SIMPLICES_NODE_H_
#define SIMPLICES_NODE_H_
/*----------------------------------------------------------------------------*/
#include<gmds/utils/CommonTypes.h>
#include<gmds/utils/BitVector.h>
/*----------------------------------------------------------------------------*/
#include<gmds/math/Point.h>
#include<gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
#include<gmds/hybridMeshAdapt/ICriterion.h>
#include<gmds/hybridMeshAdapt/CommonInfo.h>
/*----------------------------------------------------------------------------*/
#include<map>
#include<vector>
/*----------------------------------------------------------------------------*/
namespace gmds
{
  namespace hybrid
  {
      class SimplexMesh;

      namespace simplicesNode
      {
        class SimplicesNode
        {
        public:

          ////////////////////////////////
          struct nodeNeighborInfo{
            TSimplexID  cell;
            std::vector<TInt> nodes{};
          };
          ////////////////////////////////

          SimplicesNode                           (){};

          SimplicesNode                            (hybrid::SimplexMesh* simplexMesh, const TInt indexpoint);

          ~SimplicesNode                           ();

          bool checkExistance                      ();

          void buildSimplicesNode                  (/*TODO add input mesh 3D*/);

          std::vector<hybrid::TSimplexID> ballOf   (bool boundariesAccepted = false) const;

          std::vector<hybrid::TSimplexID> shell    (const SimplicesNode& simplicesNode) const;

          bool                            isInBorder() const;

          TInt getGlobalNode                       () const {return m_indexPoint;}

          math::Point    getCoords                    () const;

          math::Vector3d getNormal                    () const;

          std::vector<hybrid::TSimplexID> linksTri () const;

          std::vector<hybrid::TSimplexID> triangleUnordererdShell (const SimplicesNode& simplicesNode) const;

          //return the shell of Triangle Orderer)
          std::vector<TSimplexID>  shellTriangleOrdererdWithHole (const SimplicesNode& simplicesNode) const;

          std::vector<TInt>  commonNodeInBall(const std::vector<TInt>& nodeVec) const;

          TSimplexID directSimplex(const math::Vector3d& vector) const ;

          std::vector<TInt> neighborNodes() const;

          std::vector<TSimplexID> directSimplices(const math::Vector3d& vector) const ;

          /*return the ordered simplex that are intersected by the edge formed by this ans the argument simpliceNodeB */
          std::vector<TSimplexID> simplicesIntersectedbyRay(const SimplicesNode& simpliceNodeB) const ;

          /*return the ordered simplex that are intersected by the edge formed by *this and the argument simpliceNodeB, but return when there is more that one face in the cavity that is not visible by the node B */
          std::vector<TSimplexID> simplicesIntersectedbyRayAndCheckVisibility(std::vector<std::vector<TInt>>& facesId, const SimplicesNode& simpliceNodeB) const ;

          /*check this can see every border of this cavity*/
          bool checkStar(const std::vector<TSimplexID>& cavity, const operators::CriterionRAIS& criterion) const ;

          bool isFaceVisible(math::Vector3d& normalOfFace, math::Vector3d& vecFacePt) const;

          bool isAttachToSimplex() const ;

          std::vector<TInt> complentaryNodeShell (const SimplicesNode& simpliceNode) const;

          std::vector<TInt> directNeighboorNodeId() const ;

          /*adding some point to the mesh and check if  this point already exist*/
          void detectType(const nodeNeighborInfo& nodeInfo) const;

          SimplicesNode& operator=                         (const SimplicesNode& simpliceNodeA)
          {
            this->m_simplex_mesh = simpliceNodeA.m_simplex_mesh;
            this->m_indexPoint   = simpliceNodeA.m_indexPoint;

            return *this;
          }

          bool  operator==                         (const SimplicesNode& simpliceNodeA)
          {
            bool flag = false;
            if(simpliceNodeA.m_simplex_mesh == this->m_simplex_mesh)
            {
              if(simpliceNodeA.m_indexPoint == this->m_indexPoint)
              {
                flag = true;
              }
            }
            return flag;
          }

          friend std::ostream&  operator<<(std::ostream& os, const SimplicesNode& simpliceNode)
          {
            os << "Node = " << simpliceNode.getGlobalNode() << "| Coords : " << simpliceNode.getCoords()<< std::endl;
            return os;
          }

        private:

          hybrid::SimplexMesh* m_simplex_mesh = nullptr;

          TInt                   m_indexPoint =      -1;

        };
      }
  }
}


#endif //SIMPLEX_NODE_H_
