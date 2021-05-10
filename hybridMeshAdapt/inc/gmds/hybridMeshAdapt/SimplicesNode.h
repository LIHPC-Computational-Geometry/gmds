#ifndef SIMPLICES_NODE_H_
#define SIMPLICES_NODE_H_
/*----------------------------------------------------------------------------*/
#include<vector>
#include<gmds/utils/CommonTypes.h>
#include<gmds/math/Point.h>
/*----------------------------------------------------------------------------*/
#include<gmds/hybridMeshAdapt/CommonInfo.h>
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

          SimplicesNode                           (){};

          SimplicesNode                            (hybrid::SimplexMesh* simplexMesh, const TInt indexpoint);

          ~SimplicesNode                           ();

          void buildSimplicesNode                  (/*TODO add input mesh 3D*/);

          std::vector<hybrid::TSimplexID> ballOf   (bool boundariesAccepted = false) const;

          std::vector<hybrid::TSimplexID> shell    (const SimplicesNode& simplicesNode) const;

          bool                            isInBorder() const;

          TInt getGlobalNode                       () const {return m_indexPoint;}

          math::Point getCoords                    () const;

          std::vector<hybrid::TSimplexID> linksTri () const;

          std::vector<hybrid::TSimplexID> triangleUnordererdShell (const SimplicesNode& simplicesNode) const;

          //return the shell of Triangle Orderer)
          std::vector<TSimplexID>  shellTriangleOrdererdWithHole (const SimplicesNode& simplicesNode) const;

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
