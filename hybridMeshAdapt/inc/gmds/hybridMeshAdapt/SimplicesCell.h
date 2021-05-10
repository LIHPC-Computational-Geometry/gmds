#ifndef SIMPLICES_CELLE_H_
#define SIMPLICES_CELLE_H_
/**************************************************************/
#include<vector>
#include<gmds/utils/CommonTypes.h>
#include<gmds/math/Vector.h>
#include<gmds/math/Point.h>
#include <Eigen/Dense>
#include <iostream>
/**************************************************************/
#include<gmds/hybridMeshAdapt/CommonInfo.h>
#include<gmds/hybridMeshAdapt/SimplicesNode.h>
/**************************************************************/
namespace gmds
{
  namespace hybrid
  {
    class SimplexMesh;

    struct tetraFacesOrientation
    {
      int FacesOrientation[4][3] = {
        {1, 3, 2},
        {0, 2, 3},
        {1, 0, 3},
        {0, 1, 2},
      };
    };

    namespace simplicesCell
    {
      class SimplicesCell
      {
      public:

        SimplicesCell                        () = delete;

        SimplicesCell                         (hybrid::SimplexMesh* simplicesMesh, const TSimplexID simplexId);

        ~SimplicesCell                        ();

        SimplicesCell                         (const SimplicesCell& simplicesCell);

        SimplicesCell                         (SimplicesCell&& simplicesCell);

        std::vector<TSimplexID> getNodes      () const ;

        double getVolumeOfCell                ();

        std::vector<TSimplexID> neighborTetra (const TInt indexNodeGlobal, bool boundariesAccepted = false) const;

        std::vector<TSimplexID> neighborTri   () const;

        std::vector<TSimplexID> neighborTri   (const TInt indexNodeGlobal) const;

        std::vector<TSimplexID> neighborTri   (const simplicesNode::SimplicesNode& simplicesNode) const;

        bool containNode                      (const simplicesNode::SimplicesNode& simplicesNode, TInt& indexLocal) const;

        bool containNode                      (const simplicesNode::SimplicesNode& simplicesNode) const;

        bool containNodes                     (const std::vector<simplicesNode::SimplicesNode>& simplicesNode) const;

        bool containNodes                     (const std::vector<TInt>& simplicesNode) const;


        bool containAtLeast                   (TInt N, const std::vector<simplicesNode::SimplicesNode>& simplicesNode) const;

        /*return the opposite tetra or errorId if there is nothing*/
        TSimplexID oppositeTetraIdx           (const simplicesNode::SimplicesNode& simplicesNode) const;

        /*return the opposite tetra or errorId if there is nothing*/
        TSimplexID oppositeTetraIdx           (const TInt indexLocal) const;

        /*return the opposite tetra or errorId if there is nothing*/
        std::vector<TSimplexID> oppositeTetraVectorPrivated(const simplicesNode::SimplicesNode& simplicesNode) const;

        /*return the tetra adjacent to this tetra*/
        std::vector<TSimplexID> adjacentTetra ();

        /*return a vector of nodes conain in this and tetra (intersection nodes of both tetra)*/
        std::vector<TInt> intersectionNodes(const SimplicesCell& simplicesCell);

        /*Reorient the Tetra (the normal is out of the Tetra ) if the normal is inside the Tetra*/
        void reorientTet();

        /*return the normal of the face (node[0] node[1] node[2] )*/
        math::Vector3d normalOfFace(const std::vector<TInt>& nodes);

        /*return the signed volume of the cell p,p1,p2,p3 with p1,p2,p3 are the node of the cell \ index */
        double signedBarycentricNormalized     (const TInt index, const gmds::math::Point& pt);

        /*return true if the normal of the tet is the same from the tet formed by the tet \ index and the point*/
        double signedBarycentric               (const TInt index, const gmds::math::Point& pt);

        /*return the signed volume of the cell p,p1,p2,p3 with p1,p2,p3 are the node of the cell \ index */
        double signedBarycentric               (const gmds::math::Point& pt, const std::vector<simplicesNode::SimplicesNode>& nodes);

        /*return true if the localIndex correspond to the generalIndex*/
        bool correspondance                    (const TInt localIndex, const TInt generalIndex) const;

        TInt getLocalNode                       (const TInt generalIndex) const;

        std::vector<TInt> getOtherNodeInSimplex             (const std::vector<TInt>& v) const;


        simplicesNode::SimplicesNode  getNode  (const TInt indexLocal) const ;

        hybrid::TSimplexID simplexId           () const {return m_simplexId;}

        std::vector<unsigned int> nodes        () const;

        bool inBorder                          () const;

        std::vector<simplicesNode::SimplicesNode> removeNodeNotInSimplex            (const std::vector<simplicesNode::SimplicesNode>& nodes) const;

        simplicesNode::SimplicesNode removeNodesFromVec       (std::vector<simplicesNode::SimplicesNode>& nodes) const;

        friend std::ostream&  operator<<(std::ostream& os, const SimplicesCell& simpliceCell)
        {
          const TInt Node0 = simpliceCell.getNode(0).getGlobalNode();
          const TInt Node1 = simpliceCell.getNode(1).getGlobalNode();
          const TInt Node2 = simpliceCell.getNode(2).getGlobalNode();
          const TInt Node3 = simpliceCell.getNode(3).getGlobalNode();

          os << "cellId = " << simpliceCell.simplexId() << "| Node : " << Node0 << " " << Node1 << " " << Node2 << " " << Node3 << std::endl;
          return os;
        }

      private:

        hybrid::SimplexMesh*        m_simplex_mesh         = nullptr;

        hybrid::TSimplexID          m_simplexId            = - 1;
      };
    }
  }
}
#endif //SIMPLICES_CELLE_H_
