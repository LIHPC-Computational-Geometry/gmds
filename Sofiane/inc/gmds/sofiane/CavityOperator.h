#ifndef CAVITY_OPERATOR_H_
#define CAVITY_OPERATOR_H_
/******************************************************************/
#include "gmds/sofiane/SimplexMesh.h"
#include "gmds/sofiane/SimplicesNode.h"
#include "gmds/sofiane/SimplicesCell.h"
#include "gmds/sofiane/SimplicesTriangle.h"
/******************************************************************/
namespace gmds
{

    namespace hybrid
    {

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
            CavityIO(SimplexMesh* simplexMesh, std::vector<TSimplexID>&& cavityIn);

            CavityIO() = delete;

            /*je supprime les constructeur par copies et d'affectation non semantic pour eviter qu'on ne les utilise sans aire expres...*/
            CavityIO(const CavityIO& cavity) = delete;

            CavityIO& operator=(const CavityIO& cavityIO) = delete;

            CavityIO(CavityIO&& cavity);

            CavityIO& operator=(CavityIO&& cavityIO);

            ~CavityIO(){};


            /*return a vector of node that are in the cavity*/
            std::vector<TInt> &            nodeInCavity() ;

            std::vector<TInt> &            getNodeInCavity(){return m_nodeInCavity;}

            void addNodeInCavity           (const TInt node){m_nodeInCavity.push_back(node);}

            const std::vector<TSimplexID>  &  simplexInCavityUnmodifiable() const {return m_cavityIn;}

            std::vector<TSimplexID>  &        simplexInCavity() {return m_cavityIn;}

            /*Tableau of node to connect in order to create tetraedre */
            std::vector<std::vector<TInt>> pointToConnect();

          private:
            /*represente the cavity In*/
            std::vector<TSimplexID> m_cavityIn;

            std::vector<TInt> m_nodeInCavity;

            SimplexMesh*      m_simplex_mesh   = nullptr;

          };




          CavityOperator                     () = delete;

          CavityOperator                     (hybrid::SimplexMesh* simplexMesh);

          ~CavityOperator                    (){};

          CavityIO cavity(std::vector<TSimplexID>& initCavity, const simplicesNode::SimplicesNode& previousNode, const simplicesNode::SimplicesNode& nextNode, const CriterionRAIS& criterion, const std::vector<TSimplexID> v = std::vector<TSimplexID>{});

          CavityIO cavity(std::vector<TSimplexID>& initCavity, const simplicesNode::SimplicesNode& nextNode, const CriterionRAIS& criterion, const std::vector<TSimplexID> v = std::vector<TSimplexID>{});


        private:

          SimplexMesh*      m_simplex_mesh   = nullptr;

        };
    }
  }
}

//#include "../../../src/CavityOperator.cpp"


#endif //CAVITY_OPERATOR_H_
