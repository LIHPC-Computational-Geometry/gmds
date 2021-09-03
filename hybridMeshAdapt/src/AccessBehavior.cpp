#include "gmds/hybridMeshAdapt/AccessBehavior.h"


using namespace gmds;
using namespace accessBehavior;


SparseBehavior::SparseBehavior(std::unique_ptr<Eigen::SparseMatrix<TInt>>&& sparseIndex):
m_sparseIndex(std::move(sparseIndex))
{
  sparseIndex = nullptr;
}


SparseBehavior::SparseBehavior(SparseBehavior& sparseBehavior):
m_sparseIndex(std::move(sparseBehavior.m_sparseIndex))
{
  sparseBehavior.m_sparseIndex = nullptr;
}

SparseBehavior::SparseBehavior(SparseBehavior&& sparseBehavior):
m_sparseIndex(std::move(sparseBehavior.m_sparseIndex))
{
  sparseBehavior.m_sparseIndex = nullptr;
}

SparseBehavior::~SparseBehavior()
{

}



const TInt SparseBehavior::neighborCellOf(const TInt pointIndex, const TInt cellIndex)
{
  //TODO..

  return static_cast<TInt>(0);
}



std::vector<TInt> SparseBehavior::ballOf(const TInt pointIndex)
{
  std::vector<TInt> ballIndexCell;
  const TInt nbrIndexBall = m_sparseIndex->col(pointIndex).sum(); //!!!! IMPORTANT IL FAUT ABSOLUMENT QUE LES VALEURS DE M_SPARSEINDEX SOIT TOUT EGAUX A 1 // AFIN QUE LE RESIZONG DE BALLINDEX SOIT CORRECT!!!!!//////
  ballIndexCell.resize(nbrIndexBall);

  unsigned int indexTmp = 0;
  for (Eigen::SparseMatrix<int>::InnerIterator it(*m_sparseIndex,pointIndex); it; ++it)
  {
    ballIndexCell[indexTmp] = it.value();
    indexTmp++;
  }

  //a voir si le return utilise bien le constructeur de deplacement et non pas de copie !!
  return std::move(ballIndexCell);
}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


MapBehavior::MapBehavior()
{

}

MapBehavior::MapBehavior(MapBehavior&& mapBehavior):
m_N2C(std::move(mapBehavior.m_N2C))
{
  mapBehavior.m_N2C = nullptr;
}



MapBehavior::MapBehavior(std::unique_ptr<N2C>&& mapBehavior):
m_N2C(std::move(mapBehavior))
{
  mapBehavior = nullptr;
}


MapBehavior::~MapBehavior()
{

}



const TInt MapBehavior::neighborCellOf(const TInt pointIndex, const TInt cellIndex)
{
  //à Optimisation
  N2C::iterator iter = m_N2C->find(pointIndex);
  //premier if permert de trouver les points adjacents aux tera/triangle que l'on cherche
  if(iter != m_N2C->end())
  {
    std::vector<TInt> v(3); // a voir si les triangles ne poseront pas de soucis....
    unsigned int idxPtAdj = 0;
    for(auto const& it : *m_N2C)
    {
      if(it != *iter)//pour éviter de repasser dans le point d'ou l'on cherche le tera/triangle voisin
      {
        if(it.second.find(cellIndex) != it.second.end());
        {
          v[idxPtAdj] = it.first;
          idxPtAdj++;
        }
      }
    }



    //On retrouve maintenant l'index du triangle/tetra que l'on cherche grâce aux indices dans v trouvé précedement
    const std::set<TInt> & set1 = m_N2C->find(v[0])->second;
    const std::set<TInt> & set2 = m_N2C->find(v[1])->second;
    const std::set<TInt> & set3 = m_N2C->find(v[2])->second;
  }


}



std::vector<TInt> MapBehavior::ballOf(const TInt pointIndex)
{
  N2C::iterator it = m_N2C->find(pointIndex);
  if(it != m_N2C->end())
  {
    std::vector<TInt> v(it->second.begin(), it->second.end());
    return v;
  }

  return std::vector<TInt>({});
}



void MapBehavior::buildDatabehavior(/*const gmds::Mesh mesh*/ )
{

}
