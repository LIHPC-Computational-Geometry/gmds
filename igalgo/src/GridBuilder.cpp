/*----------------------------------------------------------------------------*/
#include <gmds/igalgo/GridBuilder.h>
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
GridBuilder::GridBuilder(Mesh* AMesh, const TInt ADim)
        :m_mesh(AMesh), m_dim(ADim)
{}
/*----------------------------------------------------------------------------*/
GridBuilder::~GridBuilder()
{}
/*----------------------------------------------------------------------------*/
bool GridBuilder::isValid() const
{
    if(m_dim==3)
        return (m_mesh->getModel()==(DIM3|R|N|R2N));
    else if(m_dim==2)
        return (m_mesh->getModel()==(DIM3|F|N|F2N) ||
                m_mesh->getModel()==(DIM2|F|N|F2N));

    //dimension error
    return false;
}
/*----------------------------------------------------------------------------*/
void GridBuilder::execute(const gmds::TInt AXNb, const gmds::TCoord AXStep, const gmds::TInt AYNb,
                          const gmds::TCoord AYStep, const gmds::TInt AZNb, const gmds::TCoord AZStep) {
    m_mesh->clear();
    if (m_dim == 2){
        build2D(AXNb, AXStep, AYNb, AYStep);
    }
    else if(m_dim==3){
        build3D(AXNb,AXStep,AYNb,AYStep,AZNb,AZStep);
    }
}
/*----------------------------------------------------------------------------*/
void GridBuilder::build2D(const gmds::TInt AXNb,
                          const gmds::TCoord AXStep,
                          const gmds::TInt AYNb,
                          const gmds::TCoord AYStep)
{
    TCellID  node_ids[AXNb][AYNb];
    for(auto x=0;x<AXNb;x++){
        for(auto y=0;y<AYNb;y++){
            node_ids[x][y]=m_mesh->newNode(x*AXStep,y*AYStep,0).id();
        }
    }
    for(auto x=0;x<AXNb-1;x++){
        for(auto y=0;y<AYNb-1;y++){
            m_mesh->newQuad(node_ids[x][y],
                            node_ids[x+1][y],
                            node_ids[x+1][y+1],
                            node_ids[x][y+1]);
        }
    }
}
/*----------------------------------------------------------------------------*/
void GridBuilder::build3D(const gmds::TInt AXNb,
                          const gmds::TCoord AXStep,
                          const gmds::TInt AYNb,
                          const gmds::TCoord AYStep,
                          const gmds::TInt AZNb,
                          const gmds::TCoord AZStep)
{
    TCellID  node_ids[AXNb][AYNb][AZNb];
    for(auto x=0;x<AXNb;x++){
        for(auto y=0;y<AYNb;y++){
            for(auto z=0;z<AZNb;z++){
                Node n = m_mesh->newNode(x*AXStep,y*AYStep,z*AZStep);
                node_ids[x][y][z]=n.id();
            }
        }
    }
    for(auto x=0;x<AXNb-1;x++){
        for(auto y=0;y<AYNb-1;y++){
            for(auto z=0;z<AZNb-1;z++){
                m_mesh->newHex(node_ids[x][y][z],
                               node_ids[x+1][y][z],
                               node_ids[x+1][y+1][z],
                               node_ids[x][y+1][z],
                               node_ids[x][y][z+1],
                               node_ids[x+1][y][z+1],
                               node_ids[x+1][y+1][z+1],
                               node_ids[x][y+1][z+1]);
            }
        }
    }
}