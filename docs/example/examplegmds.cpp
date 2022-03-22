#include <iostream>

// gmds headers
#include <gmds/ig/Mesh.h>

int main()
{
  std::cout<<"Example of linking against gmds"<<std::endl;
  gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::F2N));
  gmds::Node n0 = m.newNode(0,0,0);
  gmds::Node n1 = m.newNode(1,1,1);
  gmds::Node n2 = m.newNode(0,1,1);
  gmds::Node n3 = m.newNode(1,0,1);

  m.newTriangle(n0,n1,n3);
  m.newTriangle(n0,n3,n2);
  m.newQuad(n0,n1,n2,n3);

  std::cout<<"nbNodes "<<m.getNbNodes()<<std::endl;
  std::cout<<"nbNodes "<<m.getNbFaces()<<std::endl;
}
