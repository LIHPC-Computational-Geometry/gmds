#include <iostream>
#include <gmds/ig/Mesh.h>

using namespace gmds;

int main() {
        Mesh m(MeshModel(DIM3|F|N|F2N));

        Node n0 = m.newNode(0,0,0);
        Node n1 = m.newNode(0,10,0);
        Node n2 = m.newNode(10,10,0);
        Node n3 = m.newNode(10,0,0);

        Node n4 = m.newNode(0,0,10);
        Node n5 = m.newNode(0,10,10);
        Node n6 = m.newNode(10,10,10);
        Node n7 = m.newNode(10,0,10);

        m.newTriangle(n0,n1,n2);
        m.newTriangle(n0,n2,n3);
        m.newTriangle(n4,n6,n5);
        m.newTriangle(n4,n7,n6);
        m.newQuad(n0,n1,n2,n3);

        return 0;
}
