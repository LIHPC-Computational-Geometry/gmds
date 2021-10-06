/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testMeshNode) {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));

    gmds::Node n0 = m.newNode(0, 0, 0);
    ASSERT_EQ(n0.point(), gmds::math::Point(0, 0, 0));
    n0.setX(1);
    ASSERT_EQ(n0.X(),1);
    n0.setY(2);
    ASSERT_EQ(n0.Y(),2);
    n0.setZ(3);
    ASSERT_EQ(n0.Z(),3);
    n0.setXYZ(4,5,6);
    ASSERT_EQ(n0.point(), gmds::math::Point(4, 5, 6));
    n0.setPoint(gmds::math::Point(1,1,1));
    ASSERT_EQ(n0.point(), gmds::math::Point(1, 1, 1));

}
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testMeshNodeAndFaceIterate)
    {
        gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::F2N));

        gmds::Node n0 = m.newNode(0,0,0);
        gmds::Node n1 = m.newNode(1,1,1);
    gmds::Node n2 = m.newNode(0,1,1);
    gmds::Node n3 = m.newNode(1,0,1);

    m.newTriangle(n0,n1,n3);
    m.newTriangle(n0,n3,n2);
    m.newQuad(n0,n1,n2,n3);

    std::cout<<" ==== nodes ==== "<<std::endl;
    for(auto i: m.nodes()){
        gmds::Node ni = m.get<gmds::Node>(i);
        std::cout<<"Node "<< ni.id()<<std::endl;
    }

    std::cout<<" ==== faces ==== "<<std::endl;
    for(auto i:m.faces()){
        gmds::Face fi = m.get<gmds::Face>(i);
        std::cout<<"Face "<< fi.id()<<" of type "<< fi.type()<<" with nodes "<<std::endl;
        std::vector<gmds::TCellID > n_ids = fi.getIDs<gmds::Node>();
        for(auto j:n_ids){
            std::cout<<" "<<j;
        }
        std::cout<<std::endl;
    }
    std::cout<<" ==== nodes ==== "<<std::endl;
    for(gmds::Mesh::nodes_iterator itn = m.nodes_begin(); itn!=m.nodes_end(); ++itn){
        std::cout<<"Node "<<*itn<<std::endl;
    }
	m.deleteNode(2);


    for(auto i:m.faces()){
        gmds::Face fi = m.get<gmds::Face>(i);
        std::vector<gmds::TCellID > n_ids = fi.getIDs<gmds::Node>();
        for(auto j:n_ids){
            std::cout<<" "<<j;
        }
        std::cout<<std::endl;
    }


    ASSERT_EQ(m.getNbNodes(),3);

}


/*----------------------------------------------------------------------------*/
TEST(MeshClass, testNodeVariableAndGroups) {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));

    gmds::Node n0 = m.newNode(0, 0, 0);
    gmds::Node n1 = m.newNode(1, 1, 1);
    gmds::Node n2 = m.newNode(0, 1, 1);
    gmds::Node n3 = m.newNode(1, 0, 1);

    n1.id();

    m.newTriangle(n0, n1, n3);
    m.newTriangle(n0, n3, n2);
    m.newQuad(n0, n1, n2, n3);

    gmds::Variable<int>* var = m.newVariable<int, gmds::GMDS_NODE>( "int_var");

    for (auto i: m.nodes()) {
        ASSERT_EQ(var->value(i), 0);
        var->value(i) = i;
    }

    for (auto i: m.nodes()) {
        ASSERT_EQ((*var)[i], i);
    }
    gmds::Variable<int>* var2 = m.getVariable<int, gmds::GMDS_NODE>( "int_var");

    for (auto i: m.nodes()) {
        ASSERT_EQ((*var2)[i], i);
    }

    gmds::CellGroup<gmds::Node>* grp = m.newGroup<gmds::Node>("cloud");
    grp->add(n0);
    grp->add(n1);
    grp->add(n2.id());
    ASSERT_EQ(grp->size(),3);

    ASSERT_EQ((*grp)[0],n0.id());
    ASSERT_EQ((*grp)[1],n1.id());
    ASSERT_EQ((*grp)[2],n2.id());

    ASSERT_EQ(grp->value(0),n0.id());
    ASSERT_EQ(grp->value(1),n1.id());
    ASSERT_EQ(grp->value(2),n2.id());
}
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testAccessorAll) {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E
                                 | gmds::F2N | gmds::F2E | gmds::E2N
                                 | gmds::E2F | gmds::N2E| gmds::N2F) );

    ASSERT_EQ(3, m.getDim());

    gmds::Node n0 = m.newNode(0, 0, 1);
    gmds::Node n1 = m.newNode(1, 0, 1);
    gmds::Node n2 = m.newNode(1, 1, 1);
    gmds::Node n3 = m.newNode(0, 1, 1);
    gmds::Node n4 = m.newNode(2, 0, 1);
    gmds::Node n5 = m.newNode(2, 1, 1);

    n1.id();

    gmds::Face f1 = m.newTriangle(n0, n1, n3);
    gmds::Face f2 = m.newTriangle(n1, n3, n2);
    gmds::Face f3 = m.newQuad(n1, n2, n5, n4);


    ASSERT_EQ(f1.getAllIDs<gmds::Edge>().size(),3);
    ASSERT_EQ(f2.getAllIDs<gmds::Edge>().size(),3);
    ASSERT_EQ(f3.getAllIDs<gmds::Edge>().size(),4);
    ASSERT_EQ(f1.getIDs<gmds::Edge>().size(),0);
    ASSERT_EQ(f2.getIDs<gmds::Edge>().size(),0);
    ASSERT_EQ(f3.getIDs<gmds::Edge>().size(),0);

    ASSERT_EQ(n0.getIDs<gmds::Edge>().size(),0);
    ASSERT_EQ(n0.getAllIDs<gmds::Edge>().size(),0);
    ASSERT_EQ(n0.getIDs<gmds::Face>().size(),0);
    ASSERT_EQ(n0.getAllIDs<gmds::Face>().size(),0);

    gmds::MeshDoctor doc(&m);

    doc.buildEdgesAndX2E();

    ASSERT_EQ(f1.getIDs<gmds::Edge>().size(),3);
    ASSERT_EQ(f2.getIDs<gmds::Edge>().size(),3);
    ASSERT_EQ(f3.getIDs<gmds::Edge>().size(),4);

    gmds::Edge e0 = m.get<gmds::Edge>(*m.edges_begin());

    ASSERT_EQ(e0.getIDs<gmds::Face>().size(),0);
    ASSERT_EQ(e0.getAllIDs<gmds::Face>().size(),0);
    ASSERT_EQ(e0.getIDs<gmds::Node>().size(),2);
    ASSERT_EQ(e0.getAllIDs<gmds::Node>().size(),2);

    doc.updateUpwardConnectivity();

    ASSERT_EQ(f1.getAllIDs<gmds::Node>().size(),3);
    ASSERT_EQ(f2.getAllIDs<gmds::Node>().size(),3);
    ASSERT_EQ(f3.getAllIDs<gmds::Node>().size(),4);



    ASSERT_EQ(n0.getIDs<gmds::Edge>().size(),2);
    ASSERT_EQ(n0.getAllIDs<gmds::Edge>().size(),2);
    ASSERT_EQ(n0.getIDs<gmds::Face>().size(),1);
    ASSERT_EQ(n0.getAllIDs<gmds::Face>().size(),1);

    ASSERT_EQ(n1.getIDs<gmds::Edge>().size(),4);
    ASSERT_EQ(n1.getAllIDs<gmds::Edge>().size(),4);
    ASSERT_EQ(n1.getIDs<gmds::Face>().size(),3);
    ASSERT_EQ(n1.getAllIDs<gmds::Face>().size(),3);
}
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testTetra) {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::N | gmds::F | gmds::E |
                                 gmds::F2N | gmds::R2N| gmds::E2N |
                                 gmds::R2F | gmds::F2R | gmds::R2E ));

    gmds::Node n0 = m.newNode(0, 0, 0);
    gmds::Node n1 = m.newNode(0, 1, 0);
    gmds::Node n2 = m.newNode(1, 1, 0);
    gmds::Node n3 = m.newNode(0, 1, 0);

    gmds::Node n4 = m.newNode(0, 0, 1);
    gmds::Node n5 = m.newNode(0, 1, 1);
    gmds::Node n6 = m.newNode(1, 1, 1);
    gmds::Node n7 = m.newNode(0, 1, 1);


    gmds::Region r_null;
    ASSERT_EQ(r_null.id(),gmds::NullID);


    gmds::Region h =m.newHex(n0, n1, n2, n3, n4, n5, n6, n7);
    ASSERT_EQ(h.nbNodes(),8);
    ASSERT_EQ(h.get<gmds::Node>().size(),8);
    ASSERT_EQ(h.getIDs<gmds::Node>().size(),8);

    gmds::Region t =m.newTet(n0, n1, n2, n5);
    ASSERT_EQ(t.nbNodes(),4);
    ASSERT_EQ(t.get<gmds::Node>().size(),4);
    ASSERT_EQ(t.getIDs<gmds::Node>().size(),4);

    gmds::Region pr =m.newPrism3(n0, n1, n3, n4, n5, n7);
    ASSERT_EQ(pr.nbNodes(),6);
    ASSERT_EQ(pr.get<gmds::Node>().size(),6);
    ASSERT_EQ(pr.getIDs<gmds::Node>().size(),6);

    gmds::Region py =m.newPyramid(n0, n1, n2, n3, n4);
    ASSERT_EQ(py.nbNodes(),5);
    ASSERT_EQ(py.get<gmds::Node>().size(),5);
    ASSERT_EQ(py.getIDs<gmds::Node>().size(),5);

    ASSERT_FALSE(py==pr);
    ASSERT_TRUE(t==t);
    ASSERT_TRUE(py!=pr);
    ASSERT_FALSE(t!=t);

}
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testHexa) {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::N | gmds::F | gmds::E |
                                 gmds::F2N | gmds::R2N| gmds::E2N |
                                 gmds::R2F | gmds::F2R | gmds::R2E ));

    gmds::Node n0 = m.newNode(0, 0, 0);
    gmds::Node n1 = m.newNode(0, 1, 0);
    gmds::Node n2 = m.newNode(1, 1, 0);
    gmds::Node n3 = m.newNode(0, 1, 0);

    gmds::Node n4 = m.newNode(0, 0, 1);
    gmds::Node n5 = m.newNode(0, 1, 1);
    gmds::Node n6 = m.newNode(1, 1, 1);
    gmds::Node n7 = m.newNode(0, 1, 1);


    gmds::Region r_null;
    ASSERT_EQ(r_null.id(),gmds::NullID);


    gmds::Region h =m.newHex(n0, n1, n2, n3, n4, n5, n6, n7);
    ASSERT_EQ(h.nbNodes(),8);
    ASSERT_EQ(h.get<gmds::Node>().size(),8);
    ASSERT_EQ(h.getIDs<gmds::Node>().size(),8);

    ASSERT_EQ(h.nbEdges(),12);


    gmds::MeshDoctor doc(&m);
    doc.buildEdgesAndX2E();
    doc.buildFacesAndR2F();
    doc.updateUpwardConnectivity();

    ASSERT_EQ(h.get<gmds::Edge>().size(),12);
    ASSERT_EQ(h.getIDs<gmds::Edge>().size(),12);
    ASSERT_EQ(h.get<gmds::Face>().size(),6);
    ASSERT_EQ(h.getIDs<gmds::Face>().size(),6);
}
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testBndBox) {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::N | gmds::F | gmds::E |
                                 gmds::F2N | gmds::R2N| gmds::E2N ));

    gmds::Node n0 = m.newNode(0, 0, 0);
    gmds::Node n1 = m.newNode(0, 1, 0);
    gmds::Node n2 = m.newNode(1, 1, 0);
    gmds::Node n3 = m.newNode(0, 1, 0);

    gmds::Node n4 = m.newNode(0, 0, 1);
    gmds::Node n5 = m.newNode(0, 1, 1);
    gmds::Node n6 = m.newNode(1, 1, 1);
    gmds::Node n7 = m.newNode(0, 1, 1);



    gmds::Region h =m.newHex(n0, n1, n2, n3, n4, n5, n6, n7);

    gmds::TCoord minBB[3], maxBB[3];
    h.computeBoundingBox(minBB,maxBB);

    ASSERT_DOUBLE_EQ(minBB[0],0);
    ASSERT_DOUBLE_EQ(minBB[1],0);
    ASSERT_DOUBLE_EQ(minBB[2],0);

    ASSERT_DOUBLE_EQ(maxBB[0],1);
    ASSERT_DOUBLE_EQ(maxBB[1],1);
    ASSERT_DOUBLE_EQ(maxBB[2],1);

    gmds::Edge e = m.newEdge(n0,n5);

    e.computeBoundingBox(minBB,maxBB);

    ASSERT_DOUBLE_EQ(minBB[0],0);
    ASSERT_DOUBLE_EQ(minBB[1],0);
    ASSERT_DOUBLE_EQ(minBB[2],0);

    ASSERT_DOUBLE_EQ(maxBB[0],0);
    ASSERT_DOUBLE_EQ(maxBB[1],1);
    ASSERT_DOUBLE_EQ(maxBB[2],1);

    gmds::Face f = m.newTriangle(n2,n4,n7);

    f.computeBoundingBox(minBB,maxBB);

    ASSERT_DOUBLE_EQ(minBB[0],0);
    ASSERT_DOUBLE_EQ(minBB[1],0);
    ASSERT_DOUBLE_EQ(minBB[2],0);

    ASSERT_DOUBLE_EQ(maxBB[0],1);
    ASSERT_DOUBLE_EQ(maxBB[1],1);
    ASSERT_DOUBLE_EQ(maxBB[2],1);

    ASSERT_ANY_THROW(
            n0.computeBoundingBox(minBB,maxBB);
    );

}

/*----------------------------------------------------------------------------*/
TEST(MeshClass, testCellReplaceAndSetForNodes) {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::N | gmds::F | gmds::E |
                                 gmds::F2N | gmds::R2N| gmds::E2N ));

    gmds::Node n0 = m.newNode(0, 0, 0);
    gmds::Node n1 = m.newNode(0, 1, 0);
    gmds::Node n2 = m.newNode(1, 1, 0);
    gmds::Node n3 = m.newNode(0, 1, 0);

    gmds::Node n4 = m.newNode(0, 0, 1);
    gmds::Node n5 = m.newNode(0, 1, 1);
    gmds::Node n6 = m.newNode(1, 1, 1);
    gmds::Node n7 = m.newNode(0, 1, 1);


    gmds::Node n8 = m.newNode(0,2,2);

    gmds::Region h =m.newHex(n0, n1, n2, n3, n4, n5, n6, n7);

    h.replace(n7,n8);

    std::vector<gmds::TCellID > hn = h.getIDs<gmds::Node>();
    ASSERT_EQ(hn.size(),8);

    ASSERT_TRUE(std::find(hn.begin(), hn.end(),n7.id())==hn.end());
    ASSERT_FALSE(std::find(hn.begin(), hn.end(),n8.id())==hn.end());

    std::vector<gmds::Node> new_n;
    std::vector<gmds::TCellID > new_nid;
    new_n.push_back(n7);
    new_n.push_back(n6);
    new_n.push_back(n5);
    new_n.push_back(n4);
    new_n.push_back(n3);
    new_n.push_back(n2);
    new_n.push_back(n1);
    new_n.push_back(n0);
    h.set(new_n);
    hn = h.getIDs<gmds::Node>();
    ASSERT_EQ(hn.size(),8);

    new_nid.push_back(n0.id());
    new_nid.push_back(n1.id());
    new_nid.push_back(n2.id());
    new_nid.push_back(n3.id());
    new_nid.push_back(n4.id());
    new_nid.push_back(n5.id());
    new_nid.push_back(n6.id());
    new_nid.push_back(n7.id());


    h.set<gmds::Node>(new_nid);
    hn = h.getIDs<gmds::Node>();
    ASSERT_EQ(hn.size(),8);
    ASSERT_EQ(hn[0],n0.id());
    ASSERT_EQ(hn[2],n2.id());
    ASSERT_EQ(hn[4],n4.id());
    ASSERT_EQ(hn[5],n5.id());

    h.remove(n0);
    h.remove(n0);
    ASSERT_EQ(h.get<gmds::Node>().size(),7);


    h.add(n0);

    ASSERT_EQ(h.get<gmds::Node>().size(),8);

    ASSERT_ANY_THROW(
            h.add(n0);
    );
}

/*----------------------------------------------------------------------------*/
TEST(MeshClass, testCellReplaceAndSetForEF) {

    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::N | gmds::F | gmds::E |
                                 gmds::F2N | gmds::R2N| gmds::E2N |
                                 gmds::R2F | gmds::F2R | gmds::R2E ));

    gmds::Node n0 = m.newNode(0, 0, 0);
    gmds::Node n1 = m.newNode(0, 1, 0);
    gmds::Node n2 = m.newNode(1, 1, 0);
    gmds::Node n3 = m.newNode(0, 1, 1);

    gmds::Node n4 = m.newNode(0, 1, 2);
    gmds::Edge new_e = m.newEdge(n0,n4);
    gmds::Face new_f = m.newTriangle(n0,n1,n4);

    gmds::Region t =m.newTet(n0, n1, n2, n3);


    gmds::MeshDoctor doc(&m);
    doc.buildEdgesAndX2E();
    doc.buildFacesAndR2F();
    doc.updateUpwardConnectivity();

    std::vector<gmds::TCellID > te = t.getIDs<gmds::Edge>();
    std::vector<gmds::TCellID > tf = t.getIDs<gmds::Face>();
    ASSERT_EQ(te.size(),6);
    ASSERT_EQ(tf.size(),4);

    gmds::Edge e0 = m.get<gmds::Edge>(te[0]);
    t.replace(e0,new_e);
    ASSERT_FALSE(std::find(te.begin(), te.end(), e0.id())==te.end());
    ASSERT_TRUE (std::find(te.begin(), te.end(), new_e.id())==te.end());

    gmds::Face f0 = m.get<gmds::Face>(tf[0]);
    t.replace(f0,new_f);
    ASSERT_FALSE(std::find(tf.begin(), tf.end(), f0.id())==tf.end());
    ASSERT_TRUE (std::find(tf.begin(), tf.end(), new_f.id())==tf.end());

    std::vector<gmds::Edge> new_te;
    new_te.resize(te.size());
    new_te[0]=m.get<gmds::Edge>(te[5]);
    new_te[1]=m.get<gmds::Edge>(te[4]);
    new_te[2]=m.get<gmds::Edge>(te[3]);
    new_te[3]=m.get<gmds::Edge>(te[2]);
    new_te[4]=m.get<gmds::Edge>(te[1]);
    new_te[5]=m.get<gmds::Edge>(te[0]);

    t.set(new_te);
    new_te = t.get<gmds::Edge>();

    ASSERT_EQ(new_te[0].id(),te[5]);
    ASSERT_EQ(new_te[1].id(),te[4]);
    ASSERT_EQ(new_te[2].id(),te[3]);
    ASSERT_EQ(new_te[3].id(),te[2]);
    ASSERT_EQ(new_te[4].id(),te[1]);
    ASSERT_EQ(new_te[5].id(),te[0]);
}
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testCellCenter) {

    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::N | gmds::F | gmds::E |
                                 gmds::F2N | gmds::R2N | gmds::E2N ));

    gmds::Node n0 = m.newNode(0, 0, 0);
    gmds::Node n1 = m.newNode(0, 1, 0);
    gmds::Node n2 = m.newNode(1, 1, 0);
    gmds::Node n3 = m.newNode(1, 0, 0);

    gmds::Node n4 = m.newNode(0, 0, 1);
    gmds::Node n5 = m.newNode(0, 1, 1);
    gmds::Node n6 = m.newNode(1, 1, 1);
    gmds::Node n7 = m.newNode(1, 0, 1);

    gmds::Node n8 = m.newNode(0.5, 0.5, 1);


    gmds::math::Point p = n0.center();
    ASSERT_EQ(p, n0.point());

    gmds::Edge e = m.newEdge(n0,n1);
    ASSERT_EQ(e.center(),gmds::math::Point(0,0.5,0));

    std::vector<gmds::Node> quad_nodes;
    quad_nodes.push_back(n0);
    quad_nodes.push_back(n1);
    quad_nodes.push_back(n2);
    quad_nodes.push_back(n3);
    gmds::Face q = m.newFace(quad_nodes);
    ASSERT_EQ(q.center(),gmds::math::Point(0.5,0.5,0));

    gmds::Region h = m.newHex(n0,n1,n2,n3,n4,n5,n6,n7);
    ASSERT_EQ(h.center(),gmds::math::Point(0.5,0.5,0.5));

}
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testMeshGlobalQueries) {

    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::N | gmds::F | gmds::E |
                                 gmds::F2N | gmds::R2N | gmds::E2N ));

    gmds::Node n0 = m.newNode(0, 0, 0);
    gmds::Node n1 = m.newNode(0, 1, 0);
    gmds::Node n2 = m.newNode(1, 1, 0);
    gmds::Node n3 = m.newNode(1, 0, 0);

    gmds::Node n4 = m.newNode(0, 0, 1);
    gmds::Node n5 = m.newNode(0, 1, 1);
    gmds::Node n6 = m.newNode(1, 1, 1);
    gmds::Node n7 = m.newNode(1, 0, 1);

    gmds::Node n8 = m.newNode(0.5, 0.5, 1);


    gmds::Edge e = m.newEdge(n0,n1);
    ASSERT_EQ(e.nbNodes(),2);
    ASSERT_EQ(e.length(),1);
    gmds::Face q = m.newQuad(n0,n1,n2,n3);
    ASSERT_EQ(q.nbNodes(),4);

    gmds::Region h = m.newHex(n0,n1,n2,n3,n4,n5,n6,n7);
    ASSERT_EQ(h.nbNodes(),8);


    ASSERT_TRUE(m.has<gmds::Node>(n3.id()));
    ASSERT_TRUE(m.has<gmds::Edge>(e.id()));
    ASSERT_TRUE(m.has<gmds::Face>(q.id()));
    ASSERT_TRUE(m.has<gmds::Region>(h.id()));



    ASSERT_TRUE(m.has<gmds::Node>(n8.id()));
    m.deleteNode(n8);
    ASSERT_FALSE(m.has<gmds::Node>(n8.id()));
    ASSERT_EQ(m.getNbNodes(),8);


    ASSERT_EQ(m.getNbEdges(),1);
    m.deleteEdge(e);
    ASSERT_EQ(m.getNbEdges(),0);

    ASSERT_EQ(m.getNbFaces(),1);
    m.deleteFace(q);
    ASSERT_EQ(m.getNbFaces(),0);

    ASSERT_EQ(m.getNbRegions(),1);
    m.deleteRegion(h);
    ASSERT_EQ(m.getNbRegions(),0);

}

/*----------------------------------------------------------------------------*/
TEST(MeshClass, testFaceInfo) {

    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::F | gmds::F2N ));

    gmds::Node n0 = m.newNode(0, 0, 0);
    gmds::Node n1 = m.newNode(0, 1, 0);
    gmds::Node n2 = m.newNode(1, 1, 0);
    gmds::Node n3 = m.newNode(1, 0, 0);

    gmds::Face q = m.newQuad(n0, n1, n2, n3);
    ASSERT_EQ(q.area(), 1);
    gmds::math::Vector3d n(0,0,1);
    ASSERT_EQ(std::abs(q.normal().dot(n)),1);
    ASSERT_EQ(std::abs(q.normal(n1).dot(n)),1);

    gmds::Face t = m.newTriangle(n0, n1, n2);
    ASSERT_EQ(t.area(), 0.5);
    ASSERT_EQ(std::abs(t.normal().dot(n)),1);
    ASSERT_EQ(std::abs(t.normal(n0).dot(n)),1);

    ASSERT_EQ(q.nbNodes(),4);
    ASSERT_EQ(t.nbNodes(),3);

    gmds::math::Point p(0.2,0.2,2);
    gmds::math::Point p_on_t = t.project(p);
    gmds::math::Point p_on_q = q.project(p);

    ASSERT_EQ(p_on_q.X(),0.2);
    ASSERT_EQ(p_on_q.Y(),0.2);
    ASSERT_EQ(p_on_q.Z(),0);

    ASSERT_EQ(p_on_t.X(),0.2);
    ASSERT_EQ(p_on_t.Y(),0.2);
    ASSERT_EQ(p_on_t.Z(),0);

    gmds::Node a0,a1;
    q.getAdjacentNodes(n0,a0,a1);
    ASSERT_EQ(a0,n3);
    ASSERT_EQ(a1,n1);

    q.getAdjacentNodes(n1,a0,a1);
    ASSERT_EQ(a0,n0);
    ASSERT_EQ(a1,n2);

    q.getAdjacentNodes(n2,a0,a1);
    ASSERT_EQ(a0,n1);
    ASSERT_EQ(a1,n3);
    q.getAdjacentNodes(n3,a0,a1);
    ASSERT_EQ(a0,n2);
    ASSERT_EQ(a1,n0);

    ASSERT_EQ(std::abs(q.computeScaledJacobian2D()),1);
    ASSERT_ANY_THROW(t.computeScaledJacobian2D(););


}
