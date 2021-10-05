#ifndef GMDS_MESH_DATA_PROP_TESTSUITE_H
#define GMDS_MESH_DATA_PROP_TESTSUITE_H

#include "gtest/gtest.h"
#include <gmds/services/DataMesh.h>
#include <gmds/services/PropertyMesh.h>
#include "IdleService.h"

TEST(MeshDataPropertyTestSuite, MeshPropTest){
    gmds::MeshModel mod(gmds::DIM3|gmds::F|gmds::N|gmds::F2N);
    gmds::Mesh m(mod);
    gmds::DataMesh data(&m);

    gmds::PropertyMeshBuilder builder;

    gmds::Property* prop_empty=builder.build(gmds::PropertyMeshBuilder::IS_EMPTY);

    IdleService s;
    s.addInput(&data);

    s.addConstraint(&data,prop_empty);

    ASSERT_TRUE(s.checkInput());

}
#endif //GMDS_MESH_DATA_PROP_TESTSUITE_H
