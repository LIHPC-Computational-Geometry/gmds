#ifndef GMDS_MESH_DATA_PROP_TESTSUITE_H
#define GMDS_MESH_DATA_PROP_TESTSUITE_H

#include "gtest/gtest.h"
#include <gmds/services/DataField.h>
#include <gmds/services/PropertyField.h>
#include "IdleService.h"

TEST(FieldDataPropertyTestSuite, test1){
    gmds::MeshModel mod(gmds::DIM3|gmds::F|gmds::N|gmds::F2N);
    gmds::Mesh m(mod);
    gmds::DataField<double> d(&m);

    gmds::PropertyFieldBuilder builder;

    gmds::Property* p=builder.build<double>(&m,gmds::PropertyFieldBuilder::ON_NODE);

    IdleService s;
    s.addInput(&d);
    s.addConstraint(&d,p);

    ASSERT_FALSE(s.checkInput());

}


TEST(FieldDataPropertyTestSuite, test3){

    gmds::MeshModel mod(gmds::DIM3|gmds::F|gmds::N|gmds::F2N);
    gmds::Mesh m(mod);
    gmds::DataField<double> data(&m);

    gmds::PropertyFieldBuilder builder;

    gmds::Property* p=builder.build<double>(&m,gmds::PropertyFieldBuilder::ON_NODE);

    ASSERT_FALSE(p->isValid(&data));
}
#endif //GMDS_MESH_DATA_PROP_TESTSUITE_H
