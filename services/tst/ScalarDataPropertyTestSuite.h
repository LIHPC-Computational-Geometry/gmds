#ifndef GMDS_SCALAR_DATA_PROP_TESTSUITE_H
#define GMDS_SCALAR_DATA_PROP_TESTSUITE_H

#include "gtest/gtest.h"
#include <gmds/services/DataScalar.h>
#include <gmds/services/PropertyScalar.h>
#include <gmds/services/AbstractService.h>
#include "IdleService.h"

TEST(ScalarDataPropertyTestSuite, DataPropTest){

    gmds::DataScalar<int> data(2);

    gmds::PropertyScalarBuilder builder;
    gmds::Property* prop1 = builder.build<int>(gmds::PropertyScalarBuilder::POSITIVE_STRICTLY);
    gmds::Property* prop2 = builder.build<int>(gmds::PropertyScalarBuilder::POSITIVE);
    gmds::PropertyScalarInRange<int>* prop3 = builder.buildInRange<int>();
    prop3->setRange(1,3);
    ASSERT_EQ(data.value(),2);
    ASSERT_TRUE(prop1->isValid(&data));
    ASSERT_TRUE(prop2->isValid(&data));
    ASSERT_TRUE(prop3->isValid(&data));

    IdleService s;
    s.addInput(&data);

    s.addConstraint(&data,prop1);
    s.addConstraint(&data,prop2);
    s.addConstraint(&data,prop3);

    ASSERT_TRUE(s.checkInput());

    gmds::DataScalar<double> data2(3.5);

    gmds::Property* prop4 = builder.build<double>(gmds::PropertyScalarBuilder::NEGATIVE);


    ASSERT_FALSE(prop4->isValid(&data2));

    s.addInput(&data2);
    s.addConstraint(&data2,prop4);
    ASSERT_FALSE(s.checkInput());

    EXPECT_ANY_THROW(
            s.executeAfterChecking();
    );
}

TEST(DataPropertyTestSuite, IncompatibleDataPropTest){

    gmds::DataScalar<int> d(2);
    gmds::PropertyScalarBuilder builder;
    gmds::Property* p = builder.build<double>(gmds::PropertyScalarBuilder::POSITIVE_STRICTLY);

    EXPECT_ANY_THROW(
            p->isValid(&d);
    );
}
#endif //GMDS_SCALAR_DATA_PROP_TESTSUITE_H
