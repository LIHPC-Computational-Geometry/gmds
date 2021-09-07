#ifndef GMDS_PARAM_TESTSUITE_H
#define GMDS_PARAM_TESTSUITE_H

#include "gtest/gtest.h"

#include <unit_test_config.h>
#include "gmds/utils/Parameters.h"
using namespace gmds;

TEST(ParamTestSuite, param_1){
    Parameters p;
    p.add("section_A","e", Parameters::STRING_P);
    p.add("section_A","b",  Parameters::BOOL_P);
    p.add("section_A","str",  Parameters::STRING_P);
    p.add("section_B","solver", Parameters::STRING_P);
    p.add("section_B","nb_steps", Parameters::INT_P);
    p.add("section_B","epsilon", Parameters::DOUBLE_P);

    std::string dir(TEST_SAMPLES_DIR);
    std::string init_file = dir+"/param_sample.ini";


    p.parseIni(init_file);

    std::string val_string;
    int         val_int;
    double      val_double;
    bool        val_bool;

    ASSERT_TRUE(p.get("section_A","e", val_string));
    ASSERT_EQ(val_string,"HEX_CAV");
    ASSERT_TRUE(p.get("section_A","b", val_bool));
    ASSERT_EQ(val_bool,false);
    ASSERT_TRUE(p.get("section_A","str", val_string));
    ASSERT_EQ(val_string,"mesh/input/toto.mesh");
    ASSERT_TRUE(p.get("section_B","solver", val_string));
    ASSERT_EQ(val_string,"EIGEN");
    ASSERT_TRUE(p.get("section_B","nb_steps", val_int));
    ASSERT_EQ(val_int,100);
    ASSERT_TRUE(p.get("section_B","epsilon", val_double));
    ASSERT_EQ(val_double,1e-4);


}
#endif //GMDS_UTILS_TESTSUITE_H
