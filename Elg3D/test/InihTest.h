/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <gmds/utils/Parameters.h>

#include <iniparser.h>
/*----------------------------------------------------------------------------*/
class InihTest : public ::testing::Test
{
protected:
    InihTest()
    {
        ;
    }
    virtual ~InihTest()
    {
        ;
    }

    static void
    SetUpTestCase()
    {
        // Kokkos::Serial::initialize();
        // Kokkos::Threads::initialize();
        Kokkos::InitArguments kargs;
        kargs.num_threads = 1;
//        int num_threads = 3;
//        int use_numa = 1;
//        int use_core = 1;
//        Kokkos::OpenMP::initialize(num_threads, use_numa, use_core);
        Kokkos::initialize(kargs);
    }

    static void
    TearDownTestCase()
    {
        // Kokkos::Serial::finalize();
        // Kokkos::Threads::finalize();
//        Kokkos::OpenMP::finalize();
        Kokkos::finalize();
    }
};
/*----------------------------------------------------------------------------*/
TEST_F(InihTest, dummy)
{
    EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(InihTest, some_test)
{
    gmds::Parameters p;
//    p.add("section_A","e", gmds::Parameters::STRING_P);
//    p.add("section_A","b",  gmds::Parameters::BOOL_P);
//    p.add("section_A","str",  gmds::Parameters::STRING_P);
//    p.add("section_B","solver", gmds::Parameters::STRING_P);
//    p.add("section_B","nb_steps", gmds::Parameters::INT_P);
//    p.add("section_B","epsilon", gmds::Parameters::DOUBLE_P);
//    p.add("inputs","type", gmds::Parameters::STRING_P);
//    p.add("pixels","nbsub", gmds::Parameters::INT_P);

    p.add("inputs","type", gmds::Parameters::STRING_P);
    p.add("pixels","nbsub", gmds::Parameters::INT_P);
    p.add("mat_to_fuse","list", gmds::Parameters::STRING_P);


    p.parseIni("Elg3D/test/Samples/parameter_file.ini");



    std::string val_string;
    int         val_int;
    double      val_double;
    bool        val_bool;

//    ASSERT_TRUE(p.get("section_A","e", val_string));
//    ASSERT_EQ(val_string,"HEX_CAV");
//    ASSERT_TRUE(p.get("section_A","b", val_bool));
//    ASSERT_EQ(val_bool,false);
//    ASSERT_TRUE(p.get("section_A","str", val_string));
//    ASSERT_EQ(val_string,"mesh/input/toto.mesh");
//    ASSERT_TRUE(p.get("section_B","solver", val_string));
//    ASSERT_EQ(val_string,"EIGEN");
//    ASSERT_TRUE(p.get("section_B","nb_steps", val_int));
//    ASSERT_EQ(val_int,100);
//    ASSERT_TRUE(p.get("section_B","epsilon", val_double));
//    ASSERT_EQ(val_double,1e-4);


    ASSERT_TRUE(p.get("inputs","type", val_string));
    EXPECT_EQ(val_string, "internal");
    std::cout<<"val_string "<<val_string<<std::endl;
    ASSERT_TRUE(p.get("pixels","nbsub", val_int));
    EXPECT_EQ(val_int, 10);
//    ASSERT_TRUE(p.get("mat_to_fuse","list", val_string));
//    std::cout<<"val_string "<<val_string<<std::endl;


    dictionary  *   ini ;
    const char  *   s ;
    ini = iniparser_load("Elg3D/test/Samples/parameter_file.ini");
    if(ini == NULL) {
        std::cout<<"could not open poyop"<<std::endl;
    }

    s = iniparser_getstring(ini, "inputs:type", NULL);
    EXPECT_EQ(std::string(s), "internal");

    int value = iniparser_getint(ini, "pixels:nbsub", -1);
    EXPECT_EQ(value, 10);
    s = iniparser_getstring(ini, "pixels:nbsub", NULL);
    EXPECT_EQ(std::string(s), "10");
//    s = iniparser_getstring(ini, "mat_to_fuse:list", NULL);
//    EXPECT_EQ(std::string(s), "10");
//    std::cout<<"s "<<s<<std::endl;
}

/*----------------------------------------------------------------------------*/