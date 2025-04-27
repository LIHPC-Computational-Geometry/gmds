#ifndef GMDS_UTILS_TESTSUITE_H
#define GMDS_UTILS_TESTSUITE_H

#include "gtest/gtest.h"
#include "gmds/utils/Exception.h"
#include "gmds/utils/Log.h"
#include "gmds/utils/LogStream.h"
#include "gmds/utils/CommonTypes.h"
#include "gmds/utils/FileSystem.h"
#include "gmds/utils/BitVector.h"


TEST(UtilsTestSuite, FileSystemTest){
    EXPECT_ANY_THROW(
            throw gmds::GMDSException("test");
    );

}
TEST(UtilsTestSuite, ExceptionWhat){

    gmds::GMDSException e("test");
    ASSERT_STREQ(e.what(),"test");

}
TEST(UtilsTestSuite, LogTest){


    std::stringstream string_str;
    gmds::LogStream log_stream(&string_str);

    ASSERT_EQ(log_stream.level(),gmds::LOG_INFO);

    gmds::Log::mng().addStream(log_stream);
    ASSERT_EQ(gmds::Log::mng().reportingLevel(),gmds::LOG_INFO);
    gmds::Log::reportingLevel()=gmds::LOG_ERROR;
    ASSERT_EQ(gmds::Log::mng().reportingLevel(),gmds::LOG_ERROR);

}

TEST(UtilsTestSuite, common){
    std::vector<gmds::TCellID> v1, v2, v3;
    v1.push_back(1);
    v1.push_back(2);
    v1.push_back(3);
    v1.push_back(4);
    v2.push_back(3);
    v2.push_back(4);
    v2.push_back(2);
    v2.push_back(6);
    v3 = gmds::getCommonBut(v1,v2,4);

    ASSERT_TRUE(std::find(v3.begin(),v3.end(),2)!=v3.end());
    ASSERT_TRUE(std::find(v3.begin(),v3.end(),3)!=v3.end());
    ASSERT_TRUE(std::find(v3.begin(),v3.end(),4)==v3.end());


}
TEST(UtilsTestSuite, keepFilter){
    std::vector<gmds::TCellID> v, v2, v3;
    v.push_back(1);
    v.push_back(2);
    v.push_back(3);
    v.push_back(4);
    v.push_back(3);
    v.push_back(4);
    v.push_back(3);
    v.push_back(6);

    v2 = gmds::keepFilter(v,2);
    ASSERT_EQ(v2.size(),2);
    ASSERT_TRUE(std::find(v2.begin(),v2.end(),3)!=v2.end());
    ASSERT_TRUE(std::find(v2.begin(),v2.end(),4)!=v2.end());
    v3 = gmds::keepFilter(v,3);
    ASSERT_EQ(v3.size(),1);
    ASSERT_TRUE(std::find(v3.begin(),v3.end(),3)!=v3.end());


}
TEST(UtilsTestSuite, bitVector){
    gmds::BitVector v;

    ASSERT_EQ(v.size(),0);
    ASSERT_EQ(v.capacity(),gmds::GChunkSize);
    ASSERT_EQ(v.top(),0);

    ASSERT_EQ(v.selectNewBit(),0);
    ASSERT_EQ(v.selectNewBit(),1);
    ASSERT_EQ(v.selectNewBit(),2);
    ASSERT_EQ(v.selectNewBit(),3);
    ASSERT_EQ(v.selectNewBit(),4);

    int nb_vals=0;
    for(gmds::TInt val:v){
        nb_vals++;
    }
    ASSERT_EQ(nb_vals,5);

    v.unselect(3);
    v.unselect(1);
    ASSERT_EQ(v.size(),3);
    nb_vals=0;
    for(auto & val:v){
        nb_vals++;
    }
    ASSERT_EQ(nb_vals,3);

    std::stringstream str_in;
    v.serialize(str_in);

    gmds::BitVector v2;
    v2.unserialize(str_in);
    std::stringstream str_out;
    v2.serialize(str_out);

    ASSERT_EQ(str_in.str(),str_out.str());

}
#endif //GMDS_UTILS_TESTSUITE_H
