/*----------------------------------------------------------------------------*/
#include "../inc/gmds/utils/Log.h"
#include "../inc/gmds/utils/Exception.h"
#include "../inc/gmds/utils/Assert.h"
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
namespace {
    AssertionMode m_assert_mode = ASSERT_MODE_THROW;
}

/*----------------------------------------------------------------------------*/
void setAssertMode(AssertionMode mode) {
    m_assert_mode = mode;
}

/*----------------------------------------------------------------------------*/
AssertionMode getAssertMode() {
    return m_assert_mode;
}
/*----------------------------------------------------------------------------*/
void assertFailed(const std::string& ACondition,
                  const std::string& AFileName,
                  int ALineNum)
{
    std::ostringstream stream;
    stream << "Assertion failed - " << ACondition <<std::endl;
    stream << "\t File: " << AFileName <<std::endl;
    stream << "\t Line: " << ALineNum;
    
    if(m_assert_mode == ASSERT_MODE_THROW) {
        throw GMDSException(stream.str());
    } else {
        Log::mng()<< stream.str() << "\n";
        exit(0);
    }
}
/*----------------------------------------------------------------------------*/
void assertRangeFailed(const double& AVal,
                       const double& AMin,
                       const double& AMax,
                       const std::string& AFileName,
                       int ALineNum)
{
    std::ostringstream os;
    std::ostringstream stream;
    stream <<"Range assertion failed - "<< AVal <<" in ["
    <<AMin<<", "<<AMax<<"]"<<std::endl;
    stream <<"\t File: "<< AFileName <<std::endl;
    stream <<"\t Line: "<< ALineNum;
   
    if(m_assert_mode == ASSERT_MODE_THROW) {
        throw GMDSException(stream.str());
    } else {
        Log::mng()<< stream.str() <<"\n";
        exit(0);
    }
}
/*----------------------------------------------------------------------------*/

