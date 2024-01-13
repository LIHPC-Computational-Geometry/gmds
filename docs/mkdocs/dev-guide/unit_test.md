
# Test directory structure
Every gmds component has a *tst* directory that contains test suites to perform. All our tests are based on Google test. The classical structure of a test directory is :
```Shell
AFirstTestSuite.h
ASeconTestSuite.h
CMakeLists.txt
main.cpp
```
#Process to test a new feature
Each time you need to test a new class and methods you want to develop, the process is the following one.
1. Create a class *C* by declaring it in **inc/C.h** and provides all the method definitions in **src/C.cpp**. At this stage, method definition are either empty methods or methods that throw an exception.
2. Create in *tst* a new class for testing *C* behaviour. By convention, you will note it *CTestSuite.h*. You will then addd test method using google test. As an example below, we test some getters ans setters of the **math::Point** class. Only one test method is given but the test class (here **PointClass**) can contain as much as test method you want. **TEST** and **ASSERT_NEAR** are macros provided by the Google Test framework.

```cpp
/*----------------------------------------------------------------------------*/
#ifndef GMDS_POINT_TESTSUITE_H
#define GMDS_POINT_TESTSUITE_H
/*----------------------------------------------------------------------------*/
#include "gtest/gtest.h"
#include <gmds/math/Point.h>
/*----------------------------------------------------------------------------*/
TEST(PointClass, Setter)
{
    gmds::math::Point p(1,2,3);
    gmds::TCoord z = p.Z();
    ASSERT_NEAR(p.X(), 1, 1e-6);
    ASSERT_NEAR(p.Y(), 2, 1e-6);
    ASSERT_NEAR(z, 3, 1e-6);
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_POINT_TESTSUITE_H
/*----------------------------------------------------------------------------*/
```
3. Once created the file must be added in the CMakeLists.txt to be used in the compilation stage. Get into this file and only modify the list of header/source files used to build the executable test.

```cmake
add_executable(GMDS_MATH_TEST
        ChartTestSuite.h
        CrossTestSuite.h
        Cross2DTestSuite.h
        MathTestSuite.h
        PointTestSuite.h
        QuaternionTestSuite.h
        OrientationTestSuite.h
        ... CTestSuite.h ...
        main_test.cpp)
```
If the directory *tst* is empty and you are the first one to create a unit test suite for a component. Just copy and paste a CMakeLists.txt file of another component test directory and change the executable name and likely the list of dependencies (in the CMake macro *target_link_libraries*).

4. It remains a last step, which is to modify (or create) the **main_test.cpp** file, which looks like
```cpp
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch

#include <ChartTestSuite.h>
#include <Cross2DTestSuite.h>
#include <CrossTestSuite.h>
#include <OrientationTestSuite.h>
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/
```
You need to include the new file you just created in the list of files included  at the beginning of this file;

