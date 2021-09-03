/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 7/27/19.
//
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Point.h>
#include <gmds/math/Matrix.h>
/*----------------------------------------------------------------------------*/
#include <gmds/sofiane/SimplexMesh.h>
#include <gmds/sofiane/SimplicesCell.h>
#include <gmds/sofiane/SimplicesNode.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <bitset>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::hybrid;
using namespace gmds::math;
using namespace hybrid;
using namespace hybrid::simplicesCell;
using namespace hybrid::simplicesNode;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    std::cout << "==== SOFIANE ====" << std::endl;

    std::cout << "limite min : " << std::numeric_limits<int64_t>::min() << std::endl;
    std::cout << "limite max : " << std::numeric_limits<int32_t>::max() << std::endl;

}
