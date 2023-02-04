/*----------------------------------------------------------------------------*/
/** \file    RandomGenerator.cpp
 *  \author  legoff
 *  \date    10/21/2015
 */
/*----------------------------------------------------------------------------*/
#include <gmds/utils/RandomGenerator.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <fstream>
#ifdef _WIN32
#include <ctime>
#endif

/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
RandomGenerator::RandomGenerator()
{

}
/*----------------------------------------------------------------------------*/
RandomGenerator::~RandomGenerator()
{

}
/*----------------------------------------------------------------------------*/
void
RandomGenerator::init()
{
	unsigned int seed;

	std::string input_file_name;
	input_file_name.append("/dev/random");
	std::ifstream input_file_stream(input_file_name.c_str(),std::ios::in);
	input_file_stream >> seed;
	input_file_stream.close();

	//seed = 2;
	srand(time(nullptr) + seed);
}
/*----------------------------------------------------------------------------*/
double
RandomGenerator::value()
{
  return ((double)rand() / (double) RAND_MAX);
}
/*----------------------------------------------------------------------------*/
