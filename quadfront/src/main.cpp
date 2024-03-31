#include <iostream>
#include "gmds/quadfront/Quadfront.h"

int main() {
	std::cout<<"hello world"<<std::endl;

	gmds::quadfront::Quadfront qf;
	gmds::quadfront::Quadfront::STATUS st = qf.execute();

}