#include <iostream>
#include "gmds/milp/milp.h"

int main() {
	std::cout<<"hello world"<<std::endl;

	gmds::milp::milp qf;
	gmds::milp::milp::STATUS st = qf.execute();

}