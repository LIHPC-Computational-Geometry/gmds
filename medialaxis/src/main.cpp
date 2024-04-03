#include <iostream>
#include "gmds/medialaxis/Medialaxis.h"

int main() {
	std::cout<<"hello world"<<std::endl;

	gmds::medialaxis::Medialaxis qf;
	gmds::medialaxis::Medialaxis::STATUS st = qf.execute();

}