/*----------------------------------------------------------------------------*/
#include "gmds/blocking/Blocking.h"
/*----------------------------------------------------------------------------*/
#include <fstream>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/Exception.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace blocking {
/*----------------------------------------------------------------------------*/
Blocking::Blocking() {}
/*----------------------------------------------------------------------------*/
Blocking::~Blocking() {}
/*----------------------------------------------------------------------------*/
Blocking::STATUS
Blocking::execute()
{
	return Blocking::SUCCESS;
}
/*----------------------------------------------------------------------------*/
void Blocking::createGrid()
{
	lcc_.make_hexahedron(Point(0,0,0), Point(5,0,0),
								Point(5,5,0), Point(0,5,0),
								Point(0,5,4), Point(0,0,4),
								Point(5,0,4), Point(5,5,4));
}
/*----------------------------------------------------------------------------*/
void Blocking::writeMokaFile(std::string AFileName) const
{
	std::ofstream stream(AFileName, std::ios::out);
	if (!stream.is_open()){
		std::string s ="Impossible to create a Moka File (ASCII format): "+AFileName;
		throw gmds::GMDSException(s);
	}

	stream << "Moka file [ascii]" << std::endl;;

	// TODO investigate the magix number
	const int magic_number = 128;
	stream << magic_number <<" 7 0 0 0 0 0 0" << std::endl;

	for(auto it = lcc_.darts().begin(); it != lcc_.darts().end(); ++it) {

		stream << lcc_.darts().index(lcc_.get_alpha<0>(it)) << " ";
		stream << lcc_.darts().index(lcc_.get_alpha<1>(it)) << " ";
		stream << lcc_.darts().index(lcc_.get_alpha<2>(it)) << " ";
		stream << lcc_.darts().index(lcc_.get_alpha<3>(it)) << " ";

		stream << magic_number << " 0 0 0 ";

		// check whether the dart is associated to a point
		LCC_3::Vertex_attribute_const_handle v = lcc_.vertex_attribute(it);
		if (v->dart() == it) {
			stream << "1 " << lcc_.point_of_vertex_attribute(v) << std::endl;
		}
		else {
			stream << "0" << std::endl;
		}
	}
	stream.close();
}
/*----------------------------------------------------------------------------*/
}  // namespace blocking
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/