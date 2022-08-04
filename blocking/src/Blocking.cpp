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
void Blocking::createGrid2d()
{
	this->createGrid2d(gmds::math::Point(0,0,0), gmds::math::Point(1,1,1), 3,3);
}
/*----------------------------------------------------------------------------*/
void Blocking::createGrid3d()
{
	this->createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(1,1,1), 3,3,3);
}
/*----------------------------------------------------------------------------*/
void Blocking::createGrid2d(gmds::math::Point APmin, gmds::math::Point APmax, int ANx, int ANy)
{
	Dart_handle* dhs = new Dart_handle[ANx*ANy];

	for(int i=0; i<ANx; i++) {
		for(int j=0; j<ANy; j++) {

			double x0 = APmin.X() + ((double) i / (double) ANx) * (APmax.X() + (-1. * APmin.X()));
			double y0 = APmin.Y() + ((double) j / (double) ANy) * (APmax.Y() + (-1. * APmin.Y()));
			double x1 = APmin.X() + ((double) (i+1) / (double) ANx) * (APmax.X() + (-1. * APmin.X()));
			double y1 = APmin.Y() + ((double) (j+1) / (double) ANy) * (APmax.Y() + (-1. * APmin.Y()));

			Dart_handle dh =
				lcc_.make_quadrangle(Point(x0, y0, 0.),
											Point(x1, y0, 0.),
											Point(x1, y1, 0.),
											Point(x0, y1, 0.));

			dhs[i*ANy+j] = dh;
		}
	}


	for(int i=0; i<ANx; i++) {
		for (int j=0; j<ANy; j++) {
			Dart_handle dh = dhs[i*ANy+j];

			if(i != ANx-1) {
				Dart_handle dhi = dhs[(i+1)*ANy+j];
				lcc_.sew<2>(lcc_.alpha(dh, 0, 1), lcc_.alpha(dhi, 1));
			}
			if(j != ANy-1) {
				Dart_handle dhj = dhs[i*ANy+j+1];
				lcc_.sew<2>(lcc_.alpha(dh, 1, 0, 1), dhj);
			}
		}
	}

	delete[] dhs;

	if (!lcc_.is_valid()) {
		std::string s ="Blocking::createGrid lcc not valid";
		throw gmds::GMDSException(s);
	}
}
/*----------------------------------------------------------------------------*/
void Blocking::createGrid3d(gmds::math::Point APmin, gmds::math::Point APmax, int ANx, int ANy, int ANz)
{
	Dart_handle* dhs = new Dart_handle[ANx*ANy*ANz];

	for(int i=0; i<ANx; i++) {
		for(int j=0; j<ANy; j++) {
			for(int k=0; k<ANz; k++) {

				double x0 = APmin.X() + ((double) i / (double) ANx) * (APmax.X() + (-1. * APmin.X()));
				double y0 = APmin.Y() + ((double) j / (double) ANy) * (APmax.Y() + (-1. * APmin.Y()));
				double z0 = APmin.Z() + ((double) k / (double) ANz) * (APmax.Z() + (-1. * APmin.Z()));
				double x1 = APmin.X() + ((double) (i+1) / (double) ANx) * (APmax.X() + (-1. * APmin.X()));
				double y1 = APmin.Y() + ((double) (j+1) / (double) ANy) * (APmax.Y() + (-1. * APmin.Y()));
				double z1 = APmin.Z() + ((double) (k+1) / (double) ANz) * (APmax.Z() + (-1. * APmin.Z()));

				Dart_handle dh =
				   lcc_.make_hexahedron(Point(x0, y0, z0),
												Point(x1, y0, z0),
				                        Point(x1, y1, z0),
												Point(x0, y1, z0),
												Point(x0, y1, z1),
												Point(x0, y0, z1),
												Point(x1, y0, z1),
												Point(x1, y1, z1));

				dhs[i*(ANy*ANz)+j*ANz+k] = dh;
			}
		}
	}

	for(int i=0; i<ANx; i++) {
		for (int j=0; j<ANy; j++) {
			for (int k=0; k<ANz; k++) {
				Dart_handle dh = dhs[i*(ANy*ANz)+j*ANz+k];

				if(i != ANx-1) {
					Dart_handle dhi = dhs[(i+1)*(ANy*ANz)+j*ANz+k];
					lcc_.sew<3>(lcc_.alpha(dh, 1, 0, 1, 2), lcc_.alpha(dhi, 2));
				}
				if(j != ANy-1) {
					Dart_handle dhj = dhs[i*(ANy*ANz)+(j+1)*ANz+k];
					lcc_.sew<3>(lcc_.alpha(dh, 2, 1, 0, 1, 2), dhj);
				}
				if(k != ANz-1) {
					Dart_handle dhk = dhs[i*(ANy*ANz)+j*ANz+(k+1)];
					lcc_.sew<3>(lcc_.alpha(dh, 0, 1, 2), lcc_.alpha(dhk, 1, 2));
				}
			}
		}
	}

	delete[] dhs;

	if (!lcc_.is_valid()) {
		std::string s ="Blocking::createGrid lcc not valid";
		throw gmds::GMDSException(s);
	}
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

	// TODO investigate the magic_number
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