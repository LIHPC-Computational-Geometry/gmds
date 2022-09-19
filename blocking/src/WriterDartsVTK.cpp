/*----------------------------------------------------------------------------*/
#include "gmds/blocking/WriterDartsVTK.h"
/*----------------------------------------------------------------------------*/
#include <fstream>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/utils/Exception.h>

#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace blocking {
/*----------------------------------------------------------------------------*/
WriterDartsVTK::WriterDartsVTK() {}
/*----------------------------------------------------------------------------*/
WriterDartsVTK::~WriterDartsVTK() {}
/*----------------------------------------------------------------------------*/
WriterDartsVTK::STATUS
WriterDartsVTK::execute(std::string AFilename, int AToMark)
{
	LCC_3::size_type mark = bl()->lcc()->get_new_mark();

	// In some (all?) cases when in 2d a unique 3d cell exists and all darts self reference themselves by alpha3
	// check whether we have "real" 3d cells that we need to effectively write
	bool is2d = true;
	if(bl()->lcc()->one_dart_per_cell<3>().size() == 1) {
		int index = 0;
		for(auto it = bl()->lcc()->darts().begin(); it != bl()->lcc()->darts().end(); ++it) {
			if(index != bl()->lcc()->darts().index(bl()->lcc()->get_alpha<3>(it))) {
				is2d = false;
				break;
			}
			index++;
		}
	} else {
		is2d = false;
	}

	// edges
	std::map<Dart_handle, std::pair<gmds::math::Point, gmds::math::Point> > dart_pos;

	for (auto it = bl()->lcc()->one_dart_per_cell<1>().begin(); it != bl()->lcc()->one_dart_per_cell<1>().end(); it++) {

		Dart_handle d = it;

		LCC_3::Point barycenter = lcc()->barycenter<1>(d);
		gmds::math::Point barycenter_gmds(barycenter.x(), barycenter.y(), barycenter.z());

		for(LCC_3::Dart_of_orbit_range<0,2,3>::iterator
				 itbis(lcc()->darts_of_orbit<0,2,3>(d).begin()),
				 itend(lcc()->darts_of_orbit<0,2,3>(d).end()); itbis!=itend; ++itbis) {

			LCC_3::Dart_handle d0 = itbis;

			if(!bl()->lcc()->is_marked(d0, mark)) {
				LCC_3::Dart_handle d1 = bl()->lcc()->alpha(d0, 0);
				bl()->lcc()->mark(d0, mark);
				bl()->lcc()->mark(d1, mark);

				LCC_3::Point pt0 = bl()->lcc()->vertex_attribute(d0)->point();
				LCC_3::Point pt1 = bl()->lcc()->vertex_attribute(d1)->point();
				gmds::math::Point pt0_gmds(pt0.x(), pt0.y(), pt0.z());
				gmds::math::Point pt1_gmds(pt1.x(), pt1.y(), pt1.z());

				dart_pos.emplace(d0, std::pair<gmds::math::Point, gmds::math::Point> (pt0_gmds, pt0_gmds + ratio0_ * (barycenter_gmds - pt0_gmds)));
				dart_pos.emplace(d1, std::pair<gmds::math::Point, gmds::math::Point> (pt1_gmds, pt1_gmds + ratio0_ * (barycenter_gmds - pt1_gmds)));

				gmds::math::Point e0_0 = (1.- ratio1_) * (barycenter_gmds - dart_pos[d0].first) + dart_pos[d0].first;
				gmds::math::Point e0_1 = (1.- ratio1_) * (barycenter_gmds - dart_pos[d0].second) + dart_pos[d0].second;
				gmds::math::Point e1_0 = (1.- ratio1_) * (barycenter_gmds - dart_pos[d1].first) + dart_pos[d1].first;
				gmds::math::Point e1_1 = (1.- ratio1_) * (barycenter_gmds - dart_pos[d1].second) + dart_pos[d1].second;
				dart_pos[d0].first = e0_0;
				dart_pos[d0].second = e0_1;
				dart_pos[d1].first = e1_0;
				dart_pos[d1].second = e1_1;
			}
		}
	}

	for (auto it = bl()->lcc()->one_dart_per_cell<2>().begin(); it != bl()->lcc()->one_dart_per_cell<2>().end(); it++) {
		Dart_handle d = it;

		if(dart_pos.find(d) == dart_pos.end()) {
			std::cout<<"a dart from a 2-cell was not previously treated with the 1-cells"<<std::endl;
			return WriterDartsVTK::FAIL;
		}

		LCC_3::Point barycenter = lcc()->barycenter<2>(d);
		gmds::math::Point barycenter_gmds(barycenter.x(), barycenter.y(), barycenter.z());

		for(LCC_3::Dart_of_orbit_range<0,1,3>::iterator
				 itbis(lcc()->darts_of_orbit<0,1,3>(d).begin()),
				 itend(lcc()->darts_of_orbit<0,1,3>(d).end()); itbis!=itend; ++itbis) {

			LCC_3::Dart_handle dbis = itbis;

			if(dart_pos.find(dbis) == dart_pos.end()) {
				std::cout<<"a dart from a 2-cell was not previously treated with the 1-cells"<<std::endl;
				return WriterDartsVTK::FAIL;
			}

			gmds::math::Point e0_0 = (1.- ratio2_) * (barycenter_gmds - dart_pos[dbis].first) + dart_pos[dbis].first;
			gmds::math::Point e0_1 = (1.- ratio2_) * (barycenter_gmds - dart_pos[dbis].second) + dart_pos[dbis].second;
			dart_pos[dbis].first = e0_0;
			dart_pos[dbis].second = e0_1;
		}
	}

	if(!is2d) {

		for (auto it = bl()->lcc()->one_dart_per_cell<3>().begin(); it != bl()->lcc()->one_dart_per_cell<3>().end(); it++) {

			Dart_handle d = it;

			if(dart_pos.find(d) == dart_pos.end()) {
				std::cout<<"a dart from a 3-cell was not previously treated with the 1-cells"<<std::endl;
				return WriterDartsVTK::FAIL;
			}

			LCC_3::Point barycenter = lcc()->barycenter<3>(d);
			gmds::math::Point barycenter_gmds(barycenter.x(), barycenter.y(), barycenter.z());

			for(LCC_3::Dart_of_orbit_range<0,1,2>::iterator
					 itbis(lcc()->darts_of_orbit<0,1,2>(d).begin()),
					 itend(lcc()->darts_of_orbit<0,1,2>(d).end()); itbis!=itend; ++itbis) {

				LCC_3::Dart_handle dbis = itbis;

				if(dart_pos.find(dbis) == dart_pos.end()) {
					std::cout<<"a dart from a 3-cell was not previously treated with the 1-cells"<<std::endl;
					return WriterDartsVTK::FAIL;
				}

				gmds::math::Point e0_0 = (1.- ratio3_) * (barycenter_gmds - dart_pos[dbis].first) + dart_pos[dbis].first;
				gmds::math::Point e0_1 = (1.- ratio3_) * (barycenter_gmds - dart_pos[dbis].second) + dart_pos[dbis].second;
				dart_pos[dbis].first = e0_0;
				dart_pos[dbis].second = e0_1;
			}
		}
	}

	bl()->lcc()->free_mark(mark);

	// fill the gmds mesh
	// We create triangles because paraview seems to not be able to display edge cell types
	gmds::MeshModel model(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::F|gmds::F2N));
	gmds::Mesh m(model);

	gmds::Variable<int>* var_node_type = m.newVariable<int, gmds::GMDS_NODE>("node_type");
	gmds::Variable<int>* var_edge_type = m.newVariable<int, gmds::GMDS_FACE>("edge_type");

	gmds::Variable<int>* var_edge_marked = m.newVariable<int, gmds::GMDS_FACE>("edge_marked");

	for(auto d: dart_pos) {
		gmds::Node n0 = m.newNode(d.second.first);
		gmds::Node n1 = m.newNode(d.second.second);
		gmds::Face f = m.newTriangle(n0, n1, n1);
		if(AToMark != -1) {
			if(bl()->lcc()->is_marked(d.first, AToMark)) {
				(*var_edge_marked)[f.id()] = 1;
			}
		}
		(*var_node_type)[n0.id()] = 1;
		(*var_node_type)[n1.id()] = 0;
		(*var_edge_type)[f.id()] = 4;
	}

	// now for the alpha links
	// avoid creating the kinks twice by marking the darts
	std::set<Dart_handle> set_alpha0;
	std::set<Dart_handle> set_alpha1;
	std::set<Dart_handle> set_alpha2;
	std::set<Dart_handle> set_alpha3;

	for(auto d: dart_pos) {
		Dart_handle dh = d.first;
		Dart_handle d0 = bl()->lcc()->alpha(dh, 0);
		Dart_handle d1 = bl()->lcc()->alpha(dh, 1);
		Dart_handle d2 = bl()->lcc()->alpha(dh, 2);
		Dart_handle d3 = bl()->lcc()->alpha(dh, 3);

		if(d0 != dh) {
			if(dart_pos.find(d0) != dart_pos.end()) {
				if(set_alpha0.find(dh) == set_alpha0.end()) {
					gmds::Node n0 = m.newNode(d.second.second);
					gmds::Node n1 = m.newNode(dart_pos[d0].second);
					gmds::Face f = m.newTriangle(n0, n1, n1);
					(*var_node_type)[n0.id()] = 0;
					(*var_node_type)[n1.id()] = 0;
					(*var_edge_type)[f.id()] = 0;

					set_alpha0.insert(dh);
					set_alpha0.insert(d0);
				}
			}
		}
		if(d1 != dh) {
			if(dart_pos.find(d1) != dart_pos.end()) {
				if(set_alpha1.find(dh) == set_alpha1.end()) {
					gmds::Node n0 = m.newNode(d.second.first + display_alpha1_ * (d.second.second - d.second.first));
					gmds::Node n1 = m.newNode(dart_pos[d1].first + display_alpha1_ * (dart_pos[d1].second - dart_pos[d1].first));
					gmds::Face f = m.newTriangle(n0, n1, n1);
					(*var_node_type)[n0.id()] = 0;
					(*var_node_type)[n1.id()] = 0;
					(*var_edge_type)[f.id()] = 1;

					set_alpha1.insert(dh);
					set_alpha1.insert(d1);
				}

			}
		}
		if(d2 != dh) {
			if(dart_pos.find(d2) != dart_pos.end()) {
				if(set_alpha2.find(dh) == set_alpha2.end()) {
					gmds::Node n0 = m.newNode(d.second.first + display_alpha2_ * (d.second.second - d.second.first));
					gmds::Node n1 = m.newNode(dart_pos[d2].first + display_alpha2_ * (dart_pos[d2].second - dart_pos[d2].first));
					gmds::Face f = m.newTriangle(n0, n1, n1);
					(*var_node_type)[n0.id()] = 0;
					(*var_node_type)[n1.id()] = 0;
					(*var_edge_type)[f.id()] = 2;

					set_alpha2.insert(dh);
					set_alpha2.insert(d2);
				}

			}
		}
		if(d3 != dh) {
			if(dart_pos.find(d3) != dart_pos.end()) {
				if(set_alpha3.find(dh) == set_alpha3.end()) {
					gmds::Node n0 = m.newNode(d.second.first + display_alpha3_ * (d.second.second - d.second.first));
					gmds::Node n1 = m.newNode(dart_pos[d3].first + display_alpha3_ * (dart_pos[d3].second - dart_pos[d3].first));
					gmds::Face f = m.newTriangle(n0, n1, n1);
					(*var_node_type)[n0.id()] = 0;
					(*var_node_type)[n1.id()] = 0;
					(*var_edge_type)[f.id()] = 3;

					set_alpha3.insert(dh);
					set_alpha3.insert(d3);
				}

			}
		}
	}

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N | gmds::F);
	vtkWriter.setDataOptions(gmds::N | gmds::F);
	vtkWriter.write(AFilename);

	return WriterDartsVTK::SUCCESS;
}
/*----------------------------------------------------------------------------*/
}  // namespace blocking
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/