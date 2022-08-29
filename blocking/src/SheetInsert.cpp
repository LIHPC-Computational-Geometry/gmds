/*----------------------------------------------------------------------------*/
#include "gmds/blocking/SheetInsert.h"
/*----------------------------------------------------------------------------*/
#include <array>
/*----------------------------------------------------------------------------*/
//#include <gmds/ig/Mesh.h>
//#include <gmds/utils/Exception.h>
//
//#include <gmds/io/IGMeshIOService.h>
//#include <gmds/io/VTKReader.h>
//#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace blocking {
/*----------------------------------------------------------------------------*/
SheetInsert::SheetInsert() {}
/*----------------------------------------------------------------------------*/
SheetInsert::~SheetInsert() {}
/*----------------------------------------------------------------------------*/
SheetInsert::STATUS
SheetInsert::execute()
{
	return SheetInsert::NOT_YET_IMPLEMENTED;
}
/*----------------------------------------------------------------------------*/
SheetInsert::STATUS
SheetInsert::pillow()
{
	//TODO get the set of 3-cells in another way
	//


	// check that the implementation can handle the data
	if(!lcc()->is_valid()) {
		std::string s ="SheetInsert::pillow can be applied on a valid gmap only";
		throw gmds::GMDSException(s);
	}
	// TODO check the validity of the shrink set


	// mark the darts to extrude
	LCC_3::size_type m = lcc()->get_new_mark();
	LCC_3::size_type m_new = lcc()->get_new_mark();

//	lcc.mark(dh1, mark);
	// TODO for now we mark the darts of the first cell
//	for (auto it = getlcc()->one_dart_per_cell<3>().begin(); it != getlcc()->one_dart_per_cell<3>().end(); it++) {
//		int nb = getlcc()->darts_of_orbit<0, 1, 2>(it).size();
//		std::cout<<"nbdarts "<<nb<<std::endl;
//	}
	auto d = lcc()->one_dart_per_cell<3>().begin();
	int nb = lcc()->darts_of_orbit<0, 1, 2>(d).size();
	std::cout<<"nbdarts "<<nb<<std::endl;

	LCC_3::Point center;

	for (LCC_3::Dart_of_cell_range<3>::iterator
			  it(lcc()->darts_of_cell<3>(d).begin()),
			  itend(lcc()->darts_of_cell<3>(d).end()); it!=itend; ++it) {
		lcc()->mark(it, m);
		center = lcc()->barycenter<3>(it);
	}

//	for (LCC_3::Dart_range::iterator it(getlcc()->darts().begin()),
//			  itend(getlcc()->darts().end()); it!=itend; ++it)
//	{
//		getlcc()->mark(it, m);
//	}

	std::cout<<"nbMarked "<< lcc()->number_of_marked_darts(m)<<std::endl;

	// first mark the darts on the boundary of the shrink set


	std::map<LCC_3::Dart_handle, LCC_3::Dart_handle> old_alpha3;

	for (LCC_3::Dart_range::iterator it(lcc()->darts().begin()),
			  itend(lcc()->darts().end()); it!=itend; ++it) {

		if(lcc()->is_marked(it, m)) {
			LCC_3::Dart_handle d0 = it;
			LCC_3::Dart_handle d3 = lcc()->alpha(d0,3);

			old_alpha3.insert(std::pair<LCC_3::Dart_handle, LCC_3::Dart_handle> (d0, d3));
		}
	}
	std::cout<<"old_alpha3.size() "<< old_alpha3.size()<<std::endl;

	// create the pattern for the marked darts
	std::map<LCC_3::Dart_handle, std::array<LCC_3::Dart_handle, 10> > old_pattern;

	for(auto d: old_alpha3) {
//		LCC_3::Vertex_attribute_const_handle v = getlcc()->vertex_attribute(d.first);
		LCC_3::Vertex_attribute_handle v = lcc()->vertex_attribute(d.first);

		// TODO determine position of new vertex and which darts are assigned to it
//		LCC_3::Vertex_attribute_handle vbis = lcc()->create_vertex_attribute(LCC_3::Point ());
		double posx = center.x() + 0.8 * (v->point().x() - center.x());
		double posy = center.y() + 0.8 * (v->point().y() - center.y());
		double posz = center.z() + 0.8 * (v->point().z() - center.z());
//		gmds::math::Point newPos =
//		LCC_3::Vertex_attribute_handle vbis = lcc()->create_vertex_attribute( center + 0.8 * (center + v->point()) );
//		LCC_3::Vertex_attribute_handle vbis = lcc()->create_vertex_attribute( LCC_3::Point (v->point().x(), v->point().y(), v->point().z()));
//		lcc()->point_of_vertex_attribute(v) = LCC_3::Point (posx, posy, posz);
		LCC_3::Vertex_attribute_handle vbis = lcc()->create_vertex_attribute( LCC_3::Point (posx, posy, posz));
		lcc()->set_vertex_attribute_of_dart(d.first, vbis);

//		LCC_3::Dart_handle ddd = lcc()->create_dart(v);
//		LCC_3::Dart_handle d0 = lcc()->create_dart(v->point());
//		LCC_3::Dart_handle d1 = lcc()->create_dart(v->point());
//		LCC_3::Dart_handle d2 = lcc()->create_dart(v->point());
//		LCC_3::Dart_handle d3 = lcc()->create_dart(v->point());
//		LCC_3::Dart_handle d4 = lcc()->create_dart(v->point());
//		LCC_3::Dart_handle d5 = lcc()->create_dart(v->point());
//		LCC_3::Dart_handle d6 = lcc()->create_dart(v->point());
//		LCC_3::Dart_handle d7 = lcc()->create_dart(v->point());
//		LCC_3::Dart_handle d8 = lcc()->create_dart(v->point());
//		LCC_3::Dart_handle d9 = lcc()->create_dart(v->point());
//		LCC_3::Dart_handle d0 = lcc()->create_dart(v);
//		LCC_3::Dart_handle d1 = lcc()->create_dart(v);
//		LCC_3::Dart_handle d2 = lcc()->create_dart(v);
//		LCC_3::Dart_handle d3 = lcc()->create_dart(v);
//		LCC_3::Dart_handle d4 = lcc()->create_dart(v);
//		LCC_3::Dart_handle d5 = lcc()->create_dart(v);
//		LCC_3::Dart_handle d6 = lcc()->create_dart(v);
//		LCC_3::Dart_handle d7 = lcc()->create_dart(v);
//		LCC_3::Dart_handle d8 = lcc()->create_dart(v);
//		LCC_3::Dart_handle d9 = lcc()->create_dart(v);

		LCC_3::Dart_handle d0 = lcc()->create_dart(vbis);
		LCC_3::Dart_handle d1 = lcc()->create_dart(v);
		LCC_3::Dart_handle d2 = lcc()->create_dart(vbis);
		LCC_3::Dart_handle d3 = lcc()->create_dart(vbis);
		LCC_3::Dart_handle d4 = lcc()->create_dart(v);
		LCC_3::Dart_handle d5 = lcc()->create_dart(v);
		LCC_3::Dart_handle d6 = lcc()->create_dart(vbis);
		LCC_3::Dart_handle d7 = lcc()->create_dart(vbis);
		LCC_3::Dart_handle d8 = lcc()->create_dart(v);
		LCC_3::Dart_handle d9 = lcc()->create_dart(v);

		lcc()->mark(d0, m_new);
		lcc()->mark(d1, m_new);
		lcc()->mark(d2, m_new);
		lcc()->mark(d3, m_new);
		lcc()->mark(d4, m_new);
		lcc()->mark(d5, m_new);
		lcc()->mark(d6, m_new);
		lcc()->mark(d7, m_new);
		lcc()->mark(d8, m_new);
		lcc()->mark(d9, m_new);

		lcc()->link_alpha<2>(d0, d2);
		lcc()->link_alpha<2>(d5, d1);

		lcc()->link_alpha<1>(d2, d3);
		lcc()->link_alpha<0>(d3, d4);
		lcc()->link_alpha<1>(d4, d5);
		lcc()->link_alpha<1>(d6, d7);
		lcc()->link_alpha<0>(d7, d8);
		lcc()->link_alpha<1>(d8, d9);

		lcc()->link_alpha<3>(d2, d6);
		lcc()->link_alpha<3>(d3, d7);
		lcc()->link_alpha<3>(d4, d8);
		lcc()->link_alpha<3>(d5, d9);

		old_pattern.insert(std::pair<LCC_3::Dart_handle, std::array<LCC_3::Dart_handle, 10> > (d.first, {d0,d1,d2,d3,d4,d5,d6,d7,d8,d9}));
	}

	// link the created darts between the patterns
	// first step with the always known links
	for(auto it: old_alpha3) {
		LCC_3::Dart_handle d = it.first;
		LCC_3::Dart_handle d0 = old_pattern[d][0];
		LCC_3::Dart_handle d1 = old_pattern[d][1];
		LCC_3::Dart_handle d2 = old_pattern[d][2];
		LCC_3::Dart_handle d3 = old_pattern[d][3];
		LCC_3::Dart_handle d4 = old_pattern[d][4];
		LCC_3::Dart_handle d5 = old_pattern[d][5];
		LCC_3::Dart_handle d6 = old_pattern[d][6];
		LCC_3::Dart_handle d9 = old_pattern[d][9];

		lcc()->link_alpha<0>(d0, old_pattern[lcc()->alpha(d,0)][0]);
		lcc()->link_alpha<1>(d0, old_pattern[lcc()->alpha(d,1)][0]);

		lcc()->link_alpha<0>(d1, old_pattern[lcc()->alpha(d,0)][1]);
		lcc()->link_alpha<1>(d1, old_pattern[lcc()->alpha(d,1)][1]);

		lcc()->link_alpha<0>(d2, old_pattern[lcc()->alpha(d,0)][2]);
		lcc()->link_alpha<0>(d5, old_pattern[lcc()->alpha(d,0)][5]);

		lcc()->link_alpha<0>(d6, old_pattern[lcc()->alpha(d,0)][6]);
		lcc()->link_alpha<0>(d9, old_pattern[lcc()->alpha(d,0)][9]);

		// alpha2 for d3 d4
		lcc()->link_alpha<2>(d3, old_pattern[lcc()->alpha(d,1)][3]);
		lcc()->link_alpha<2>(d4, old_pattern[lcc()->alpha(d,1)][4]);
	}

	// link the created darts between the patterns
	// second step with the links that depend on the orbits
	for(auto it: old_alpha3) {
		LCC_3::Dart_handle d = it.first;
		LCC_3::Dart_handle d6 = old_pattern[d][6];
		LCC_3::Dart_handle d7 = old_pattern[d][7];
		LCC_3::Dart_handle d8 = old_pattern[d][8];
		LCC_3::Dart_handle d9 = old_pattern[d][9];

		lcc()->darts_of_orbit<2,3>(d);
		std::cout<<"orbit.size() "<< lcc()->darts_of_orbit<2,3>(d).size()<<std::endl;

		// find the marked dart
		// TODO we assume that there is only one at the moment; which is false with self-intersecting/touching sheets
		bool found_opp = false;

		for(LCC_3::Dart_of_orbit_range<2,3>::iterator
				  itbis(lcc()->darts_of_orbit<2,3>(d).begin()),
				  itend(lcc()->darts_of_orbit<2,3>(d).end()); itbis!=itend; ++itbis) {

			LCC_3::Dart_handle dbis = itbis;

			if((d != dbis) && (lcc()->is_marked(dbis, m))) {

				found_opp = true;

				LCC_3::Dart_handle dbis6 = old_pattern[dbis][6];
				LCC_3::Dart_handle dbis7 = old_pattern[dbis][7];
				LCC_3::Dart_handle dbis8 = old_pattern[dbis][8];
				LCC_3::Dart_handle dbis9 = old_pattern[dbis][9];

				lcc()->link_alpha<2>(d6, dbis6);
				lcc()->link_alpha<2>(d7, dbis7);
				lcc()->link_alpha<2>(d8, dbis8);
				lcc()->link_alpha<2>(d9, dbis9);

				break;
			}
		}

//		if(!found_opp) {
//
//			lcc()->erase_dart(d6);
//			lcc()->erase_dart(d7);
//			lcc()->erase_dart(d8);
//			lcc()->erase_dart(d9);
//		}

//		for(auto dbis: getlcc()->darts_of_orbit<2,3>(d)) {
//			if(getlcc()->is_marked(dbis, m)) {
//
//			}
//		}



	}

	std::cout<<"nbBlocks "<<bl()->nbBlocks()<<std::endl;
//	lcc()->correct_invalid_attributes();

	// alpha3 for the marked darts
	for(auto it: old_alpha3) {
		LCC_3::Dart_handle d = it.first;
		LCC_3::Dart_handle d0 = old_pattern[d][0];
//		LCC_3::Dart_handle d1 = old_pattern[d][1];
//
//		Dart_handle dopp = lcc()->alpha(d, 3);
//		if(d != dopp) {
//			lcc()->link_alpha<3>(d1, dopp);
//		}

		LCC_3::Dart_handle dlink = lcc()->alpha(d, 3);
		if(dlink != d0) {
			lcc()->sew<3>(d, d0);
		}

//		lcc()->link_alpha<3>(d, d0);

	}

	std::cout<<"nbVertices "<<bl()->nbVertices()<<std::endl;
	std::cout<<"nbBlocks "<<bl()->nbBlocks()<<std::endl;
	// TODO clear the orphaned darts


	// free the marks
	lcc()->free_mark(m);
	lcc()->free_mark(m_new);

	return SheetInsert::NOT_YET_IMPLEMENTED;
}
/*----------------------------------------------------------------------------*/
}  // namespace blocking
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/