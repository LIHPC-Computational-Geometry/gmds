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

#include "gmds/blocking/WriterDartsVTK.h"
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
SheetInsert::execute(LCC_3::size_type AMark)
{
	// check that the implementation can handle the data
	if(!lcc()->is_valid()) {
		std::string s ="SheetInsert::pillow can be applied on a valid gmap only";
		throw gmds::GMDSException(s);
	}
	// TODO check the validity of the shrink set


//	std::cout<<"SheetInsert::execute nbMarked "<< lcc()->number_of_marked_darts(AMark)<<std::endl;

	// TODO check validity of the marked darts
//	// check that on orbit<2,3> two and only two darts are marked
//	for (LCC_3::Dart_range::iterator it(lcc()->darts().begin()),
//	     itend(lcc()->darts().end()); it!=itend; ++it) {
//
//		int nbfound = 1;
//
//		if (lcc()->is_marked(it, m)) {
//			for (LCC_3::Dart_of_orbit_range<2, 3>::iterator itbis(lcc()->darts_of_orbit<2, 3>(it).begin()), itend(lcc()->darts_of_orbit<2, 3>(it).end());
//			     itbis != itend; ++itbis) {
//
//				LCC_3::Dart_handle dbis = itbis;
//
//				if ((it != dbis) && (lcc()->is_marked(dbis, m))) {
//					nbfound++;
//				}
//			}
//
//			if (nbfound != 2) {
//				std::string s ="wrong number of marked darts on orbit<2,3> " + std::to_string(nbfound);
//				throw gmds::GMDSException(s);
//			}
//		}
//	}

	// store current alpha links that will be invalidated by the unsew
	std::map<LCC_3::Dart_handle, LCC_3::Dart_handle> old_alpha3;

	for (LCC_3::Dart_range::iterator it(lcc()->darts().begin()),
	     itend(lcc()->darts().end()); it!=itend; ++it) {

		if(lcc()->is_marked(it, AMark)) {
			LCC_3::Dart_handle d0 = it;
			LCC_3::Dart_handle d3 = lcc()->alpha(d0,3);

			old_alpha3.emplace(d0, d3);
		}
	}
//	std::cout<<"old_alpha3.size() "<< old_alpha3.size()<<std::endl;

	// unsew the alpha3
	for(auto d: old_alpha3) {
		if(lcc()->alpha(d.first, 3) != d.first) {
			lcc()->unsew<3>(d.first);
		}
	}
	lcc()->is_valid();

	// create the pattern for the marked darts
	LCC_3::size_type m_new = lcc()->get_new_mark();
	std::map<LCC_3::Dart_handle, std::array<LCC_3::Dart_handle, 12> > old_pattern;

	std::map<LCC_3::Vertex_attribute_handle, LCC_3::Vertex_attribute_handle> old2new_vertices;

	LCC_3::size_type m_cells_center = lcc()->get_new_mark();

	for(auto d: old_alpha3) {


		LCC_3::Vertex_attribute_handle v = lcc()->vertex_attribute(d.first);
		if(old2new_vertices.find(v) == old2new_vertices.end()) {

			// TODO determine position of new vertex and which darts are assigned to it
			double posx = 0.;
			double posy = 0.;
			double posz = 0.;
			int nb_marked_Cells = 0;

			lcc()->darts_of_orbit<1,2,3>(d.first);
			std::set<LCC_3::Dart_handle> darts;
			for (LCC_3::Dart_of_orbit_range<1, 2, 3>::iterator it(lcc()->darts_of_orbit<1, 2, 3>(d.first).begin()), itend(lcc()->darts_of_orbit<1, 2, 3>(d.first).end());
			     it != itend; ++it) {

				LCC_3::Dart_handle dbis = it;

				if(!lcc()->is_marked(it, m_cells_center)) {
					lcc()->mark_cell<3>(it, m_cells_center);
					nb_marked_Cells++;
					LCC_3::Point center_tmp = lcc()->barycenter<3>(it);

					posx += center_tmp.x();
					posy += center_tmp.y();
					posz += center_tmp.z();
				}
			}
			posx /= nb_marked_Cells;
			posy /= nb_marked_Cells;
			posz /= nb_marked_Cells;
			lcc()->unmark_all(m_cells_center);

			posx += ratio_vertex_insert_ * (v->point().x() - posx);
			posy += ratio_vertex_insert_ * (v->point().y() - posy);
			posz += ratio_vertex_insert_ * (v->point().z() - posz);


			if(screenshot_) {
//				posx = v->point().x();
//				posy = v->point().y();
//				posz = 0.5;
				if ((posx > 3.24) && (posx < 3.27) && (posz > 1.24) && (posz < 1.27)) {
					posz = 1.1;
				}

				if ((posx > 2.73) && (posx < 2.75) && (posz > 0.69) && (posz < 0.72)) {
					posx = 2.9;
				}

				if ((posx > 3.24) && (posx < 3.3) && (posz > 0.21) && (posz < 0.26)) {
					posz = 0.;
				}

				if ((posx > 3.79) && (posx < 3.9) && (posz > 0.71) && (posz < 0.76)) {
					posx = 4.;
				}
			}

			LCC_3::Point newpt(v->point().x(), v->point().y(), v->point().z());
			v->point() = LCC_3::Point (posx, posy, posz);
			LCC_3::Vertex_attribute_handle vbis = lcc()->create_vertex_attribute( newpt);
			old2new_vertices.emplace(v, vbis);
		}
		LCC_3::Vertex_attribute_handle vbis = old2new_vertices[v];

		LCC_3::Dart_handle d0 = lcc()->create_dart(v);
		LCC_3::Dart_handle d1 = lcc()->create_dart(vbis);
		LCC_3::Dart_handle d2 = lcc()->create_dart(v);
		LCC_3::Dart_handle d3 = lcc()->create_dart(v);
		LCC_3::Dart_handle d4 = lcc()->create_dart(vbis);
		LCC_3::Dart_handle d5 = lcc()->create_dart(vbis);
		LCC_3::Dart_handle d6 = lcc()->create_dart(v);
		LCC_3::Dart_handle d7 = lcc()->create_dart(v);
		LCC_3::Dart_handle d8 = lcc()->create_dart(vbis);
		LCC_3::Dart_handle d9 = lcc()->create_dart(vbis);
		LCC_3::Dart_handle d10 = lcc()->create_dart(v);
		LCC_3::Dart_handle d11 = lcc()->create_dart(vbis);

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
		lcc()->mark(d10, m_new);
		lcc()->mark(d11, m_new);

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

		lcc()->link_alpha<0>(d10, d11);
		lcc()->link_alpha<2>(d7, d10);
		lcc()->link_alpha<2>(d8, d11);

		old_pattern.emplace(d.first, std::array<LCC_3::Dart_handle, 12> ({d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11}));
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

	gmds::blocking::WriterDartsVTK writer;
	writer.setBl(bl());
	gmds::blocking::WriterDartsVTK::STATUS st = writer.execute("darts_first_link.vtk");

	// link the created darts between the patterns
	// second step with the links that depend on the orbits
	for(auto it: old_alpha3) {
		LCC_3::Dart_handle d = it.first;
		LCC_3::Dart_handle d6 = old_pattern[d][6];
		LCC_3::Dart_handle d7 = old_pattern[d][7];
		LCC_3::Dart_handle d8 = old_pattern[d][8];
		LCC_3::Dart_handle d9 = old_pattern[d][9];
		LCC_3::Dart_handle d10 = old_pattern[d][10];
		LCC_3::Dart_handle d11 = old_pattern[d][11];

		lcc()->darts_of_orbit<2, 3>(d);
		//		std::cout<<"orbit.size() "<< lcc()->darts_of_orbit<2,3>(d).size()<<std::endl;

		// find the marked dart
		bool found_other = false;

		for (LCC_3::Dart_of_orbit_range<2, 3>::iterator itbis(lcc()->darts_of_orbit<2, 3>(d).begin()), itend(lcc()->darts_of_orbit<2, 3>(d).end());
		     itbis != itend; ++itbis) {

			LCC_3::Dart_handle dbis = itbis;

			if ((d != dbis) && (lcc()->is_marked(dbis, AMark))) {

				found_other = true;

				LCC_3::Dart_handle dbis6 = old_pattern[dbis][6];
				LCC_3::Dart_handle dbis7 = old_pattern[dbis][7];
				LCC_3::Dart_handle dbis8 = old_pattern[dbis][8];
				LCC_3::Dart_handle dbis9 = old_pattern[dbis][9];
				LCC_3::Dart_handle dbis10 = old_pattern[dbis][10];
				LCC_3::Dart_handle dbis11 = old_pattern[dbis][11];

				lcc()->link_alpha<2>(d6, dbis6);
				lcc()->link_alpha<1>(d10, dbis10);

				break;
			}
		}

		// check if we are in an autointersect situation
		LCC_3::Dart_handle d_opp = it.second;
		bool found_opp = false;

		for (LCC_3::Dart_of_orbit_range<2, 3>::iterator itbis(lcc()->darts_of_orbit<2, 3>(d_opp).begin()), itend(lcc()->darts_of_orbit<2, 3>(d_opp).end());
		     itbis != itend; ++itbis) {

			LCC_3::Dart_handle dopp = itbis;

			if (lcc()->is_marked(dopp, AMark)) {
				found_opp = true;

				LCC_3::Dart_handle dopp6 = old_pattern[dopp][6];
				LCC_3::Dart_handle dopp10 = old_pattern[dopp][10];


				// beware of the link order : the opposite should carry the good vertex attribute
//				lcc()->link_alpha<2>(d9, dopp6);
//				lcc()->link_alpha<1>(d11, dopp10);
				lcc()->link_alpha<2>(dopp6, d9);
				lcc()->link_alpha<1>(dopp10, d11);
				break;
			}
		}

		// this is the regular case
		if (found_other && (!found_opp)) {

			for (LCC_3::Dart_of_orbit_range<2, 3>::iterator itbis(lcc()->darts_of_orbit<2, 3>(d).begin()), itend(lcc()->darts_of_orbit<2, 3>(d).end());
			     itbis != itend; ++itbis) {

				LCC_3::Dart_handle dbis = itbis;

				if ((d != dbis) && (lcc()->is_marked(dbis, AMark))) {

					LCC_3::Dart_handle dbis9 = old_pattern[dbis][9];
					LCC_3::Dart_handle dbis11 = old_pattern[dbis][11];

					lcc()->link_alpha<2>(d9, dbis9);
					lcc()->link_alpha<1>(d11, dbis11);
				}
			}
		}
	}
	lcc()->correct_invalid_attributes();

	writer.execute("darts_second_link.vtk");

	LCC_3::size_type mark_chord = lcc()->get_new_mark();

	// last round to close the inserted chord
	for(auto it: old_alpha3) {
		LCC_3::Dart_handle d = it.first;
		LCC_3::Dart_handle d6 = old_pattern[d][6];
		LCC_3::Dart_handle d9 = old_pattern[d][9];
		LCC_3::Dart_handle d10 = old_pattern[d][10];
		LCC_3::Dart_handle d11 = old_pattern[d][11];

		if((lcc()->alpha(d11, 1) == d11) && (lcc()->alpha(d10, 1) != d10)) {
			lcc()->link_alpha<1>(d11, lcc()->alpha(d11, 0, 1, 0, 1, 0, 1, 0));
			lcc()->link_alpha<2>(d9, lcc()->alpha(d9, 1, 0, 1, 2, 1, 0, 1, 2, 1, 0, 1, 2, 1, 0, 1));

			lcc()->mark(d11, mark_chord);
			lcc()->mark(lcc()->alpha(d11, 0, 1, 0, 1, 0, 1, 0), mark_chord);
			lcc()->mark(lcc()->alpha(d11, 0, 1, 0, 1, 0, 1), mark_chord);
			lcc()->mark(lcc()->alpha(d11, 0, 1, 0, 1, 0), mark_chord);
			lcc()->mark(lcc()->alpha(d11, 0, 1, 0, 1), mark_chord);
			lcc()->mark(lcc()->alpha(d11, 0, 1, 0), mark_chord);
			lcc()->mark(lcc()->alpha(d11, 0, 1), mark_chord);
			lcc()->mark(lcc()->alpha(d11, 0), mark_chord);
		}
	}
	lcc()->correct_invalid_attributes();

	// remove unused patterns
	for(auto it: old_alpha3) {
		LCC_3::Dart_handle d = it.first;
		LCC_3::Dart_handle d6 = old_pattern[d][6];
		LCC_3::Dart_handle d7 = old_pattern[d][7];
		LCC_3::Dart_handle d8 = old_pattern[d][8];
		LCC_3::Dart_handle d9 = old_pattern[d][9];
		LCC_3::Dart_handle d10 = old_pattern[d][10];
		LCC_3::Dart_handle d11 = old_pattern[d][11];

		if (lcc()->alpha(d6, 2) == d6) {

			lcc()->unlink_alpha(d6, 3);
			lcc()->unlink_alpha(d7, 3);
			lcc()->unlink_alpha(d8, 3);
			lcc()->unlink_alpha(d9, 3);

			lcc()->erase_dart(d6);
			lcc()->erase_dart(d7);
			lcc()->erase_dart(d8);
			lcc()->erase_dart(d9);
			lcc()->erase_dart(d10);
			lcc()->erase_dart(d11);
		}
	}

	// alpha3 for the marked darts
	for(auto it: old_alpha3) {
		LCC_3::Dart_handle d = it.first;
		LCC_3::Dart_handle d0 = old_pattern[d][0];
		LCC_3::Dart_handle d1 = old_pattern[d][1];

		Dart_handle dopp = it.second;
		if(d == dopp) {
			lcc()->link_alpha<3>(d, d0);
		} else {
			lcc()->link_alpha<3>(d, d0);
			lcc()->link_alpha<3>(dopp, d1);
		}
	}
	//	TODO try to get rid of this automatic correction
	lcc()->correct_invalid_attributes();

	writer.execute("darts_third_link.vtk");

	// TODO clear the orphaned darts and flat cells
	int nbRemoved = cleanup_flat_cells(AMark);
//	std::cout<<"insertsheet nbRemoved "<<nbRemoved<<std::endl;

	lcc()->correct_invalid_attributes();

	// very last round for alpha-3 of the chords
	for(auto it: old_alpha3) {
		LCC_3::Dart_handle d = it.first;
		LCC_3::Dart_handle d10 = old_pattern[d][10];
		LCC_3::Dart_handle d11 = old_pattern[d][11];

		if (bl()->lcc()->is_marked(d10,mark_chord)) {

			if (bl()->lcc()->is_free(d10, 3)) {
				// consider only darts from hex
				if (bl()->lcc()->darts_of_orbit<0, 1, 2>(d10).size() == 48) {
				//	std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAA " << lcc()->darts_of_orbit<2, 3>(d10).size() << std::endl;

					for (LCC_3::Dart_of_orbit_range<2, 3>::iterator itbis(lcc()->darts_of_orbit<2, 3>(d10).begin()), itend(lcc()->darts_of_orbit<2, 3>(d10).end());
					     itbis != itend; ++itbis) {

						if (bl()->lcc()->is_marked(itbis,mark_chord)) {
							if (bl()->lcc()->alpha(itbis, 3) == itbis) {
						//		std::cout << "FREE" << std::endl;
								if ((d10 != itbis)) {
						//			std::cout << "DIFFERENT " << bl()->lcc()->darts_of_orbit<0, 1, 2>(itbis).size() << std::endl;
									// check that itbis is a dart from an hex
									if (bl()->lcc()->darts_of_orbit<0, 1, 2>(itbis).size() == 48) {

										//						bl()->lcc()->link_alpha<3>(d10, itbis);
										//						bl()->lcc()->link_alpha<3>(d11, bl()->lcc()->alpha(itbis, 0));
										bl()->lcc()->sew<3>(d10, itbis);
								//		std::cout << "QQQQQQQQQQQQQQQQQQQQQQQQQQ" << std::endl;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	lcc()->correct_invalid_attributes();

	writer.execute("darts_chord_link.vtk");

	// free the marks
	lcc()->free_mark(m_new);
	lcc()->free_mark(m_cells_center);

	return SheetInsert::SUCCESS;
}
/*----------------------------------------------------------------------------*/
SheetInsert::STATUS
SheetInsert::pillow(LCC_3::size_type AMark)
{
	//TODO get the set of 3-cells in another way


	// check that the implementation can handle the data
	if(!lcc()->is_valid()) {
		std::string s ="SheetInsert::pillow can be applied on a valid gmap only";
		throw gmds::GMDSException(s);
	}
	// TODO check the validity of the shrink set


//	std::cout<<"SheetInsert::pillow nbMarked "<< lcc()->number_of_marked_darts(AMark)<<std::endl;

	// check that on orbit<2,3> two and only two darts are marked
	for (LCC_3::Dart_range::iterator it(lcc()->darts().begin()),
			  itend(lcc()->darts().end()); it!=itend; ++it) {

		int nbfound = 1;

		if (lcc()->is_marked(it, AMark)) {
			for (LCC_3::Dart_of_orbit_range<2, 3>::iterator itbis(lcc()->darts_of_orbit<2, 3>(it).begin()), itend(lcc()->darts_of_orbit<2, 3>(it).end());
			     itbis != itend; ++itbis) {

				LCC_3::Dart_handle dbis = itbis;

				if ((it != dbis) && (lcc()->is_marked(dbis, AMark))) {
					nbfound++;
				}
			}

			if (nbfound != 2) {
				std::string s ="wrong number of marked darts on orbit<2,3> " + std::to_string(nbfound);
				throw gmds::GMDSException(s);
			}
		}
	}

	// store current alpha links that will be invalidated by the unsew
	std::map<LCC_3::Dart_handle, LCC_3::Dart_handle> old_alpha3;

	for (LCC_3::Dart_range::iterator it(lcc()->darts().begin()),
			  itend(lcc()->darts().end()); it!=itend; ++it) {

		if(lcc()->is_marked(it, AMark)) {
			LCC_3::Dart_handle d0 = it;
			LCC_3::Dart_handle d3 = lcc()->alpha(d0,3);

			old_alpha3.emplace(d0, d3);
		}
	}
	//std::cout<<"old_alpha3.size() "<< old_alpha3.size()<<std::endl;

	// unsew the alpha3
	for(auto d: old_alpha3) {
		if(lcc()->alpha(d.first, 3) != d.first) {
			lcc()->unsew<3>(d.first);
		}
	}
	lcc()->is_valid();
//	lcc()->correct_invalid_attributes();

	// create the pattern for the marked darts
	LCC_3::size_type m_new = lcc()->get_new_mark();
	std::map<LCC_3::Dart_handle, std::array<LCC_3::Dart_handle, 12> > old_pattern;

	std::map<LCC_3::Vertex_attribute_handle, LCC_3::Vertex_attribute_handle> old2new_vertices;
	LCC_3::size_type m_cells_center = lcc()->get_new_mark();

	for(auto d: old_alpha3) {
		LCC_3::Vertex_attribute_handle v = lcc()->vertex_attribute(d.first);
		if(old2new_vertices.find(v) == old2new_vertices.end()) {
			// TODO determine position of new vertex and which darts are assigned to it=
			double posx = 0.;
			double posy = 0.;
			double posz = 0.;
			int nb_marked_Cells = 0;

			lcc()->darts_of_orbit<1,2,3>(d.first);
			std::set<LCC_3::Dart_handle> darts;
			for (LCC_3::Dart_of_orbit_range<1, 2, 3>::iterator it(lcc()->darts_of_orbit<1, 2, 3>(d.first).begin()), itend(lcc()->darts_of_orbit<1, 2, 3>(d.first).end());
			     it != itend; ++it) {

				LCC_3::Dart_handle dbis = it;

				if(!lcc()->is_marked(it, m_cells_center)) {
					lcc()->mark_cell<3>(it, m_cells_center);
					nb_marked_Cells++;
					LCC_3::Point center_tmp = lcc()->barycenter<3>(it);

					posx += center_tmp.x();
					posy += center_tmp.y();
					posz += center_tmp.z();
				}
			}
			posx /= nb_marked_Cells;
			posy /= nb_marked_Cells;
			posz /= nb_marked_Cells;
			lcc()->unmark_all(m_cells_center);

			posx += ratio_vertex_insert_ * (v->point().x() - posx);
			posy += ratio_vertex_insert_ * (v->point().y() - posy);
			posz += ratio_vertex_insert_ * (v->point().z() - posz);

			LCC_3::Point newpt(v->point().x(), v->point().y(), v->point().z());
			v->point() = LCC_3::Point (posx, posy, posz);
			LCC_3::Vertex_attribute_handle vbis = lcc()->create_vertex_attribute( newpt);
			old2new_vertices.emplace(v, vbis);
		}
		LCC_3::Vertex_attribute_handle vbis = old2new_vertices[v];

		LCC_3::Dart_handle d0 = lcc()->create_dart(v);
		LCC_3::Dart_handle d1 = lcc()->create_dart(vbis);
		LCC_3::Dart_handle d2 = lcc()->create_dart(v);
		LCC_3::Dart_handle d3 = lcc()->create_dart(v);
		LCC_3::Dart_handle d4 = lcc()->create_dart(vbis);
		LCC_3::Dart_handle d5 = lcc()->create_dart(vbis);
		LCC_3::Dart_handle d6 = lcc()->create_dart(v);
		LCC_3::Dart_handle d7 = lcc()->create_dart(v);
		LCC_3::Dart_handle d8 = lcc()->create_dart(vbis);
		LCC_3::Dart_handle d9 = lcc()->create_dart(vbis);
		LCC_3::Dart_handle d10 = lcc()->create_dart(v);
		LCC_3::Dart_handle d11 = lcc()->create_dart(vbis);

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
		lcc()->mark(d10, m_new);
		lcc()->mark(d11, m_new);

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

		lcc()->link_alpha<0>(d10, d11);
		lcc()->link_alpha<2>(d7, d10);
		lcc()->link_alpha<2>(d8, d11);

		old_pattern.emplace(d.first, std::array<LCC_3::Dart_handle, 12> ({d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11}));
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
		LCC_3::Dart_handle d10 = old_pattern[d][10];
		LCC_3::Dart_handle d11 = old_pattern[d][11];

		lcc()->darts_of_orbit<2,3>(d);
//		std::cout<<"orbit.size() "<< lcc()->darts_of_orbit<2,3>(d).size()<<std::endl;

		// find the marked dart
		// TODO we assume that there is only one at the moment; which is false with self-intersecting/touching sheets
		bool found_opp = false;

		for(LCC_3::Dart_of_orbit_range<2,3>::iterator
				  itbis(lcc()->darts_of_orbit<2,3>(d).begin()),
				  itend(lcc()->darts_of_orbit<2,3>(d).end()); itbis!=itend; ++itbis) {

			LCC_3::Dart_handle dbis = itbis;

			if((d != dbis) && (lcc()->is_marked(dbis, AMark))) {

				found_opp = true;

				LCC_3::Dart_handle dbis6 = old_pattern[dbis][6];
//				LCC_3::Dart_handle dbis7 = old_pattern[dbis][7];
//				LCC_3::Dart_handle dbis8 = old_pattern[dbis][8];
				LCC_3::Dart_handle dbis9 = old_pattern[dbis][9];
				LCC_3::Dart_handle dbis10 = old_pattern[dbis][10];
				LCC_3::Dart_handle dbis11 = old_pattern[dbis][11];

				lcc()->link_alpha<2>(d6, dbis6);
//				lcc()->link_alpha<2>(d7, dbis7);
//				lcc()->link_alpha<2>(d8, dbis8);
				lcc()->link_alpha<2>(d9, dbis9);
				lcc()->link_alpha<1>(d10, dbis10);
				lcc()->link_alpha<1>(d11, dbis11);

				break;
			}
		}

		if(!found_opp) {
			std::string s = "SheetInsert::pillow cannot have only one dart marked on orbit 2,3";
			throw gmds::GMDSException(s);
		}
	}

//	std::cout<<"nbBlocks "<<bl()->nbBlocks()<<std::endl;
//	lcc()->correct_invalid_attributes();

	// alpha3 for the marked darts
	for(auto it: old_alpha3) {
		LCC_3::Dart_handle d = it.first;
		LCC_3::Dart_handle d0 = old_pattern[d][0];
		LCC_3::Dart_handle d1 = old_pattern[d][1];

		Dart_handle dopp = it.second;
		if(d == dopp) {
			lcc()->link_alpha<3>(d, d0);
		} else {
			lcc()->link_alpha<3>(d, d0);
			lcc()->link_alpha<3>(d1, dopp);
		}
	}
	// TODO try to get rid of this automatic correction
	lcc()->correct_invalid_attributes();

	//std::cout<<"nbVertices "<<bl()->nbVertices()<<std::endl;
//	std::cout<<"nbBlocks "<<bl()->nbBlocks()<<std::endl;

	// TODO clear the orphaned darts and flat cells
	int nbRemoved = cleanup_flat_cells(AMark);
	//std::cout<<"pillow nbRemoved "<<nbRemoved<<std::endl;

	// free the marks
	lcc()->free_mark(m_new);
	lcc()->free_mark(m_cells_center);

	return SheetInsert::SUCCESS;
}
/*----------------------------------------------------------------------------*/
SheetInsert::STATUS
SheetInsert::buildCADfromGrid(gmds::math::Point APmin, gmds::math::Point APmax, int ANx, int ANy, int ANz)
{
	gmds::Mesh& mesh = cad_.getMeshView();
	gmds::Node corners[8];
	corners[0] = mesh.newNode(APmin.X(), APmin.Y(), APmin.Z());
	corners[1] = mesh.newNode(APmax.X(), APmin.Y(), APmin.Z());
	corners[2] = mesh.newNode(APmax.X(), APmax.Y(), APmin.Z());
	corners[3] = mesh.newNode(APmin.X(), APmax.Y(), APmin.Z());
	corners[4] = mesh.newNode(APmin.X(), APmin.Y(), APmax.Z());
	corners[5] = mesh.newNode(APmax.X(), APmin.Y(), APmax.Z());
	corners[6] = mesh.newNode(APmax.X(), APmax.Y(), APmax.Z());
	corners[7] = mesh.newNode(APmin.X(), APmax.Y(), APmax.Z());

	gmds::Face triangles[12];
	triangles[ 0] = mesh.newTriangle(corners[0], corners[1], corners[5]); // front
	triangles[ 1] = mesh.newTriangle(corners[0], corners[5], corners[4]);
	triangles[ 2] = mesh.newTriangle(corners[2], corners[3], corners[7]); // back
	triangles[ 3] = mesh.newTriangle(corners[2], corners[7], corners[6]);
	triangles[ 4] = mesh.newTriangle(corners[3], corners[0], corners[4]); // left
	triangles[ 5] = mesh.newTriangle(corners[3], corners[4], corners[7]);
	triangles[ 6] = mesh.newTriangle(corners[1], corners[2], corners[6]); // right
	triangles[ 7] = mesh.newTriangle(corners[1], corners[6], corners[5]);
	triangles[ 8] = mesh.newTriangle(corners[0], corners[3], corners[2]); // bottom
	triangles[ 9] = mesh.newTriangle(corners[0], corners[2], corners[1]);
	triangles[10] = mesh.newTriangle(corners[4], corners[5], corners[6]); // top
	triangles[11] = mesh.newTriangle(corners[4], corners[6], corners[7]);

	gmds::CellGroup<gmds::Node>* corner0 = mesh.newGroup<gmds::Node>("corner_0");
	corner0->add(corners[0]);
	gmds::CellGroup<gmds::Node>* corner1 = mesh.newGroup<gmds::Node>("corner_1");
	corner1->add(corners[1]);
	gmds::CellGroup<gmds::Node>* corner2 = mesh.newGroup<gmds::Node>("corner_2");
	corner2->add(corners[2]);
	gmds::CellGroup<gmds::Node>* corner3 = mesh.newGroup<gmds::Node>("corner_3");
	corner3->add(corners[3]);
	gmds::CellGroup<gmds::Node>* corner4 = mesh.newGroup<gmds::Node>("corner_4");
	corner4->add(corners[4]);
	gmds::CellGroup<gmds::Node>* corner5 = mesh.newGroup<gmds::Node>("corner_5");
	corner5->add(corners[5]);
	gmds::CellGroup<gmds::Node>* corner6 = mesh.newGroup<gmds::Node>("corner_6");
	corner6->add(corners[6]);
	gmds::CellGroup<gmds::Node>* corner7 = mesh.newGroup<gmds::Node>("corner_7");
	corner7->add(corners[7]);

	gmds::CellGroup<gmds::Node>* curve_0 = mesh.newGroup<gmds::Node>("curve_0");
	curve_0->add(corners[0]);
	curve_0->add(corners[1]);
	gmds::CellGroup<gmds::Node>* curve_1 = mesh.newGroup<gmds::Node>("curve_1");
	curve_1->add(corners[1]);
	curve_1->add(corners[2]);
	gmds::CellGroup<gmds::Node>* curve_2 = mesh.newGroup<gmds::Node>("curve_2");
	curve_2->add(corners[2]);
	curve_2->add(corners[3]);
	gmds::CellGroup<gmds::Node>* curve_3 = mesh.newGroup<gmds::Node>("curve_3");
	curve_3->add(corners[3]);
	curve_3->add(corners[0]);
	gmds::CellGroup<gmds::Node>* curve_4 = mesh.newGroup<gmds::Node>("curve_4");
	curve_4->add(corners[4]);
	curve_4->add(corners[5]);
	gmds::CellGroup<gmds::Node>* curve_5 = mesh.newGroup<gmds::Node>("curve_5");
	curve_5->add(corners[5]);
	curve_5->add(corners[6]);
	gmds::CellGroup<gmds::Node>* curve_6 = mesh.newGroup<gmds::Node>("curve_6");
	curve_6->add(corners[6]);
	curve_6->add(corners[7]);
	gmds::CellGroup<gmds::Node>* curve_7 = mesh.newGroup<gmds::Node>("curve_7");
	curve_7->add(corners[7]);
	curve_7->add(corners[4]);
	gmds::CellGroup<gmds::Node>* curve_8 = mesh.newGroup<gmds::Node>("curve_8");
	curve_8->add(corners[0]);
	curve_8->add(corners[4]);
	gmds::CellGroup<gmds::Node>* curve_9 = mesh.newGroup<gmds::Node>("curve_9");
	curve_9->add(corners[1]);
	curve_9->add(corners[5]);
	gmds::CellGroup<gmds::Node>* curve_10 = mesh.newGroup<gmds::Node>("curve_10");
	curve_10->add(corners[2]);
	curve_10->add(corners[6]);
	gmds::CellGroup<gmds::Node>* curve_11 = mesh.newGroup<gmds::Node>("curve_11");
	curve_11->add(corners[3]);
	curve_11->add(corners[7]);

	gmds::CellGroup<gmds::Face>* surf_front = mesh.newGroup<gmds::Face>("surf_front");
	surf_front->add(triangles[0]);
	surf_front->add(triangles[1]);
	gmds::CellGroup<gmds::Face>* surf_back = mesh.newGroup<gmds::Face>("surf_back");
	surf_back->add(triangles[2]);
	surf_back->add(triangles[3]);
	gmds::CellGroup<gmds::Face>* surf_left = mesh.newGroup<gmds::Face>("surf_left");
	surf_left->add(triangles[4]);
	surf_left->add(triangles[5]);
	gmds::CellGroup<gmds::Face>* surf_right = mesh.newGroup<gmds::Face>("surf_right");
	surf_right->add(triangles[6]);
	surf_right->add(triangles[7]);
	gmds::CellGroup<gmds::Face>* surf_bottom = mesh.newGroup<gmds::Face>("surf_bottom");
	surf_bottom->add(triangles[8]);
	surf_bottom->add(triangles[9]);
	gmds::CellGroup<gmds::Face>* surf_top = mesh.newGroup<gmds::Face>("surf_top");
	surf_top->add(triangles[10]);
	surf_top->add(triangles[11]);

	cad_.updateFromMesh();

	std::vector<gmds::cad::GeomPoint*> vertices;
	std::vector<gmds::cad::GeomCurve*> curves;
	std::vector<gmds::cad::GeomSurface*> surfaces;
	cad_.getPoints(vertices);
	cad_.getCurves(curves);
	cad_.getSurfaces(surfaces);
//	std::map<std::string, std::uintptr_t> name2Points;
//	std::map<std::string, std::uintptr_t> name2Curves;
//	std::map<std::string, std::uintptr_t> name2Surfaces;
//	for(auto v: vertices) {
//		name2Points[v->name()] = reinterpret_cast<std::uintptr_t> (v);
//	}
//	for(auto c: curves) {
//		name2Curves[c->name()] = reinterpret_cast<std::uintptr_t> (c);
//	}
//	for(auto s: surfaces) {
//		name2Surfaces[s->name()] = reinterpret_cast<std::uintptr_t> (s);
//	}
	std::map<std::string, gmds::cad::GeomPoint*> name2Points;
	std::map<std::string, gmds::cad::GeomCurve*> name2Curves;
	std::map<std::string, gmds::cad::GeomSurface*> name2Surfaces;
	for(auto v: vertices) {
		name2Points[v->name()] = v;
	}
	for(auto c: curves) {
		name2Curves[c->name()] = c;
	}
	for(auto s: surfaces) {
		name2Surfaces[s->name()] = s;
	}

	// TODO this is valid only because we know how the grid was built
	Dart_handle d_first = lcc()->one_dart_per_cell<3>().begin();

	// associate corners
	Dart_handle d = d_first;
	gmds::cad::FACPoint* pt = dynamic_cast<gmds::cad::FACPoint*> (name2Points["corner_0"]);
	d2p_.emplace(d, pt);

	d = d_first;
	for(int i=0; i<ANx; i++) {
		d = lcc()->alpha(d, 1, 0, 1, 2, 3, 2);
	}
	d = lcc()->alpha(d, 1, 0);
	pt = dynamic_cast<gmds::cad::FACPoint*> (name2Points["corner_1"]);
	d2p_.emplace(d, pt);

	d = d_first;
	for(int i=0; i<ANx; i++) {
		d = lcc()->alpha(d, 1, 0, 1, 2, 3, 2);
	}
	d = lcc()->alpha(d, 1, 0, 1, 2);
	for(int j=0; j<ANy; j++) {
		d = lcc()->alpha(d, 1, 0, 1, 2, 3, 2);
	}
	d = lcc()->alpha(d, 1, 0);
	pt = dynamic_cast<gmds::cad::FACPoint*> (name2Points["corner_2"]);
	d2p_.emplace(d, pt);


//	gmds::cad::FACPoint* pt = reinterpret_cast<gmds::cad::FACPoint*> (name2Points["corner_0"]);
//	d2p_.emplace(d, dynamic_cast<gmds::cad::FACPoint*> (name2Points["corner_0"]));


	return SheetInsert::NOT_YET_IMPLEMENTED;
}
/*----------------------------------------------------------------------------*/
int
SheetInsert::cleanup_flat_cells(LCC_3::size_type AMark)
{
	int nbRemoved = 0;

	for (auto it = bl()->lcc()->one_dart_per_cell<3>().begin(); it != bl()->lcc()->one_dart_per_cell<3>().end(); it++) {
		int nb = bl()->lcc()->darts_of_orbit<0, 1, 2>(it).size();
		switch (nb) {
		case 24: {
			LCC_3::Dart_handle d = it;

			bool removed = false;

			// find a dart with an alpha3 and its
			for (LCC_3::Dart_of_orbit_range<1, 2, 3>::iterator it(lcc()->darts_of_orbit<1, 2, 3>(d).begin()), itend(lcc()->darts_of_orbit<1, 2, 3>(d).end());
			     it != itend; ++it) {

				Dart_handle alpha3_0 = bl()->lcc()->alpha(d, 3);
				if(alpha3_0 != d) {

					Dart_handle alpha2_0 = bl()->lcc()->alpha(d, 2);
					Dart_handle alpha3_1 = bl()->lcc()->alpha(d, 2, 3);
					if(alpha3_1 != alpha2_0) {
						bl()->lcc()->remove_cell<3>(d);
						bl()->lcc()->sew<3> (alpha3_0, alpha3_1);
						removed = true;
						break;
					}
				}
			}
			if(removed) {
				nbRemoved++;
			} else {
				std::string s = "SheetInsert::cleanup_flat_cells could not remove a 24 darts cell";
				throw gmds::GMDSException(s);
			}
			break;
		}
		default:
			break;
		}
	}

	return nbRemoved;
}
/*----------------------------------------------------------------------------*/
}  // namespace blocking
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/