/*----------------------------------------------------------------------------*/
#include "gmds/blocking/InputMarkedDarts.h"
/*----------------------------------------------------------------------------*/
#include <set>
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
InputMarkedDarts::InputMarkedDarts() {}
/*----------------------------------------------------------------------------*/
InputMarkedDarts::~InputMarkedDarts() {}
/*----------------------------------------------------------------------------*/
int
InputMarkedDarts::insertsheet_mark_first_cells3d(LCC_3* ALcc, LCC_3::size_type AMark, int ANbCells)
{
	// mark the darts to extrude
	LCC_3::size_type m = ALcc->get_new_mark();

	// TODO for now we mark the darts of the n first cell
	const int nbCellsToMark = 5;
	int nbCellsMarked = 0;
	for (auto it = ALcc->one_dart_per_cell<3>().begin(); it != ALcc->one_dart_per_cell<3>().end(); it++) {

		for (LCC_3::Dart_of_cell_range<3>::iterator
		        itbis(ALcc->darts_of_cell<3>(it).begin()),
		     itend(ALcc->darts_of_cell<3>(it).end()); itbis!=itend; ++itbis) {
			ALcc->mark(itbis, m);
		}

		nbCellsMarked++;
		if(nbCellsMarked>=nbCellsToMark) {
			break;
		}
	}

	for (LCC_3::Dart_range::iterator it(ALcc->darts().begin()),
	     itend(ALcc->darts().end()); it!=itend; ++it)
	{
		if(ALcc->is_marked(it, m)) {
			// Do not mark the boundary darts
			if(ALcc->alpha(it, 3) != it) {
				if(!ALcc->is_marked(ALcc->alpha(it, 3), m)) {
					ALcc->mark(it, AMark);
				}
			}
		}
	}

	const int nbMarkedDarts = ALcc->number_of_marked_darts(AMark);
	std::cout<<"InputMarkedDarts::insertsheet_mark_first_cells3d nbMarkedDarts "<< nbMarkedDarts<<std::endl;

	// free the marks
	ALcc->free_mark(m);

	return nbMarkedDarts;
}
/*----------------------------------------------------------------------------*/
int
InputMarkedDarts::pillow_mark_first_cells3d(LCC_3 *ALcc, LCC_3::size_type AMark, int ANbCells)
{
	// mark the darts to extrude
	LCC_3::size_type m = ALcc->get_new_mark();

	// TODO for now we mark the darts of the n first cell
	int nbCellsMarked = 0;
	for (auto it = ALcc->one_dart_per_cell<3>().begin(); it != ALcc->one_dart_per_cell<3>().end(); it++) {

		for (LCC_3::Dart_of_cell_range<3>::iterator
		        itbis(ALcc->darts_of_cell<3>(it).begin()),
		     itend(ALcc->darts_of_cell<3>(it).end()); itbis!=itend; ++itbis) {
			ALcc->mark(itbis, m);
		}

		nbCellsMarked++;
		if(nbCellsMarked>=ANbCells) {
			break;
		}
	}

	for (LCC_3::Dart_range::iterator it(ALcc->darts().begin()),
	     itend(ALcc->darts().end()); it!=itend; ++it)
	{
		if(ALcc->is_marked(it, m)) {
			if(ALcc->alpha(it, 3) == it) {
				ALcc->mark(it, AMark);
			} else {
				if(!ALcc->is_marked(ALcc->alpha(it, 3), m)) {
					ALcc->mark(it, AMark);
				}
			}
		}
	}

	const int nbMarkedDarts = ALcc->number_of_marked_darts(AMark);
	std::cout<<"InputMarkedDarts::pillow_mark_first_cells3d nbMarkedDarts "<< nbMarkedDarts<<std::endl;

	// free the marks
	ALcc->free_mark(m);

	return nbMarkedDarts;
}
/*----------------------------------------------------------------------------*/
int
InputMarkedDarts::insertsheet_mark_intersect_3d(LCC_3 *ALcc, LCC_3::size_type AMark, int ANi, int ANj, int ANk)
{
	std::set<int> cells_to_mark_top;
	std::set<int> cells_to_mark_bottom;
	std::set<int> cells_to_mark_left;
	std::set<int> cells_to_mark_right;
	std::set<int> cells_to_mark_front;
	std::set<int> cells_to_mark_back;

	int index = 0;
	for(int i=0; i<ANi; i++) {
		for(int j=0; j<ANj; j++) {
			for(int k=0; k<ANk; k++) {

				if((k==0) && (i!=0)) {
					cells_to_mark_top.insert(index);
				}
				if((k==ANk-1) && (i!=0) && (i!=ANi-1)) {
					cells_to_mark_bottom.insert(index);
				}
				if((i==ANi-1) && (k!=ANk-1)) {
					cells_to_mark_left.insert(index);
				}
				if((i==0) && (k!=0) && (k!=ANk-1)) {
					cells_to_mark_right.insert(index);
				}
				index++;
			}
		}
	}
	
	// mark the darts to extrude
	index = 0;
	for (auto it = ALcc->one_dart_per_cell<3>().begin(); it != ALcc->one_dart_per_cell<3>().end(); it++) {

		if(cells_to_mark_top.find(index) != cells_to_mark_top.end()) {
			Dart_handle d = ALcc->darts_of_cell<3>(it).begin();
			d = ALcc->alpha(d, 0, 1, 2);

			for (LCC_3::Dart_of_orbit_range<0,1>::iterator
					  it(ALcc->darts_of_orbit<0,1>(d).begin()),
					  itend(ALcc->darts_of_orbit<0,1>(d).end()); it!=itend; ++it) {
				ALcc->mark(it, AMark);
			}
		}
		if(cells_to_mark_bottom.find(index) != cells_to_mark_bottom.end()) {
			Dart_handle d = ALcc->darts_of_cell<3>(it).begin();
			d = ALcc->alpha(d, 1, 2);

			for (LCC_3::Dart_of_orbit_range<0,1>::iterator
			        it(ALcc->darts_of_orbit<0,1>(d).begin()),
			     itend(ALcc->darts_of_orbit<0,1>(d).end()); it!=itend; ++it) {
				ALcc->mark(it, AMark);
			}
	   }
		if(cells_to_mark_left.find(index) != cells_to_mark_left.end()) {
			Dart_handle d = ALcc->darts_of_cell<3>(it).begin();
			d = ALcc->alpha(d, 2);

			for (LCC_3::Dart_of_orbit_range<0,1>::iterator
			        it(ALcc->darts_of_orbit<0,1>(d).begin()),
			     itend(ALcc->darts_of_orbit<0,1>(d).end()); it!=itend; ++it) {
				ALcc->mark(it, AMark);
			}
		}
		if(cells_to_mark_right.find(index) != cells_to_mark_right.end()) {
			Dart_handle d = ALcc->darts_of_cell<3>(it).begin();
			d = ALcc->alpha(d, 1, 0, 1, 2);

			for (LCC_3::Dart_of_orbit_range<0,1>::iterator
			        it(ALcc->darts_of_orbit<0,1>(d).begin()),
			     itend(ALcc->darts_of_orbit<0,1>(d).end()); it!=itend; ++it) {
				ALcc->mark(it, AMark);
			}
		}
		if(cells_to_mark_front.find(index) != cells_to_mark_front.end()) {
			Dart_handle d = ALcc->darts_of_cell<3>(it).begin();
//			d = ALcc->alpha(d, 1, 2);

			for (LCC_3::Dart_of_orbit_range<0,1>::iterator
			        it(ALcc->darts_of_orbit<0,1>(d).begin()),
			     itend(ALcc->darts_of_orbit<0,1>(d).end()); it!=itend; ++it) {
				ALcc->mark(it, AMark);
			}
		}
		if(cells_to_mark_back.find(index) != cells_to_mark_back.end()) {
			Dart_handle d = ALcc->darts_of_cell<3>(it).begin();
			d = ALcc->alpha(d, 2, 1, 0, 1, 2);

			for (LCC_3::Dart_of_orbit_range<0,1>::iterator
			        it(ALcc->darts_of_orbit<0,1>(d).begin()),
			     itend(ALcc->darts_of_orbit<0,1>(d).end()); it!=itend; ++it) {
				ALcc->mark(it, AMark);
			}
		}
		index++;
	}

	return ALcc->number_of_marked_darts(AMark);
}
/*----------------------------------------------------------------------------*/
}  // namespace blocking
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/