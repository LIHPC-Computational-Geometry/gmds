/*----------------------------------------------------------------------------*/
#include "gmds/blocking/InputMarkedDarts.h"
/*----------------------------------------------------------------------------*/
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
}  // namespace blocking
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/