/*----------------------------------------------------------------------------*/
#include "gmds/blocking/SheetCollapse.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace blocking {
/*----------------------------------------------------------------------------*/
SheetCollapse::SheetCollapse() {}
/*----------------------------------------------------------------------------*/
SheetCollapse::~SheetCollapse() {}
/*----------------------------------------------------------------------------*/
SheetCollapse::STATUS SheetCollapse::execute(){

	LCC_3::size_type dart_to_collapse = lcc()->get_new_mark();

	std::vector<Dart_handle> sheet = bl_->getSheet(m_origin);
	LCC_3 ::size_type mark = lcc()->get_new_mark();

	for(auto d : sheet){
		lcc()->mark_cell<3>(d,dart_to_collapse);

		lcc()->mark(lcc()->alpha<2,1,0,1,2, 1,0,1,0>(d),mark);
		lcc()->mark(lcc()->alpha<2,1,0,1,2, 1,0,1>(d),mark);
		lcc()->mark(lcc()->alpha<2,1,0,1,2, 1,0>(d),mark);
		lcc()->mark(lcc()->alpha<2,1,0,1,2, 1>(d),mark);
		lcc()->mark(lcc()->alpha<2,1,0,1,2, 0,1,0>(d),mark);
		lcc()->mark(lcc()->alpha<2,1,0,1,2, 0,1>(d),mark);
		lcc()->mark(lcc()->alpha<2,1,0,1,2, 0>(d),mark);
		lcc()->mark(lcc()->alpha<2,1,0,1,2  >(d),mark);

		lcc()->mark(lcc()->alpha<1,0,1,2, 1,0,1,0>(d),mark);
		lcc()->mark(lcc()->alpha<1,0,1,2, 1,0,1>(d),mark);
		lcc()->mark(lcc()->alpha<1,0,1,2, 1,0>(d),mark);
		lcc()->mark(lcc()->alpha<1,0,1,2, 1>(d),mark);
		lcc()->mark(lcc()->alpha<1,0,1,2, 0,1,0>(d),mark);
		lcc()->mark(lcc()->alpha<1,0,1,2, 0,1>(d),mark);
		lcc()->mark(lcc()->alpha<1,0,1,2, 0>(d),mark);
		lcc()->mark(lcc()->alpha<1,0,1,2  >(d),mark);

		lcc()->mark(lcc()->alpha<2, 1,0,1,0>(d),mark);
		lcc()->mark(lcc()->alpha<2, 1,0,1>(d),mark);
		lcc()->mark(lcc()->alpha<2, 1,0>(d),mark);
		lcc()->mark(lcc()->alpha<2, 1>(d),mark);
		lcc()->mark(lcc()->alpha<2, 0,1,0>(d),mark);
		lcc()->mark(lcc()->alpha<2, 0,1>(d),mark);
		lcc()->mark(lcc()->alpha<2, 0>(d),mark);
		lcc()->mark(lcc()->alpha<2  >(d),mark);

		lcc()->mark(lcc()->alpha<1,0,1,0>(d),mark);
		lcc()->mark(lcc()->alpha<1,0,1>(d),mark);
		lcc()->mark(lcc()->alpha<1,0>(d),mark);
		lcc()->mark(lcc()->alpha<1>(d),mark);
		lcc()->mark(lcc()->alpha<0,1,0>(d),mark);
		lcc()->mark(lcc()->alpha<0,1>(d),mark);
		lcc()->mark(lcc()->alpha<0>(d),mark);
		lcc()->mark(d,mark);
	}

	std::map<Dart_handle,Dart_handle> darts_opp;
	//We build the link between the two sides of the sheet to collapse
	//For each hex in the sheet we get the opposite dart of each side
	for(auto d : sheet){
		Dart_handle dart_left;
		dart_left = lcc()->alpha<1,2>(d);

		Dart_handle dart_right;
		dart_right = lcc()->alpha<0,1,2>(d);
		if(!lcc()->is_marked(dart_left,mark) && !lcc()->is_marked(dart_right,mark))
			darts_opp[dart_left] = dart_right;

	}

	//Erasing darts representing parallel edges of the sheet
	erase(mark);

	std::map<Dart_handle,Dart_handle> darts_to_sew;


	for(auto d : darts_opp){

		if(lcc()->is_free<3>(d.second)){
			lcc()->sew<2>(d.second,d.first);
			lcc()->sew<2>(lcc()->alpha<0,1>(d.second),lcc()->alpha<0,1>(d.first));
			lcc()->sew<2>(lcc()->alpha<0,1,0,1>(d.second),lcc()->alpha<0,1,0,1>(d.first));
			lcc()->sew<2>(lcc()->alpha<1,0>(d.second),lcc()->alpha<1,0>(d.first));


		}else{
			lcc()->sew<2>(d.first,d.second);
			lcc()->sew<2>(lcc()->alpha<0,1>(d.first),lcc()->alpha<0,1>(d.second));
			lcc()->sew<2>(lcc()->alpha<0,1,0,1>(d.first),lcc()->alpha<0,1,0,1>(d.second));
			lcc()->sew<2>(lcc()->alpha<1,0>(d.first),lcc()->alpha<1,0>(d.second));
		}
	}
	lcc()->unmark_all(mark);

	for(auto d : darts_opp){
		if(!lcc()->is_marked(d.first,mark)) {
			lcc()->mark(d.first, mark);
			if (!lcc()->is_free<3>(d.first) && !lcc()->is_free<3>(d.second)) {
				Dart_handle dart1 = lcc()->alpha<3>(d.first);
				Dart_handle dart2 = lcc()->alpha<3>(d.second);

				if (lcc()->is_marked(dart1, dart_to_collapse)) {
					lcc()->mark_cell<3>(dart1, mark);
					dart1 = lcc()->alpha<3, 2, 3>(d.first); // a<3,2,3> from d.first because i dont want to iterate from variable we are modifying
																						// -> maybe useless and a<2,3>(dart1) better ?
				}
				if (lcc()->is_marked(dart2, dart_to_collapse)) {
					lcc()->mark_cell<3>(dart2, mark);
					dart2 = lcc()->alpha<3, 2, 3>(d.second);
				}
				darts_to_sew[dart1] = dart2;
			}
		}
	}

	for(auto d : darts_opp){
		lcc()->remove_cell<3>(d.first);
	}

	for(auto d : darts_to_sew){
		lcc()->sew<3>(d.first,d.second);
	}

	lcc()->unmark_all(mark);
	lcc()->free_mark(mark);
	lcc()->unmark_all(dart_to_collapse);
	lcc()->free_mark(dart_to_collapse);

	return SheetCollapse::SUCCESS;
}
/*----------------------------------------------------------------------------*/
unsigned int SheetCollapse::erase(LCC_3::size_type amark)
{
	unsigned int res = 0;
	Dart_handle d;
	for ( typename LCC_3::Dart_range::iterator it(lcc()->darts().begin()),
	     itend(lcc()->darts().end()); it!=itend; )
	{
		d = it++;
		if (lcc()->is_marked(d, amark))
		{
			if (!lcc()->is_free(d, 0)){
				if(lcc()->alpha<0>(lcc()->alpha<0>(d)) == d) {
					if(lcc()->dart_of_attribute<0>(lcc()->vertex_attribute(d)) == d){
						//trying to handle attribute correctly
						lcc()->set_dart_attribute<0>(lcc()->alpha<0>(d),lcc()->vertex_attribute(d));
						if(lcc()->dart_of_attribute<0>(lcc()->vertex_attribute(d)) != d){
						}
					}
					lcc()->unlink_alpha<0>(lcc()->alpha<0>(d));
				}
			}
			if (!lcc()->is_free(d, 1)){
				//lcc()->topo_unsew<0>(d);
				if(lcc()->alpha<1>(lcc()->alpha<1>(d)) == d)
					lcc()->unlink_alpha<1>(lcc()->alpha<1>(d));
			}
			if (!lcc()->is_free(d, 2)){
				//lcc()->topo_unsew<0>(d);
				if(lcc()->alpha<2>(lcc()->alpha<2>(d)) == d)
					lcc()->unlink_alpha<2>(lcc()->alpha<2>(d));
			}
			if (!lcc()->is_free(d, 3)){
				//lcc()->topo_unsew<0>(d);
				if(lcc()->alpha<3>(lcc()->alpha<3>(d)) == d)
					lcc()->unlink_alpha<3>(lcc()->alpha<3>(d));
			}

			lcc()->erase_dart(d); ++res;
		}
	}
	return res;
}
/*----------------------------------------------------------------------------*/
}  // namespace blocking
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/