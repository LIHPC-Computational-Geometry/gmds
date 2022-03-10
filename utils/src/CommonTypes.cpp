/*----------------------------------------------------------------------------*/
/*
 * CommonTypes.cpp
 *
 *  Created on: 6 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include "../inc/gmds/utils/CommonTypes.h"
#include <algorithm>
#include <set>
/*----------------------------------------------------------------------------*/
namespace gmds{
 const int MeshModel::Full2D =  DIM3 | F | E | N | F2N | E2N | E2F | F2E | N2F;

    /*----------------------------------------------------------------------------*/
std::vector<TCellID>
getCommonBut(const std::vector<TCellID>& AS1,const std::vector<TCellID>& AS2, const TCellID ABut)
{
	std::vector<TCellID> result;
	unsigned int size1 = AS1.size();
	unsigned int size2 = AS2.size();
	for (unsigned int i1=0; i1<size1; i1++){
		TCellID id1= AS1[i1];
		if(id1!=ABut)
		{
			bool found = false;
			for (unsigned int i2=0; i2<size2 && !found; i2++){
				TCellID id2= AS2[i2];
				if(id1==id2){
					result.push_back(id1);
					found =true;
				}
			}
		}
	}
	return result;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID>
keepFilter(const std::vector<TCellID>& ASet, const TInt ANb)
{
	std::vector<TCellID> result;
	std::set<TCellID> filter;
	unsigned int s = ASet.size();

	for (unsigned int i=0; i<s; i++){
		TCellID current_id = ASet[i];
		int nb_occ=1;
		for (unsigned int j=i+1; j<s; j++){
			TCellID other_id = ASet[j];
			if(other_id==current_id)
				nb_occ++;
		}

		if(nb_occ>=ANb && std::find(result.begin(),result.end(),current_id)==result.end())
			filter.insert(current_id);
	}
	for(auto x:filter)
        result.push_back(x);
	return result;
}
/*----------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& stream, const MeshModel& model) {
	int m = model.getDef();
	stream << "|";
	if(m&N)
		stream<<"N|";
	if(m&E)
		stream<<"E|";
	if(m&F)
		stream<<"F|";
	if(m&R)
		stream<<"R|";
	if(m&N2E)
		stream<<"N2E|";
	if(m&N2F)
		stream<<"N2F|";
	if(m&N2R)
		stream<<"N2R|";
	if(m&E2N)
		stream<<"E2N|";
	if(m&E2F)
		stream<<"E2F|";
	if(m&E2R)
		stream<<"E2R|";
	if(m&F2N)
		stream<<"F2N|";
	if(m&F2E)
		stream<<"F2E|";
	if(m&F2F)
		stream<<"F2F|";
	if(m&F2R)
		stream<<"F2R|";
	if(m&R2N)
		stream<<"R2N|";
	if(m&R2E)
		stream<<"R2E|";
	if(m&R2F)
		stream<<"R2F|";
	if(m&R2R)
		stream<<"R2R|";

	return stream;
}
template<int N>
std::ostream & operator << (std::ostream & AStream, const TabCellID<N> & ATab) {
	for (int i = 0; i < ATab.size(); i++)
		AStream << ATab[i] << " ";
	return AStream;
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/

