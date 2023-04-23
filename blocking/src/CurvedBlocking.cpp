/*----------------------------------------------------------------------------*/
#include "gmds/blocking/CurvedBlocking.h"
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::blocking;
/*----------------------------------------------------------------------------*/
CurvedBlocking::CurvedBlocking()
{}
/*----------------------------------------------------------------------------*/
CurvedBlocking::~CurvedBlocking()
{}
/*----------------------------------------------------------------------------*/
DartHandler CurvedBlocking::createHex()
{
	DartHandler dh = m_gmap.make_combinatorial_hexahedron();
	// Initialize attribute for the created hexahedron
	ClassificationInfo ci;
	ci.id=1;
	ci.dim=2;
	for (auto it=m_gmap.one_dart_per_incident_cell<0,3>(dh).begin(),
	          itend=m_gmap.one_dart_per_incident_cell<0,3>(dh).end(); it!=itend; ++it)
	{m_gmap.set_attribute<0>(it, m_gmap.create_attribute<0>(ci));}
	for (auto it=m_gmap.one_dart_per_incident_cell<1,3>(dh).begin(),
	          itend=m_gmap.one_dart_per_incident_cell<1,3>(dh).end(); it!=itend; ++it)
	{m_gmap.set_attribute<1>(it, m_gmap.create_attribute<1>(ci));}
	for (auto it=m_gmap.one_dart_per_incident_cell<2,3>(dh).begin(),
	          itend=m_gmap.one_dart_per_incident_cell<2,3>(dh).end(); it!=itend; ++it)
	{m_gmap.set_attribute<2>(it, m_gmap.create_attribute<2>(ci));}
	for (auto it=m_gmap.one_dart_per_incident_cell<3,3>(dh).begin(),
	          itend=m_gmap.one_dart_per_incident_cell<3,3>(dh).end(); it!=itend; ++it)
	{m_gmap.set_attribute<3>(it, m_gmap.create_attribute<3>(ci));}

	return dh;
}
/*----------------------------------------------------------------------------*/
std::string CurvedBlocking::info() const{
	std::ostringstream mess;
	mess<<"Blocking Info: "<<std::endl;
	m_gmap.display_characteristics(mess);
	mess<<", validity="<<(isValidTopology()?"true":"false");
	mess<<"\nOd-attributes: ";
	for (auto it=m_gmap.attributes<0>().begin(),
	          itend=m_gmap.attributes<0>().end(); it!=itend; ++it)
	{ mess<<m_gmap.info_of_attribute<0>(it).id<<"; "; }
	mess<<"\n1d-attributes: ";
	for (auto it=m_gmap.attributes<1>().begin(),
	          itend=m_gmap.attributes<1>().end(); it!=itend; ++it)
	{ mess<<m_gmap.info_of_attribute<1>(it).id<<"; "; }
	mess<<"\n2d-attributes: ";
	for (auto it=m_gmap.attributes<2>().begin(),
	          itend=m_gmap.attributes<2>().end(); it!=itend; ++it)
	{ mess<<m_gmap.info_of_attribute<2>(it).id<<"; "; }
	mess<<"\n3d-attributes: ";
	for (auto it=m_gmap.attributes<3>().begin(),
	          itend=m_gmap.attributes<3>().end(); it!=itend; ++it)
	{ mess<<m_gmap.info_of_attribute<3>(it).id<<"; "; }
	return mess.str();
}
/*----------------------------------------------------------------------------*/
