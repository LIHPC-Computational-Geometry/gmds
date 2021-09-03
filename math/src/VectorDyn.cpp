/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * The GMDS library is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/*
 * VectorDyn.cpp
 *
 *  Created on: 15 juin 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gmds/utils/Exception.h>
#include <gmds/math/VectorDyn.h>
/*----------------------------------------------------------------------------*/
#include<cmath>
#include <string.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
VectorDyn::VectorDyn()
:size_(0)
{}
/*----------------------------------------------------------------------------*/
VectorDyn::VectorDyn(const TCoord& v1, const TCoord& v2)
{
	size_ = 2;
	tab_.resize(2);
	tab_[0]=v1;
	tab_[1]=v2;
}
/*----------------------------------------------------------------------------*/
VectorDyn::VectorDyn(const TCoord& v1, const TCoord& v2, const TCoord& v3)
{
	size_ = 3;
	tab_.resize(3);
	tab_[0]=v1;
	tab_[1]=v2;
	tab_[2]=v3;
}
/*----------------------------------------------------------------------------*/
VectorDyn::VectorDyn(const TInt& AN)
{
	size_ = AN;
	tab_.resize(AN);
}
/*----------------------------------------------------------------------------*/
VectorDyn::VectorDyn(const VectorDyn& AVec)
{
	size_ = AVec.size_;
	tab_.resize(size_);
	memcpy(&tab_[0],&AVec.tab_[0],size_*sizeof(TCoord));
}
/*----------------------------------------------------------------------------*/
VectorDyn::VectorDyn(const Vector& AVec)
{
	size_ = 3;
	tab_.resize(3);
	tab_[0] = AVec[0];
	tab_[1] = AVec[1];
	tab_[2] = AVec[2];
}
/*----------------------------------------------------------------------------*/
VectorDyn& VectorDyn::operator=(const VectorDyn& vec)
{
	if(vec==*this)
		return *this;

	size_ = vec.size_;
	tab_.resize(size_);
	memcpy(&tab_[0],&vec.tab_[0],size_*sizeof(TCoord));
	return *this;
}
/*----------------------------------------------------------------------------*/
VectorDyn::~VectorDyn()
{}
/*----------------------------------------------------------------------------*/
TInt VectorDyn::size() const
{
	return size_;
}
/*----------------------------------------------------------------------------*/
bool VectorDyn::
operator==(const VectorDyn& vec) const
{
	if(size_!=vec.size_)
		GMDSException("2 vectors with different dimension cannot be compared");

	bool result = true;
	for(unsigned int i=0;i<size_;i++)
		if (tab_[i] != vec.tab_[i])
			return false;

	return result;
}
/*----------------------------------------------------------------------------*/
TCoord VectorDyn::norm2() const
{
	TCoord s=0;
 	for(int i=0;i<size_;i++)
		s=s+tab_[i]*tab_[i];
	return s;
}
/*----------------------------------------------------------------------------*/
TCoord VectorDyn::norm() const
{
	return sqrt(norm2());
}
/*----------------------------------------------------------------------------*/
TCoord VectorDyn::normL2() const
{
	return norm();
}
/*----------------------------------------------------------------------------*/
TCoord VectorDyn::normL1() const
{
	TCoord s=0;
 	for(int i=0;i<size_;i++)
		s=s+fabs(tab_[i]);
	return s;
}
/*----------------------------------------------------------------------------*/
TCoord VectorDyn::normLinf() const
{
	TCoord s=fabs(tab_[0]);
 	for(int i=1;i<size_;i++)
		if(s<fabs(tab_[i]))
			s = fabs(tab_[i]);

	return s;
}
/*----------------------------------------------------------------------------*/
TCoord VectorDyn::normLp(const TInt& AP) const
{
	TCoord s=0;
 	for(int i=0;i<size_;i++)
		s=s+pow(tab_[i],AP);
 	double inv_p = 1.0/(double)AP;
	return pow(s,inv_p);
}
/*----------------------------------------------------------------------------*/
TCoord VectorDyn::sumComponents() const
{
	TCoord s=0;
 	for(int i=0;i<size_;i++)
		s=s+tab_[i];
	return s;
}

/*----------------------------------------------------------------------------*/
void VectorDyn::normalize()
{
	TCoord n = norm();
	if(n!=0.0){
		for(int i=0;i<size_;i++)
			tab_[i] = tab_[i]/n;
	}
}
/*----------------------------------------------------------------------------*/
void VectorDyn::set(const TCoord* ATab, const TInt ANb)
{

	size_ = ANb;
	tab_.resize(ANb);
	memcpy(&tab_[0],&ATab[0],size_*sizeof(TCoord));
}
/*----------------------------------------------------------------------------*/
TCoord VectorDyn::dot(const VectorDyn& vec) const
{
	if(size_!=vec.size_)
		GMDSException("A dot product cannot be done between 2 vectors with different dimensions");

	TCoord r=0.0;
	for(int i=0;i<size_;i++)
		r=r+tab_[i]*vec.tab_[i];

	return r;
}
/*----------------------------------------------------------------------------*/
VectorDyn VectorDyn::cross(const VectorDyn& vec) const
{
	if(size_!=vec.size_)
		GMDSException("A cross product cannot be done between 2 vectors with different dimensions");
	if(size_>3){
		throw GMDSException("Cross product not yet implemented for 4 and higher dimension vector");
	}
	VectorDyn v;
	if(size_==2){
		TCoord c[2] ={tab_[0] * vec[1], - tab_[1] * v[0]};
		v.set(c,2);
	}
	else if (size_==3){
		TCoord c[3]={tab_[1] * vec[2] - tab_[2] * vec[1],
				   -tab_[0] * vec[2] + tab_[2] * vec[0],
				    tab_[0] * vec[1] - tab_[1] * vec[0]};
		v.set(c,3);
	}
	return v;
}
/*----------------------------------------------------------------------------*/
bool VectorDyn::isColinear(const VectorDyn& vec) const
{
	if(size_!=vec.size_)
		GMDSException("We cannot check the colinearity of 2 vectors with different dimensions");

	VectorDyn v(*this), v2=vec;
	v.normalize();
	v2.normalize();
	TCoord result =v.dot(v2);
	return (result==1e+0);
}
/*----------------------------------------------------------------------------*/
bool VectorDyn::isOrthogonal(const VectorDyn& vec) const
{
	if(size_!=vec.size_)
		GMDSException("We cannot check the orthogonality of 2 vectors with different dimensions");
	return (dot(vec)==0.0);
}
/*------------------------------------------------------------------------*/
VectorDyn VectorDyn::operator^(const TInt n)
{
	VectorDyn R(size_);
	for( int i=0;i<size_;i++)
		R[i]=pow(tab_[i],n);

	return R;
}


/*----------------------------------------------------------------------------*/
VectorDyn operator-(const VectorDyn& a, const VectorDyn& b)
{
	if(a.size_!=b.size_)
		GMDSException("We cannot make a difference betwen 2 vectors with different dimensions");
	VectorDyn v(a.size_);
	for(int i=0;i<a.size_;i++)
		v.tab_[i] = a.tab_[i]-b.tab_[i];
	return v;
}
/*----------------------------------------------------------------------------*/
VectorDyn operator+(const VectorDyn& a, const VectorDyn& b)
{
	if(a.size_!=b.size_)
		GMDSException("We cannot sum 2 vectors with different dimensions");
	VectorDyn v(a.size_);
	for(int i=0;i<a.size_;i++)
		v.tab_[i] = a.tab_[i]+b.tab_[i];
	return v;
}
/*----------------------------------------------------------------------------*/
VectorDyn operator*(const VectorDyn& a, const TCoord& k)
{
	VectorDyn v(a.size_);
	for(int i=0;i<a.size_;i++)
		v.tab_[i] = k*a.tab_[i];
	return v;
}
/*----------------------------------------------------------------------------*/
VectorDyn operator*(const TCoord& k, const VectorDyn& a)
{
	VectorDyn v(a.size_);
	for(int i=0;i<a.size_;i++)
		v.tab_[i] = k*a.tab_[i];
	return v;
}
/*----------------------------------------------------------------------------*/
VectorDyn operator/(const VectorDyn& a, const  TCoord& k)
{
	VectorDyn v(a.size_);
	for(int i=0;i<a.size_;i++)
		v.tab_[i] = k/a.tab_[i];
	return v;
}
/*----------------------------------------------------------------------------*/
template<typename TCoord>
std::ostream& operator<<(std::ostream& stream, const VectorDyn& vec)
{
	stream<<"[";
	for(int i=0;i<vec.size()-1;i++)
		stream<<vec[i]<<", ";
	stream<<vec[vec.size()-1]<<"]";
	return stream;
}
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
