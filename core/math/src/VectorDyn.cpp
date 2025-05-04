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
		throw GMDSException("2 vectors with different dimensions cannot be compared");

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
		throw GMDSException("A dot product cannot be done between 2 vectors with different dimensions");

	TCoord r=0.0;
	for(int i=0;i<size_;i++)
		r=r+tab_[i]*vec.tab_[i];

	return r;
}
/*----------------------------------------------------------------------------*/
VectorDyn VectorDyn::cross(const VectorDyn& vec) const
{
	if(size_!=vec.size_)
		throw GMDSException("A cross product cannot be done between 2 vectors with different dimensions");
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
		throw GMDSException("We cannot check the colinearity of 2 vectors with different dimensions");

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
		throw GMDSException("We cannot check the orthogonality of 2 vectors with different dimensions");
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
		throw GMDSException("We cannot make a difference betwen 2 vectors with different dimensions");
	VectorDyn v(a.size_);
	for(int i=0;i<a.size_;i++)
		v.tab_[i] = a.tab_[i]-b.tab_[i];
	return v;
}
/*----------------------------------------------------------------------------*/
VectorDyn operator+(const VectorDyn& a, const VectorDyn& b)
{
	if(a.size_!=b.size_)
		throw GMDSException("We cannot sum 2 vectors with different dimensions");
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
