//
// Created by rochec on 01/12/23.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/TransfiniteInterpolation_3D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
bool TransfiniteInterpolation_3D::
   computeHex(Array3D<math::Point> &AG) {
	//We first check that the dimensions of AG are correct
	TInt nb_i = AG.nbElements(0);
	TInt nb_j = AG.nbElements(1);
	TInt nb_k = AG.nbElements(2);
	if (nb_i == 0 || nb_j == 0 || nb_k == 0)
	{
		return false;
	}

	double di = 1.0/(double(nb_i)-1.0) ;
	double dj = 1.0/(double(nb_j)-1.0) ;
	double dk = 1.0/(double(nb_k)-1.0) ;

	for (auto i=1;i<nb_i-1;i++)
	{
		for (auto j=1;j<nb_j-1;j++)
		{
			for (auto k=1;k<nb_k-1;k++)
			{
				math::Point U = (1.0-i*di)*AG(0,j,k) + i*di*AG(nb_i-1,j,k) ;
				math::Point V = (1.0-j*dj)*AG(i,0,k) + j*dj*AG(i,nb_j-1,k) ;
				math::Point W = (1.0-k*dk)*AG(i,j,0) + k*dk*AG(i,j,nb_k-1) ;

				math::Point UW = (1.0-i*di)*(1.0-k*dk)*AG(0,j,0) + (1-i*di)*k*dk*AG(0,j,nb_k-1)
				   + i*di*(1.0-k*dk)*AG(nb_i-1,j,0) + i*di*k*dk*AG(nb_i-1,j,nb_k-1);

				math::Point UV = (1.0-i*di)*(1.0-j*dj)*AG(0,0,k) + i*di*j*dj*AG(nb_i-1,nb_j-1,k)
				   + i*di*(1.0-j*dj)*AG(nb_i-1,0,k) + (1.0-i*di)*j*dj*AG(0,nb_j-1,k) ;

				math::Point VW = (1.0-j*dj)*(1.0-k*dk)*AG(i,0,0) + (1.0-j*dj)*k*dk*AG(i,0,nb_k-1)
					+ j*dj*(1.0-k*dk)*AG(i,nb_j-1,0) + j*dj*k*dk*AG(i,nb_j-1,nb_k-1) ;

				math::Point UVW = (1.0-i*di)*(1.0-j*dj)*(1.0-k*dk)*AG(0,0,0) + (1.0-i*di)*(1.0-j*dj)*k*dk*AG(0,0,nb_k-1)
				   + (1.0-i*di)*j*dj*(1.0-k*dk)*AG(0,nb_j-1,0) + i*di*(1.0-j*dj)*k*dk*AG(nb_i-1,0,nb_k-1)
				   + i*di*j*dj*(1.0-k*dk)*AG(nb_i-1,nb_j-1,0) + i*di*j*dj*k*dk*AG(nb_i-1,nb_j-1,nb_k-1)
				   + (1.0-i*di)*j*dj*k*dk*AG(0,nb_j-1,nb_k-1) + i*di*(1.0-j*dj)*(1.0-k*dk)*AG(nb_i-1,0,0);

				AG(i,j,k) = U + V + W + (-1.0)*UW + (-1.0)*UV + (-1.0)*VW + UVW ;
			}
		}
	}

	return true;

}
/*------------------------------------------------------------------------*/