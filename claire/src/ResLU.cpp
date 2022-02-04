//
// Created by rochec on 04/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/ResLU.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

ResLU::ResLU(Eigen::SparseMatrix<double> A_A, Eigen::VectorXd A_b) {
	m_A = A_A ;
	m_b = A_b;
	m_n = m_b.size();
	m_LU = A_A ;
	Eigen::VectorXd m_x(m_n);
}


/*------------------------------------------------------------------------*/
ResLU::STATUS ResLU::execute()
{
	computeLU();
	solveLU();

	return ResLU::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void ResLU::computeLU(){

	for(int i=0;i<m_n-1;i++){
		for(int j=i+1;j<m_n;j++){
			m_LU.coeffRef(j,i) = (m_LU.coeffRef(j,i) / m_LU.coeffRef(i,i)) ;
			for(int k=i+1;k<m_n;k++){
				m_LU.coeffRef(j,k) -= m_LU.coeffRef(i,k)*m_LU.coeffRef(j,i) ;
			}
		}
	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void ResLU::solveLU(){
	Eigen::SparseMatrix<double> A_1(m_n,m_n);

	for(int l=0;l<m_n;l++){
		A_1.coeffRef(l,l) = 1.0 ;
		for(int i=0;i<m_n;i++){
			for(int j=0;j<i;j++){
				A_1.coeffRef(i,l) -= m_LU.coeffRef(i,j)*A_1.coeffRef(j,l) ;
			}
		}
		for(int i=m_n-1; i>=0;i--){
			for (int j=m_n-1;j>i;j--){
				A_1.coeffRef(i,l) -= m_LU.coeffRef(i,j)*A_1.coeffRef(j,l);
			}
			A_1.coeffRef(i,l) /= m_LU.coeffRef(i,i) ;
		}
	}

	m_x = A_1*m_b;

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
Eigen::VectorXd ResLU::getSolution(){
	return m_x;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
Eigen::SparseMatrix<double> ResLU::getLU(){
	return m_LU;
}
/*------------------------------------------------------------------------*/