
/*----------------------------------------------------------------------------*/
/*
 * cross.cpp
 *
 *  Created on: Sep 05, 2014
 *      Author: franck Ledoux
 */
/*----------------------------------------------------------------------------*/
#include "gmds/math/Cross2D.h"
#include "gmds/math/Quaternion.h"
#include "gmds/math/Vector.h"
#include "gmds/math/Constants.h"
#include "gmds/math/Numerics.h"
/*-----------------------------------------------------------------*/
#include <iostream>
#include <math.h>
namespace gmds {
  /*-----------------------------------------------------------------*/
  namespace math{
    /*-----------------------------------------------------------------*/
    Cross2D::Cross2D() : m_angle(0),m_known_vectors(false){
        // m_vectors is then emppty
    }
      /*-----------------------------------------------------------------*/
      Cross2D::Cross2D(Vector3d& AV1, Vector3d& AV2)
      {
          if (!math::isZero(AV1.dot(AV2), 1e-10)) {
              std::cout<<AV1<<std::endl;
              std::cout<<AV2<<std::endl;

              std::cout<<AV1.dot(AV2)<<std::endl;
              throw GMDSException("A Cross2D object can only be built from 2 orthogonal vectors");
          }
          Vector3d v(1,0,0);
          TCoord a = v.orientedAngle(AV1);
      /* //a must be converted from degree to radians
          a = a*Constants::PIDIV180;
          if(a<0){
              a = Constants::PI2+a;
          }
        */
          m_angle = modulo2PI(4*a);
          m_known_vectors = false;
      }
      
      /*-----------------------------------------------------------------*/
      Cross2D::
      Cross2D(const Vector3d& ARef)
      {
          Vector3d v(1,0,0);
          TCoord a = v.orientedAngle(ARef);
          if(a<0){
              a = Constants::PI2+a;
          }
          
          m_angle = a;
          m_known_vectors = false;
      }
      /*-----------------------------------------------------------------*/
      Cross2D::
      Cross2D(const TCoord& ARefAngle)
      {
          
      m_angle = modulo2PI(ARefAngle);
    
      m_known_vectors = false;
    }
    /*-----------------------------------------------------------------*/
    Cross2D::Cross2D(const Cross2D& AC)
    {
      m_angle = AC.m_angle;
      m_known_vectors = AC.m_known_vectors;
      if(m_known_vectors)
	m_vectors = AC.m_vectors;
    }
    /*-----------------------------------------------------------------*/
    Vector3d Cross2D::closestComponentVector(const Vector3d& AN) const
    {
      std::vector<Vector3d> vecs = (m_known_vectors)?m_vectors:componentVectors();

        Vector3d n = AN;
      n.normalize();

        Vector3d v = vecs[0];
      TCoord val = n.dot(v);

      for(unsigned int i=1;i<4;i++){
          Vector3d v_i = vecs[i];
	TCoord val_i=n.dot(v_i);
	if(val_i> val){
	  val = val_i;
	  v = v_i;
	}
      }

      return v;
    } 
      /*-----------------------------------------------------------------*/
	 
	unsigned int Cross2D::closestComponentVectorAsIndex(const Vector3d& AN) const
    {
      	std::vector<Vector3d> vecs = (m_known_vectors)?m_vectors:componentVectors();

       	Vector3d n = AN;
      	n.normalize();

        	Vector3d v = vecs[0];
      	TCoord val = n.dot(v);
	 	unsigned int closestCompIdx = 0;

      	for(unsigned int i=1;i<4;i++){
         		Vector3d v_i = vecs[i];
			TCoord val_i=n.dot(v_i);
			if(val_i> val){
	 			val = val_i;
	  			closestCompIdx = i;
			}
      	}

      	return closestCompIdx;
    } 
    
    
    /*-----------------------------------------------------------------*/
    bool Cross2D::hasStoredComponentVectors() const{
      return m_known_vectors;
    }
      /*-----------------------------------------------------------------*/
      std::vector<Vector3d> Cross2D::componentVectors() const
      {
          std::vector<Vector3d> v;
          if(m_known_vectors)
              v = m_vectors;
          else{
              v.resize(4);
              TCoord a = m_angle/4.0;
              for(unsigned int i=0;i<4;i++){
                  TCoord a_i = a +i*Constants::PIDIV2;
                  v[i] = Vector3d(cos(a_i),sin(a_i),0.0);
              }
          }
          return v;
      }
      /*-----------------------------------------------------------------*/
    void Cross2D::computeComponentVectors() 
    {
      if(!m_known_vectors){
	m_known_vectors = true;
	m_vectors.resize(4);
	TCoord a = m_angle/4.0;
	for(unsigned int i=0;i<4;i++){
	  TCoord a_i = a +i*Constants::PIDIV2;
	  m_vectors[i] = Vector3d(cos(a_i),sin(a_i),0.0);
	}
      }
    }
    /*-----------------------------------------------------------------*/
    Vector3d Cross2D::referenceVector() const
    {
      TCoord a = m_angle;
      return Vector3d(cos(a),sin(a),0.0);
    }
    /*-----------------------------------------------------------------*/
    Cross2D  Cross2D::mean(const vector<Cross2D> & ACrosses, 
			   const vector<TCoord> & AWeights,
			   const TInt ANbSteps)
    {   
      if (ACrosses.empty()){
	throw GMDSException("The mean of zero 2D crosses is undefined.");
      }

      if (AWeights.size() != ACrosses.size())	{
	throw GMDSException("There must exactly the same number of crosses and weights.");
      }

      TCoord ref_angle =  ACrosses[0].referenceAngle(); 
      TCoord prev_ref_angle =  0;
        Vector3d ref_vector =  ACrosses[0].referenceVector();
      TInt step_index=0;
      do{
	step_index++;
	prev_ref_angle =  ref_angle;
	TCoord pen_angle=0;
	for (unsigned int i = 0; i < ACrosses.size(); i++)
	  {
          Vector3d vec_i = ACrosses[i].referenceVector();
	    const TCoord pen_i = ref_vector.orientedAngle(vec_i);
	    pen_angle += pen_i;//*AWeights[i];
	  }
	pen_angle /=ACrosses.size();
	
	ref_angle  = modulo2PI(ref_angle+pen_angle);
	ref_vector = Cross2D(ref_angle).referenceVector();
      }  
      while( (step_index<ANbSteps) && (ref_angle!=prev_ref_angle));

      TCoord new_angle = ref_angle;
      return Cross2D(new_angle);
      
    }
    
        /*-----------------------------------------------------------------*/
    Cross2D  Cross2D::meanNotMedian(const vector<Cross2D> & ACrosses, 
			   const vector<TCoord> & AWeights,
			   const TInt ANbSteps)
    {   
      if (ACrosses.empty()){
	throw GMDSException("The mean of zero 2D crosses is undefined.");
      }

      if (AWeights.size() != ACrosses.size())	{
	throw GMDSException("There must exactly the same number of crosses and weights.");
      }

      TCoord ref_angle =  ACrosses[0].referenceAngle(); 
      TCoord prev_ref_angle =  0.0;
        Vector3d ref_vector =  ACrosses[0].referenceVector();
      TInt step_index=0;
      do{
	step_index++;
	prev_ref_angle =  ref_angle;
	TCoord pen_angle=0.0;
	for (unsigned int i = 0; i < ACrosses.size(); i++)
	  {
          Vector3d vec_i = ACrosses[i].referenceVector();
	    const TCoord pen_i = ACrosses[i].referenceAngle(); 
	    pen_angle += pen_i*AWeights[i];
	  }
	
	
	ref_angle  = modulo2PI(pen_angle);
	ref_vector = Cross2D(ref_angle).referenceVector();
      }  
      while( (step_index<ANbSteps) && (ref_angle!=prev_ref_angle));

      TCoord new_angle = ref_angle;
	 //////////////////////////////////////
	  vector<math::Vector3d> c_vectors0 = ACrosses[0].componentVectors();
	 vector<math::Vector3d> c_vectors1 = ACrosses[1].componentVectors();
	 vector<math::Vector3d> c_vectors2 = ACrosses[2].componentVectors();
	 	
	 gmds::math::Vector3d first_vect = c_vectors0[0];
	 gmds::math::Vector3d second_closest_vect =  ACrosses[1].closestComponentVector(first_vect);
	 gmds::math::Vector3d third_closest_vect =  ACrosses[2].closestComponentVector(first_vect);
	 gmds::math::Vector3d center_vect = first_vect * AWeights[0] + second_closest_vect * AWeights[1] + third_closest_vect * AWeights[2];
	 gmds::math::Vector3d center_vect_second = center_vect.getOneOrtho();
     
      return Cross2D(center_vect,center_vect_second);
      
    }
    /*-----------------------------------------------------------------*/
    Cross2D  Cross2D::mean(const Cross2D& AC1, const TCoord& AW1,
			   const Cross2D& AC2, const TCoord& AW2)
    {
        Vector3d v1 = AC1.componentVectors()[0];
        Vector3d v2 = AC2.closestComponentVector(v1);
        Vector3d v = AW1*v1 + AW2*v2;
        Vector3d v_ortho = v.cross(math::Vector3d(0,0,1));

      return Cross2D(v,v_ortho);
    }

    /*-----------------------------------------------------------------*/
    TCoord Cross2D::angle(const Cross2D& AC) const
    {
        Vector3d v1 = referenceVector();
        Vector3d v2 = AC.referenceVector();

      return v1.angle(v2);
    }
    /*-----------------------------------------------------------------*/
    int Cross2D::index(const Cross2D& AC1,
		       const Cross2D& AC2,
		       const Cross2D& AC3)
    {
        Vector3d vi = AC1.referenceVector();
        Vector3d vj = AC2.referenceVector();
        Vector3d vk = AC3.referenceVector();

      
      double wij = vi.orientedAngle(vj);
      double wjk = vj.orientedAngle(vk);
      double wki = vk.orientedAngle(vi);
      
	 return std::round((wij+wjk+wki)/Constants::PI2);
     
    }
    /*-----------------------------------------------------------------*/
    ostream & operator << (ostream & AStr, const Cross2D & AC)
    {
      AStr << "Cross2D (" << AC.referenceAngle() << ")";
      return AStr;
    }
    /*-----------------------------------------------------------------*/
    Cross2D operator+(const Cross2D& AC1, const Cross2D& AC2){
      TCoord a1 =  AC1.referenceAngle();
        Vector3d v1 =  AC1.referenceVector();
        Vector3d v2 =  AC2.referenceVector();
      TCoord penalty =v1.orientedAngle(v2);
      TCoord new_angle = modulo2PI(a1 + penalty/2.0);      
      return Cross2D(new_angle);
    }
  }
  /*-----------------------------------------------------------------*/
}
/*-----------------------------------------------------------------*/
