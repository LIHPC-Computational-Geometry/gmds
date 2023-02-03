
/*-----------------------------------------------------------------*/
/*
 * cross.cpp
 *
 *  Created on: Sep 05, 2014
 *      Author: franck Ledoux
 */
/*-----------------------------------------------------------------*/
#include "gmds/math/Cross.h"
#include "gmds/math/Quaternion.h"
#include "gmds/math/Chart.h"
#include "gmds/math/Vector.h"
/*-----------------------------------------------------------------*/
#include <iostream>
#include <cmath>
namespace gmds {
/*-----------------------------------------------------------------*/
namespace math {
/*-----------------------------------------------------------------*/
Cross::Cross()
{
	m_x = math::Vector3d({1., 0., 0.});
	m_y = math::Vector3d({0., 1., 0.});
}

/*-----------------------------------------------------------------*/
Cross::Cross(math::Vector3d &AV1, math::Vector3d &AV2)
{
	if (AV1.dot(AV2) != 0.0) {
		throw GMDSException("A CROSS object can only be built from 2 orthogonal vectors");
	}
	m_x = AV1;
	m_y = AV2;
	m_x.normalize();
	m_y.normalize();
}

/*-----------------------------------------------------------------*/
Cross::Cross(const Quaternion &AQ, const math::Vector3d &AN)
{
	Quaternion q = AQ;
	// we rotate q to be aligned with AN
	q.alignWith(AN);
	Chart tq(q);
	if (fabs(AN.dot(tq.X())) > 0.5) {
		m_x = tq.Y();
		m_y = tq.Z();
	}
	else if (fabs(AN.dot(tq.Y())) > 0.5) {
		m_x = tq.X();
		m_y = tq.Z();
	}
	else {
		m_x = tq.X();
		m_y = tq.Y();
	}
	m_x.normalize();
	m_y.normalize();
	// std::cout << "normal: " << AN << std::endl;
	// std::cout << " -> Cross X: " << m_x << std::endl;
	// std::cout << " -> Cross Y: " << m_y << std::endl;
}
/*-----------------------------------------------------------------*/
math::Vector3d
Cross::closestVector(const math::Vector3d &AN)
{
	math::Vector3d v = AN;
	v.normalize();
	math::Vector3d x_opp({-m_x.X(), -m_x.Y(), -m_x.Z()});
	math::Vector3d y_opp({-m_y.X(), -m_y.Y(), -m_y.Z()});

	math::Vector3d result = m_x;
	double valX = v.dot(m_x);
	double valXOpp = v.dot(x_opp);
	double valY = v.dot(m_y);
	double valYOpp = v.dot(y_opp);

	double val = valX;
	if (valXOpp > val) {
		val = valXOpp;
		result = x_opp;
	}
	if (valY > val) {
		val = valY;
		result = m_y;
	}
	if (valYOpp > val) {
		result = y_opp;
	}

	return result;
}
/*-----------------------------------------------------------------*/
math::Chart
Cross::chart() const
{
	return {m_x,m_y,m_x.cross(m_y)};
}
/*-----------------------------------------------------------------*/

ostream &
operator<<(ostream &str, const Cross &c)
{
	str << "Cross (" << c.X() << ", " << c.Y() << ")";
	return str;
}
}
  /*-----------------------------------------------------------------*/
}
/*-----------------------------------------------------------------*/
