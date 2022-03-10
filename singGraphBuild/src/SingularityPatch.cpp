/*----------------------------------------------------------------------------*/
/*
 * SingularityPatch.cpp
 *
 *  Created on: June 4, 2015
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#include <gmds/singGraphBuild/SingularityLine.h>
#include <gmds/singGraphBuild/SingularityPatch.h>
#include <gmds/singGraphBuild/SingularityPoint.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
SingularityPatch::SingularityPatch() : m_valid_flag(true) {}
/*----------------------------------------------------------------------------*/
SingularityPatch::~SingularityPatch() {}
/*----------------------------------------------------------------------------*/
void
SingularityPatch::addPoint(SingularityPoint *APoint)
{
	if (APoint == 0) throw GMDSException("SingularityPatch - Try to add a null point");
	m_points.push_back(APoint);
	m_valid_flag = false;
}
/*----------------------------------------------------------------------------*/
void
SingularityPatch::addLine(SingularityLine *ALine)
{
	if (ALine == 0) throw GMDSException("SingularityPatch - Try to add a null line");
	m_lines.push_back(ALine);
	m_valid_flag = false;
}
/*----------------------------------------------------------------------------*/
bool
SingularityPatch::getPoints(std::vector<SingularityPoint *> &APoints)
{
	if (!m_valid_flag) m_valid_flag = checkValidityAnReorder();

	APoints = m_points;

	return m_valid_flag;
}
/*----------------------------------------------------------------------------*/
bool
SingularityPatch::getLines(std::vector<SingularityLine *> &ALines)
{
	if (!m_valid_flag) m_valid_flag = checkValidityAnReorder();

	ALines = m_lines;

	return m_valid_flag;
}
/*----------------------------------------------------------------------------*/
bool
SingularityPatch::checkValidityAnReorder()
{
	//==========================================================================
	// We check that every point is connected to  one or two lines in
	// the patch
	//==========================================================================
	for (unsigned int i = 0; i < m_points.size(); i++) {

		SingularityPoint *pi = m_points[i];
		std::vector<SingularityLine *> pi_lines = pi->getLines();
		int nb_connected = 0;

		for (unsigned int j = 0; j < pi_lines.size(); j++) {

			SingularityLine *lj = pi_lines[j];

			for (unsigned int k = 0; k < m_lines.size(); k++) {
				SingularityLine *lk = m_lines[k];
				if (lk == lj) nb_connected++;
			}

		}     // for(unsigned int j=0; j < pi_lines.size(); j++)

		if (nb_connected < 1 || nb_connected > 2) return false;

	}     // for(unsigned int i=0; i < m_points.size(); i++)

	//==========================================================================
	// We check that every line is connected to points available in the patch
	//==========================================================================
	for (unsigned int i = 0; i < m_lines.size(); i++) {

		SingularityLine *li = m_lines[i];
		std::vector<SingularityPoint *> li_points = li->getEndPoints();

		for (unsigned int j = 0; j < li_points.size(); j++) {

			SingularityPoint *pj = li_points[j];

			bool found_pj = false;
			for (unsigned int k = 0; k < m_points.size(); k++) {
				SingularityPoint *pk = m_points[k];
				if (pk == pj) found_pj = true;
			}
			// Point pj is missing in the patch
			if (!found_pj) return false;

		}     // for(unsigned int j=0; j<li_points.size();j++)

	}     // for(unsigned int i=0; i < m_lines.size(); i++)

	//==========================================================================
	// We check that lines and points form a single loop and we reorder on the
	// fly
	//==========================================================================
	std::vector<SingularityPoint *> new_points;
	std::vector<SingularityLine *> new_lines;
	SingularityPoint *first_pnt = m_points[0];
	SingularityPoint *current_pnt = first_pnt;
	std::vector<SingularityLine *> adj_lines = current_pnt->getLines();
	std::vector<SingularityLine *> candidate_lines;
	for (unsigned int i = 0; i < adj_lines.size(); i++)
		for (unsigned int j = 0; j < m_lines.size(); j++)
			if (adj_lines[i] == m_lines[j]) candidate_lines.push_back(adj_lines[i]);

	SingularityLine *current_line = candidate_lines[0];

	do {
		new_points.push_back(current_pnt);
		new_lines.push_back(current_line);

		std::vector<SingularityPoint *> end_line_points = current_line->getEndPoints();

		if (current_pnt == end_line_points[0])
			current_pnt = end_line_points[1];
		else
			current_pnt = end_line_points[0];

		std::vector<SingularityLine *> adj_lines = current_pnt->getLines();
		std::vector<SingularityLine *> candidate_lines;

		for (unsigned int i = 0; i < adj_lines.size(); i++)
			for (unsigned int j = 0; j < m_lines.size(); j++)
				if (adj_lines[i] == m_lines[j]) candidate_lines.push_back(adj_lines[i]);

		if (current_line == candidate_lines[0])
			current_line = candidate_lines[1];
		else
			current_line = candidate_lines[0];

	} while (current_pnt != first_pnt);

	m_points = new_points;
	m_lines = new_lines;
	return true;
}
/*----------------------------------------------------------------------------*/
SingularityLine *
SingularityPatch::computeOpposedLine(SingularityLine *Aline)
{
	for (int i = 0; i < m_lines.size(); i++) {
		if (m_lines[i] == Aline) {
			unsigned int opposedLineID = (i + 2) % m_lines.size();
			return m_lines[opposedLineID];
		}
	}
	throw GMDSException("unable to find the requested line inside the patch");
}
