//
// Created by chenyt on 19/04/24.
//

#ifndef GMDS_MEDIALAXISMATH_H
#define GMDS_MEDIALAXISMATH_H

#include <gmds/math/Triangle.h>
#include <gmds/math/Point.h>
#include <gmds/math/Tetrahedron.h>
#include <gmds/math/Matrix.h>
#include <gmds/math/Vector.h>
#include <gmds/ig/Mesh.h>
using namespace gmds;

// Computes the maximum distance to a sphere 1 of a point of a sphere 2 which doesn't belong to the sphere 1
double maxRange(const math::Point& ACenter1, const double& ARadius1, const math::Point& ACenter2, const double& ARadius2);

// Compute the oriented angle between two 2D vectors, the result is in [-pi,pi[;
double oriented_angle(const math::Vector &AV1, const math::Vector &AV2);

// Return the vector corresponding to the given edge with first node the given node
math::Vector edge2vec(const gmds::Edge &AE, const gmds::Node &AN);

// Return the common node, if it exists, of twe edges
Node getCommonNode(const Edge &AE1, const Edge &AE2);

// Takes a node and a list of edges containing this node. Return a list of the sorted edges, minimizing the
// rotation between two consecutive edges
std::vector<Edge> sortEdges(Node &AN, std::vector<Edge> &AV);

// Useful functions for oriented graphs represented by a mesh
std::vector<Node> getNextNodes(Node &AN);
void propagateValue(Node &AN, Variable<double> &AVar);
bool isASource(Node &AN);
bool isAWell(Node &AN);
int NbWells(Mesh* AMesh);

// Orientate edges of a face in the direct sense
std::vector<Edge> orientateEdges(Face &AF);

// Regroup aligned conformal edges of a face
std::vector<std::vector<Edge>> groupsOfAlignedEdges(Face &AF);
#endif     // GMDS_MEDIALAXISMATH_H
