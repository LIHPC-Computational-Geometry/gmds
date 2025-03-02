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

// Check if the first point is on the segment formed by the two last points.
bool isOnSegment(math::Point AP0, math::Point AP1, math::Point AP2);

// Insert a given point in the given set of points.
std::vector<TCellID> insertPoint(Node AN, std::vector<Node> AV);

// Merge the twp given vectors of points.
std::vector<math::Point> merge(std::vector<math::Point> AV1, std::vector<math::Point> AV2);

// Takes a vector of edges and a vector of the same size containing quantities assotiated to each edge. Returns the 
// vector of the edges sorted with respect to this quantity.
std::vector<Edge> order(std::vector<Edge> AVE, std::vector<double> AVX);

// Returns, if it exists, the edge corresponding to the two given nodes.
Edge getEdge(Node &AN1, Node &AN2);

// Returns true if the node is internal.
bool isInterior(Node &AN);

// Returns the area delimited by the edges containted in the input vector. These edges must sourround a convex region.
double delimitedArea(std::vector<Edge> AV);

// Returns true if the edge touches the boundary.
bool touchesBoundary(Edge &AE);

// Returns the common face of the two input edges.
Face getCommonFace(Edge &AE1, Edge &AE2);

// Returns the projection of the second input vector on the first input vector.
math::Vector projection(math::Vector &AV1, math::Vector &AV2);

// Takes as input two nodes of the input minimal triangulation, and returns the nodes forming the shortest path to go from the first to the second.
std::vector<Node> shortestPathAlongBoundaryOrConstraints(Node &AN1, Node &AN2, Mesh &AMesh);

// Takes a quad and an edge belonging to the face, and return the opposite face
Edge opp(Face AFace, Edge AEdge);


#endif     // GMDS_MEDIALAXISMATH_H
