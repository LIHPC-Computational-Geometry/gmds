/*----------------------------------------------------------------------------*/
/*
 * SingGraphBuilder2DShortestPath.h
 *
 *  Created on: April 10 2014
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#pragma once
/*----------------------------------------------------------------------------*/
#include <unordered_map>
/*----------------------------------------------------------------------------*/
#include <glpk.h>
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_SINGGRAPHBUILD_export.h"
#include <gmds/singGraphBuild/SingularityGraph.h>
#include <gmds/singGraphBuild/SingularityGraphBuilder2D.h>
/*----------------------------------------------------------------------------*/

namespace gmds {

class IllegalLineCrossingFinder;
class LineIntersectionsDetector;

class LIB_GMDS_SINGGRAPHBUILD_API SingGraphBuilder2DShortestPath : public SingularityGraphBuilder2D
{
 public:
	SingGraphBuilder2DShortestPath(Mesh *AMesh, Variable<math::Cross2D> *AField, const bool ABuildGeomSing = true);

	/*  Set robustness level, range from 0 to 3. The less the level, the less time it takes to compute, but also the less guarantee about the final mesh quality.
	   If solution is not found or not optimal, consider upgrading this level.
	   Warning: with a high level and a complex case, glpk time limit can be reached easily and exeption will be trown. */
	void setRobustness(int level);

	/*  Same as setRobustness(int level), but with more complex inputs. Do not use if not familiar. */
	void setRobustness(unsigned int nBoundaryTarget, unsigned int nSurfaceTarget, unsigned int glpkTimeLimit = 30000);

	/*  set the time limit in millisecond to solve the system with glpk. This is NOT the max limit of the entire algorithm. */
	void setGLPKTimeLimit(int glpkTimeLimit);

 private:
	bool m_visualizeCost = false;
	bool m_visualizeShortestPaths = false;

	// @ { Both of those parameters allows to reduce the number of potential path found in dijkstra algo, and thus reduce the difficulty for the GLPK solver.
	//     GLPK will struggle due to a very high number of constraints between potential paths (illegal line crossing)
	//      -> the less path, the less constrain, the less time
	//		-> with too low parameters, it may however hinder GLPK from finding the best solution, or even only one.
	unsigned int m_maxNBoundaryTarget = 6;         // dijkstra algorithm will keep the the N best paths ending on boundaries
	unsigned int m_maxNValidSurfaceTarget = 5;     // dijkstra algorithm will stop after having found N potential VALID surfaces target (valid = without penalty)
	// @ }
	unsigned int m_glpkTimeLimit = 30000;     // time limit to solve with glpk, in millisecond

	unsigned int m_singPointNo;                // number of singularity points in the initial graph
	unsigned int m_totalNumberOfSlots;         // = 5 * m_singPointNo
	unsigned int m_totalNumberOfVariables;     // = m_totalNumberOfSlots + m_maxNBoundaryTarget

	vector<math::Vector3d> m_bdry_edge_normals;
	vector<vector<Face>> m_face2Face_neighbours_by_verts;     // face neighbours (the adjacent faces to its 3 nodes )
	vector<bool> m_is_bdry_face;                              // vector of boolean values indicating whether a triangle is situated along a boundary
	vector<bool> m_singOrGeomFaces;                           // indicates if a face is in a singularity

	vector<vector<vector<TCellID>>> m_finalPaths;      // all cells indexes by which lines will be discretize
	vector<SingularityPoint::Slot *> m_targets;        // all slots in a flat vector
	vector<vector<double>> m_distances;                // costs of the detected paths
	vector<TCellID> m_slotFaces;                       // initial face index of each line

	// @ { aliases to specify slots IDs
	using SourceID = unsigned int;
	using TargetID = unsigned int;
	using Slot = SingularityPoint::Slot;
	// @ }

	std::vector<std::pair<SourceID, TargetID>> m_solutions;     // pair of slot connected to the final singularity lines computed in the solver

	struct BoundaryPathEndParam     // data about the last discret point of a shortest path ending on a boundary
	{
		unsigned int dim = 0;
		TCellID cellId = 0;
		math::Point point;
		TCellID finalFace = 0;
	};

	std::unordered_map<SourceID, std::vector<BoundaryPathEndParam>> m_bdryPathEndParam;     // data about all shortest path to ending on a boundary

	std::vector<std::pair<unsigned int, unsigned int>>
	   m_illegalOverlappingPaths;     // pairs of mutually excluding shortest paths (that shouldn't co-exist in the final graph)

	/*----------------------------------------------------------------------------------------------------*/
	/*----------------------    Graph building functions     ---------------------------------------------*/
	/*----------------------------------------------------------------------------------------------------*/
	void createSingularityLines() override;
	void createSolvedSingularityLinesInGraph();

	/*----------------------------------------------------------------------------------------------------*/
	/*----------------------    initialization/util functions     ----------------------------------------*/
	/*----------------------------------------------------------------------------------------------------*/
	void computeFaceNeighboursInfo();
	void initializeFieldsValue();
	bool targetIsBoundary(const TargetID &contTarget) const;
	bool findBoundary(const math::Ray &ray, const TCellID currentFace, BoundaryPathEndParam &bdryParam);
	void computeDijkstraStartingFaces();

	/*----------------------------------------------------------------------------------------------------*/
	/*----------------------    Shortest path system functions     ---------------------------------------*/
	/*----------------------------------------------------------------------------------------------------*/
	void computeDijkstra();

	/** \brief get the shortest path between a singular triangle and all other singular triangles
	   \param[in] targets : all the other singular triangles
	   \param[in] contSource : value indicating the "order" of the source -> contSource = 5*i+j, where i represents
	                     the singularity point number and j represents the number of the slot (j<=5)
	*/
	void getShortestPathBtwFacesOptimized(const vector<TCellID> &targets, const SourceID &contSource);

	std::vector<std::pair<SourceID, TargetID>> glpkSolve();

	/*----------------------------------------------------------------------------------------------------*/
	/*----------------------    illegal crossing functions     -------------------------------------------*/
	/*----------------------------------------------------------------------------------------------------*/
	void computeIllegalOverlappingPaths();

	// Build pairs of mutually excluding shortest paths, that shouldn't co-exist in the final graph.
	class IllegalLineCrossingFinder
	{
	 public:
		IllegalLineCrossingFinder(SingGraphBuilder2DShortestPath *graphBuilder) : m_graphBuilder(graphBuilder) {};
		std::vector<std::pair<unsigned int, unsigned int>> getIllegalOverlappingPaths();
		void registerAllLineSegmentsOnCells();

	 private:
		SingGraphBuilder2DShortestPath *m_graphBuilder;
		struct LineSegmentOnCell     // data stored in each cell crossed by a line
		{
			const unsigned int contSource;
			const unsigned int contTarget;
			const math::Segment segment;
			LineSegmentOnCell(unsigned int contSource, unsigned int contTarget, math::Segment segment) :
			  contSource(contSource), contTarget(contTarget), segment(segment)
			{
			}
		};
		std::unordered_map<TCellID, vector<LineSegmentOnCell>> m_traversedFacesByLine;     // store line segments on each cells
		void registerOneLineSegment(const math::Point &startPnt,
		                            const math::Point &endPnt,
		                            const TCellID &startFace,
		                            const TCellID &endFace,
		                            const SourceID &contSource,
		                            const TargetID &contTarget);
		void registerLineSegmentOnCommonNode(
		   const TCellID &f1Id, const TCellID &f2Id, const math::Point &p1, const math::Point &p2, const SourceID &contSource, const TargetID &contTarget);
	};
	/*----------------------------------------------------------------------------------------------------*/
	/*----------------------    line intersections functions     -----------------------------------------*/
	/*----------------------------------------------------------------------------------------------------*/

	void detectLineIntersections() override;

	class LineIntersectionDetector
	{
	 public:
		LineIntersectionDetector(SingGraphBuilder2DShortestPath *graphBuilder) : m_graphBuilder(graphBuilder), m_nFaces(graphBuilder->m_original_faces_number) {};
		void registerAllSolutionsOnCells();
		void detectAndCreateIntersections();

	 private:
		SingGraphBuilder2DShortestPath *m_graphBuilder;
		unsigned int m_nFaces;
		vector<bool>
		   m_singOrGeomFacesNeighbors;     // indicates if a face is near a singularity (use m_face2Face_neighbours_by_verts, but F2F connection should be better)
		struct TraversedCell               // singularity line data stored in each cell traversed by it
		{
			unsigned int prevPnt;
			SurfaceSingularityLine *singLine = nullptr;
			TraversedCell(unsigned int prevPnt, SurfaceSingularityLine *singLine) : prevPnt(prevPnt), singLine(singLine) {}
		};
		std::unordered_map<TCellID, vector<TraversedCell>> m_traversedCellsByLine;

		void registerOneSingularityLine(
		   SurfaceSingularityLine *singline, const TCellID &startFace, const TCellID &endFace, const vector<TCellID> &facePath, const bool isBoundary);
		void registerLineBetweenTwoFaces(SurfaceSingularityLine *singLine, const TCellID &f1Id, const TCellID &f2Id, const unsigned int &prevDiscId);
		void registerLineOnNodesNearSlots(SurfaceSingularityLine *singLine,
		                                  const TCellID &slotFace,
		                                  const TCellID &finalPathFace,
		                                  const unsigned int &slotDiscPointID,
		                                  const unsigned int &singDiscPointID,
		                                  const bool isBoundary);
		void registerLineOnFacesNearSlots(SurfaceSingularityLine *singLine,
		                                  const unsigned int &slotDiscPointID,
		                                  const unsigned int &singDiscPointID,
		                                  const bool isBoundary);
	};
};

}     // namespace gmds
