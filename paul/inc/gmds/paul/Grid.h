//
// Created by Paul Bourmaud on 23/03/2022.
//

#ifndef GMDS_GRID_H
#define GMDS_GRID_H
/*----------------------------------------------------------------------------*/
#include "gmds/ig/Mesh.h"
#include "gmds/igalgo/BoundaryExtractor2D.h"
#include "GMDSIgAlgo_export.h"
/*----------------------------------------------------------------------------*/

namespace gmds {

	class Grid
   {
 	public:
		Grid(const int AX = 3, const int AY = 3);
		int getX() const;
		int getY() const;

 	private:
		int m_X;
		int m_Y;
	};
   class GridBuilderAround
   {
	 public:
	   GridBuilderAround(Mesh* AMesh, const TInt ADim=2);
	   const Mesh getMesh() const;
	   int getDim() const;

	   virtual ~GridBuilderAround();

	   /*------------------------------------------------------------------------*/
	   /** @brief Check if the mesh fits algorithm requirements, which are:
         *         - 3D mesh with R, N and R2N
         *         - 2D mesh with F, N and F2N
	    */
	   bool isValid() const;

	   /** @brief  Performs the grid building
         * @param ANb      NB nodes in X  and Y directions
         */
	   void executeGrid2D(const TInt  ANb);

	   /* a mesh*/
	   gmds::Mesh m_mesh;


	   /* add variable activate -> Face*/
	   Variable<int> *activate = m_mesh.newVariable<int,GMDS_FACE>("activate");

	   /** Grid dimension*/
	   TInt m_dim;

	   /** @brief grid Build 2D around the target mesh
	    *
	    * @param AXNb number of nodes on X
	    * @param AXStep range on X
	    * @param AYNb number of nodes on Y
	    * @param AYStep range on Y
	    * @param Xmin  coord Xmin
	    * @param Ymin coord Ymin
	    */
   	void gridBuild2D(const TInt  AXNb, const TCoord AXStep,
      	          const TInt  AYNb, const TCoord AYStep, const TCoord Xmin,const TCoord Ymin
   	);

	   /** \brief get the value au the variable activate for a specify Face
	    *
	    * @param AFace the face
	    * @return the int value of activate (0 -> not activate, 1 -> activate
	    */
	   int getActivate(gmds::Face AFace);

	   /** \brief change the value of activate for a 0 (the false)
	    *
	    * @param faceID id face
	    */

	   void flipActivate(const int faceID);
   };
}
#endif     // GMDS_GRID_H
