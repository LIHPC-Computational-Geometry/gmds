//
// Created by Paul Bourmaud on 23/03/2022.
//

#ifndef GMDS_GRID_H
#define GMDS_GRID_H
/*----------------------------------------------------------------------------*/
#include "gmds/ig/Mesh.h"
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

	   /** @brief  Performs the grid building
         * @param AXNb      NB cells in X direction
         * @param AXStep    Length of each step in in X direction
         * @param AYNb      NB cells in Y direction
         * @param AYStep    Length of each step in in Y direction
         */
	   void executeGrid2D(const TInt  AXNb, const TCoord AXStep,
	                    const TInt  AYNb, const TCoord AYStep);

	   /** a mesh*/
	   gmds::Mesh m_mesh;

	   /** add variable activate -> Face*/
	   Variable<int> *activate = m_mesh.newVariable<int,GMDS_FACE>("exist");

	   /** Grid dimension*/
	   TInt m_dim;

   	void gridBuild2D(const TInt  AXNb, const TCoord AXStep,
      	          const TInt  AYNb, const TCoord AYStep
   	);

	   void flipActivate(const int faceID);
   };
}
#endif     // GMDS_GRID_H
