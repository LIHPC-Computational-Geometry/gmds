//
// Created by calderans on 02/03/23.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MDLREADER_H
#define GMDS_MDLREADER_H

/*----------------------------------------------------------------------------*/
#include <fstream>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
	class MdlReader{
	 public:


	   /** @brief Constructor
         *
         * @param AMesh the mesh in which we want to copy the content of a Lima
			* file.
	    */
	   explicit MdlReader(Mesh& AMesh, std::string AString = "");

	   /*------------------------------------------------------------------------*/
	   /** \brief  Destructor.	*/
	   virtual ~MdlReader();

	   /*------------------------------------------------------------------------*/
	   /** \brief  Read the content of the file named outputName_ and write it in
     *   		mesh_.
	    */
	   void read(const std::string &AFileName);

	   /*------------------------------------------------------------------------*/
	   /** \brief Get the result mesh obtained after read
	    *
	    * @return a Mesh
	    */
	   Mesh *getMesh();

	   /*------------------------------------------------------------------------*/
	   /** \brief Create a mesh variable for each cells group in the mesh,
	    * mainly used in order to write the mesh
	    */
	   void createVariablesFromGroups();

	 private:
	   /**@brief Move the stream pointer onto the first occurrence of @AString
         *
         * @param AString the word we look for in the file
         *
         * @return true if we find @AString, false otherwise
	    */
	   bool moveStreamOntoSupport(std::string &AString);



	 private:
	   //the input mesh we work with
	   Mesh& m_mesh;
	   /* length unit */
	   double m_lengthUnit;
	   /** the file stream */
	   std::ifstream* m_stream;

	   std::string m_name2find;
   };
   /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif     // GMDS_MDLREADER_H
           /*----------------------------------------------------------------------------*/