/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    FracPres.h
 *  \author  legoff
 *  \date    20/11/2017
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_FRACPRES_H_
#define ELG3D_FRACPRES_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <map>
#include <string>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <KM/Utils/KTypes.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/
const double FracPres_TRACE_AMOUNT = 10e-10;

/*----------------------------------------------------------------------------*/
/** \class FracPres
 *  \brief This is a dummy class.
 */
/*----------------------------------------------------------------------------*/
    class FracPres
    {
    public:
        /*------------------------------------------------------------------------*/
        /** \brief  Constructor
         */
        FracPres();

        /*------------------------------------------------------------------------*/
        /** \brief Copy constructor
         */
        FracPres(const FracPres&) =delete;

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor
         */
        ~FracPres() =default;

        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator=
         */
        FracPres& operator=(const FracPres&);

        /*------------------------------------------------------------------------*/
        /** \brief  Empties the frac pres container
         *
         */
        void clear();

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the number of materials
         *
         *  \return the number of materials
         */
        int getNbMaterials() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the id of a material
         *
         *  \param[in]  AName the name of the material
         *
         *  \return the id of the material
         */
        int getMaterialID(std::string AName) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the name of a material
         *
         *  \param[in]  AId the id of the material
         *
         *  \return the name of the material
         */
        std::string getMaterialName(int AId) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the names of the materials
         *
         *  \return the name of the material
         */
        std::map<int, std::string> getMaterialList() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Sets the materials, from a list
         *
         *  \param[in] the list of the materials (map between their id and name)
         */
         void setMaterialList(const std::map<int, std::string> AMatList);

        /*------------------------------------------------------------------------*/
        /** \brief  Checks whether a material exists
         *
         *  \param[in]  AName the name of the material
         *
         *  \return true if the material exists, false otherwise
         */
        bool materialExists(std::string AName) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Adds a material
         *
         *  \param[in]  AName the name of the material
         *
         *  \return the id of the added material
         */
        int createMaterial(std::string AName);

        /*------------------------------------------------------------------------*/
        /** \brief  Sets the fracpres of a material for a cell
         *
         *  \param[in]  AMat the id of the material
         *  \param[in]  AId the cell id
         *  \param[in]  AFracPres the fracpres to set
         */
// TODO thread safety
        void setFracPres(int AMat, kmds::TCellID AId, double AFracPres);

        /*------------------------------------------------------------------------*/
        /** \brief  Gets the fracpres of a material for a cell
         *
         *  \param[in]  AMat the id of the material
         *  \param[in]  AId the cell id
         *
         *  \return  the fracpres of material AMat in cell AId
         */
        double getFracPres(int AMat, kmds::TCellID AId) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns whether the cell is a mixed cell or not.
         *
         *  \param[in]  AId the cell id
         *
         *  \return  whether the cell is mixed (meaning it contains several materials)
         */
        bool isMixedCell(kmds::TCellID AId) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Gets the majority material for a cell
         *
         *  \param[in]  AId the cell id
         *
         *  \return  the material with majority fracpres in cell AId
         */
        int getMaxMatFracPresIndex(kmds::TCellID AId) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Gets the majority material frac pres for a cell
         *
         *  \param[in]  AId the cell id
         *
         *  \return  the frac pres of the material with majority fracpres in cell AId
         */
        double getMaxMatFracPresValue(kmds::TCellID AId) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Gets the fracpres of a material, for at least the cells that contain some of it. Some of the cells can carry a fracpres of zero.
         *
         *  \param[in]  AMat the id of the material
         *
         *  \return  a mapping between the cells IDs and the fracpres for AMat material
         */
        std::map<kmds::TCellID, double> getMatFracPres(int AMat) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Gets the fracpres of a material, for at least the cells that contain some of it. Some of the cells can carry a fracpres of zero.
         *
         *  \param[in]  AMat the id of the material
         *  \param[in]  AMatFracPres the frac pres for a material
         *
         */
        void setMatFracPres(int AMat, const std::map<kmds::TCellID, double> AMatFracPres);

        /*------------------------------------------------------------------------*/
        /** \brief  Normalize the frac pres for a cell
         *
         *  \param[in]  AId the cell id
         *
         */
        void normalize(kmds::TCellID AId);

        /*------------------------------------------------------------------------*/
        /** \brief  Normalize the frac pres for all cells if needed
         *
         *  \return  Whether the fp was modified for at least one cell
         *
         */
        bool normalize();

        /*------------------------------------------------------------------------*/
        /** \brief  Check whether the frac pres is normalized for all cells
         *
         *  \return  Whether the fp if normalized
         *
         */
        bool checkNormalized() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Remove traces of materials in cells. Should produce more pure cells
         *
         *  \return  Whether the fp was modified for at least one cell
         *
         */
        bool removeTraceAmounts();


    private:

        /* mapping between material ids and their names */
        std::map<int, std::string> m_matId2Name;

        /* mapping between material names and their ids */
        std::map<std::string, int> m_matName2Id;

        /* mapping between cells and their fracpres, for each material */
        std::map<int, std::map<kmds::TCellID, double> > m_mat2CellFracPres;

    };
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_FRACPRES_H_ */
