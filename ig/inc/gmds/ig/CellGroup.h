/*----------------------------------------------------------------------------*/
/** \file    CellGroup.h
 *  \author  N. Le Goff
 *  \date    16/09/2013
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_CELLGROUP_H_
#define GMDS_CELLGROUP_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <utility>
#include <vector>
#include <string>
#include <cmath>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/Cell.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
    /** \class CellGroup
     *  \brief Define a group of cells having the same dimension (nodes, edges,
     *  	   faces, regions). These groups are used as containers for node clouds,
     *  	   lines, surfaces by the mesh class.
     *
     *   	   WARNING: if a cell contained in the cell group is removed by the owner
     *  	            mesh class, the cell group keeps pointing on it (
     *  	            can induce an invalid pointer).
     */
	template<typename TCellType>
	class CellGroup{
	public:

		/*--------------------------------------------------------------------*/
		/** \brief  Constructor.
         */
		explicit CellGroup(Mesh* AMesh, std::string  AName="")
				: m_mesh(AMesh), m_name(std::move(AName)){}

		/*------------------------------------------------------------------------*/
		/** \brief  Copy constructor.
         */
		CellGroup(const CellGroup<TCellType>& ACGroup)
				: m_mesh(ACGroup.m_mesh), m_name(ACGroup.m_name), m_cells(ACGroup.m_cells){}


		bool operator==(const CellGroup<TCellType>& ACGroup)
		{
			return  (m_mesh=ACGroup.m_mesh) &&
					(m_name==ACGroup.m_name) &&
					(m_cells==ACGroup.m_cells);
		}

		/*------------------------------------------------------------------------*/
		/** \brief  Destructor.
         */
		virtual ~CellGroup(){}

		/*------------------------------------------------------------------------*/
		/** \brief  Add a cell into the group. Warning, A cell can be twice inside
         * 			a group.
         *
         *  \param ACell the cell to add
         */
		void add(const TCellType& ACell){m_cells.push_back(ACell.id());}
		/*------------------------------------------------------------------------*/
		/** \brief  Add a cell into the group using the ID
		 * 			WARNING, A cell can be twice inside a group.
         *
         *  \param AId id of the cell to add
         */
		void add(const TCellID AId){m_cells.push_back(AId);}


		/*------------------------------------------------------------------------*/
		/** \brief  Del a cell from the group by his id. Warning, A cell can be
		 *          twice inside a group.
         *
         *  \param AId if of the cell to add
         */
		void del(const TCellID AId){
			const unsigned int size = m_cells.size();
			bool notfound = true;
			unsigned int k=0;
			for(unsigned int i=0; i<size && notfound;i++)
			{
				if(m_cells[i]==AId)
				{
					notfound = false;
					k=i;
				}
			}
			if(!notfound){
				m_cells[k]=m_cells[size-1];
				m_cells.pop_back();
			}
		}

		/*------------------------------------------------------------------------*/
		/** \brief  Del a cell from the group. Warning, A cell can be twice inside
         * 			a group.
         *
         *  \param ACell the cell to add
         */
		void del(const TCellType& ACell) { del(ACell.id()); }

		/*------------------------------------------------------------------------*/
		/** \brief  Indicates if ACell belongs to this.
         *
         *  \param ACell the cell to find
         */
		bool has(TCellType ACell) const
		{
			const unsigned int size = m_cells.size();
			for(unsigned int i=0; i<size;i++)
				if(m_cells[i]==ACell.id())
					return true;

			return false;
		}


		/*------------------------------------------------------------------------*/
		/** \brief Indicates if the group is empty.
         *
         *  \return a boolean whose value is true if the the group is empty
         */
		bool empty() const {return m_cells.empty();}

		/*------------------------------------------------------------------------*/
		/** \brief  Get cells.
         *
         *  \return the cells contained in this group
         */
		std::vector<TCellID> cells() {
			return m_cells;
		}

		/*------------------------------------------------------------------------*/
		/** \brief  Get group name.
         *
         *  \return a STL string which is the group name
         */
		std::string name() const {return m_name;}

		/*------------------------------------------------------------------------*/
		/** \brief  Set the group name.
         *
         *  \param AName the new name
         */
		void setName(const std::string& AName) {m_name=AName;}

		/*------------------------------------------------------------------------*/
		/** \brief  Gives the number of cells composing this.
         *
         *  \return the number of cells
         */
		int size() const {return m_cells.size();}

		/*------------------------------------------------------------------------*/
		/** \brief  Provides Overload the bracket operator
         *
         *  \param i the index of the elt you want to get
         *
         *  \return parameter in location i
         */
		TCellID operator[] (const int AIndex)
		{
			return m_cells[AIndex];
		}

		TCellID  value(const int AIndex){
		    return m_cells[AIndex];
		}

		/*------------------------------------------------------------------------*/
		/** \brief serialize (*this) in AStr
         *
         * \param AStr an output streammap
         */
		void serialize(std::ostream& AStr){
			const int name_size= m_name.size();
			const int nb_values = m_cells.size();
			const int total_size = sizeof(int) 				/* 1) total size */
								   + sizeof(int)  				/* 2) name size*/
								   +name_size*sizeof(char) 	/* 3) name */
								   +sizeof(int)				/* 4) nb values*/
								   + nb_values*sizeof(TCellID);/* 5) values */


			AStr.write((char*)&total_size,sizeof(int)); 			/* fill (1) */
			AStr.write((char*)&name_size,sizeof(int));				/* fill (2) */
			AStr.write(m_name.c_str(),sizeof(char)*m_name.size());	/* fill (3) */
			AStr.write((char*)&nb_values,sizeof(int));				/* fill (4) */
			AStr.write((char*)&m_cells[0],nb_values*sizeof(TCellID));	/* fill (5) */

		}

		/*------------------------------------------------------------------------*/
		/** \brief unserialize (*this) from AStr
         *
         * \param AStr an input stream
         */
		void unserialize(std::istream& AStr){
				int total_size = 0;
				int name_size  = 0;
				int nb_values = 0;
				AStr.read((char*)&total_size,sizeof(int)); /* read (1) */
				AStr.read((char*)&name_size,sizeof(int));  /* read (2) */
				char *n = new char[name_size];
				AStr.read(n,sizeof(char)*name_size);   	 /* read (3) */
				m_name.clear();
				m_name.assign(n,n+name_size);
				delete[] n;


				AStr.read((char*)&nb_values,sizeof(int));	/* read (4) */
				m_cells.clear();
				m_cells.resize(nb_values);

				for(int i=0;i<nb_values;i++){
					TCellID elt;
					AStr.read((char*)&elt  ,sizeof(TCellID));
					m_cells[i]=elt;

				}
		}

	private:

		/** mesh containing the cells of (*this) */
		Mesh* m_mesh;

		/** name of the group */
		std::string m_name;

		/** ids of the cells contained in this group*/
		std::vector<TCellID> m_cells;
	};
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_CELLGROUP_H_ */
/*----------------------------------------------------------------------------*/
