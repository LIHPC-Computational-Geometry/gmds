/*----------------------------------------------------------------------------*/
/* Variable.t.h
 *
 *  Created on: 3 aout 2010
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_VARIABLE_H_
#define GMDS_VARIABLE_H_
/*----------------------------------------------------------------------------*/
#include "CommonTypes.h"
#include "SmartVector.h"
#include "GMDSUtils_export.h"
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds{

/*----------------------------------------------------------------------------*/
/** \class VariableItf
 *  \brief Defines the API of a mesh variable.
 */
    class GMDSUtils_API VariableItf{
    public:

        VariableItf() {}
        virtual ~VariableItf() {}

        /*------------------------------------------------------------------------*/
        /** \brief  Accessor to the variable name
         */
        virtual std::string getName() const = 0;
        /*------------------------------------------------------------------------*/
        /** \brief  Add a new entry for this variable
         */
        virtual void addEntry(const TCellID& i) = 0;
        /*------------------------------------------------------------------------*/
        /** \brief  Remove an entry for this variable
         */
        virtual void removeEntry(const TCellID& i) = 0;
        /*------------------------------------------------------------------------*/
        /** \brief Set the domain entry to [0,i]
         */
        virtual void setDomain(const TCellID& i) = 0;
        /*------------------------------------------------------------------------*/
        /** \brief Get the domain size
         */
        virtual int getDomainSize() const = 0;
        /*------------------------------------------------------------------------*/
        /** \brief Set the domain entry to [0,i] with some default values
         */
        virtual void setDomainWithDefault(const TCellID& i, const std::vector<TInt>& def) = 0;
        /*------------------------------------------------------------------------*/
        /** \brief Clear all the domain
         */
        virtual void clear() = 0;
        /*------------------------------------------------------------------------*/
        /** \brief Serialize the variable into stream str. Warning this method does
         * 		   not support typename T where pointers would be present.
         */
        virtual void serialize(std::ostream& stream) = 0;
        /*------------------------------------------------------------------------*/
        /** \brief Unserialize the variable from stream str. Warning this method
         * 		   does not support typename T where pointers would be present.
         */
        virtual void unserialize(std::istream& stream) = 0;
        /*------------------------------------------------------------------------*/
        /** \brief compact the variable data
         */
        virtual void compact() = 0;
    };
/*----------------------------------------------------------------------------*/
/** \class Variable
 *  \brief Defines what a mesh variable is
 *
 *  \param T  type of the variable
 */
    template <typename T>
    class Variable : public VariableItf {
    public:
/*------------------------------------------------------------------------*/
/** \brief  Constructor. Every variable has a name
 */
        Variable(const std::string& AName = "no name")
                : m_name(AName) {}

/*------------------------------------------------------------------------*/
/** \brief  Destructor.
 */
        virtual ~Variable() {}

/*------------------------------------------------------------------------*/
/** \brief  overload operator[] const
 */
        T const& operator[](const TCellID& i) const { return m_data[i]; }
/*------------------------------------------------------------------------*/
/** \brief  overload operator[]. It allows to modify a variable value
 */
        T& operator[](const TCellID& i) { return m_data[i]; }
        /*------------------------------------------------------------------------*/
/** \brief  overload operator[] const
 */
        T const& value(const TCellID &i) const { return m_data[i]; }
/*------------------------------------------------------------------------*/
/** \brief  overload operator[]. It allows to modify a variable value
 */
        T& value(const TCellID &i) { return m_data[i]; }


/*------------------------------------------------------------------------*/
/** \brief Update all the elements to value AVal
 */
        void setValuesTo(const T& AVal)
        {
            typename SmartVector<T>::iterator it = m_data.begin();

            for (; it!=m_data.end(); it.operator++()) {
                (*it)= AVal;
            }
        }

/*------------------------------------------------------------------------*/
/** \brief  Accessor to the variable name
 */
        std::string getName() const {return m_name; }

/*------------------------------------------------------------------------*/
/** \brief  Get the number of values the variable is defined on
 */
        int getNbValues() const { return m_data.size(); }

/*------------------------------------------------------------------------*/
/** \brief  Add a new entry for this variable
 */
        virtual void
        addEntry(const TCellID& i) {
            if (m_data.isOutOfContainer(i))
                m_data.resize(2 * i);
            // This choice is maybe too expensive
        }

/*------------------------------------------------------------------------*/
/** \brief  Remove an entry for the variable
 */
        virtual void
        removeEntry(const TCellID& i) { m_data.remove(i); }
/*------------------------------------------------------------------------*/
/** \brief Set entry i to value val
 */
        void set(const TCellID& i, const T& val) { m_data.assign(val, i); }

/*------------------------------------------------------------------------*/
/** \brief Set the domain entry to [0,i]
 */
        void setDomain(const TCellID& i) { m_data.resize(i); }

/*------------------------------------------------------------------------*/
/** \brief Get the domain size
 */
        virtual int getDomainSize() const {return m_data.capacity(); }

/*------------------------------------------------------------------------*/
/** \brief Set the domain entry to [0,i] with some default values
 */
        void
        setDomainWithDefault(const TCellID& i, const std::vector<TInt>& def)
        {
            m_data.resize(i);
            if (def.size() > (unsigned int)i)
                return;

            m_data.mark(def);

            // to make iterators consistent
            m_data.update();
        }
/*------------------------------------------------------------------------*/
/** \brief Clear all the domain
 */
        void clear() { m_data.clear(); }

/*------------------------------------------------------------------------*/
/** \brief Serialize the variable into stream str. Warning this method does
 * 		   not support typename T where pointers would be present.
 */
        void serialize(std::ostream& stream) {
            // necessary to use iterators
            m_data.update();

            const size_t name_size = m_name.size();
            const size_t nb_values = m_data.top();

            const size_t total_size = sizeof(int)            /* total size */
                                      + sizeof(int)            /* name size*/
                                      + name_size * sizeof(char) /* name */
                                      + sizeof(int)            /*nb values*/
                                      + nb_values * sizeof(T);
//							(sizeof(id)	+		/* id per value */
//							 sizeof(T));		/* one value*/

            stream.write((char *) &total_size, sizeof(int));
            stream.write((char *) &name_size, sizeof(int));
            stream.write(m_name.c_str(), sizeof(char) * m_name.size());

            m_data.serialize(stream);
        }
/*------------------------------------------------------------------------*/
/** \brief Unserialize the variable from stream str. Warning this method
 * 		   does not support typename T where pointers would be present.
 */
        void unserialize(std::istream& stream) {
            int total_size = 0;
            int name_size  = 0;

            stream.read((char*)&total_size,sizeof(int));
            stream.read((char*)&name_size,sizeof(int));
            char *n = new char[name_size];
            stream.read(n,sizeof(char)*name_size);
            m_name.clear();
            m_name.assign(n,n+name_size);
            delete[] n;
            std::cout<<"VAR UNSERIALIZATION: "<<name_size<<" - "<<m_name<<std::endl;
            m_data.unserialize(stream);
        }


/*------------------------------------------------------------------------*/
/** \brief compact the variable data
 */
        void compact(){ m_data.compact(); }

    private:
/* variable name*/
        std::string m_name;

/* data*/
        SmartVector<T> m_data;
    };

/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_VARIABLE_H_ */
/*----------------------------------------------------------------------------*/
