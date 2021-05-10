/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/* Variable_def.h
 *
 *  Created on: 03/10/2017
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_VARIABLE_DEF_H_
#define KMDS_VARIABLE_DEF_H_
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// KMDS headers
#include "KM/DS/CellHandle.h"
#include "KM/Utils/Exception.h"
#include "KM/Utils/KTypes.h"
/*----------------------------------------------------------------------------*/
// STL headers
#include <iostream>
#include <string>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/
/** \class VariableItf
 *  \brief Defines the API of a mesh variable.
 */
class EXPORT_KMDS VariableItf
{
 public:

    enum export_type{
        var_int,
        var_double,
        var_double_vec,
        var_unknown
    };


        /*------------------------------------------------------------------------*/
        /** \brief  Accessor to the variable name
         */
        virtual std::string getName() const = 0;

        /*------------------------------------------------------------------------*/
        /** \brief  Gives access to the type of variable if it has been specified
         */
        export_type getType();

        /*------------------------------------------------------------------------*/
        /** \brief Set the domain entry to \p ASize
         */
        virtual void resize(const int ASize) = 0;
        /*------------------------------------------------------------------------*/
        /** \brief Get the domain size
         */
        virtual int getSize() const = 0;
        /*------------------------------------------------------------------------*/
        /** \brief initialize the value at index \p AI
         */
        virtual void initialize(const TCellID AI) = 0;
        /*------------------------------------------------------------------------*/
        /** \brief initialize the values from index \p AI to \p AI+\p ANb
         */
        virtual void initialize(const TCellID AI, const TInt32 ANb) = 0;
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
class EXPORT_KMDS Variable : public VariableItf
{
 public:
        /*------------------------------------------------------------------------*/
        /** \brief  Constructor. Every variable has a name
         */
        Variable(const T ADefaultValue, const std::string& AName = "unnamed", const TInt32& ASize = DefaultContainerSize)
                : m_defaultval(ADefaultValue)
                , m_name(AName)
                , m_data(AName, ASize)
        {
                Kokkos::parallel_for(ASize, KOKKOS_LAMBDA(const int i) { initialize(i); });
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor.
         */
        virtual ~Variable()
        {
                ;
        }

        /*------------------------------------------------------------------------*/
        /** \brief  overload operator[] const
         */
        T const& operator[](const Cell& AC) const
        {
                return m_data(AC.id);
        }

        /*------------------------------------------------------------------------*/
        /** \brief  overload operator[]. It allows to modify a variable value
         */
        T& operator[](const Cell& AC)
        {
                return m_data(AC.id);
        }

        /*------------------------------------------------------------------------*/
        /** \brief  overload operator[] const
         */
        T const& operator[](const int& i) const
        {
                return m_data(i);
        }

        /*------------------------------------------------------------------------*/
        /** \brief  overload operator[]. It allows to modify a variable value
         */
        T& operator[](const int& i)
        {
                return m_data(i);
        }

        /*------------------------------------------------------------------------*/
        /** \brief Update all the elements to value AVal. Cannot be used in a
         * parallel kernel
         */
        void
        setValuesTo(const T& AVal)
        {
                int s = m_data.extent(0);
                Kokkos::parallel_for(s, KOKKOS_LAMBDA(const int i) { m_data(i) = AVal; });
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Accessor to the variable name
         */
        std::string
        getName() const
        {
                return m_name;
        }

        /*------------------------------------------------------------------------*/
        /** \brief Set entry \p AI to value \p AVal
         */
        void
        set(const int AI, const T& AVal)
        {
                m_data(AI) = AVal;
        }

        /*------------------------------------------------------------------------*/
        /** \brief Set the domain entry to \p ASize
         */
        void
        resize(const int ASize)
        {
                size_t oldsize = m_data.size();
                Kokkos::resize(m_data, ASize);

                if(ASize > oldsize) {
                        Kokkos::parallel_for(ASize-oldsize, KOKKOS_LAMBDA(const int i) { initialize(oldsize + i); });
                }
        }
        /*------------------------------------------------------------------------*/
        /** \brief initialize the value at index \p AI
         */
         
        void
        initialize(const TCellID AI)
        {
                m_data(AI) = m_defaultval;
        };

        /*------------------------------------------------------------------------*/
        /** \brief initialize the values from index \p AI to \p AI+\p ANb
         */
        void
        initialize(const TCellID AI, const TInt32 ANb)
        {
                for (auto i = AI; i < AI + ANb; i++)
                        initialize(i);
        }

        /*------------------------------------------------------------------------*/
        /** \brief initialize all the elements. Cannot be used in a
         * parallel kernel
         */
        void
        initialize()
        {
            int s = m_data.extent(0);
            Kokkos::parallel_for(s, KOKKOS_LAMBDA(const int i) { m_data(i) = m_defaultval; });
        };

        /*------------------------------------------------------------------------*/
        /** \brief Get the domain size
         */
        int
        getSize() const
        {
                return m_data.extent(0);
        }

        /*------------------------------------------------------------------------*/
        /** \brief Serialize the variable into stream str. Warning this method does
         * 		   not support typename T where pointers would be present.
         */
        void serialize(std::ostream& stream);

        /*------------------------------------------------------------------------*/
        /** \brief Unserialize the variable from stream str. Warning this method
         * 		   does not support typename T where pointers would be present.
         */
        void unserialize(std::istream& stream);

        /*------------------------------------------------------------------------*/
        /** \brief compact the variable data
         */
        void compact();

 private:
        /* variable name*/
        std::string m_name;

        /* data*/
        Kokkos::View<T*> m_data;

        /* default value */
        T m_defaultval;
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* KMDS_VARIABLE_DEF_H_ */
/*----------------------------------------------------------------------------*/