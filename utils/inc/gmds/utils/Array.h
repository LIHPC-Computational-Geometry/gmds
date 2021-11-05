/*----------------------------------------------------------------------------*/
#ifndef GMDS_ARRAY_H
#define GMDS_ARRAY_H
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <iostream>
#include <vector>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <GMDSUtils_export.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    /** \class Array2D
     *  \brief class implementing a 2D array made of T type objects
     *
     */
    /*----------------------------------------------------------------------------*/
    template<typename  T>
    class GMDSUtils_API Array2D {
    public:
        /*------------------------------------------------------------------------*/
        /** @brief  Default constructor with @p AI lines and @p AJ columns. A
         *          AN exception is thrown if @p AI or @p AJ is <=0.
         *
         * @param AI number of lines
         * @param AJ number of columns
         */
        Array2D(const int AI, const int AJ) {
            m_i = AI; m_j = AJ;
            m_tab.resize(m_i);
            for(auto i=0;i<m_i; i++){
                m_tab[i].resize(m_j);
            }
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Copy constructor
         */
        Array2D(const Array2D& A) {
            m_i = A.m_i;
            m_j = A.m_j;
            m_tab.resize(m_i);
            for (auto i = 0; i < m_i; i++) {
                m_tab[i].resize(m_j);
                for (auto j = 0; j < m_j; j++) {
                    m_tab[i][j] = A[i][j];
                }
            }
        }
        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator =
        */
        Array2D& operator=(const Array2D& A){
            m_i = A.m_i;
            m_j = A.m_j;
            m_tab.resize(m_i);
            for (auto i = 0; i < m_i; i++) {
                m_tab[i].resize(m_j);
                for (auto j = 0; j < m_j; j++) {
                    m_tab[i][j] = A[i][j];
                }
            }
        }
        /*------------------------------------------------------------------------*/
        /** \brief Destructor
        */
        virtual ~Array2D(){;}

        /*------------------------------------------------------------------------*/
        /** \brief give the number of lines
        */
        TInt nbLines() const {return m_i;}
        /*------------------------------------------------------------------------*/
        /** \brief give the number of columns
        */
        TInt nbColumns() const {return m_j;}

        /*------------------------------------------------------------------------*/
        /** \brief  Access to  @p AI th line of the array
        *       \return the line contain
        *       \exception if i<0 or i>size()
        */
        inline T const& operator()(const TInt AI ,const TInt AJ ) const{
            return m_tab[AI][AJ];
        }

        inline T& operator()(const TInt AI ,const TInt AJ ) {
            return m_tab[AI][AJ];
        }

    private:
        TInt m_i;
        TInt m_j;
        std::vector<std::vector<T> > m_tab;

    };
    /** \class Array2D
     *  \brief class implementing a 2D array made of T type objects
     *
     */
    /*----------------------------------------------------------------------------*/
    template<typename  T>
    class GMDSUtils_API Array3D {
    public:
        /*------------------------------------------------------------------------*/
        /** @brief  Default constructor with @p AI, @p AJ and @p AK dimension. A
         *          AN exception is thrown if @p AI, @p AJ or @p AK is <=0.
         *
         * @param AI dimension 1
         * @param AJ dimension 2
         * @param AK dimension 3
         */
        Array3D(const int AI, const int AJ, const int AK) {
            m_dim[0] = AI;
            m_dim[1] = AJ;
            m_dim[2] = AK;
            m_tab.resize(m_dim[0]);
            for(auto i=0;i<m_dim[0]; i++){
                m_tab[i].resize(m_dim[1]);
                for(auto j=0;j<m_dim[1]; j++){
                    m_tab[i][j].resize(m_dim[2]);
                }
            }
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Copy constructor
         */
        Array3D(const Array3D& A) {
            m_dim[0] = A.m_dim[0];
            m_dim[1] = A.m_dim[1];
            m_dim[2] = A.m_dim[2];
            m_tab.resize(m_dim[0]);
            for(auto i=0;i<m_dim[0]; i++){
                m_tab[i].resize(m_dim[1]);
                for(auto j=0;j<m_dim[1]; j++){
                    m_tab[i][j].resize(m_dim[2]);
                    for(auto k=0;k<m_dim[2]; k++){
                        m_tab[i][j][k]=A(i,j,k);
                    }
                }
            }
        }
        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator =
        */
        Array3D& operator=(const Array3D& A){
            m_dim[0] = A.m_dim[0];
            m_dim[1] = A.m_dim[1];
            m_dim[2] = A.m_dim[2];
            m_tab.resize(m_dim[0]);
            for(auto i=0;i<m_dim[0]; i++){
                m_tab[i].resize(m_dim[1]);
                for(auto j=0;j<m_dim[1]; j++){
                    m_tab[i][j].resize(m_dim[2]);
                    for(auto k=0;k<m_dim[2]; k++){
                        m_tab[i][j][k]=A(i,j,k);
                    }
                }
            }
        }
        /*------------------------------------------------------------------------*/
        /** \brief Destructor
        */
        virtual ~Array3D(){;}

        /*------------------------------------------------------------------------*/
        /** \brief give the number of elements for dim @p ADim
        */
        TInt nbElements(const int ADim) const {return m_dim[ADim];}

        /*------------------------------------------------------------------------*/
        /** \brief  Access to  @p AI th line of the array
        *       \return the line contain
        *       \exception if i<0 or i>size()
        */
        inline T const& operator()(const TInt AI ,const TInt AJ ,
                const TInt AK) const{
            return m_tab[AI][AJ][AK];
        }

        inline T& operator()(const TInt AI ,const TInt AJ,const TInt AK ) {
            return m_tab[AI][AJ][AK];
        }

    private:
        TInt m_dim[3];
        std::vector<std::vector<std::vector<T> > >m_tab;

    };

/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_ARRAY_H
/*----------------------------------------------------------------------------*/
