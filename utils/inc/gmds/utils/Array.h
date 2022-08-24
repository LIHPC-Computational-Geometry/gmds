/*----------------------------------------------------------------------------*/
#ifndef GMDS_ARRAY_H
#define GMDS_ARRAY_H
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <iostream>
#include <vector>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    /** \class Array2D
     *  \brief class implementing a 2D array made of T type objects
     *
     */
    /*----------------------------------------------------------------------------*/
    template<typename  T>
    class Array2D {
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
                    m_tab[i][j] = A(i,j);
                }
            }
        }
        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator =
        */
        Array2D&  operator=(const Array2D& A){
            m_i = A.m_i;
            m_j = A.m_j;
            m_tab.resize(m_i);
            for (auto i = 0; i < m_i; i++) {
                m_tab[i].resize(m_j);
                for (auto j = 0; j < m_j; j++) {
                    m_tab[i][j] = A(i,j);
                }
            }
		      return *this;
        }
        /*------------------------------------------------------------------------*/
        /** \brief Destructor
        */
        virtual  ~Array2D(){;}

        /*------------------------------------------------------------------------*/
        /** \brief give the number of lines
        */
        TInt  nbLines() const {return m_i;}
        /*------------------------------------------------------------------------*/
        /** \brief give the number of columns
        */
        TInt  nbColumns() const {return m_j;}

        /*------------------------------------------------------------------------*/
        /** \brief  Access to  @p AI th line of the array
        *       \return the line contain
        *       \exception if i<0 or i>size()
        */
        inline  T const& operator()(const TInt AI ,const TInt AJ ) const{
            return m_tab[AI][AJ];
        }

        inline  T& operator()(const TInt AI ,const TInt AJ ) {
            return m_tab[AI][AJ];
        }

    private:
        TInt m_i;
        TInt m_j;
        std::vector<std::vector<T> > m_tab;

    };
/*----------------------------------------------------------------------------*/
    /** \class TriArray
     *  \brief class implementing a triangular array, that is a kind of array,
     *         that store data in a N-side triangles:
     *         1 2 3 ......... N
     *         1 2 3 ..... N -1
     *         .
     *         .
     *         1 2 3
     *         1 2
     *         1
     *
     */
    /*----------------------------------------------------------------------------*/
    template<typename  T>
    class  TriArray {
    public:
        /*------------------------------------------------------------------------*/
        /** @brief  Default constructor with @p AN the triangular size
         *
         * @param AN triangular side size
         *
         */
         TriArray(const int AN) {
            m_N=AN;
            /* the value are stored in a 1-dim array containing the N elements of
             * the 1st line, then the N-1 of the second line then ... the single
             * element of the Nth line */
            TInt s = m_N;
            TInt tab_size =0;
            while(s>0){
                tab_size +=s;
                s-=1;
            }
            m_tab.resize(tab_size);
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Copy constructor
         */
         TriArray(const TriArray& A): m_N(A.m_N) {
            m_tab.resize(A.m_tab.size());
            for(auto i=0; i<A.m_tab.size();i++){
                m_tab[i]=A.m_tab[i];
            }
        }
        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator =
        */
        TriArray& operator=(const TriArray& A){
            m_tab.resize(A.m_tab.size());
            for(auto i=0; i<A.m_tab.size();i++){
                m_tab[i]=A.m_tab[i];
            }
        }
        /*------------------------------------------------------------------------*/
        /** \brief Destructor
        */
        virtual  ~TriArray(){;}

        /*------------------------------------------------------------------------*/
        /** \brief give the number of lines
        */
        TInt  size() const {return m_N;}

        TInt  nbElements() const {return m_tab.size();}
        /*------------------------------------------------------------------------*/
        /** \brief  Access to  @p AI th line of the array
        *       \return the line contain
        *       \exception if i<0 or i>size()
        */
        inline  T const& operator()(const TInt AI ,const TInt AJ ) const{
            return m_tab[compute1DIndexFromIJ(AI,AJ)];
        }

        inline  T& operator()(const TInt AI ,const TInt AJ) {
            return m_tab[compute1DIndexFromIJ(AI,AJ)];
        }
        inline  T const& operator()(const TInt AI ,const TInt AJ, const TInt AK) const{
            return m_tab[compute1DIndexFromIJK(AI,AJ,AK)];
        }

        inline  T& operator()(const TInt AI ,const TInt AJ, const TInt AK) {
            return m_tab[compute1DIndexFromIJK(AI,AJ,AK)];
        }

        void  print(){
            for(auto i:m_tab){
                std::cout<<i<<" ";
            }
            std::cout<<std::endl;
        }
    private:
        int compute1DIndexFromIJ(const TInt AI ,const TInt AJ) {
            if(AI>=m_N || AI+AJ>=m_N)
                throw GMDSMathException("invalid values");
            int index_line = 0;
            int i_tmp = 1;
            while (i_tmp<=AI){
                index_line += m_N-i_tmp+1;
                i_tmp++;
            }
            return  index_line+AJ;
        }
        int compute1DIndexFromIJK(const TInt AI ,const TInt AJ, const int AK) {
            if(AI<0 || AJ<0 ||AK<0 || AI+AJ+AK!=m_N-1)
                throw GMDSMathException("invalid (i,j,k) barycentric values: ("
                                        + std::to_string(AI)+", "
                                        + std::to_string(AJ)+", "
                                        + std::to_string(AJ)+")");

            //(N,0,0) is the first tab element
            //(O,N,0) is the last element of the first line
            //(0,0,N) is the last tabular element
            int index_line = 0;
            int i_tmp = 1;
            while (i_tmp<=AK){
                index_line += m_N-i_tmp+1;
                i_tmp++;
            }
            //we've got the line. No we traverse it according to AJ value
            return  index_line+AJ;
        }
    private:
        TInt m_N;
        std::vector<T> m_tab;

    };
    /** \class Array3D
     *  \brief class implementing a 2D array made of T type objects
     *
     */
    /*----------------------------------------------------------------------------*/
    template<typename  T>
    class Array3D {
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
		      return *this;
        }
        /*------------------------------------------------------------------------*/
        /** \brief Destructor
        */
        virtual  ~Array3D(){;}

        /*------------------------------------------------------------------------*/
        /** \brief give the number of elements for dim @p ADim
        */
        TInt  nbElements(const int ADim) const {return m_dim[ADim];}

        /*------------------------------------------------------------------------*/
        /** \brief  Access to  @p AI th line of the array
        *       \return the line contain
        *       \exception if i<0 or i>size()
        */
        inline  T const& operator()(const TInt AI ,const TInt AJ ,
                const TInt AK) const{
            return m_tab[AI][AJ][AK];
        }

        inline  T& operator()(const TInt AI ,const TInt AJ,const TInt AK ) {
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
