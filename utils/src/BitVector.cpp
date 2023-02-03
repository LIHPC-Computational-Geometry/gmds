/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/3/19.
//
/*----------------------------------------------------------------------------*/
#include "../inc/gmds/utils/BitVector.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
    /*------------------------------------------------------------------------*/
    BitVector::BitVector(const int capacity): m_top(0),m_size(0) {
        m_bits.resize(capacity,false);
        m_capacity = capacity;
    }
    /*------------------------------------------------------------------------*/
    BitVector::BitVector(const BitVector& v)
            : m_top(v.m_top),m_size(v.m_size), m_capacity(v.m_capacity),
              m_free_stack(v.m_free_stack), m_bits(v.m_bits)
    {}
    /*------------------------------------------------------------------------*/
    BitVector::~BitVector() {}
    /*------------------------------------------------------------------------*/
    BitVector::size_type BitVector::size() const  {return m_size;}
    /*------------------------------------------------------------------------*/
    BitVector::size_type BitVector::top() const {return m_top;}
    /*------------------------------------------------------------------------*/
    BitVector::size_type BitVector::capacity() const {return m_bits.size();}
    /*------------------------------------------------------------------------*/
    bool BitVector::empty() const { return (m_size==0); }
    /*------------------------------------------------------------------------*/
    void BitVector::clear()
    {
        m_free_stack.clear();

        m_bits.clear();
        m_bits.resize(GChunkSize,false);

        m_top=0;
        m_size=0;
        m_capacity=GChunkSize;
    }
    /*------------------------------------------------------------------------*/
    void BitVector::resize(const BitVector::size_type & ASize) {
        m_bits.resize(ASize, false);
        m_capacity = ASize;
        update();

    }

    /*----------------------------------------------------------------------------*/
    void BitVector::fillAll() {
        TInt cap = m_capacity;
        for(int i=0;i<cap;i++){m_bits[i]=1;}
        update();
    }
    /*----------------------------------------------------------------------------*/
    void BitVector::update()
    {
        /* the free stack is made empty */
        m_free_stack.clear();

        TInt free=0;
        /* we traverse the m_bits tabular to know which items are free and we link
         * them together in the m_bits tabular.
         */
        m_size=0;

        int i =m_bits.size()-1;
        while(i>=0 && m_bits[i]==0){
            i--;
        }

        /* in this case, all the items are free */
        if(i==-1){
            m_top=0;
            return;
        }
        /* we are on the greatest used item */
        m_top=i+1;
        free = i;

        for(;i>=0;i--){
            if(m_bits[i]==0)
                m_free_stack.push_back(i);
            else
                m_size++;
        }

    }
    /*----------------------------------------------------------------------------*/
    bool BitVector::operator[]  (const size_type AIndex) const
    {
#ifdef __DEBUG__
        if(AIndex>=m_capacity || m_bits[AIndex]==0)
		throw GMDSException("Bad index in SmartBitVector::operator[]");
#endif //__DEBUG__
        return m_bits[AIndex];
        //return val;
    }
    /*----------------------------------------------------------------------------*/
    BitVector::size_type BitVector::selectNewBit()
    {
        TInt index = getFreeIndex();
        assign(index);
        return index;
    }
/*----------------------------------------------------------------------------*/
    void BitVector::unselect(const TInt& AIndex)
    {
#ifdef __DEBUG__
        if(AIndex>=m_capacity || AIndex<0)
		throw GMDSException("Bad index in SmartBitVector::unselect()");
#endif //__DEBUG__
        /* the item AIndex is already free */
        if(m_bits[AIndex]==false)
            return;

        m_bits[AIndex]=false;
        m_size--;

        if(AIndex==m_top)
            m_top--;
        else
            m_free_stack.push_back(AIndex);

    }
/*----------------------------------------------------------------------------*/
    bool BitVector::isAvailable(const TInt& AIndex) const
    {

        if(AIndex>=m_capacity)
            return false;
        else
            return !(m_bits[AIndex]);
    }
/*----------------------------------------------------------------------------*/
    bool BitVector::isOutOfContainer(const TInt& AIndex) const
    {
        return (AIndex>=m_capacity);
    }
/*----------------------------------------------------------------------------*/
    void BitVector::serialize(std::ostream& stream)
    {
        int container_size = m_top;
        stream.write((char*)&container_size,sizeof(int));

        for(int i=0;i<m_top;i++)
        {
            int m = m_bits[i];
            stream.write((char*)&m,sizeof(int));
        }
    }
/*----------------------------------------------------------------------------*/
    void BitVector::unserialize(std::istream& stream)
    {
        int nb_items=0;

        stream.read((char*)&nb_items,sizeof(int));
        clear();
        resize(nb_items);
        m_top=nb_items;
        for(int i=0;i<nb_items;i++){
            int m;
            stream.read((char*)&m  ,sizeof(int));
            m_bits[i]=m;
            if(m)
                m_size++;
        }
        update();
    }
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
