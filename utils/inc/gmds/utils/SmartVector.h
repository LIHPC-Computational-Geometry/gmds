/*----------------------------------------------------------------------------*/
/** \file    SmartVector.t.h
 *  \author  F. LEDOUX
 *  \date    04/06/2009
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_SMARTVECTOR_H_
#define GMDS_SMARTVECTOR_H_
/*----------------------------------------------------------------------------*/
#include <vector>
/*----------------------------------------------------------------------------*/
#include "CommonTypes.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    template<typename T> class SmartVector{

    public:
        typedef int size_type;

        class iterator {
        public:

            friend class SmartVector<T>;

            using self_type = iterator;
            using iterator_category = std::forward_iterator_tag;
            using value_type = T;
            using difference_type = int;
            using pointer = T*;
            using reference = T&;

            iterator(SmartVector<T>* AContainer):container_(AContainer),current_(0){
                // we value the location of the first real item
                while(current_!=container_->top_ && container_->mark_[current_]==0)
                    current_++;
            }

            iterator(SmartVector<T>* AContainer, size_type ACurrent):container_(AContainer),current_(ACurrent){}
            iterator(const iterator& AIt):container_(AIt.container_),current_(AIt.current_){}


            self_type operator++() {
                // move to the first real successor of index_
                do{ current_++;
                }while(current_!=container_->top_ && container_->mark_[current_]==0);
                return *this;
            }

            reference operator*()  { return container_->vec_[current_]; }
            pointer operator->() { return &(container_->vec_[current_]); }
            bool operator==(const self_type& rhs) { return (container_== rhs.container_ && current_==rhs.current_);}
            bool operator!=(const self_type& rhs) { return (container_!= rhs.container_ || current_!=rhs.current_);}
        private:
            SmartVector<T>* container_;
            size_type current_;
        };

        class const_iterator
        {
        public:
            using self_type= const_iterator;
            using iterator_category = std::forward_iterator_tag;
            using value_type = T;
            using difference_type = int;
            using pointer = T*;
            using reference = T&;


            const_iterator(const SmartVector<T>* AContainer):container_(AContainer),current_(0){
                // we value the location of the first real item
                while(current_!=container_->top_ && container_->mark_[current_]==0)
                    current_++;
            }

            const_iterator(const SmartVector<T>* AContainer, size_type ACurrent):container_(AContainer),current_(ACurrent){}
            const_iterator(const const_iterator& AIt):container_(AIt.container_),current_(AIt.current_){}


            self_type operator++() {
                // move to the first real successor of index_
                do{ current_++;
                }while(current_!=container_->top_ && container_->mark_[current_]==0);
            }
            reference  operator*() { return container_->vec_[current_]; }
            pointer operator->() { return &(container_->vec_[current_]); }
            bool operator==(const self_type& rhs) { return (container_== rhs.container_ && current_==rhs.current_);}
            bool operator!=(const self_type& rhs) { return (container_!= rhs.container_ || current_!=rhs.current_);}
        private:
            const SmartVector<T>* container_;
            size_type current_;
        };

        /*------------------------------------------------------------------------*/
        /** \brief  Default Constructor
         */
        SmartVector(const int capacity=GChunkSize): top_(0),size_(0)
        {
            vec_.resize(capacity);
            mark_.resize(capacity,false);
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Default Constructor
         */
        SmartVector(std::vector<bool> AMarks): mark_(AMarks)
        {
            vec_.resize(AMarks.size());
            update();
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Copy constructor. Note the algorithm is linear in its size.
         */
        SmartVector(const SmartVector<T>& AVec)
                :top_(AVec.top_),size_(AVec.size_),
                 free_stack_(AVec.free_stack_), vec_(AVec.vec_), mark_(AVec.mark_)
        {}

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor. Note the algorithm is linear in the number of holes
         * 			in the container.
         */
        ~SmartVector(){}

        /*------------------------------------------------------------------------*/
        /** \brief  Gives the number of elements stored in the container.
         */
        TInt size() const {return size_;}

        /*------------------------------------------------------------------------*/
        /** \brief  Gives the next right index after the last used item.
         */
        TInt top() const {return top_;}


        /*------------------------------------------------------------------------*/
        /** \brief  Gives the size of the container, i.e, the effective number of
         * 			items that are available in the container. It corresponds to
         * 			the number of stored elements + the number of free spaces.
         */
        TInt capacity() const {return vec_.size();}


        /*------------------------------------------------------------------------*/
        /** \brief  Indicates if the container is empty.
         */
        bool empty() const {return (size_==0);}

        /*------------------------------------------------------------------------*/
        /** \brief  Clear the container and resizes it to 2.
         */
        void clear() {
            if(!free_stack_.empty())
                free_stack_.clear();

            if(!vec_.empty()){
                vec_.clear();
                vec_.resize(GChunkSize);
            }
            if(!mark_.empty()){
                mark_.clear();
                mark_.resize(GChunkSize,false);
            }
            top_=0;
            size_=0;
        }
        /*------------------------------------------------------------------------*/
        /** \brief  resizes the container to ASize. size and top are  both
         * 			equal to ASize. If ASize>current(size) it is just extended,
         * 			new items are marked free. Otherwise, all the items between
         * 			ASize and current(size) are lost.
         */
        void resize(const TInt& ASize){
            vec_.resize(ASize);
            mark_.resize(ASize,false);
            update();
        }


        /*------------------------------------------------------------------------*/
        /** \brief  make the container full of elements putting all the marks to 1
         * 		    and top_ equals to size_;
         */
        void fillAll() {
            TInt cap = capacity();
            for(int i=0;i<cap;i++)
                mark_[i]=true;
            update();

        }

        /*------------------------------------------------------------------------*/
        /** \brief  This method is necessary when you want to regularize the
         * 			container after several assignements. Indeed assignement method
         * 			has the advantage to check nothing before insertion. The problem
         * 			is then that the container is no more coherent. This method fix
         * 			the container.
         */
        void update(){
            /* the free stack is made empty */
            free_stack_.clear();

            /* we traverse the marks_ tabular to know which items are free and we link
             * them together in the vec_ tabular.
             */
            size_=0;

            int i =mark_.size()-1;
            while(i>=0 && mark_[i]==false){
                i--;
		    }

            /* in this case, all the items are free */
            if(i==-1){
                top_=0;
                return;
            }
            /* we are on the greatest used item */
            top_=i+1;

            for(;i>=0;i--){
                if(mark_[i]==false)
                    free_stack_.push_back(i);
                else
                    size_++;
            }

        }

        /*------------------------------------------------------------------------*/
        /** \brief  Give the value stored in the AIndex item but do not allow to
         * 			modify it.
         *
         * 			Warning, in release mode, no test is performed to ensure the
         * 			choice of AIndex.
         */
        inline T const& operator[](const TInt& AIndex) const{
#ifdef __DEBUG__
            if(AIndex>=capacity() || mark_[AIndex]==0)
			throw GMDSException("Bad index in SmartVector<T>::operator[]");
#endif //__DEBUG__

            return vec_[AIndex];

        }
        /*------------------------------------------------------------------------*/
        /** \brief  Give the value stored in the AIndex item and allow to
         * 			modify it.
         *
         * 			Warning, in release mode, no test is performed to ensure the
         * 			choice of AIndex. Moreover, in order to improve performances,
         * 			the container structure is partially corrupted. Thus an update
         * 			is necessary at the end (especially for the iterators)
         */
        inline T& operator[](const TInt& AIndex){
#ifdef __DEBUG__
            if(AIndex>=capacity() || mark_[AIndex]==0)
			throw GMDSException("Bad index in SmartVector<T>::operator[]");
#endif //__DEBUG__
            if(!mark_[AIndex])
            {
                size_++;
                mark_[AIndex]=1;

            }
            return vec_[AIndex];

        }


        /*------------------------------------------------------------------------*/
        /** \brief  Provides the next index where an item can be added. It also
                    reserve the necessary space for it and update some internal
                    data. It must be used with the assign method.
         */

        TInt  selectNewIndex() {
            TInt free;
            TInt cap=capacity();

            if(free_stack_.empty()){
                free= top_;
                top_++;
                if(top_>=cap){
                    vec_.resize((cap+1)*2);
                    mark_.resize((cap+1)*2,false);
                }
            }
            else{
                free = free_stack_.back();
                free_stack_.pop_back();
            }
            return free;
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Assign AElt to index AIndex. If this index is still occupied,
         * 			it contents is replaced. If AIndex is out of the vector bondary
         * 			nothing is specified to avoid the insertion.
         */
        void assign(T AElt, const TInt& AIndex){
            vec_[AIndex]=AElt;
            if(mark_[AIndex]==false)
            {
                mark_[AIndex]=true;
                size_++;
            }
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Add AElt. The container is responsible of finding the index
         * 			where to add the element. This operation is SAFER than assign.
         */
        void add(T AElt) {assign(AElt,selectNewIndex()); }

        /*------------------------------------------------------------------------*/
        /** \brief  find if an AElt is inside the vector. It returns true if it is,
         * 			false otherwise.
         */
        bool find(T AElt) const{

            for(T t: *this) {if (t == AElt){ return true; } }
            return false;
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Remove the element stored in AIndex.
         *
         * 			Warning - In release mode no test is performed on AIndex. In
         * 			debug mode, AIndex must be in [0, size[
         */
        void remove(const TInt& AIndex)
        {
#ifdef __DEBUG__
            if(AIndex>=capacity() || AIndex<0)
		throw GMDSException("Bad index in SmartVector<T>::remove()");
#endif //__DEBUG__
            /* the item AIndex is already free */
            if(mark_[AIndex]==false)
                return;
            T t = T();
            vec_[AIndex]=t;
            mark_[AIndex]=false;
            size_--;

            if(AIndex==top_)
                top_--;
            else
                free_stack_.push_back(AIndex);

        }

        /*------------------------------------------------------------------------*/
        /** \brief  Remove the element AElt
         */
        void removeElement(const T& AElt)
        {
            /* the item AIndex is already free */
            iterator it = this->begin();
            iterator ite = this->end();
            for(T t: *this) {
                if (t == AElt){
                    remove(*it);
                    return;
                }
            }

        }


        /*------------------------------------------------------------------------*/
        /** \brief  Indicates if the AIndex item is available.
         */
        bool isAvailable(const TInt& AIndex) const {
            if(AIndex>=capacity())
                return false;
            else if (mark_[AIndex]==true)
                return false;

            return true;
        }
        /*------------------------------------------------------------------------*/
        /** \brief  Indicates if the AIndex item is out of the container.
         */
        bool isOutOfContainer(const TInt& AIndex) const{
            return (AIndex>=capacity());
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Display the content of the container
         */
        void display() {
            for(unsigned int i=0;i<vec_.size();i++)
                std::cout<<vec_[i]<<" ";
            std::cout<<std::endl;
            for(unsigned int i=0;i<mark_.size();i++)
                std::cout<<mark_[i]<<" ";
            std::cout<<std::endl;
            std::cout<<"top: "<<top_
                     <<" - size: "<<size_
                     <<" - capacity:"<<capacity()<<std::endl;
        }


        void mark(const std::vector<TInt>& ARef){
            for(unsigned int i=0; i<ARef.size();i++)
                mark_[ARef[i]]=true;
            update();
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Provides an iterator on the first element of the container.
         */
        iterator begin()
        {
            return iterator(this);
        }

        iterator end()
        {
            return iterator(this, top_);
        }

        const_iterator begin() const
        {
            return const_iterator(this);
        }

        const_iterator end() const
        {
            return const_iterator(this,top_);
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Compact the container by removing all free items.
         */
        void compact() {
            TInt first_not_free=0;
            TInt last_used = top_;

            while(last_used>first_not_free){
                /* increase first free*/
                while(mark_[first_not_free]==true)
                    first_not_free++;

                /* decrease last used*/
                while(mark_[last_used]==false)
                    last_used--;

                if(last_used>first_not_free){
                    vec_[first_not_free] = vec_[last_used];
                    mark_[first_not_free]=true;
                    mark_[last_used]=false;
                }
                //else, we have a compact collection
            }

            /* new dimension of the container */
            top_=first_not_free;
            std::vector<T> new_vec;
            new_vec.resize(top_);
            std::vector<bool> new_mark;
            new_mark.resize(top_,false);

            new_vec.assign(vec_.begin(),vec_.end());

            for(int i=0;i<top_;i++)
                new_mark[i] = mark_[i];

            vec_.swap(new_vec);
            mark_.swap(new_mark);

            /* the free stack is made empty */
            free_stack_.clear();
        }
        /*------------------------------------------------------------------------*/
        /** \brief  Compact the container by removing all free items and provide
         * 			the movement of elements in AMove. Elt in location i is moved
         * 			to location AMove[i].
         *
         */
        void compactWithMemory(std::vector<int>& AMove){
            TInt first_not_free=0;
            TInt last_used = top_;

            AMove.clear();
            /* +1 is used in case of size_=0*/
            AMove.resize(top_+1);

            while(last_used>first_not_free){
                /* increase first free*/
                while(mark_[first_not_free]==true){
                    AMove[first_not_free]= first_not_free;
                    first_not_free++;
                }
                /* decrease last used*/
                while(mark_[last_used]==false){
                    AMove[last_used]=last_used;
                    last_used--;
                }

                if(last_used>first_not_free){
                    vec_[first_not_free] = vec_[last_used];
                    mark_[first_not_free]=true;
                    mark_[last_used]=false;
                    AMove[last_used]= first_not_free;
                }
                //else, we have a compact collection

            }

            /* new dimension of the container */
            top_=first_not_free;
            std::vector<T> new_vec;
            new_vec.resize(top_);
            std::vector<bool> new_mark;
            new_mark.resize(top_,false);

            new_vec.assign(&vec_[0],&vec_[top_]);

            for(int i=0;i<top_;i++)
                new_mark[i] = mark_[i];

            vec_.swap(new_vec);
            mark_.swap(new_mark);

            /* the free stack is made empty */
            free_stack_.clear();
        }


        /*------------------------------------------------------------------------*/
        /** \brief Serialize the variable into stream str. Warning this method does
         * 		   not support typename T where pointers would be present.
         */
        inline void serialize(std::ostream& stream)
        {
            int container_size = top_;
            stream.write((char*)&container_size,sizeof(int));
            for(int i=0;i<top_;i++)
            {
                T elt = vec_[i];
                int m = mark_[i];
                stream.write((char*)&elt,sizeof(T));
                stream.write((char*)&m,sizeof(int));
            }
        }

        /*------------------------------------------------------------------------*/
        /** \brief Unserialize the variable from stream str. Warning this method
         * 		   does not support typename T where pointers would be present.
         */
        inline void unserialize(std::istream& stream)
        {
            int nb_items=0;

            stream.read((char*)&nb_items,sizeof(int));
            clear();
            resize(nb_items);
            top_=nb_items;
            for(int i=0;i<nb_items;i++){
                T elt;
                int m;
                stream.read((char*)&elt  ,sizeof(T  ));
                stream.read((char*)&m  ,sizeof(int));
                vec_[i]=elt;
                mark_[i]=m;
                if(m)
                    size_++;
            }
            update();
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Return the vector of the boolean mark.
         *
         * \return the vector of the boolean mark
         */
        std::vector<bool> getMarks() {
            return mark_;
        }

    protected:
        /* top_ is the index of the the right next item after the last used)*/
        TInt top_;
        TInt size_;
        std::vector<TInt> free_stack_;
        std::vector<T> vec_;
        /* mark_[i] = 0 indicates that vec_[i] is free, otherwise it is not
         * available */
        std::vector<bool> mark_;
    };

/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_SMARTVECTOR_H_ */
/*----------------------------------------------------------------------------*/
