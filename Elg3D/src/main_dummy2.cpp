/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    main_elg2d.cpp
 *  \author  legoff
 *  \date    07/04/2018
 */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
struct fillgv
{
    kmds::GrowingView<kmds::TCellID>* gv;

    fillgv(
            kmds::GrowingView<kmds::TCellID>* gv_)
            : gv(gv_)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void
    operator()(int i) const
    {
//        if(i%2 == 0) {
            gv->push_back(i);
        //}
    }
};


int main(int argc, char* argv[])
{
    int num_threads = -1;

    std::istringstream iss( argv[1] );
    int val;

    if (iss >> val)
    {
        num_threads = val;
    } else {
        std::cerr<<"could not convert number of threads argument."<<std::endl;
        exit(-1);
    }

    Kokkos::InitArguments kargs;
    kargs.num_threads = num_threads;
    Kokkos::initialize(kargs);



    Kokkos::Timer timer;
    timer.reset();

    const int nbelems = 100000000;


    std::cout<<"timer kmap "<<timer.seconds()<<std::endl;

    std::cout<<"kmap.size() "<<kmap.size()<<std::endl;

    Kokkos::finalize();
    return 0;
}
*/

#include <Kokkos_Core.hpp>
#include <cstdio>
#include <cstdlib>
#include <cmath>

// Type of a one-dimensional length-N array of int.
typedef Kokkos::View<int*> view_type;
typedef view_type::HostMirror host_view_type;
// This is a "zero-dimensional" View, that is, a View of a single
// value (an int, in this case).  Access the value using operator()
// with no arguments: e.g., 'count()'.
//
// Zero-dimensional Views are useful for reduction results that stay
// resident in device memory, as well as for irregularly updated
// shared state.  We use it for the latter in this example.
typedef Kokkos::View<int> count_type;
typedef count_type::HostMirror host_count_type;



// Functor for finding a list of primes in a given set of numbers.  If
// run in parallel, the order of results is nondeterministic, because
// hardware atomic updates do not guarantee an order of execution.
struct findprimes {
    view_type data;
    view_type result;
    count_type count;

    findprimes (view_type data_, view_type result_, count_type count_) :
    data (data_), result (result_), count (count_)
    {}

    // Test if data(i) is prime.  If it is, increment the count of
    // primes (stored in the zero-dimensional View 'count') and add the
    // value to the current list of primes 'result'.
    KOKKOS_INLINE_FUNCTION
    void operator() (const int i) const {

        const int number = data(i); // the current number
//
//        // Test all numbers from 3 to ceiling(sqrt(data(i))), to see if
//        // they are factors of data(i).  It's not the most efficient prime
//        // test, but it works.
//        const int upper_bound = std::sqrt(1.0*number)+1;
//        bool is_prime = !(number%2 == 0);
//        int k = 3;
//        while (k < upper_bound && is_prime) {
//            is_prime = !(number%k == 0);
//            k += 2; // don't have to test even numbers
//        }

//        if (is_prime) {


            // Use an atomic update both to update the current count of
            // primes, and to find a place in the current list of primes for
            // the new result.
            //
            // atomic_fetch_add results the _current_ count, but increments
            // it (by 1 in this case).  The current count of primes indexes
            // into the first unoccupied position of the 'result' array.
            const int idx = Kokkos::atomic_fetch_add (&count(), 1);
            result(idx) = number;

//        }
    }

};

struct findprimes_chunk {
    view_type data;
    view_type result;
    count_type count;
    int nnumbers;
    int chunk_size;


    findprimes_chunk (view_type data_, view_type result_, count_type count_, int nnumbers_, int chunk_size_) :
            data (data_), result (result_), count (count_), nnumbers(nnumbers_), chunk_size(chunk_size_)
    {}

    // Test if data(i) is prime.  If it is, increment the count of
    // primes (stored in the zero-dimensional View 'count') and add the
    // value to the current list of primes 'result'.
    KOKKOS_INLINE_FUNCTION
    void operator() (const int i) const {

        int stored[chunk_size];
        int nbstored = 0;

        int startindex = i*chunk_size;
        int endindex = std::min(startindex, nnumbers);
        for(int ii=startindex; ii<endindex; ii++) {

            const int number = data(ii); // the current number

            stored[nbstored] = number;
            nbstored++;
        }

//        // Test all numbers from 3 to ceiling(sqrt(data(i))), to see if
//        // they are factors of data(i).  It's not the most efficient prime
//        // test, but it works.
//        const int upper_bound = std::sqrt(1.0*number)+1;
//        bool is_prime = !(number%2 == 0);
//        int k = 3;
//        while (k < upper_bound && is_prime) {
//            is_prime = !(number%k == 0);
//            k += 2; // don't have to test even numbers
//        }

//        if (is_prime) {


        // Use an atomic update both to update the current count of
        // primes, and to find a place in the current list of primes for
        // the new result.
        //
        // atomic_fetch_add results the _current_ count, but increments
        // it (by 1 in this case).  The current count of primes indexes
        // into the first unoccupied position of the 'result' array.
        const int idx = Kokkos::atomic_fetch_add (&count(), nbstored);
        for(int ii=0; ii<nbstored; ii++) {
            result(idx+ii) = stored[ii];
        }

//        }
    }

};

int main (int argc, char* argv[])
{
    int num_threads = -1;

    std::istringstream iss( argv[1] );
    int val;

    if (iss >> val)
    {
        num_threads = val;
    } else {
        std::cerr<<"could not convert number of threads argument."<<std::endl;
        exit(-1);
    }

    Kokkos::InitArguments kargs;
    kargs.num_threads = num_threads;
    Kokkos::initialize(kargs);

    srand (61391); // Set the random seed

    Kokkos::Timer timer;
    timer.reset();

//    int nnumbers = 100000000;
    const int nnumbers = 100000000;
    view_type data ("RND", nnumbers);
    view_type result ("Prime", nnumbers);
    count_type count ("Count");

    host_view_type h_data = Kokkos::create_mirror_view (data);
    host_view_type h_result = Kokkos::create_mirror_view (result);
    host_count_type h_count = Kokkos::create_mirror_view (count);

    typedef view_type::size_type size_type;
    // Fill the 'data' array on the host with random numbers.  We assume
    // that they come from some process which is only implemented on the
    // host, via some library.  (That's true in this case.)
    for (size_type i = 0; i < data.dimension_0 (); ++i) {
        h_data(i) = rand () % nnumbers;
    }
    Kokkos::deep_copy (data, h_data); // copy from host to device


    timer.reset();
//    Kokkos::parallel_for (nnumbers, findprimes (data, result, count));
//    Kokkos::deep_copy (h_count, count); // copy from device to host
//    std::cout<<"timer poyop "<<timer.seconds()<<std::endl;

    const int chunk_size = 10000;

    timer.reset();
    const int niter = (double) nnumbers / (double) chunk_size;
    Kokkos::parallel_for (niter, findprimes_chunk (data, result, count, nnumbers, chunk_size));
    std::cout<<"timer chunk "<<timer.seconds()<<std::endl;

    printf ("Found %i prime numbers in %i random numbers\n", h_count(), nnumbers);
    Kokkos::finalize ();
}