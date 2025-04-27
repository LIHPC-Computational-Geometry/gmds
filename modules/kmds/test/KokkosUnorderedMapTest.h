/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <Kokkos_Core.hpp>
#include <Kokkos_UnorderedMap.hpp>
/*----------------------------------------------------------------------------*/
#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class KokkosUnorderedMapTest : public ::testing::Test
{
 protected:
    KokkosUnorderedMapTest()
        {
                ;
        }
        virtual ~KokkosUnorderedMapTest()
        {
                ;
        }

    static void
    SetUpTestCase()
    {
    }

    static void
    TearDownTestCase()
    {
    }
};
/*----------------------------------------------------------------------------*/
TEST_F(KokkosUnorderedMapTest, KokkosMap_0)
{
        Kokkos::UnorderedMap<kmds::TCellID, double> kmap(10);

        Kokkos::UnorderedMapInsertResult res;
        res = kmap.insert(0, 0.);
        bool success = res.success();
        bool fail = res.failed();

        EXPECT_EQ(true, success);
        EXPECT_EQ(false, fail);
        EXPECT_EQ(1, kmap.size());
}
/*----------------------------------------------------------------------------*/
TEST_F(KokkosUnorderedMapTest, KokkosMap_1)
{
        Kokkos::UnorderedMap<kmds::TCellID, double> kmap(3);

        bool success;
        bool fail;

        Kokkos::UnorderedMapInsertResult res;
        for(int i=0; i<3; i++) {
                 res = kmap.insert(i, 0.7);
                 success = res.success();
                 fail = res.failed();
                 EXPECT_EQ(true, success);
                 EXPECT_EQ(false, fail);
                 EXPECT_EQ(i+1, kmap.size());
        }

}
/*----------------------------------------------------------------------------*/
TEST_F(KokkosUnorderedMapTest, KokkosMap_2)
{
    Kokkos::UnorderedMap<kmds::TCellID, double> kmap(100);

    Kokkos::parallel_for(100, KOKKOS_LAMBDA(const int i) {
        Kokkos::UnorderedMapInsertResult res = kmap.insert(i, i);
        bool success = res.success();
        bool fail = res.failed();

        EXPECT_EQ(true, success);
        EXPECT_EQ(false, fail);
    });

    EXPECT_EQ(100, kmap.size());

}
/*----------------------------------------------------------------------------*/
TEST_F(KokkosUnorderedMapTest, KokkosMap_3)
{
    Kokkos::UnorderedMap<kmds::TCellID, double> kmap(100);

    Kokkos::parallel_for(100, KOKKOS_LAMBDA(const int i) {

        int toinsert = i;

        if(i == 17) {
            toinsert = 2;
        }
        if(i == 27) {
            toinsert = 2;
        }
        if(i == 30) {
            toinsert = 45;
        }
        if(i == 80) {
            toinsert = 98;
        }
        if(i == 99) {
            toinsert = 67;
        }

        Kokkos::UnorderedMapInsertResult res = kmap.insert(toinsert, i);
        bool success = res.success();
        bool exist = res.existing();
        bool fail = res.failed();


        EXPECT_EQ(true, success || exist);
        EXPECT_EQ(false, fail);
    });

    EXPECT_EQ(95, kmap.size());
}
/*----------------------------------------------------------------------------*/
TEST_F(KokkosUnorderedMapTest, KokkosMap_4)
{
    Kokkos::UnorderedMap<kmds::TCellID, double> kmap(100);

    EXPECT_FALSE(kmap.is_set);

    int result = 0;

    Kokkos::parallel_reduce(100, KOKKOS_LAMBDA(const int i, int &sum) {

                                int toinsert = i;

                                if(i == 17) {
                                    toinsert = 2;
                                }
                                if(i == 27) {
                                    toinsert = 2;
                                }
                                if(i == 30) {
                                    toinsert = 45;
                                }
                                if(i == 80) {
                                    toinsert = 98;
                                }
                                if(i == 99) {
                                    toinsert = 67;
                                }

                                Kokkos::UnorderedMapInsertResult res = kmap.insert(toinsert, i);
                                bool success = res.success();
                                bool exist = res.existing();
                                bool fail = res.failed();


                                EXPECT_EQ(true, success || exist);
                                EXPECT_EQ(false, fail);

                                if(exist) {
                                    sum++;
                                }
                            },
                            result
    );

    EXPECT_EQ(95, kmap.size());
    EXPECT_EQ(5, result);


    int keys_sum = 0;
    // This does not work, accessing the Kokkos::UnorderedMap is index-based up to capacity ?
//    Kokkos::parallel_reduce(kmap.size(),
//                            KOKKOS_LAMBDA(const int i, int& ls) {
//                                kmds::TCellID id = kmap.key_at(i);
//                                ls += id;
//                                        }
//                                ,
//                                result_bis);
    Kokkos::parallel_reduce(kmap.capacity(),
                            KOKKOS_LAMBDA(const int i, int& ls) {
                                if(kmap.valid_at(i)) {
                                    kmds::TCellID id = kmap.key_at(i);
                                    ls += id;
                                }
                            }
            ,
                            keys_sum);

    EXPECT_EQ(4900+50-17-27-30-80-99, keys_sum);
}
/*----------------------------------------------------------------------------*/
TEST_F(KokkosUnorderedMapTest, KokkosMap_5)
{
    Kokkos::UnorderedMap<kmds::TCellID, void> kmap(100);

    EXPECT_TRUE((bool) kmap.is_set);

    int result = 0;

    Kokkos::parallel_reduce(100, KOKKOS_LAMBDA(const int i, int &sum) {

                                int toinsert = i;

                                if(i == 17) {
                                    toinsert = 2;
                                }
                                if(i == 27) {
                                    toinsert = 2;
                                }
                                if(i == 30) {
                                    toinsert = 45;
                                }
                                if(i == 80) {
                                    toinsert = 98;
                                }
                                if(i == 99) {
                                    toinsert = 67;
                                }

                                Kokkos::UnorderedMapInsertResult res = kmap.insert(toinsert, i);
                                bool success = res.success();
                                bool exist = res.existing();
                                bool fail = res.failed();


                                EXPECT_EQ(true, success || exist);
                                EXPECT_EQ(false, fail);

                                if(exist) {
                                    sum++;
                                }
                            },
                            result
    );

    EXPECT_EQ(95, kmap.size());
    EXPECT_EQ(5, result);

    int result_bis = 0;
    // This does not work, accessing the Kokkos::UnorderedMap is index-based up to capacity ?
//    Kokkos::parallel_reduce(kmap.size(),
//                            KOKKOS_LAMBDA(const int i, int& ls) {
//                                kmds::TCellID id = kmap.key_at(i);
//                                ls += id;
//                                        }
//                                ,
//                                result_bis);
    Kokkos::parallel_reduce(kmap.capacity(),
                            KOKKOS_LAMBDA(const int i, int& ls) {
                                if(kmap.valid_at(i)) {
                                    kmds::TCellID id = kmap.key_at(i);
                                    ls += id;
                                }
                            }
            ,
                            result_bis);

    EXPECT_EQ(4900+50-17-27-30-80-99, result_bis);
}
/*----------------------------------------------------------------------------*/
TEST_F(KokkosUnorderedMapTest, KokkosMap_6) {
    Kokkos::UnorderedMap<kmds::TCellID, double> kmap(100);

    EXPECT_FALSE(kmap.is_set);

    int result = 0;

    Kokkos::parallel_reduce(100, KOKKOS_LAMBDA(const int i, int &sum) {

                                Kokkos::UnorderedMapInsertResult res = kmap.insert(i, i + 1);
                                bool success = res.success();
                                bool exist = res.existing();
                                bool fail = res.failed();


                                EXPECT_EQ(true, success || exist);
                                EXPECT_EQ(false, fail);

                                if (exist) {
                                    sum++;
                                }
                            },
                            result
    );

    int keys_sum = 0;
    Kokkos::parallel_reduce(kmap.capacity(),
                            KOKKOS_LAMBDA(const int i, int &ls) {
                                if (kmap.valid_at(i)) {
                                    kmds::TCellID id = kmap.key_at(i);
                                    ls += id;
                                }
                            },
                            keys_sum);

    EXPECT_EQ(4900 + 50, keys_sum);

    int values_sum = 0;
    Kokkos::parallel_reduce(kmap.capacity(),
                            KOKKOS_LAMBDA(const int i, int& ls) {
                                if(kmap.valid_at(i)) {
                                    ls += kmap.value_at(i);
                                    kmap.value_at(i) -= 2;
                                }
                            }
            ,
                            values_sum);

    EXPECT_EQ(4900+50 + 100, values_sum);

    int values_newsum = 0;
    Kokkos::parallel_reduce(kmap.capacity(),
                            KOKKOS_LAMBDA(const int i, int& ls) {
                                if(kmap.valid_at(i)) {
                                    ls += kmap.value_at(i);
                                }
                            }
            ,
                            values_newsum);

    EXPECT_EQ(4900+50 + 100 - 200, values_newsum);
}
/*----------------------------------------------------------------------------*/
TEST_F(KokkosUnorderedMapTest, KokkosMap_7) {
    Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> kmap(200000);

    EXPECT_FALSE(kmap.is_set);

    Kokkos::parallel_for(200000, KOKKOS_LAMBDA(const int i) {

        int toinsert = i;

        if(i == 17) {
            toinsert = 2;
        }
        if(i == 27) {
            toinsert = 2;
        }
        if(i == 30) {
            toinsert = 45;
        }
        if(i == 80) {
            toinsert = 98;
        }
        if(i == 99) {
            toinsert = 67;
        }

        Kokkos::UnorderedMapInsertResult res = kmap.insert(toinsert, i);
        bool success = res.success();
        bool exist = res.existing();
        bool fail = res.failed();


        EXPECT_EQ(true, success || exist);
        EXPECT_EQ(false, fail);
    });

    EXPECT_EQ(200000 - 5, kmap.size());
}
/*----------------------------------------------------------------------------*/