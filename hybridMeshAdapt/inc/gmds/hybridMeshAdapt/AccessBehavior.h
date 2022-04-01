#ifndef ACCESS_BEHAVIOR_H_
#define ACCESS_BEHAVIOR_H_
/*----------------------------------------------------------------------------*/
#include<Eigen/Sparse>
#include<gmds/utils/CommonTypes.h>
#include<memory>
#include<set>
#include<map>
/*----------------------------------------------------------------------------*/
namespace gmds
{
  namespace accessBehavior
  {
    class IBehavior
    {
    public:
      IBehavior                                 (){}

      ~IBehavior                                (){}

      virtual const TInt neighborCellOf         (const TInt pointIndex, const TInt cellIndex) = 0;

      virtual std::vector<TInt> ballOf          (const TInt pointIndex) = 0;

      virtual void buildDatabehavior                 (/*const gmds::Mesh mesh*/ ) = 0;

    };




    class SparseBehavior : public IBehavior
    {
    public:
      SparseBehavior              (std::unique_ptr<Eigen::SparseMatrix<TInt>>&& sparseIndex);

      SparseBehavior              ( SparseBehavior& sparseBehavior);

      SparseBehavior              (SparseBehavior&& sparseBehavior);

      ~SparseBehavior             ();

      virtual const TInt neighborCellOf   (const TInt pointIndex, const TInt cellIndex);

      virtual std::vector<TInt> ballOf    (const TInt pointIndex);

      virtual void buildDatabehavior        (/*const gmds::Mesh mesh*/ ){}

    private:
      //a voir pour mettre juste Eigen::SparseMatrix<TInt> sans le unique_ptr
      std::unique_ptr<Eigen::SparseMatrix<TInt>> m_sparseIndex;
    };





    class MapBehavior : public IBehavior
    {

    public:

      using N2C = std::map<TInt, std::set<TInt>>;

      MapBehavior              ();
      //for the debug in order to give directly the useful data (the map) whereas building it directly with build.../
      MapBehavior              (MapBehavior&& mapBehavior);

      MapBehavior              (const MapBehavior& mapBehavior);

      MapBehavior              (std::unique_ptr<N2C>&& mapBehavior);

      ~MapBehavior             ();



      virtual const TInt neighborCellOf   (const TInt pointIndex, const TInt cellIndex);

      virtual std::vector<TInt> ballOf    (const TInt pointIndex);

      virtual void buildDatabehavior      (/*const gmds::Mesh mesh*/ );



    private:

      //function that research the common index in few set of the map without the previous cell
      TInt commonIdxInSet();
      // a voir pour le remplacer par un non pointer...
      std::unique_ptr<N2C> m_N2C;

      //sans unique_ptr
      std::vector<TInt> m_base;

    };



  }
}




#endif //ACCESS_BEHAVIOR_H_
