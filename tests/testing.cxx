#include "gtest/gtest.h"

#include "needle.h"
#include "sparsediagonalmatrix.h"

typedef double Real;

typedef BendingNeedleModel<Real> BendingNeedleModelD;

TEST(NeedleTest, SimulationComparable)
{
  BendingNeedleModelD needle( 150.0, 30, 5000.0 );
  int last = needle.getX().size()-1;
  double s = needle.getSegmentLength();
  needle.setSpring( last, Vector( last * s, 10, 0 ), 0.5 );
  needle.setSpring(0,Vector( 0, 0, 0 ), 0.5);

  double error = 0;
  for(int i = 0; i < 1000; i++)
  {
    error = needle.simulateImplicitStatic(0.001);
  }
  ASSERT_NEAR( error, 0, 0.1 );
}

TEST(NeedleTest, MinimalWorkingExample)
{
  BendingNeedleModelD needle( 150.0, 31, 1000.0 );
  for(int i = 0; i < 1000; i++)
  {
    needle.simulateImplicitDynamic(0.001);
  }
}

TEST(NeedleTest, MinimalWorkingExampleXDynamic)
{
  BendingNeedleModelD needle( 150.0, 31, 1000.0 );
  double error = 0;
  for(int i = 0; i < 1000; i++)
  {
    error = needle.simulateImplicitDynamic(0.001);
  }
  ASSERT_LE( error, 10 );
}


TEST(NeedleTest, MinimalWorkingExampleZDynamic)
{
  BendingNeedleModelD needle( 150.0, 31, 1000.0 );
  needle.setBaseDirection(Vector(0,0,1));
  double error = 0;
  for(int i = 0; i < 1000; i++)
  {
    error = needle.simulateImplicitDynamic(0.001);
  }
  ASSERT_LE( error, 10 );
}

TEST(NeedleTest, MinimalWorkingExampleXStatic)
{
  BendingNeedleModelD needle( 150.0, 31, 1000.0 );
  double error = 0;
  for(int i = 0; i < 1000; i++)
  {
    error = needle.simulateImplicitStatic(0.001);
  }
  ASSERT_LE( error, 10 );
}


TEST(NeedleTest, MinimalWorkingExampleZStatic)
{
  BendingNeedleModelD needle( 150.0, 31, 1000.0 );
  needle.setBaseDirection(Vector(0,0,1));
  double error = 0;
  for(int i = 0; i < 1000; i++)
  {
    error = needle.simulateImplicitStatic(0.001);
  }
  ASSERT_LE( error, 10 );
}



TEST(NeedleTest, Stiff)
{
  BendingNeedleModelD needle( 150.0, 31, 10000.0 );
  for(int i = 0; i < 1000; i++)
  {
    needle.simulateImplicitDynamic(0.001);
  }
}
TEST(NeedleTest, ManyNodesDynamic)
{
  BendingNeedleModelD needle( 150.0, 51, 1000.0 );
  for(int i = 0; i < 1000; i++)
  {
    needle.simulateImplicitDynamic(0.001);
  }
}
TEST(NeedleTest, AllDynamic)
{
  BendingNeedleModelD needle( 150.0, 51, 10000.0 );
  int last = needle.getX().size()-1;
  double s = needle.getSegmentLength();
  needle.setSpring( last, Vector( last * s, 0, 0 ), 0.1 );
  needle.setSpring( last-1, Vector( (last-1) * s, 0, 0 ), 0.1 );
  needle.setSpring( last-2, Vector( (last-2) * s, 0, 0 ), 0.1 );

  for(int i = 0; i < 1000; i++)
  {
    needle.simulateImplicitDynamic(0.001);
  }
}

TEST(NeedleTest, ManyNodesStatic)
{
  BendingNeedleModelD needle( 150.0, 51, 1000.0 );
  for(int i = 0; i < 1000; i++)
  {
    needle.simulateImplicitStatic(0.001);
  }
}
TEST(NeedleTest, AllStatic)
{
  BendingNeedleModelD needle( 150.0, 51, 10000.0 );
  int last = needle.getX().size()-1;
  double s = needle.getSegmentLength();
  needle.setSpring( last, Vector( last * s, 0, 0 ), 0.1 );
  needle.setSpring( last-1, Vector( (last-1) * s, 0, 0 ), 0.1 );
  needle.setSpring( last-2, Vector( (last-2) * s, 0, 0 ), 0.1 );

  for(int i = 0; i < 1000; i++)
  {
    needle.simulateImplicitStatic(0.001);
  }
}

