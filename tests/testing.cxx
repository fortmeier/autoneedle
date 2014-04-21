#include "gtest/gtest.h"

#include "needle.h"

TEST(NeedleTest, MinimalWorkingExample)
{
  BendingNeedleModel needle( 150.0, 31, 1000.0 );
  for(int i = 0; i < 1000; i++)
  {
    needle.simulateImplicitDynamic(0.001);
  }
}
TEST(NeedleTest, Stiff)
{
  BendingNeedleModel needle( 150.0, 31, 10000.0 );
  for(int i = 0; i < 1000; i++)
  {
    needle.simulateImplicitDynamic(0.001);
  }
}
TEST(NeedleTest, ManyNodesDynamic)
{
  BendingNeedleModel needle( 150.0, 51, 1000.0 );
  for(int i = 0; i < 1000; i++)
  {
    needle.simulateImplicitDynamic(0.001);
  }
}
TEST(NeedleTest, AllDynamic)
{
  BendingNeedleModel needle( 150.0, 51, 10000.0 );
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
  BendingNeedleModel needle( 150.0, 51, 1000.0 );
  for(int i = 0; i < 1000; i++)
  {
    needle.simulateImplicitStatic(0.001);
  }
}
TEST(NeedleTest, AllStatic)
{
  BendingNeedleModel needle( 150.0, 51, 10000.0 );
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
