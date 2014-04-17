#include "gtest/gtest.h"

#include "needle.h"

TEST(NeedleTest, MinimalWorkingExample)
{
  BendingNeedleModel needle;
  for(int i = 0; i < 10000; i++)
  {
    needle.simulateImplicitChentanez(0.01);
  }
}
