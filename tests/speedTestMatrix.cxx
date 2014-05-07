#include "gtest/gtest.h"

#include "needle.h"
#include "sparsediagonalmatrixOpt.h"
#include "sparsediagonalmatrix.h"


void speedTest( BandMatrixInterface& m, int reps, int mode )
{
  int size = m.getSize();
  cml::vectord x(size);

  for(int i = 0; i < size; i++) 
  {
    x[i] = 1; 
    m(i,i) = i;
  }
  
  cml::vectord r(size);

  for(int i = 0; i < reps; i++)
  {
    switch(mode)
    {
      case 0:
        r = m * x;
        break;
      case 1:
        m.multiplyWith( x, r );
        break;
    }
  }
}

TEST(SpeedTest, SparseDiagonalMatrixOpt)
{
  SparseDiagonalMatrixOpt m(150,19);
  speedTest( m, 100000, 0 );
}

TEST(SpeedTest, SparseDiagonalMatrixOptIntrinsics)
{
  SparseDiagonalMatrixOpt m(150,19);
  speedTest( m, 100000, 1 );
}

TEST(SpeedTest, SparseDiagonalMatrix)
{
  SparseDiagonalMatrix m(150,19);
  speedTest( m, 100000, 0 );
}
