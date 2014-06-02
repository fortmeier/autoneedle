#include "gtest/gtest.h"

#include "needle.h"
#include "sparsediagonalmatrixOpt.h"
#include "sparsediagonalmatrix.h"

typedef double Real;

typedef BendingNeedleModel<Real> BendingNeedleModelD;
typedef SparseDiagonalMatrix<Real> SparseDiagonalMatrixD;

typedef SparseDiagonalMatrixOpt<Real> SparseDiagonalMatrixDOpt;
typedef BandMatrixInterface<Real> BandMatrixInterfaceD;

void speedTest( BandMatrixInterfaceD& m, int reps, int mode )
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

TEST(SpeedTest, SparseDiagonalMatrixDOpt)
{
  SparseDiagonalMatrixDOpt m(150,19);
  speedTest( m, 100000, 0 );
}

TEST(SpeedTest, SparseDiagonalMatrixDOptIntrinsics)
{
  SparseDiagonalMatrixDOpt m(150,19);
  speedTest( m, 100000, 1 );
}

TEST(SpeedTest, SparseDiagonalMatrixD)
{
  SparseDiagonalMatrixD m(150,19);
  speedTest( m, 100000, 0 );
}




void speedTest2( int size, int reps, int mode )
{
  SparseDiagonalMatrixDOpt m(size,19);

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

TEST(SparseDiagonalMatrixDOptTest, SpeedTest1)
{
  speedTest2( 150, 100000, 0);
}

TEST(SparseDiagonalMatrixDOptTest, SpeedTest2)
{
  speedTest2( 150, 100000, 1);
}