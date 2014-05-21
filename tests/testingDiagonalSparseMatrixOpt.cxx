#include "gtest/gtest.h"

#include "needle.h"
#include "sparsediagonalmatrixOpt.h"

extern int sum( cml::vectord b );

using namespace std;

TEST(SparseDiagonalMatrixOptTest, Test15)
{
  SparseDiagonalMatrixOpt m(30,15);

  cml::vectord x(30);

  int j = 0;
  
  for(int i = 0; i < 30; i++) 
  {
    x[i] = 1; 
    m(i,i) = i;
    if(i>6) m(i-7,i) = i;
    if(i<30-7) m(i+7,i) = -i;
  }

  cout<<m<<endl;
  cout<<x<<endl;
  cml::vectord b=m*x;
  cout<<b<<endl;
  cout<<sum(b)<<endl;
  ASSERT_EQ(sum(b), 596);
}

TEST(SparseDiagonalMatrixOptTest, Test15b)
{
  SparseDiagonalMatrixOpt m(31,15);

  cml::vectord x(31);

  int j = 0;
  
  for(int i = 0; i < 31; i++) 
  {
    x[i] = 1; 
    m(i,i) = i;
    if(i>6) m(i-7,i) = i;
    if(i<30-7) m(i+7,i) = -i;
  }

  cout<<m<<endl;
  cout<<x<<endl;
  cml::vectord b=m*x;
  cout<<b<<endl;
  cout<<sum(b)<<endl;
  ASSERT_EQ(sum(b), 656);
}

TEST(SparseDiagonalMatrixOptTest, Test17)
{
  SparseDiagonalMatrixOpt m(30,17);

  cml::vectord x(30);

  int j = 0;
  
  for(int i = 0; i < 30; i++) 
  {
    x[i] = 1; 
    m(i,i) = i;
    if(i>6) m(i-7,i) = i;
    if(i<30-7) m(i+7,i) = -i;
  }

  cout<<m<<endl;
  cout<<x<<endl;
  cml::vectord b=m*x;
  cout<<b<<endl;
  //cout<<sum(b)<<endl;
  ASSERT_EQ(sum(b), 596);
}

TEST(SparseDiagonalMatrixOptTest, Test17b)
{
  SparseDiagonalMatrixOpt m(29,17);

  cml::vectord x(29);

  int j = 0;
  
  for(int i = 0; i < 29; i++) 
  {
    x[i] = 1; 
    m(i,i) = i;
    if(i>6) m(i-7,i) = i;
    if(i<29-7) m(i+7,i) = -i;
  }

  cout<<m<<endl;
  cout<<x<<endl;
  cml::vectord b=m*x;
  cout<<b<<endl;
  cout<<sum(b)<<endl;
  ASSERT_EQ(sum(b), 560);
}

// TEST(SparseDiagonalMatrixOptTest, Test19)
// {
//   SparseDiagonalMatrixOpt m(30,19);

//   cml::vectord x(30);

//   int j = 0;
  
//   for(int i = 0; i < 30; i++) 
//   {
//     x[i] = 1; 
//     m(i,i) = i;
//     if(i>6) m(i-7,i) = i;
//     if(i<30-7) m(i+7,i) = -i;
//   }

//   cout<<m<<endl;
//   //cout<<x<<endl;
//   cml::vectord b=m*x;
//   //cout<<b<<endl;
//   //cout<<sum(b)<<endl;
//   ASSERT_EQ(sum(b), 596);
// }

// TEST(SparseDiagonalMatrixOptTest, Test21)
// {
//   SparseDiagonalMatrixOpt m(30,21);

//   cml::vectord x(30);

//   int j = 0;
  
//   for(int i = 0; i < 30; i++) 
//   {
//     x[i] = 1; 
//     m(i,i) = i;
//     if(i>6) m(i-7,i) = i;
//     if(i<30-7) m(i+7,i) = -i;
//   }

//   cout<<m<<endl;
//   //cout<<x<<endl;
//   cml::vectord b=m*x;
//   //cout<<b<<endl;
//   //cout<<sum(b)<<endl;
//   ASSERT_EQ(sum(b), 596);
// }

// TEST(SparseDiagonalMatrixOptTest, Test19b)
// {
//   SparseDiagonalMatrixOpt m(30,19);

//   cml::vectord x(30);

  
//   int summ = 0;
//   for(int i = 0; i < 30; i++) 
//   {
//     x[i] = i;
//     for(int j = 0; j < 19; j++)
//     {
//       m._at(j,i) = j;
//       summ += j * i;
//     }
//   }

//   cout<<m<<endl;
//   //cout<<x<<endl;
//   cml::vectord b=m*x;
//   cout<<b<<endl;
//   //cout<<sum(b)<<endl;
//   ASSERT_EQ(sum(b), summ);
// }


// TEST(SparseDiagonalMatrixOptTest, MinimalWorkingExample)
// {
//   SparseDiagonalMatrixOpt m(9,5);
//   m(8,8) = 1;
//   cml::vectord x(9);
//   for(int i = 0; i < 9; i++) 
//   { 
//     m(i,i) = i+1;
//     if(i>0) m(i-1,i) = i+1;
//     if(i<8) m(i+1,i) = i+1;
//   }
//   //cout<<m<<endl;
  
//   x[0] = 1;
//   cml::vectord b=m*x;
//   //cout<<b<<endl;
//   ASSERT_EQ( sum(b), 3 );
//   x[0] = 0;

//   x[1] = 1;
//   b=m*x;
//   //cout<<b<<endl;
//   ASSERT_EQ( sum(b), 6 );
//   x[1] = 0;
  
//   x[2] = 1;
//   b=m*x;
//   //cout<<b<<endl;
//   ASSERT_EQ( sum(b), 9 );
//   x[2] = 0;

//   x[3] = 1;
//   b=m*x;
//   //cout<<b<<endl;
//   ASSERT_EQ( sum(b), 12 );
//   x[3] = 0;

//   x[4] = 1;
//   b=m*x;
//   //cout<<b<<endl;
//   ASSERT_EQ( sum(b), 15 );
//   x[4] = 0;

//   x[5] = 1;
//   b=m*x;
//   //cout<<b<<endl;
//   ASSERT_EQ( sum(b), 18 );
//   x[5] = 0;

//   x[6] = 1;
//   b=m*x;
//   //cout<<b<<endl;
//   ASSERT_EQ( sum(b), 21 );
//   x[6] = 0;

//   x[7] = 1;
//   b=m*x;
//   //cout<<b<<endl;
//   ASSERT_EQ( sum(b), 24 );
//   x[7] = 0;

//   x[8] = 1;
//   b=m*x;
//   //cout<<b<<endl;
//   ASSERT_EQ( sum(b), 17 );
//   x[8] = 0;

//   for(int i = 0; i < 9; i++) 
//   { 
//     x[i] = i+1;
//   }


//   m.zero();
//   m(8,8) = 1;
//   b=m*x;
//   ASSERT_EQ(b[8], 9.0);
//   m(8,8) = 0;
//   m(7,8) = 1;
//   b=m*x;
//   ASSERT_EQ(b[8], 8.0);

// }

// TEST(SparseDiagonalMatrixOptTest, SecondTest)
// {
//   SparseDiagonalMatrixOpt m(9,5);

//   cml::vectord x(9);

//   for(int i = 0; i < 9; i++) 
//   {
//     x[i] = 1; 
//     m(i,i) = 1;
//     if(i>0) m(i-1,i) = 1;
//     if(i<8) m(i+1,i) = 1;
//     if(i>1) m(i-2,i) = 1;
//     if(i<7) m(i+2,i) = 1;
//   }
//   //cout<<m<<endl;
//   //cout<<x<<endl;
//   //cml::vectord b=m*x;
//   //cout<<b<<endl;
// }

TEST(SparseDiagonalMatrixOptTest, TestMultiplicationEquality)
{
  int size = 150;
  SparseDiagonalMatrixOpt m(size,19);

  cml::vectord x(size);

  for(int i = 0; i < size; i++) 
  {
    x[i] = 1; 
    m(i,i) = i;
  }
  
  cml::vectord r1(size);
  cml::vectord r2(size);

  r1 = m * x;

  m.multiplyWith( x, r2 );

  ASSERT_EQ(r1, r2);
}

TEST(SparseDiagonalMatrixOptTest, TestMultiplicationEquality2)
{
  int size = 150;
  SparseDiagonalMatrixOpt m(size,15);

  cml::vectord x(size);

  for(int i = 0; i < size; i++) 
  {
    x[i] = 1; 
    m(i,i) = i;
  }
  
  cml::vectord r1(size);
  cml::vectord r2(size);

  r1 = m * x;

  m.multiplyWith( x, r2 );

  ASSERT_EQ(r1, r2);
}

TEST(Bandmatrices, TestSumEquality1)
{
  int size = 30;
  SparseDiagonalMatrixOpt m1(size,15);
  SparseDiagonalMatrix m2(size,15);

  for(int i = 0; i < size; i++) 
  {
    m1(i,i) = i;
    m2(i,i) = i;
  }
  ASSERT_EQ(m1.sum(), m2.sum());

}


void subMatrixAdd3x3( BandMatrixInterface& m, int i, int j, int s = 0 )
{
  for(int x = 0; x < 3; x++ )
  {
    for( int y = 0; y < 3; y++ )
    {
      m(i+x,j+y) = ++s;
    }
  }
}

void setMatrix( BandMatrixInterface& m, int n )
{
  for(int i = 0; i < n; i++)
  {
    // pF_pxi
    if(i>=2) subMatrixAdd3x3(m, i*3, i*3 ); 
    if(i>=1 && i < n-1) subMatrixAdd3x3(m, i*3, i*3);
    if(i>=0 && i< n-2) subMatrixAdd3x3(m, i*3, i*3);

    // pF_pxi+1
    if(i>=1 && i < n-1) subMatrixAdd3x3(m, i*3, i*3+3);
    if(i>=0 && i < n-2) subMatrixAdd3x3(m, i*3, i*3+3);

    // pF_pxi+2
    if(i>=0 && i < n-2) subMatrixAdd3x3(m, i*3, i*3+6);

    // pF_pxi-1
    if(i>=1 && i < n-1) subMatrixAdd3x3(m, i*3, i*3-3);
    if(i>=2 && i < n) subMatrixAdd3x3(m, i*3, i*3-3);

    // pF_pxi+2
    if(i>=2 && i < n) subMatrixAdd3x3(m, i*3, i*3-6);
  }
}

TEST(Bandmatrices, TestSumEquality2)
{
  int nodes = 10;
  int size = nodes * 3;
  SparseDiagonalMatrixOpt m1(size,19);
  SparseDiagonalMatrix m2(size,19);

  setMatrix( m1, nodes );
  setMatrix( m2, nodes );
  cout << m1 << endl;
  cout << m2 << endl;

  ASSERT_EQ(m1.sum(), m2.sum());


  cml::vectord x(size);

  for(int i = 0; i < size; i++) 
  {
    x[i] = i+1 + 0.33; 
  }

  cml::vectord b1 = m1 * x;
  cml::vectord b2 = m2 * x;

  cout << b1 << endl;
  cout << b2 << endl;

  ASSERT_EQ(b1, b2);

  cml::vectord b3(x.size());
  cml::vectord b4(x.size());
  m1.multiplyWith(x, b3);
  m2.multiplyWith(x, b4);

  cout << b3 << endl;
  cout << b4 << endl;

  ASSERT_NEAR(sum(b3), sum(b4), 0.00001);
}

TEST(Bandmatrices, TestSumEquality3)
{
  int nodes = 30;
  int size = nodes * 3;
  SparseDiagonalMatrixOpt m1(size,19);
  SparseDiagonalMatrix m2(size,19);

  setMatrix( m1, nodes );
  setMatrix( m2, nodes );
  cout << m1 << endl;
  cout << m2 << endl;

  ASSERT_EQ(m1.sum(), m2.sum());


  cml::vectord x(size);

  for(int i = 0; i < size; i++) 
  {
    x[i] = i+1 + 0.33; 
  }

  cml::vectord b1 = m1 * x;
  cml::vectord b2 = m2 * x;

  cout << b1 << endl;
  cout << b2 << endl;
  ASSERT_EQ(b1, b2);

  cml::vectord b3(x.size());
  cml::vectord b4(x.size());
  m1.multiplyWith(x, b3);
  m2.multiplyWith(x, b4);

  cout << b3 << endl;
  cout << b4 << endl;

  ASSERT_NEAR(sum(b3), sum(b4), 0.00001);
}

// TEST(SparseDiagonalMatrixOptTest, Test3)
// {
//   SparseDiagonalMatrixOpt m(10,3);


//   for(int i = 0; i < 10; i++) 
//   {
//     m(i,i) = i+1;
//   }

//   cout<<m<<endl;
// }

// TEST(SparseDiagonalMatrixOptTest, Test5)
// {
//   SparseDiagonalMatrixOpt m(10,5);


//   for(int i = 0; i < 10; i++) 
//   {
//     m(i,i) = i+1;
//   }

//   cout<<m<<endl;
// }

// TEST(SparseDiagonalMatrixOptTest, Test7)
// {
//   SparseDiagonalMatrixOpt m(10,7);


//   for(int i = 0; i < 10; i++) 
//   {
//     m(i,i) = i+1;
//   }

//   cout<<m<<endl;
// }
