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
  //cout<<sum(b)<<endl;
  ASSERT_EQ(sum(b), 596);
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



void speedTest( int size, int reps, int mode )
{
  SparseDiagonalMatrixOpt m(size,19);

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

TEST(SparseDiagonalMatrixOptTest, SpeedTest1)
{
  speedTest( 150, 100000, 0);
}

TEST(SparseDiagonalMatrixOptTest, SpeedTest2)
{
  speedTest( 150, 100000, 1);
}

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
  SparseDiagonalMatrixOpt m(size,17);

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
