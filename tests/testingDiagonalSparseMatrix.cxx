#include "gtest/gtest.h"

#include "needle.h"
#include "sparsediagonalmatrix.h"

int sum( cml::vectord b )
{
  int s = 0;
  for(int i = 0; i < b.size(); i++)
  {
    s += b[i];
  }
  return s;
}

TEST(SparseDiagonalMatrixTest, Test15)
{
  SparseDiagonalMatrix m(30,15);

  cml::vectord x(30);

  int j = 0;
  
  for(int i = 0; i < 30; i++) 
  {
    x[i] = 1; 
    m(i,i) = i;
    if(i>6) m(i-7,i) = i;
    if(i<30-7) m(i+7,i) = -i;
  }

  std::cout<<m<<std::endl;
  //std::cout<<x<<std::endl;
  cml::vectord b=m*x;
  //std::cout<<b<<std::endl;
  //std::cout<<sum(b)<<std::endl;
  ASSERT_EQ(sum(b), 596);
}

TEST(SparseDiagonalMatrixTest, Test17)
{
  SparseDiagonalMatrix m(30,17);

  cml::vectord x(30);

  int j = 0;
  
  for(int i = 0; i < 30; i++) 
  {
    x[i] = 1; 
    m(i,i) = i;
    if(i>6) m(i-7,i) = i;
    if(i<30-7) m(i+7,i) = -i;
  }

  std::cout<<m<<std::endl;
  //std::cout<<x<<std::endl;
  cml::vectord b=m*x;
  //std::cout<<b<<std::endl;
  //std::cout<<sum(b)<<std::endl;
  ASSERT_EQ(sum(b), 596);
}

TEST(SparseDiagonalMatrixTest, Test19)
{
  SparseDiagonalMatrix m(30,19);

  cml::vectord x(30);

  int j = 0;
  
  for(int i = 0; i < 30; i++) 
  {
    x[i] = 1; 
    m(i,i) = i;
    if(i>6) m(i-7,i) = i;
    if(i<30-7) m(i+7,i) = -i;
  }

  std::cout<<m<<std::endl;
  //std::cout<<x<<std::endl;
  cml::vectord b=m*x;
  //std::cout<<b<<std::endl;
  //std::cout<<sum(b)<<std::endl;
  ASSERT_EQ(sum(b), 596);
}


TEST(SparseDiagonalMatrixTest, MinimalWorkingExample)
{
  SparseDiagonalMatrix m(9,5);
  m(8,8) = 1;
  cml::vectord x(9);
  for(int i = 0; i < 9; i++) 
  { 
    m(i,i) = i+1;
    if(i>0) m(i-1,i) = i+1;
    if(i<8) m(i+1,i) = i+1;
  }
  //std::cout<<m<<std::endl;
  
  x[0] = 1;
  cml::vectord b=m*x;
  //std::cout<<b<<std::endl;
  ASSERT_EQ( sum(b), 3 );
  x[0] = 0;

  x[1] = 1;
  b=m*x;
  //std::cout<<b<<std::endl;
  ASSERT_EQ( sum(b), 6 );
  x[1] = 0;
  
  x[2] = 1;
  b=m*x;
  //std::cout<<b<<std::endl;
  ASSERT_EQ( sum(b), 9 );
  x[2] = 0;

  x[3] = 1;
  b=m*x;
  //std::cout<<b<<std::endl;
  ASSERT_EQ( sum(b), 12 );
  x[3] = 0;

  x[4] = 1;
  b=m*x;
  //std::cout<<b<<std::endl;
  ASSERT_EQ( sum(b), 15 );
  x[4] = 0;

  x[5] = 1;
  b=m*x;
  //std::cout<<b<<std::endl;
  ASSERT_EQ( sum(b), 18 );
  x[5] = 0;

  x[6] = 1;
  b=m*x;
  //std::cout<<b<<std::endl;
  ASSERT_EQ( sum(b), 21 );
  x[6] = 0;

  x[7] = 1;
  b=m*x;
  //std::cout<<b<<std::endl;
  ASSERT_EQ( sum(b), 24 );
  x[7] = 0;

  x[8] = 1;
  b=m*x;
  //std::cout<<b<<std::endl;
  ASSERT_EQ( sum(b), 17 );
  x[8] = 0;

  for(int i = 0; i < 9; i++) 
  { 
    x[i] = i+1;
  }


  m.zero();
  m(8,8) = 1;
  b=m*x;
  ASSERT_EQ(b[8], 9.0);
  m(8,8) = 0;
  m(7,8) = 1;
  b=m*x;
  ASSERT_EQ(b[8], 8.0);

}

TEST(SparseDiagonalMatrixTest, SecondTest)
{
  SparseDiagonalMatrix m(9,5);

  cml::vectord x(9);

  for(int i = 0; i < 9; i++) 
  {
    x[i] = 1; 
    m(i,i) = 1;
    if(i>0) m(i-1,i) = 1;
    if(i<8) m(i+1,i) = 1;
    if(i>1) m(i-2,i) = 1;
    if(i<7) m(i+2,i) = 1;
  }
  //std::cout<<m<<std::endl;
  //std::cout<<x<<std::endl;
  //cml::vectord b=m*x;
  //std::cout<<b<<std::endl;
}

