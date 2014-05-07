/**
 * Copyright (C) 2014 Dirk Fortmeier
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE
 *
 */

#include <iomanip>
#include <exception>

#include <emmintrin.h>
#include <pmmintrin.h>

#include "sparsediagonalmatrixOpt.h"

SparseDiagonalMatrixOpt::SparseDiagonalMatrixOpt( int m, int b ) :
  size(m),
  bandwidth(b),
  bandwidth_2(b/2)
{

  if( b < 3 || b % 2 != 1 ) throw std::runtime_error("bandwitdh is smaller than 3 or odd");

  values = new double[m*(bandwidth+1)];
  zero();

}

SparseDiagonalMatrixOpt::~SparseDiagonalMatrixOpt()
{
	delete[] values;
}

int SparseDiagonalMatrixOpt::getOffset( int y ) const
{
  int t1 =  (int)floor( bandwidth / 2.0 ) % 2;
  int t2 =  (int)ceil( bandwidth / 2.0 ) % 2;
  int a = floor( ( y + t1) / 2.0 ) * 2;
  int b = ceil( ( bandwidth - t2 ) / 2.0 );
  return a - b; 
}


void SparseDiagonalMatrixOpt::zero()
{
  for( int i = 0; i < bandwidth+1; i++ )
  {
    for( int j = 0; j < size; j++ )
    {
      _at(i,j) = 0;
    }
  }
}

std::ostream& operator<< ( std::ostream &out, const SparseDiagonalMatrixOpt &matrix )
{
  for( int j = 0; j < matrix.size; j++ )
  {
    int o = matrix.getOffset(j);
    out << std::setw( 4 );
    out << o << "[ ";
    for( int i = o; i > 0; i-- )
    {
      out << "   .";
    }
    int s = 0;
    if( o < 0) s = -o;
    for( int i = s; i < matrix.bandwidth+1; i++ )
    {
      out << std::setw( 4 ) << matrix._at(i,j);
    }

    out << " ]\n";
  }
  return out;
}

double& SparseDiagonalMatrixOpt::_at( int i, int j ) const
{
  return values[ i + (bandwidth+1)*j ];
}

double& SparseDiagonalMatrixOpt::operator() ( int i, int j ) const
{
  return _at(i-getOffset(j),j);
}

cml::vectord SparseDiagonalMatrixOpt::operator* (const cml::vectord& x) const
{
  cml::vectord r(x.size());
  for(int j = 0; j < size; j++)
  {
    r[j] = 0;
    int o = getOffset(j);
    for( int i = 0; i < bandwidth + 1; i++) 
    {
      int X = i + o;
      if( X >= 0 && X < size )
      {
        r[j] += _at(i,j) * x[X];
        //std::cout<<_at(i,j) * x[X]<<std::endl;
      }
    }
    //std::cout<<"="<<r[j]<<std::endl;
  }

  // for(int j = 0; j < bandwidth_2+1; j++ )
  // {
  //   r[j] = 0;
  //   for( int i = 0; i < bandwidth; i++ )
  //   {
  //     r[j] += (*this)(i,j) * x[i];
  //   }

  // }

  // for(int j = bandwidth_2+1; j < size - bandwidth_2-1; j++ )
  // {
  //   r[j] = 0;
  //   for( int i = 0; i < bandwidth; i++ )
  //   {
  //     r[j] += _at(i,j) * x[i+j-bandwidth_2];
  //   }

  // }

  // for(int j = size - bandwidth_2-1; j < size; j++ )
  // {
  //   r[j] = 0;
  //   for( int i = 0; i < bandwidth; i++ )
  //   {
  //     r[j] += _at(i,j) * x[i+size-bandwidth];
  //   }

  // }

  return r;
}

int SparseDiagonalMatrixOpt::getSize()
{
  return size;
}


int SparseDiagonalMatrixOpt::getBandwidth()
{
  return bandwidth;
}

void SparseDiagonalMatrixOpt::multiplyWith( const cml::vectord& x, cml::vectord& r ) const
{
  // first, check if size of matrix and vectors are right:
  if( size != x.size() )
    throw std::runtime_error("size of x must be the same as matrix size");

  if( size != r.size() )
    throw std::runtime_error("size of r must be the same as matrix size");


  if ((this->bandwidth+1) % 4 != 0 ) 
    throw std::runtime_error("fast multiplication only works with matrix with bandwidth is multiple of 4");
  
  // now, do the multiplication
  int o = getOffset(0);
  int k = o;

  for(int j = -o; j < size+o; j++)
  {
    r[j] = 0;
    for( int i = 0; i < bandwidth + 1; i+=4) 
    {
      int X = i + k;
      //r[j] += _at(i,j) * x[X];
      //std::cout<<_at(i,j) * x[X]<<std::endl;
            //int j2 = j+1;
      //std::cout<< " try " << i << " xoffset "<< (i+0+j-bandwidth_2) <<" check "<< (j-bandwidth_2) % 2 << std::endl;
      __m128d matrixPart1 = _mm_load_pd( &_at(i,j) );
      __m128d vectorPart1 = _mm_load_pd( &x[X] );
      //std::cout << "---" << std::endl;

      __m128d prod1 = _mm_mul_pd( matrixPart1, vectorPart1 );

      __m128d matrixPart2 = _mm_load_pd( &_at(i+2,j) );
      __m128d vectorPart2 = _mm_load_pd( &x[X+2] );
      //std::cout << "---" << std::endl;

      __m128d prod2 = _mm_mul_pd( matrixPart2, vectorPart2 );


      __m128d sum = _mm_add_pd( prod1, prod2 );

      __m128d red1 = _mm_hadd_pd( sum, sum );

      double res;
      _mm_store_sd( &res, red1 );

      r[j] += res;
      
    }
    if(j%2 == 0) k+=2;
    //std::cout<<"="<<r[j]<<std::endl;
  }

/*
  // do the start non-optimized
  for(int j = 0; j < bandwidth_2+1; j++ )
  {
    r[j] = 0;
    for( int i = 0; i < bandwidth; i++ )
    {
      r[j] += (*this)(i,j) * x[i];
    }

  }

  // do the middle part fully optimized
  for(int j = bandwidth_2+1; j < size - bandwidth_2-1; j+=2 )
  {
    r[j] = 0;

    //if((j-bandwidth_2) % 2 == 1)
    for( int i = 0; i < bandwidth; i+=4 )
    {
      //r[j] += _at(i,j) * x[i+j-bandwidth_2];
       r[j] += _at(i+0,j) * x[i+0+j-bandwidth_2];
       r[j] += _at(i+1,j) * x[i+1+j-bandwidth_2];
       r[j] += _at(i+2,j) * x[i+2+j-bandwidth_2];
       r[j] += _at(i+3,j) * x[i+3+j-bandwidth_2];
    }
    //else
    for( int i = 0; i < bandwidth; i+=4 )
    {
      int j2 = j+1;
      //std::cout<< " try " << i << " xoffset "<< (i+0+j-bandwidth_2) <<" check "<< (j-bandwidth_2) % 2 << std::endl;
      __m128d matrixPart1 = _mm_load_pd( &_at(i+0,j2) );
      __m128d vectorPart1 = _mm_load_pd( &x[i+0+j2-bandwidth_2] );
      //std::cout << "---" << std::endl;

      __m128d prod1 = _mm_mul_pd( matrixPart1, vectorPart1 );

      //std::cout<< " try " << i << " + 2" << std::endl;
      __m128d matrixPart2 = _mm_load_pd( &_at(i+2,j2) );
      __m128d vectorPart2 = _mm_load_pd( &x[i+2+j2-bandwidth_2] );


      __m128d prod2 = _mm_mul_pd( matrixPart2, vectorPart2 );

      __m128d sum = _mm_add_pd( prod1, prod2 );

      __m128d red1 = _mm_hadd_pd( sum, sum );

      //__m128d red2 = _mm_hadd_pd( red1, red1 );

      double res;
      _mm_store_sd( &res, red1 );

      r[j2] += res;
      // r[j] += _at(i+0,j) * x[i+0+j-bandwidth_2];
      // r[j] += _at(i+1,j) * x[i+1+j-bandwidth_2];
      // r[j] += _at(i+2,j) * x[i+2+j-bandwidth_2];
      // r[j] += _at(i+3,j) * x[i+3+j-bandwidth_2];
    }
  }

  // do the end part non-optimized as well
  for(int j = size - bandwidth_2-1; j < size; j++ )
  {
    r[j] = 0;
    for( int i = 0; i < bandwidth; i++ )
    {
      r[j] += _at(i,j) * x[i+size-bandwidth];
    }

  }
*/
}