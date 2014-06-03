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

template<>
SparseDiagonalMatrixOpt<double>::SparseDiagonalMatrixOpt( int m, int b ) :
  size(m),
  bandwidth(b),
  bandwidth_2(b/2)
{

  if( b < 3 || b % 2 != 1 ) throw std::runtime_error("bandwidth is smaller than 3 or even");

  bandwidth_full = bandwidth + 1;
  values = new double[m*(bandwidth_full)];
  zero();

}

template<>
SparseDiagonalMatrixOpt<float>::SparseDiagonalMatrixOpt( int m, int b ) :
  size(m),
  bandwidth(b),
  bandwidth_2(b/2)
{

  if( b < 3 || b % 2 != 1 ) throw std::runtime_error("bandwidth is smaller than 3 or even");

  if( b+1 % 4 == 0 )
  {
    bandwidth_full = bandwidth + 1;
  }
  else
  {
    bandwidth_full = bandwidth + 3;
  }

  values = new float[m*(bandwidth_full)];

  zero();

}


template<typename Real>
SparseDiagonalMatrixOpt<Real>::~SparseDiagonalMatrixOpt()
{
	delete[] values;
}

template<typename Real>
int SparseDiagonalMatrixOpt<Real>::getOffset( int y ) const
{
  int t1 =  (int)floor( bandwidth / 2.0 ) % 2;
  int t2 =  (int)ceil( bandwidth / 2.0 ) % 2;
  int a = floor( ( y + t1) / 2.0 ) * 2;
  int b = ceil( ( bandwidth - t2 ) / 2.0 );
  return a - b; 
}

template<typename Real>
void SparseDiagonalMatrixOpt<Real>::zero()
{
  for( int i = 0; i < bandwidth_full; i++ )
  {
    for( int j = 0; j < size; j++ )
    {
      _at(i,j) = 0;
    }
  }
}

template<typename Real>
void SparseDiagonalMatrixOpt<Real>::zeroRow( int j )
{
  for( int i = 0; i < bandwidth_full; i++ )
  {
    _at(i,j) = 0;
  }
}

template<typename Real>
std::ostream& SparseDiagonalMatrixOpt<Real>::print ( std::ostream &out ) const
{
  out << std::endl;
  for( int j = 0; j < size; j++ )
  {
    int o = getOffset(j);
    out << std::setw( 5 );
    out << o << "[ ";
    for( int i = o; i > 0; i-- )
    {
      out << "    .";
    }
    int s = 0;
    if( o < 0) s = -o;
    for( int i = s; i < bandwidth_full; i++ )
    {
      out << std::setw( 5 ) << _at(i,j);
    }

    out << " ]\n";
  }
  return out;
}

template<typename Real>
Real& SparseDiagonalMatrixOpt<Real>::_at( int i, int j ) const
{
  return values[ i + (bandwidth_full)*j ];
}

template<typename Real>
Real& SparseDiagonalMatrixOpt<Real>::operator() ( int i, int j ) const
{
  return _at(i-getOffset(j),j);
}

template<typename Real>
VectorDyn SparseDiagonalMatrixOpt<Real>::sumRows() const
{
  VectorDyn r(getSize());
  for( int j = 0; j < getSize(); j++ )
  {
    for( int i = 0; i < bandwidth_full; i++)
    {
      r[i] += _at(i, j);
    }
  }
  return r;
}

template<typename Real>
VectorDyn SparseDiagonalMatrixOpt<Real>::operator* (const VectorDyn& x) const
{
  VectorDyn r(x.size());
  for(int j = 0; j < size; j++)
  {
    r[j] = 0;
    int o = getOffset(j);
    for( int i = 0; i < bandwidth_full; i++) 
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
  return r;

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

template<typename Real>
int SparseDiagonalMatrixOpt<Real>::getSize() const
{
  return size;
}

template<typename Real>
int SparseDiagonalMatrixOpt<Real>::getBandwidth() const
{
  return bandwidth;
}

template<>
void SparseDiagonalMatrixOpt<double>::multiplyWith( const VectorDouble& x, VectorDouble& r ) const
{
  // first, check if size of matrix and vectors are right:
  if( size != x.size() )
    throw std::runtime_error("size of x must be the same as matrix size");

  if( size != r.size() )
    throw std::runtime_error("size of r must be the same as matrix size");

  //if( size%2 != 0 )
  //  throw std::runtime_error("size of matrix must be multiple of 2");


  if ((this->bandwidth_full) % 4 != 0 ) 
    throw std::runtime_error("fast multiplication only works with matrix with bandwidth_full is multiple of 4");
  
  // now, do the multiplication
  const int o = getOffset(0);
  int k = o;

  int mox = 0;

  // first rows
  for(int j = 0; j < -o; j++)
  {
    r[j] = 0;
    for( int i = -k ; i < bandwidth_full; i+=4) 
    {
      int X = i + k;
      __m128d matrixPart1 = _mm_load_pd( &_at(i,j) );
      __m128d vectorPart1 = _mm_load_pd( &x[X] );
      __m128d prod1 = _mm_mul_pd( matrixPart1, vectorPart1 );
      __m128d matrixPart2 = _mm_load_pd( &_at(i+2,j) );
      __m128d vectorPart2 = _mm_load_pd( &x[X+2] );
      __m128d prod2 = _mm_mul_pd( matrixPart2, vectorPart2 );
      __m128d sum = _mm_add_pd( prod1, prod2 );
      __m128d red1 = _mm_hadd_pd( sum, sum );

      double res;
      _mm_store_sd( &res, red1 );

      r[j] += res;
      
    }
    if(j%2 == mox) k+=2;
  }

  #define USE_SINGLE_ROW_PROCESSING
  #ifdef USE_SINGLE_ROW_PROCESSING

  // middle part
  for(int j = -o; j < size+o; j+=1)
  {
    r[j] = 0;
    for( int i = 0; i < bandwidth_full; i+=4) 
    {
      int X = i + k;
      //r[j] += _at(i,j) * x[X];
      //std::cout<<_at(i,j) * x[X]<<std::endl;
            //int j2 = j+1;
      //std::cout<< " try " << i << " xoffset "<< (i+0+j-bandwidth_2) <<" check "<< (j-bandwidth_2) % 2 << std::endl;
      
      // first row
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
    if(j%2 == mox) k+=2;
    //std::cout<<"="<<r[j]<<std::endl;
  }

#else  
  // middle part
  for(int j = -o; j < size+o; j+=2)
  {
    r[j] = 0;
    r[j+1] = 0;
    for( int i = 0; i < bandwidth_full; i+=4) 
    {
      int X = i + k;
      //r[j] += _at(i,j) * x[X];
      //std::cout<<_at(i,j) * x[X]<<std::endl;
            //int j2 = j+1;
      //std::cout<< " try " << i << " xoffset "<< (i+0+j-bandwidth_2) <<" check "<< (j-bandwidth_2) % 2 << std::endl;
      
      // first row
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

      Real res;
      _mm_store_sd( &res, red1 );

      r[j] += res;

      // second row
      matrixPart1 = _mm_load_pd( &_at(i,j+1) );
      prod1 = _mm_mul_pd( matrixPart1, vectorPart1 );
      matrixPart2 = _mm_load_pd( &_at(i+2,j+1) );
      prod2 = _mm_mul_pd( matrixPart2, vectorPart2 );
      sum = _mm_add_pd( prod1, prod2 );
      red1 = _mm_hadd_pd( sum, sum );
      _mm_store_sd( &res, red1 );
      r[j+1] += res;
      
    }
    if(j%2 == mox) k+=2;
    //std::cout<<"="<<r[j]<<std::endl;
  }
#endif

  // end rows
  for(int j = size+o; j < size; j++)
  {
    r[j] = 0;
    for( int i = 0 ; i < size-k; i+=4) 
    {
      int X = i + k;
      __m128d matrixPart1 = _mm_load_pd( &_at(i,j) );
      __m128d vectorPart1 = _mm_load_pd( &x[X] );
      __m128d prod1 = _mm_mul_pd( matrixPart1, vectorPart1 );
      __m128d matrixPart2 = _mm_load_pd( &_at(i+2,j) );
      __m128d vectorPart2 = _mm_load_pd( &x[X+2] );
      __m128d prod2 = _mm_mul_pd( matrixPart2, vectorPart2 );
      __m128d sum = _mm_add_pd( prod1, prod2 );
      __m128d red1 = _mm_hadd_pd( sum, sum );

      double res;
      _mm_store_sd( &res, red1 );
      //std::cout<<"res  "<<i<<":  "<<res<<std::endl;

      r[j] += res;
      
    }
    if(j%2 == mox) k+=2;
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

      Real res;
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

template<>
void SparseDiagonalMatrixOpt<float>::multiplyWith( const VectorFloat& x, VectorFloat& r ) const
{

  r = (*this) * x;
  return;

  // first, check if size of matrix and vectors are right:
  if( size != x.size() )
    throw std::runtime_error("size of x must be the same as matrix size");

  if( size != r.size() )
    throw std::runtime_error("size of r must be the same as matrix size");

  //if( size%2 != 0 )
  //  throw std::runtime_error("size of matrix must be multiple of 2");


  if ((this->bandwidth_full) % 4 != 0 ) 
    throw std::runtime_error("fast multiplication only works with matrix with bandwidth_full is multiple of 4");
  
  // now, do the multiplication
  const int o = getOffset(0);
  int k = o;

  int mox = 0;

  // first rows
  std::cout<<"first rows"<<std::endl;
  for(int j = 0; j < -o; j++)
  {
    r[j] = 0;
    for( int i = -k ; i < bandwidth_full; i+=4) 
    {
      int X = i + k;
      std::cout << i << " / " << X << std::endl;
      __m128 matrixPart1 = _mm_load_ps( &_at(i,j) );
      std::cout << i << " / " << X << std::endl;
      __m128 vectorPart1 = _mm_load_ps( &x[X] );

      __m128 r1 = _mm_mul_ps(matrixPart1, vectorPart1);
      __m128 r2 = _mm_hadd_ps(r1, r1);
      __m128 r3 = _mm_hadd_ps(r2, r2);

      float res;
      _mm_store_ss(&res, r3);

      r[j] += res;
      
    }
    if(j%2 == mox) k+=2;
  }

  #define USE_SINGLE_ROW_PROCESSING
  #ifdef USE_SINGLE_ROW_PROCESSING

  std::cout<<"middle part"<<std::endl;
  // middle part
  for(int j = -o; j < size+o; j+=1)
  {
    r[j] = 0;
    for( int i = 0; i < bandwidth_full; i+=4) 
    {
      int X = i + k;
      __m128 matrixPart1 = _mm_load_ps( &_at(i,j) );
      __m128 vectorPart1 = _mm_load_ps( &x[X] );

      __m128 r1 = _mm_mul_ps(matrixPart1, vectorPart1);
      __m128 r2 = _mm_hadd_ps(r1, r1);
      __m128 r3 = _mm_hadd_ps(r2, r2);

      float res;
      _mm_store_ss(&res, r3);

      r[j] += res;
      
    }
    if(j%2 == mox) k+=2;
    //std::cout<<"="<<r[j]<<std::endl;
  }

#else  
  // middle part
  for(int j = -o; j < size+o; j+=2)
  {
    r[j] = 0;
    r[j+1] = 0;
    for( int i = 0; i < bandwidth_full; i+=4) 
    {
      int X = i + k;
      __m128 matrixPart1 = _mm_load_ps( &_at(i,j) );
      __m128 vectorPart1 = _mm_load_ps( &x[X] );

      __m128 r1 = _mm_mul_ps(matrixPart1, vectorPart1);
      __m128 r2 = _mm_hadd_ps(r1, r1);
      __m128 r3 = _mm_hadd_ps(r2, r2);

      float res;
      _mm_store_ss(&res, r3);

      r[j] += res;

      // second row
      __m128 matrixPart2 = _mm_load_ps( &_at(i,j+1) );
      __m128 t1 = _mm_mul_ps(matrixPart2, vectorPart1);
      __m128 t2 = _mm_hadd_ps(t1, t1);
      __m128 t3 = _mm_hadd_ps(t2, t2);

      float res2;
      _mm_store_ss(&res2, t3);

      r[j+1] += res2;

    }
    if(j%2 == mox) k+=2;
    //std::cout<<"="<<r[j]<<std::endl;
  }
#endif

  // end rows
  for(int j = size+o; j < size; j++)
  {
    r[j] = 0;
    for( int i = 0 ; i < size-k; i+=4) 
    {
      int X = i + k;
      __m128 matrixPart1 = _mm_load_ps( &_at(i,j) );
      __m128 vectorPart1 = _mm_load_ps( &x[X] );

      __m128 r1 = _mm_mul_ps(matrixPart1, vectorPart1);
      __m128 r2 = _mm_hadd_ps(r1, r1);
      __m128 r3 = _mm_hadd_ps(r2, r2);

      float res;
      _mm_store_ss(&res, r3);

      r[j] += res;
      
    }
    if(j%2 == mox) k+=2;
  }

}

template class SparseDiagonalMatrixOpt<double>;
template class SparseDiagonalMatrixOpt<float>;