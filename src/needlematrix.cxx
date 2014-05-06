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

#include <iostream>
#include <cassert>
#include <cstring>

#include <stdexcept>

#include "needlematrix.h"


NeedleMatrix::NeedleMatrix( int nNodes ) :
  numNodes( nNodes ),
  A( numNodes*3, 19 )
{

}

cml::vectord NeedleMatrix::operator* (const cml::vectord& x) const
{
  int numLagrangeModifiers = modifiers.size();
  //std::cout<<x.size()<<" vs. "<<numNodes * 3 + numLagrangeModifiers<<std::endl;

  assert( "length of vector is not okay!" && x.size() == numNodes * 3 + numLagrangeModifiers );
  //if( r_lambda.size() != numLagrangeModifiers ) assert ("not implemented yet!" && false );

  //std::cout<<"x: "<<x<<std::endl;

  // first compute result of sparse diagonal matrix;
#if 1  
  cml::vectord r = A * x;
#else
  // use faster method
  cml::vectord r(x.size());
  A.multiplyWith(x, r);
#endif
  //std::cout<<"r1: "<<r<<std::endl;
  r.resize( numNodes * 3 + numLagrangeModifiers );
  //std::cout<<"r2: "<<r<<std::endl;
  int n = numNodes;
  // add multiplication with lagrange modifiers
  for( Modifiers::const_iterator iter = modifiers.begin(); iter != modifiers.end(); iter++ )
  {
    int i = iter->first;
    
    // add multipliciation from elements in rows of matrix A
    r[i*3+0] += iter->second[0] * x[n*3+i];
    r[i*3+1] += iter->second[1] * x[n*3+i];
    r[i*3+2] += iter->second[2] * x[n*3+i];

    // add sum of multiplications of lagrange modifiers and x for last rows
    r[n*3+i]  = iter->second[0] * x[i*3+0];
    r[n*3+i] += iter->second[1] * x[i*3+1];
    r[n*3+i] += iter->second[2] * x[i*3+2];

    /*A(n*3+i, i*3+0) = N[0];
    A(n*3+i, i*3+1) = N[1];
    A(n*3+i, i*3+2) = N[2];

    A(i*3+0, n*3+i) = N[0];
    A(i*3+1, n*3+i) = N[1];
    A(i*3+2, n*3+i) = N[2];*/
  }

//  std::cout<<x_needle<<std::endl;
  //std::cout<<"r: "<<r<<std::endl;
  
/* todo readd lagrange multipliers
    // add lagrange multipliers
  for(int i = 0; i < n; i++) 
  {
    Vector N(1,0,0);
    if(i==0) N = Vector(0,1,0);
    if(i==1) N = Vector(0,1,0);

    if(i>1) N = (x[i]-x[i-1]).normalize();
#if 1
    A(n*3+i, i*3+0) = N[0];
    A(n*3+i, i*3+1) = N[1];
    A(n*3+i, i*3+2) = N[2];
    A(i*3+0, n*3+i) = N[0];
    A(i*3+1, n*3+i) = N[1];
    A(i*3+2, n*3+i) = N[2];
#endif

    //A(n*3+i, n*3+i) = 1;
    //b[n*3+i] = 0; //-4*cml::dot(N, v[i])*dt - dt*dt*cml::dot(N,a[i]);

    //if(i==0 || i==1) b[n*3+i] += 10.1*sin(totaltime);

  } */

  return r;
/*  size_t size_r_b = numNodes * 3 * sizeof(double); 
  memcpy( r.data(), r_b.data(), size_r_b );

  size_t size_r_lambda = numLagrangeModifiers * sizeof(double);
  std::memcpy( r.data() + size_r_b, r_lambda.data(), size_r_lambda);
  return x;*/
}

SparseDiagonalMatrix& NeedleMatrix::getSystemMatrix()
{
  return A;
}

NeedleMatrix::Modifiers& NeedleMatrix::getLagrangeModifiers()
{
  return modifiers;
}


