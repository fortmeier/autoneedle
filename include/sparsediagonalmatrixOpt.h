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

#pragma once

#include <iostream>

#include "mathheader.h"

#include "interface_bandmatrix.h"

class SparseDiagonalMatrixOpt : public BandMatrixInterface
{
private:
  double* values;
  int size;
  int bandwidth;
  int bandwidth_2;
  int offset;

  int getOffset( int y ) const;


public:
  SparseDiagonalMatrixOpt( int m, int b );
  ~SparseDiagonalMatrixOpt();

  /**
   * set the matrix to zero
   */
  void zero();

  /**
   * access the matrix a position i,j
   */
  double& operator() (int i, int j) const;
  
  /**
   * access the b x m matrix at position i,j
   * this is slightly different from () and should only be called when actually needed
   * or a time costly operation is performed.
   */
  double& _at(int i, int j) const;

  cml::vectord operator* (const cml::vectord& x) const;

  //friend std::ostream& operator<< ( std::ostream &out, const SparseDiagonalMatrixOpt &matrix );

  virtual std::ostream& print( std::ostream &out ) const;

  int getSize();
  int getBandwidth();

  /**
   * fast multiplication using intrincisc
   */
  void multiplyWith( const cml::vectord& x, cml::vectord& r ) const;

};
