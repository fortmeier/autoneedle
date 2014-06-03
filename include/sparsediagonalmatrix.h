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

/*
 * A sparse diagonal matrix only has entries close to the diagonal of the matrix:
 * | * * 0 0 0 |
 * | * * * 0 0 |
 * | 0 * * * 0 |
 * | 0 0 * * * |
 * | 0 0 0 * * |
 * This example is a 5x5 diagonal with a bandwidth of 3 ( only the 3 entries
 * closest to the diagonal are different from 0).
 * Here, the matrix is stored as a bandwithd * rows matrix:
 * | * * * |
 * | * * * |
 * | * * * |
 * | * * * |
 * | * * * |
 * This way, there are additional values that can be different from zero in a
 * few of the first and last rows:
 * | * * * 0 0 |
 * | * * * 0 0 |
 * | 0 * * * 0 |
 * | 0 0 * * * |
 * | 0 0 * * * |
 */

template<typename Real>
class SparseDiagonalMatrix : public BandMatrixInterface<Real>
{
private:
  Real* values;
  int size;
  int bandwidth;
  int bandwidth_2;
  int offset;


public:
  SparseDiagonalMatrix( int m, int b );
  ~SparseDiagonalMatrix();

  /**
   * set the matrix to zero
   */
  virtual void zero();

  virtual void zeroRow( int j );

  /**
   * access the matrix a position i,j
   */
  Real& operator() (int i, int j) const;

  virtual VectorDyn sumRows() const;
  
  /**
   * access the b x m matrix at position i,j
   * this is slightly different from () and should only be called when actually needed
   * or a time costly operation is performed.
   */
  Real& _at(int i, int j) const;

  cml::vectord operator* (const cml::vectord& x) const;

  //friend std::ostream& operator<< ( std::ostream &out, const SparseDiagonalMatrix &matrix );

  int getSize() const;
  int getBandwidth() const;

    /**
   * fast multiplication using intrincisc
   */
  void multiplyWith( const cml::vectord& x, cml::vectord& r ) const;

  virtual std::ostream& print ( std::ostream &out ) const;


};
