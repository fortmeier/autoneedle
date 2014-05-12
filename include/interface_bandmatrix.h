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

class BandMatrixInterface
{
public:
  /**
   * set the matrix to zero
   */
  virtual void zero() = 0;

  /**
   * access the matrix a position i,j
   */
  virtual double& operator() (int i, int j) const = 0;
  
  /**
   * access the b x m matrix at position i,j
   * this is slightly different from () and should only be called when actually needed
   * or a time costly operation is performed.
   */
  virtual double& _at(int i, int j) const = 0;

  virtual cml::vectord operator* (const cml::vectord& x) const = 0;

  virtual int getSize() = 0;
  virtual int getBandwidth() = 0;

    /**
   * fast multiplication using intrincisc
   */
  virtual void multiplyWith( const cml::vectord& x, cml::vectord& r ) const = 0;

  virtual std::ostream& print ( std::ostream &out ) const = 0;

  friend std::ostream& operator<< ( std::ostream &out, const BandMatrixInterface &matrix );

};

