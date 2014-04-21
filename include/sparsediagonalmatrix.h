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

#include "mathheader.h"

class SparseDiagonalMatrix
{
private:
  double* values;
  int size;
  int bandwidth;
  int bandwidth_2;
  int offset;


public:
  SparseDiagonalMatrix( int m, int b );
  ~SparseDiagonalMatrix();

  void zero();

  double& operator() (int i, int j) const;
  
  double& _at(int i, int j) const;

  cml::vectord operator* (const cml::vectord& x) const;

  friend std::ostream& operator<< ( std::ostream &out, const SparseDiagonalMatrix &matrix );

  int getSize();
  int getBandwidth();

};
