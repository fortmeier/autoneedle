
#include <iostream>

#include "needle.h"

SparseDiagonalMatrix::SparseDiagonalMatrix( int m, int n ) :
  columns(m),
  rows(n)
{

  values = new double[m*n];
  zero();

}

SparseDiagonalMatrix::~SparseDiagonalMatrix()
{
	delete[] values;
}

void SparseDiagonalMatrix::zero()
{
  for( int i = 0; i < columns; i++ )
  {
    for( int j = 0; j < rows; j++ )
    {
      (*this)(i,j) = 0;
    }
  }
}

std::ostream& operator<< ( std::ostream &out, SparseDiagonalMatrix &matrix )
{
  for( int j = 0; j < matrix.rows; j++ )
  {
    out << "[ ";
    for( int k = 0; k < j; k++) 
    {
      out << ". ";
    }

    for( int i = 0; i < matrix.columns; i++ )
    {
      out << matrix(i,j) << " ";
    }

    for( int k = 0; k < matrix.columns - j -1; k++) 
    {
      out << ". ";
    }

    out << " ]\n";
  }
  return out;
}

double& SparseDiagonalMatrix::operator() ( int i, int j )
{
  return values[ i + columns*j ];
}

BendingNeedleModel::BendingNeedleModel() :
  numberNodes(10),
  dF_dx_new( numberNodes, 5 * 3 )
{
  exit(0);
}