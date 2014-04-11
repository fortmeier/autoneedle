/**
 *
 */


#include <iostream>

#include "mathheader.h"

class SparseDiagonalMatrix
{
private:
  double* values;
  int columns;
  int rows;

public:
  SparseDiagonalMatrix( int m, int n );
  ~SparseDiagonalMatrix();

  void zero();

  double& operator() (int i, int j);

  friend std::ostream& operator<< ( std::ostream &out, SparseDiagonalMatrix &matrix );

};

// todo add constructor for diagonal matrix

class BendingNeedleModel
{
private:
  int numberNodes;
  SparseDiagonalMatrix dF_dx;
public:
  BendingNeedleModel();
};