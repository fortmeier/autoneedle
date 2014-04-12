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

  double& _at(int i, int j);

public:
  SparseDiagonalMatrix( int m, int b );
  ~SparseDiagonalMatrix();

  void zero();

  double& operator() (int i, int j);

  cml::vectord operator* (cml::vectord x);

  friend std::ostream& operator<< ( std::ostream &out, SparseDiagonalMatrix &matrix );

};
