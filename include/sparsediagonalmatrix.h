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
