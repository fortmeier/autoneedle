/**
 *
 */


#include <iostream>
#include <map>

#include "mathheader.h"
#include "sparsediagonalmatrix.h"


class NeedleMatrix
{
public:
  typedef std::map<int, Vector > Modifiers;

private:
  int numNodes;
  //int numLagrangeModifiers;
  SparseDiagonalMatrix A;

  Modifiers modifiers;

  void updateRow( int j );
public:

  NeedleMatrix( int numNodes );
  cml::vectord operator* (const cml::vectord& x) const;
  SparseDiagonalMatrix& getSystemMatrix();
  Modifiers& getLagrangeModifiers();

};

