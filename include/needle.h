/**
 *
 */


#include <iostream>

#include "mathheader.h"
#include "sparsediagonalmatrix.h"

// todo add constructor for diagonal matrix

class BendingNeedleModel
{
private:
  int numNodes;
  SparseDiagonalMatrix dF_dx_new;

  std::vector<Vector> x;
  std::vector<Vector> f;
  std::vector<Vector> v;
  std::vector<Vector> a;

  Vector calcF(int i, double k);
  Vector calcFNext(int i, double k);
  Vector calcFPrev(int i, double k);
  Vector calcSpring(Vector a, Vector b, double k);
  void cg(const cml::matrixd_r& A, const cml::vectord& b, cml::vectord& x );

  double totaltime;

  bool debugOut;

  double kSpring;
  double kNeedle;


public:
  BendingNeedleModel();
  void simulateImplicitChentanez( double dt );
  const std::vector<Vector>& getX() const;
};