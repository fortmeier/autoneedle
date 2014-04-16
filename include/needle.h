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

class BendingNeedleModel
{
private:
  typedef cml::vectord NeedleVector;
  int numNodes;
  NeedleMatrix A;
  NeedleVector b;

  SparseDiagonalMatrix dF_dx;
  SparseDiagonalMatrix dF_dv;

  std::vector<Vector> nodes;

  cml::vectord x; // positions from last step
  cml::vectord v; // velocities from last step
  cml::vectord ao; // accelerations from last step
  cml::vectord ap; // accelerations from last step

  cml::vectord m; // accelerations from last step

  Vector calcF(int i, double k);
  Vector calcFNext(int i, double k);
  Vector calcFPrev(int i, double k);
  Vector calcSpring(Vector a, Vector b, double k);
  void cg( );

  void updateJacobianForce();
  void updateJacobianVelocity();
  void updateSystemMatrix_A();
  void updateResultVector_b();
  double updateStep();

  double totaltime;
  double dt;

  bool debugOut;

  double kSpring;
  double kNeedle;


public:
  BendingNeedleModel();
  void simulateImplicitChentanez( double dt );
  void addLagrangeModifier( int nodeIndex, Vector N );
  const std::vector<Vector>& getX() const;

};
