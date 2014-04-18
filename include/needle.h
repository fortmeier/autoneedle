/**
 *
 */


#include <iostream>
#include <map>

#include "mathheader.h"
#include "needlematrix.h"

class Spring
{
public:
  Vector x;
  double k;

  Spring( Vector x = Vector(0,0,0), double k = 0);
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

  typedef std::map<int, Spring> SpringMap;
  SpringMap springs;

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

  double kNeedle;

  double segmentLength;

  Vector baseDirection;
  Vector basePosition;


public:
  BendingNeedleModel( double length, int nodes, double k );
  double simulateImplicitChentanez( double dt );
  void addLagrangeModifier( int nodeIndex, Vector N );

  void setSpring( int nodeIndex, Vector pos, double k );

  const std::vector<Vector>& getX() const;
  double getSegmentLength( ) const;

  void setBasePosition( const Vector& pos );
  void setBaseDirection( const Vector& dir );


};
