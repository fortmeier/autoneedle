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
#include <map>

#include "mathheader.h"
#include "needlematrix.h"

/**
 * auxilliary class that repersents a spring which is
 * attachted to a node of a needle
 */
class Spring
{
public:
  /**
   * position of spring end
   */
  Vector x;

  /**
   * spring stiffness
   */
  double k;

  Spring( Vector x = Vector(0,0,0), double k = 0);
};

/**
 * A bendable needle.
 * Can compute deformations and uses concepts from [1-4].
 *
 * There are two different modes available. The first is a dynamic mode, which uses
 * the integration scheme of [1], including inertia.
 * The second one is a quasi-static mode as used in [2]. ATM, this is not implemented.
 * Both methods need the computation of a Jacobian. Here, this is supported by
 * automatic differentiation of the energy terms. See \see GenerateCode.py.
 *
 * [1] Chentanez N, Alterovitz R, Ritchie D, Cho L, Hauser KK, Goldberg K, et al. 
 * Interactive simulation of surgical needle insertion and steering. 
 * ACM Trans Graph;28(3):1-10. Available from: http://portal.acm.org/citation.cfm?id=1531394
 *
 * [2] Goksel et al. ...
 *
 * [3] 
 *
 * [4]
 */
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

  /**
   * use the conjugate gradient method to solve a system of
   * linear equations (Ax=b).
   */
  void cg( );

  void updateJacobianForce();
  void updateJacobianVelocity();
  void updateSystemMatrix_A();
  void updateResultVector_b();
  void addForcesToB();
  double updateStep();

  double totaltime;
  double dt;

  bool debugOut;

  double kNeedle;

  double segmentLength;

  /**
   * position of the needle base
   */
  Vector baseDirection;

  /**
   * direction of the needle base
   */
  Vector basePosition;

  /**
   * gravity force vector;
   */
  Vector G;

public:
  BendingNeedleModel( double length, int nodes, double k );

  /**
   * Use the dynamic simulation.
   * This first solves A * a+ = b with
   * A = ( M - dF_dx * dt**2 * beta - dF_dv* dt * gamma)
   * and
   * b  = F + dt * dF_dx * ( v + dt * ( 0.5 - beta) * a ) + dF_dv * dt * (1-gamma) * a
   * and the uses Newmark's method for the velocity and position update:
   * See [1].
   */
  double simulateImplicitDynamic( double dt );

  /**
   * Use the static simulation.
   * This solves K * u = f
   * See [2].
   */
  double simulateImplicitStatic( double dt );

  void addLagrangeModifier( int nodeIndex, Vector N );

  void setSpring( int nodeIndex, Vector pos, double k );

  const std::vector<Vector>& getX() const;
  double getSegmentLength( ) const;

  void setBasePosition( const Vector& pos );
  void setBaseDirection( const Vector& dir );

  void setGravity( const Vector& g );

  double getTotalLength();

  void reset();



};
