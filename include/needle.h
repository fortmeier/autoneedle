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
 * A needle is a set of n connected nodes with a base and a tip:
 *
 * (0)-----(1)-----(2)-----(3)-----(4)- . . . -(n-1)
 *
 * The 0-th node is called the base node, which usually is connected
 * to a handle.
 * The (n-1)-th node is called the needle tip node and is interacting
 * with tissue.
 *
 * There are two different modes for the computation of deformations
 * available. The first is a dynamic mode, which uses
 * the Newton-Raphson integration scheme of [1], including inertia.
 * The second one is a quasi-static mode as used in [2].
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

  Vector calcF(int i, double k) const;
  Vector calcFNext(int i, double k) const;
  Vector calcFPrev(int i, double k) const;
  Vector calcSpring(Vector a, Vector b, double k) const;

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

  /**
   * indicates if the needle should output large
   * amounts of debugging information
   * This should only be enabled by developers for debugging
   * purposes.
   */
  bool debugOut;

  /**
   * needle Stiffness
   */
  double kNeedle;

  /**
   * distance between two nodes
   */
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

  /**
   * Adds a lagrange modifier to the node with index nodeIndex
   * N is the direction in which the movement of the node is
   * restricted. To restrict the movement along a line, use
   * two Lagrange-Modifiers perpendicular to the line.
   * One Lagrange-Modifier restricts the movement to the plane
   * perpendicular to the normal N. Three linearly independent
   * Modifiers fix a point completely at a position.
   * Currently, this is not supported.
   */
  void addLagrangeModifier( int nodeIndex, Vector N );

  /**
   * connects the node with index nodeIndex to a position pos
   * k is the spring stiffness constant
   */ 
  void setSpring( int nodeIndex, Vector pos, double k );

  /**
   * return a vector of the positions of the nodes of the
   * needle
   */
  const std::vector<Vector>& getX() const;
  
  /**
   * return the length of a single segment, that is is
   * distance between two adjacent nodes of the undeformed
   * needle
   */
  double getSegmentLength( ) const;

  /**
   * set the position of the base node of the needle
   */
  void setBasePosition( const Vector& pos );

  /**
   * set the direction of the needle at the base
   */ 
  void setBaseDirection( const Vector& dir );

  /**
   * set the force vector that acts on all nodes of the needle.
   * Per default, this is (0, -9.81m/s^2, 0)
   */
  void setGravity( const Vector& g );

  /**
   * get the torque that is acting on the base of a needle node
   */
  Vector getBaseTorque() const;

  /**
   * calculate the force that is acting on the base node of the
   * needle
   */
  Vector getBaseForce() const;

  /**
   * calculate the total length of all segments of the needle
   */
  double getTotalLength();

  /**
   * reset the needle to an undeformed configuration
   * with position and orientation as defined by
   * basePosition and baseDirection
   */
  void reset();



};
