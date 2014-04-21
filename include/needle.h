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

  double getTotalLength();


};
