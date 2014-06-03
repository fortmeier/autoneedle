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
#include <cassert>
#include <cstring>

#include <stdexcept>

#include "needle.h"
#include "generatedCode.h"

template<typename Real>
Spring<Real>::Spring( Vector _x, Real _k ) :
  x( _x ),
  k( _k ),
  needsUpdate( true )
{

}


template<typename Real>
BendingNeedleModel<Real>::BendingNeedleModel( Real length, int nNum, Real k ) :
  numNodes( nNum ),
  A( numNodes ),
  b( numNodes * 3),
  dF_dx( numNodes*3, 19 ),
  dF_dv( numNodes*3, 19 ),
  nodes( numNodes ),
  nodesLastUpdate( numNodes ), // position of node at last time of update of the needle matrix
  normals( numNodes ),
  x(numNodes*3), // positions from last step
  v(numNodes*3), // velocities from last step
  ao(numNodes*3), // accelerations from last step
  ap(numNodes*3), // accelerations from last step
  m(numNodes*3), // mass of node in x y z direction
  totaltime( 0 ),
  debugOut( false ),
  kNeedle(k),
  kBase(2000),
  segmentLength(length/(Real)(numNodes-1)),
  baseDirection(1,0,0),
  basePosition(0,0,0),
  G(0, -9.81, 0)
{

  //testMatrix();
  for(int i = 0; i < numNodes; i++)
  {
     m[i*3+0] = 0.0001;
     m[i*3+1] = 0.0001;
     m[i*3+2] = 0.0001;
  }
  reset();

  for(int i = 0; i < numNodes; i++)
  {
    //addLagrangeModifier(i, Vector(1,0,0));
  }

  setSpring( 0, Vector( 0, 0, 0 ), 200 );
  setSpring( 1, Vector( segmentLength, 0, 0 ), 200 );

}

template<typename Real>
Vector BendingNeedleModel<Real>::calcF(int i, Real k) const
{
  const std::vector<Vector>& x = nodes;
  Vector f = Vector(
    needle_Fx(x[i-1],x[i],x[i+1], k, segmentLength),
    needle_Fy(x[i-1],x[i],x[i+1], k, segmentLength),
    needle_Fz(x[i-1],x[i],x[i+1], k, segmentLength)
  );
  return  f;
}

template<typename Real>
Vector BendingNeedleModel<Real>::calcFNext(int i, Real k) const
{
  const std::vector<Vector>& x = nodes;
  Vector f = Vector(
    needle_Fx_next(x[i],x[i+1],x[i+2], k, segmentLength),
    needle_Fy_next(x[i],x[i+1],x[i+2], k, segmentLength),
    needle_Fz_next(x[i],x[i+1],x[i+2], k, segmentLength)
  );
  return  f;
}

template<typename Real>
Vector BendingNeedleModel<Real>::calcFPrev(int i, Real k) const
{
  const std::vector<Vector>& x = nodes;
  Vector f = Vector(
    needle_Fx_next(x[i],x[i-1],x[i-2], k, segmentLength),
    needle_Fy_next(x[i],x[i-1],x[i-2], k, segmentLength),
    needle_Fz_next(x[i],x[i-1],x[i-2], k, segmentLength)
  );
  return  f;
}

template<typename Real>
Vector BendingNeedleModel<Real>::calcSpring(Vector a, Vector b, Real k) const
{
  Vector f = Vector(
    spring_FSx(a,b,k),
    spring_FSy(a,b,k),
    spring_FSz(a,b,k)
  );
  return  f;
}

template<typename Real>
Vector BendingNeedleModel<Real>::calcSpringNoTangential(Vector a, Vector b, Vector n, Real k) const
{
  Vector f = Vector(
    springNoTangential_FSx(a,b,n,k),
    springNoTangential_FSy(a,b,n,k),
    springNoTangential_FSz(a,b,n,k)
  );
  return  f;
}


template<typename Real>
void BendingNeedleModel<Real>::cg()
{
  if(debugOut) std::cout<<"Starting CG"<<std::endl;
  cml::vectord x = ap;
  if(debugOut) std::cout<<"b:"<<b<<std::endl;
  x.resize(b.size());

  cml::vectord r = b - A * x;

  if(debugOut) std::cout<<"r: "<<r<<std::endl;

  cml::vectord p = r;

  Real rsold=cml::dot(r,r);
  Real rsnew;

  Real eps = 0.000001;
  int i = 0;
  do {
    cml::vectord Ap = A * p;
    if(debugOut) std::cout<<"1:"<<p<<std::endl;
    if(debugOut) std::cout<<"2:"<<Ap<<std::endl;
    Real alpha = rsold/cml::dot(p,Ap);
    x = x + alpha * p;
    r = r-alpha*Ap;
    rsnew = cml::dot(r,r);
    p=r+p*(rsnew/rsold);
    rsold=rsnew; 
    i++;
    //std::cout<<i<<": "<<rsnew<<std::endl;
  } while (sqrt(rsnew) > eps && i < 10000);
  if(debugOut) std::cout<<"cg needed "<<i<<" iterations"<<std::endl;

  //x.resize(ap.size());
  for(int i = 0; i < ap.size(); i++){
    ap[i] = x[i];
  }
}

template<typename Real>
void BendingNeedleModel<Real>::updateJacobianForce()
{
  dF_dx.zero();
  int n = numNodes;
  for(int i = 0; i < numNodes; i++)
  {
    std::vector<Vector>& x = nodes;

    // pF_pxi
    if(i>=2) needle_JacobianCCMatrixAdd( i*3, i*3, dF_dx, x[i-2], x[i-1], x[i], kNeedle, segmentLength);
    if(i>=1 && i < n-1) needle_JacobianBBMatrixAdd( i*3, i*3, dF_dx, x[i-1], x[i], x[i+1], kNeedle, segmentLength);
    if(i>=0 && i< n-2) needle_JacobianAAMatrixAdd( i*3, i*3, dF_dx, x[i], x[i+1], x[i+2], kNeedle, segmentLength);

    // pF_pxi+1
    if(i>=1 && i < n-1) needle_JacobianBCMatrixAdd( i*3, i*3+3, dF_dx, x[i-1], x[i], x[i+1], kNeedle, segmentLength);
    if(i>=0 && i < n-2) needle_JacobianABMatrixAdd( i*3, i*3+3, dF_dx, x[i], x[i+1], x[i+2], kNeedle, segmentLength);

    // pF_pxi+2
    if(i>=0 && i < n-2) needle_JacobianACMatrixAdd( i*3, i*3+6, dF_dx, x[i], x[i+1], x[i+2], kNeedle, segmentLength);

    // pF_pxi-1
    if(i>=1 && i < n-1) needle_JacobianBAMatrixAdd( i*3, i*3-3, dF_dx, x[i-1], x[i], x[i+1], kNeedle, segmentLength);
    if(i>=2 && i < n) needle_JacobianCBMatrixAdd( i*3, i*3-3, dF_dx, x[i-2], x[i-1], x[i], kNeedle, segmentLength);

    // pF_pxi+2
    if(i>=2 && i < n) needle_JacobianCAMatrixAdd( i*3, i*3-6, dF_dx, x[i-2], x[i-1], x[i], kNeedle, segmentLength);

/*    // testcopy
    if(i>=2) needle_JacobianCCMatrixAdd( i*3, i*3, dF_dx_new, x[i-2], x[i-1], x[i], kNeedle);
    if(i>=1 && i < n-1) needle_JacobianBBMatrixAdd( i*3, i*3, dF_dx_new, x[i-1], x[i], x[i+1], kNeedle);
    if(i>=0 && i< n-2) needle_JacobianAAMatrixAdd( i*3, i*3, dF_dx_new, x[i], x[i+1], x[i+2], kNeedle);

    // pF_pxi+1
    if(i>=1 && i < n-1) needle_JacobianBCMatrixAdd( i*3, i*3+3, dF_dx_new, x[i-1], x[i], x[i+1], kNeedle);
    if(i>=0 && i < n-2) needle_JacobianABMatrixAdd( i*3, i*3+3, dF_dx_new, x[i], x[i+1], x[i+2], kNeedle);

    // pF_pxi+2
    if(i>=0 && i < n-2) needle_JacobianACMatrixAdd( i*3, i*3+6, dF_dx_new, x[i], x[i+1], x[i+2], kNeedle);

    // pF_pxi-1
    if(i>=1 && i < n-1) needle_JacobianBAMatrixAdd( i*3, i*3-3, dF_dx_new, x[i-1], x[i], x[i+1], kNeedle);
    if(i>=2 && i < n) needle_JacobianCBMatrixAdd( i*3, i*3-3, dF_dx_new, x[i-2], x[i-1], x[i], kNeedle);

    // pF_pxi+2
    if(i>=2 && i < n) needle_JacobianCCMatrixAdd( i*3, i*3-6, dF_dx_new, x[i-2], x[i-1], x[i], kNeedle);
*/
    /*if(i>0 && i < n-1) {


      needle_JacobianMatrixAdd( i*3, i*3, dF_dx, x[i-1],x[i],x[i+1]);
      needle_JacobianprevMatrixAdd( i*3, i*3-3, dF_dx, x[i+1], x[i], x[i-1]);
      needle_JacobianprevMatrixAdd( i*3, i*3+3, dF_dx, x[i-1], x[i], x[i+1]);
    }*/
  } 

  //todo readd: spring_JacobianMatrixAdd( (n-1)*3, (n-1)*3, dF_dx, x[n-1], sX, kSpring ); 
  for( typename SpringMap::iterator iter = springs.begin(); iter != springs.end(); iter++) 
  {
    int i = iter->first;
    if( i == 0 )
      spring_JacobianMatrixAdd( i*3, i*3, dF_dx, nodes[i], iter->second.x, iter->second.k );
    else
      springNoTangential_JacobianMatrixAdd( i*3, i*3, dF_dx, nodes[i], iter->second.x, normals[i] , iter->second.k);
  }
  //dF_dx *= 1.0;
  //dF_dx.transpose();
  if(debugOut) std::cout<<"dF_dx: "<<std::endl<<dF_dx<<std::endl;

}

template<typename Real>
void BendingNeedleModel<Real>::updateJacobianForceLazy()
{
  //dF_dx.zero();
  int n = numNodes;
  for(int i = 0; i < numNodes; i++)
  {
    std::vector<Vector>& x = nodes;
    std::vector<Vector>& lx = nodesLastUpdate;
    typename SpringMap::iterator spring = springs.find(i);

    bool needsUpdate = false;

    Real diff = (x[i]-lx[i]).length_squared();
    //std::cout << i <<":"<< x[i] << " - " << lx[i] << " = " << diff << std::endl;
    assert ( x[i] == x[i] );
    if( diff > 0.001  ) needsUpdate = true;


    if( spring != springs.end() && spring->second.needsUpdate ) needsUpdate = true;

    if( needsUpdate == false ) continue;
    //std::cout << "not next" << std::endl;

    lx[i] = x[i];

    dF_dx.zeroRow( i*3 + 0);
    dF_dx.zeroRow( i*3 + 1);
    dF_dx.zeroRow( i*3 + 2);
    // pF_pxi
    if(i>=2) needle_JacobianCCMatrixAdd(  i*3, i*3, dF_dx, x[i-2], x[i-1], x[i], kNeedle, segmentLength);
    if(i>=1 && i < n-1) needle_JacobianBBMatrixAdd(  i*3, i*3, dF_dx, x[i-1], x[i], x[i+1], kNeedle, segmentLength);
    if(i>=0 && i< n-2) needle_JacobianAAMatrixAdd(  i*3, i*3, dF_dx, x[i], x[i+1], x[i+2], kNeedle, segmentLength);

    // pF_pxi+1
    if(i>=1 && i < n-1) needle_JacobianBCMatrixAdd(  i*3+3, i*3, dF_dx, x[i-1], x[i], x[i+1], kNeedle, segmentLength);
    if(i>=0 && i < n-2) needle_JacobianABMatrixAdd(  i*3+3, i*3, dF_dx, x[i], x[i+1], x[i+2], kNeedle, segmentLength);

    // pF_pxi+2
    if(i>=0 && i < n-2) needle_JacobianACMatrixAdd(  i*3+6, i*3, dF_dx, x[i], x[i+1], x[i+2], kNeedle, segmentLength);

    // pF_pxi-1
    if(i>=1 && i < n-1) needle_JacobianBAMatrixAdd(  i*3-3, i*3, dF_dx, x[i-1], x[i], x[i+1], kNeedle, segmentLength);
    if(i>=2 && i < n) needle_JacobianCBMatrixAdd(  i*3-3, i*3, dF_dx, x[i-2], x[i-1], x[i], kNeedle, segmentLength);

    // pF_pxi+2
    if(i>=2 && i < n) needle_JacobianCAMatrixAdd(  i*3-6, i*3, dF_dx, x[i-2], x[i-1], x[i], kNeedle, segmentLength);

    if( spring != springs.end())
    {
      if( i == 0 )
        spring_JacobianMatrixAdd( i*3, i*3, dF_dx, nodes[i], spring->second.x, spring->second.k );
      else
        springNoTangential_JacobianMatrixAdd( i*3, i*3, dF_dx, nodes[i], spring->second.x, normals[i] , spring->second.k);
    }
  }
  //dF_dx *= 1.0;
  //dF_dx.transpose();
  if(debugOut) std::cout<<"dF_dx: "<<std::endl<<dF_dx<<std::endl;

}

template<typename Real>
void BendingNeedleModel<Real>::updateJacobianVelocity()
{
  //dF_dv = alphaC * M + betaC * dF_dx;
  // OK dF_dv = -50.5 * M + 0.005 * dF_dx;
  dF_dv.zero();

  Real alphaC = -0.0;
  Real betaC = 0.0;

  for(int j = 0; j < dF_dv.getSize(); j++)
  {
    for(int i = 0; i < dF_dv.getBandwidth(); i++)
    {
      dF_dv._at(i,j) = dF_dx._at(i,j) * betaC;
    }
    dF_dv(j,j) += alphaC * m[j];
  }

  if(debugOut) std::cout<<"dF_dv: "<<std::endl<<dF_dv<<std::endl;

}

template<typename Real>
void BendingNeedleModel<Real>::updateSystemMatrix_A()
{
  // fill A: 
  //    = M - dF_dx * dt * dt * beta - dF_dv* dt * gamma;
  for(int j = 0; j < A.getSystemMatrix().getSize(); j++)
  {
    for(int i = 0; i < A.getSystemMatrix().getBandwidth(); i++)
    {
      Real a = -(dF_dx._at(i,j)* dt * dt * 0.25) - (dF_dv._at(i,j) * dt * 0.5);
      A.getSystemMatrix()._at(i,j) = a;
    }
    A.getSystemMatrix()(j,j) += m[j];
  }
  if(debugOut) std::cout<<"A: "<<std::endl<<A.getSystemMatrix()<<std::endl;

}

template<typename Real>
void dOut( const cml::vectord& d, bool shuffle = false )
{
  for(int i = 0; i < d.size() / 3; i++ )
  {
    std::cout<<d[i*3+0]<<std::endl;
    if( shuffle == false )
    {
      std::cout<<d[i*3+1]<<std::endl;
      std::cout<<d[i*3+2]<<std::endl;
    }
    else
    {
      std::cout<<d[i*3+2]<<std::endl;
      std::cout<<d[i*3+1]<<std::endl;
    }
    std::cout<<""<<std::endl;
  }
}

template<typename Real>
void BendingNeedleModel<Real>::updateResultVector_b()
{
  // fill b: F + dt * dF_dx * ( v + dt * ( 0.5 - beta) * a ) + dF_dv * dt * (1-gamma) * a
  cml::vectord tmp = (( v + dt * ( 0.25 ) * ao ) * dt);
  b  = dF_dx * tmp;
  //dOut(b);
  tmp = (ao * dt * 0.5);
  b += dF_dv * tmp;

  addForcesToB();

}

template<typename Real>
void BendingNeedleModel<Real>::addForcesToB()
{
  // std::cout<<b<<std::endl;
  int n = numNodes;

  for(int i = 0; i < n; i++)
  {
    Vector f(0,0,0);

    if( i > 0 && i < n-1 ) f += calcF(i, kNeedle);
    if( i > 1 ) f += calcFPrev(i, kNeedle);
    if( i < n-2 ) f += calcFNext(i, kNeedle);

    f[0] +=  G[0] * m[i*3+0];
    f[1] +=  G[1] * m[i*3+1];
    f[2] +=  G[2] * m[i*3+2];

    /*if(i==n-1) {

      Vector FP = calcSpring(x[i], sX, kSpring);
      if(true) std::cout<<"offset: "<<sX<<std::endl;
      f += FP;
    }*/
    // std::cout<<"calcF("<<i<<") = "<<calcF(i, kNeedle)<<std::endl;
    // std::cout<<"calcFPrev("<<i<<") = "<<calcFPrev(i, kNeedle)<<std::endl;
    // std::cout<<"calcFNext("<<i<<") = "<<calcFNext(i, kNeedle)<<std::endl;
    // std::cout<<"f("<<i<<") = "<<f<<std::endl;

    b[i*3 + 0] += f[0]; 
    b[i*3 + 1] += f[1]; 
    b[i*3 + 2] += f[2];

  }
  for(int i = 0; i < A.getLagrangeModifiers().size(); i++)
  {
    b[numNodes*3+i]=0;
  }

  for( typename SpringMap::iterator iter = springs.begin(); iter != springs.end(); iter++) 
  {
    int i = iter->first;
    Vector fSpring;
    if( i == 0)
      fSpring = calcSpring(nodes[i], iter->second.x, iter->second.k);
    else
      fSpring = calcSpringNoTangential(nodes[i], iter->second.x, normals[i], iter->second.k);

    b[i*3+0] += fSpring[0];
    b[i*3+1] += fSpring[1];
    b[i*3+2] += fSpring[2];
    // std::cout<<"fspring!"<<fSpring<<std::endl;
  } 
  if(debugOut) std::cout<<"b: "<<std::endl<<b<<std::endl;
}

template<typename Real>
Real BendingNeedleModel<Real>::updateStep()
{
  // length enforcement
  bool useEnforcement = false;

  Real error = 0;
  int n = numNodes;
  
  if(useEnforcement) 
  {
    /*ap[0] = 0;
    ap[1] = 0;
    ap[2] = 0;
    ap[3] = 0;
    ap[4] = 0;
    ap[5] = 0;*/
    for(int i = 1; i < n; i++)
    {
      Vector N = nodes[i]-nodes[i-1];
      Vector aa(ap[i*3+0], ap[i*3+1], ap[i*3+2]); 
      N.normalize();
      Vector p = cml::dot(N,aa) * N;
      ap[i*3 + 0] -= p[0];
      ap[i*3 + 1] -= p[1];
      ap[i*3 + 2] -= p[2];
    }
  }  

  // update x and v for next iteration
  /////////////////////////////////////////////////////////

  for(int i = 0; i < n; i++)
  {
    int ix = i * 3 + 0;
    int iy = ix + 1;
    int iz = ix + 2;

    // update std::vector x (giving x_plus)
    Real ux = dt * v[ix] + dt*dt*(0.25*ao[ix] + 0.25*ap[ix]);
    Real uy = dt * v[iy] + dt*dt*(0.25*ao[iy] + 0.25*ap[iy]);
    Real uz = dt * v[iz] + dt*dt*(0.25*ao[iz] + 0.25*ap[iz]);

    x[ix] += ux;
    x[iy] += uy;
    x[iz] += uz;

    error += ux*ux + uy*uy + uz*uz;

    // update std::vector v (giving v_plus)
    v[ix] += dt*(0.5*ao[ix] + 0.5*ap[ix]);
    v[iy] += dt*(0.5*ao[iy] + 0.5*ap[iy]);
    v[iz] += dt*(0.5*ao[iz] + 0.5*ap[iz]);
  }
  // length enforcement

  if(useEnforcement) 
  {
    for(int i = 1; i < n; i++)
    {
      int ix = i * 3 + 0;
      int iy = ix + 1;
      int iz = ix + 2;

      Vector N = nodes[i]-nodes[i-1];
      N.normalize();
      if(i>1)
      {
        // change direction of lagrange modifier normal
        //A.getLagrangeModifiers()[i] = N;
      }

      Vector V = Vector(v[ix], v[iy], v[iz]);
      Vector dV = cml::dot(N,V) * N;
      v[ix] -= dV[0];
      v[iy] -= dV[1];
      v[iz] -= dV[2];

      Real l = segmentLength;
      Vector X0 = Vector(x[ix], x[iy], x[iz]);
      Vector X1 = Vector(x[ix-3], x[iy-3], x[iz-3]);
      Vector dX = N * ((X1-X0).length() - l);
      x[ix] -= dX[0];
      x[iy] -= dX[1];
      x[iz] -= dX[2];

    
    }
  }

  // additional damping (to be sure...)
  for(int i = 0; i < n*3; i++)
  {
    ap[i] *= 0.99;
    v[i] *= 0.99;
  }

  // update node positions based on x
  for(int i = 0; i < n; i++)
  {
    int ix = i * 3 + 0;
    int iy = ix + 1;
    int iz = ix + 2;

    nodes[i][0] = x[ix];
    nodes[i][1] = x[iy];
    nodes[i][2] = x[iz];

    if( i < numNodes - 1)
      normals[i] = (nodes[i]-nodes[i-1]).normalize();
    else
      normals[i] = normals[i-1];
  }

  return error;
}

template<typename Real>
Real BendingNeedleModel<Real>::simulateImplicitDynamic( Real _dt )
{

  dt = _dt;
  totaltime += dt;
  int n = numNodes;

  /*
      Real offset = (sin(totaltime-3.141/2.0)+1.0)*2.5;
      Vector sX = Vector(10,5,0);
      kSpring = offset / 10.0;
      */

  //cml::vectord F(n*3+n); // sum of forces
  //cml::vectord b(n*3+n); // resulting vector for LSE

  //Real m = 0.001;

  ao = ap;
  ap.zero();

  // fill the Jacobian dF_dx
  updateJacobianForce();

  // fill the Jacobian dF_dv
  //dF_dv = alphaC * M + betaC * dF_dx;
  updateJacobianVelocity();


  // fill A: A = ( M - dF_dx * dt**2 * beta - dF_dv* dt * gamma)
  updateSystemMatrix_A();

  updateResultVector_b();

  cg();

  if(debugOut) std::cout<<"A^-1 * b = ap = "<<ap<<std::endl;

  Real error = updateStep();

  if(debugOut) std::cout<<"x: "<<std::endl<<x<<std::endl;
  if(debugOut) std::cout<<"v: "<<std::endl<<v<<std::endl;
  if(debugOut) std::cout<<"ap: "<<std::endl<<ap<<std::endl;


  if(debugOut) std::cout<<"time: "<<totaltime<<" Error: "<<error<<std::endl; 

  if( error != error ) throw std::runtime_error("error is NaN!");

  return error;
}

template<typename Real>
Real BendingNeedleModel<Real>::simulateImplicitStatic( Real _dt )
{
  Real error = 0;

  dt = _dt;
  totaltime += dt;

  b.zero();
  A.getSystemMatrix().zero();
  // Set Matrix A
  //updateJacobianForce();
  updateJacobianForceLazy();

  for(int j = 0; j < dF_dx.getSize(); j++)
  {
    for(int i = 0; i < dF_dx.getBandwidth(); i++)
    {
      A.getSystemMatrix()._at(i,j) = dF_dx._at(i,j);
    }
  }

  if(debugOut) std::cout<<"A: "<<A.getSystemMatrix()<<std::endl;

  // set vector b
  addForcesToB();
  b*=-1.0;
  if(debugOut) std::cout<<"b: "<<std::endl<<b<<std::endl;

  // perform conjugate gradient
  cg();

  if(debugOut) std::cout<<"ap: "<<ap<<std::endl;

  //x = ap;
  ap *= dt;
  if(debugOut) std::cout<<"ap: "<<ap<<std::endl;
  // update node positions based on x
  for(int i = 0; i < numNodes; i++)
  {
    int ix = i * 3 + 0;
    int iy = ix + 1;
    int iz = ix + 2;

    Vector xnew( x[ix]+ap[ix], x[iy]+ap[iy], x[iz]+ap[iz] );
    error += (xnew-nodes[i]).length();
    nodes[i] = xnew;
    if( i < numNodes - 1 && i > 0)
      normals[i] = (nodes[i]-nodes[i-1]).normalize();
    if( i == 0 )
      normals[i] = normals[i+1];
    else
      normals[i] = normals[i-1];
  }
  x+=ap;
  if(debugOut) std::cout<<"x: "<<x<<std::endl;


  return error;
}

template<typename Real>
void BendingNeedleModel<Real>::addLagrangeModifier( int nodeIndex, Vector N )
{
  A.getLagrangeModifiers()[nodeIndex] = N;
  b.resize(numNodes * 3 + A.getLagrangeModifiers().size() );
  b.zero();
  if(debugOut) std::cout<<"b size 2 "<<b.size()<<std::endl;

  

}

template<typename Real>
void BendingNeedleModel<Real>::setSpring( int nodeIndex, Vector pos, Real k )
{
  springs[nodeIndex] = Spring<Real>( pos, k );
}

template<typename Real>
const std::vector<Vector>& BendingNeedleModel<Real>::getX() const
{
  return nodes;
}

template<typename Real>
Real BendingNeedleModel<Real>::getSegmentLength( ) const
{
  return segmentLength;
}

template<typename Real>
void BendingNeedleModel<Real>::setBasePosition( const Vector& pos )
{
  basePosition = pos;
  //x[0] = basePosition[0];
  //x[1] = basePosition[1];
  //x[2] = basePosition[2];


  Vector secondPosition = basePosition + baseDirection * segmentLength;
  //x[3] = secondPosition[0];
  //x[4] = secondPosition[1];
  //x[5] = secondPosition[2];

  setSpring(0, basePosition, kBase );
  setSpring(1, secondPosition, kBase );

}

template<typename Real>
void BendingNeedleModel<Real>::setBaseDirection( const Vector& dir )
{
  baseDirection = dir;
  baseDirection.normalize();

  Vector secondPosition = basePosition + baseDirection * segmentLength;
  //x[3] = secondPosition[0];
  //x[4] = secondPosition[1];
  //x[5] = secondPosition[2];
  setSpring(1, secondPosition, kBase );
  //A.getLagrangeModifiers()[0] = baseDirection;
  //A.getLagrangeModifiers()[1] = baseDirection;
}

template<typename Real>
void BendingNeedleModel<Real>::setGravity( const Vector& g )
{
  G = g;
}

template<typename Real>
Vector BendingNeedleModel<Real>::getBaseTorque() const
{
  if( springs.find(0) != springs.end() && springs.find(1) != springs.end() )
  {
    const Spring<Real>& first = springs.find(0)->second;
    const Spring<Real>& second = springs.find(1)->second;
    Vector lever =  first.x - second.x;
    Vector force = calcSpring(nodes[0], first.x, first.k);
    Vector torque = cml::cross(lever, force);
    Vector force2 = calcSpring(nodes[1], first.x, first.k);
    torque += cml::cross(-lever, force2);
    return torque;
  }
  return Vector( 0, 0 ,0 );
}

template<typename Real>
Vector BendingNeedleModel<Real>::getBaseForce() const
{
  return calcFNext(0,kNeedle) + calcF(1,kNeedle);
}


template<typename Real>
Real BendingNeedleModel<Real>::getTotalLength()
{
  Real l = 0;
  for(int i = 0; i < numNodes-1; i++)
  {
    l += (nodes[i]-nodes[i+1]).length();
  }
  return l;
}

template<typename Real>
void BendingNeedleModel<Real>::reset()
{
  for( int i = 0; i < numNodes; i++ )
  {
    nodes[i] = (Real)i * segmentLength * baseDirection + basePosition;
    x[i*3+0] = nodes[i][0];
    x[i*3+1] = nodes[i][1];
    x[i*3+2] = nodes[i][2];

    springs.clear();

    normals[i] = baseDirection.normalize() * -1.0;
    nodesLastUpdate[i] = Vector(10000, 10000, 10000);
  }
}


template class BendingNeedleModel<double>;
template class BendingNeedleModel<float>;