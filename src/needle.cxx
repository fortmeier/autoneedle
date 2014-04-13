
#include <iostream>
#include <cassert>
#include <cstring>

#include "needle.h"
#include "generatedCode.h"


NeedleMatrix::NeedleMatrix( int nNodes ) :
  numNodes( nNodes ),
  numLagrangeModifiers( 0 ), // no lagrange modifiers at initialization 
  A( numNodes*3, 5*3 ),
  r_b( nNodes * 3 ),
  r_lambda( 0 ), // no lagrange modifiers at initialization
  r( nNodes * 3 )
{

}

cml::vectord NeedleMatrix::operator* (const cml::vectord& x) const
{
  assert( "length of vector is not okay!" && x.size() == numNodes * 3 + numLagrangeModifiers );
  if( r_lambda.size() != numLagrangeModifiers ) assert ("not implemented yet!" && false );

  cml::vectord x_needle = x;
//  std::cout<<x_needle<<std::endl;
  x_needle.resize( numNodes * 3);
//  std::cout<<"x_needle: "<<x_needle<<std::endl;
  r_b = A * x_needle;
//  std::cout<<x_needle<<std::endl;
//  std::cout<<"r_b: "<<r_b<<std::endl;

/* todo readd lagrange multipliers
    // add lagrange multipliers
  for(int i = 0; i < n; i++) 
  {
    Vector N(1,0,0);
    if(i==0) N = Vector(0,1,0);
    if(i==1) N = Vector(0,1,0);

    if(i>1) N = (x[i]-x[i-1]).normalize();
#if 1
    A(n*3+i, i*3+0) = N[0];
    A(n*3+i, i*3+1) = N[1];
    A(n*3+i, i*3+2) = N[2];
    A(i*3+0, n*3+i) = N[0];
    A(i*3+1, n*3+i) = N[1];
    A(i*3+2, n*3+i) = N[2];
#endif

    //A(n*3+i, n*3+i) = 1;
    //b[n*3+i] = 0; //-4*cml::dot(N, v[i])*dt - dt*dt*cml::dot(N,a[i]);

    //if(i==0 || i==1) b[n*3+i] += 10.1*sin(totaltime);

  } */

  size_t size_r_b = numNodes * 3 * sizeof(double); 
  memcpy( r.data(), r_b.data(), size_r_b );

  size_t size_r_lambda = numLagrangeModifiers * sizeof(double);
  std::memcpy( r.data() + size_r_b, r_lambda.data(), size_r_lambda);
  return x;
}

SparseDiagonalMatrix& NeedleMatrix::getSystemMatrix()
{
  return A;
}


BendingNeedleModel::BendingNeedleModel() :
  numNodes( 6 ),
  A( numNodes ),
  b( numNodes ),
  dF_dx( numNodes*3, 5*3 ),
  dF_dv( numNodes*3, 5*3 ),
  nodes( numNodes ),
  x(numNodes*3), // positions from last step
  v(numNodes*3), // velocities from last step
  ao(numNodes*3), // accelerations from last step
  ap(numNodes*3), // accelerations from last step
  m(numNodes*3), // mass of node in x y z direction
  totaltime( 0 ),
  debugOut( true ),
  kSpring(1.0),
  kNeedle(1000.0)
{

  //testMatrix();
  for(int i = 0; i < numNodes; i++)
  {
     nodes[i] = Vector(i, 0.0f,0);
     x[i*3+0] = nodes[i][0];
     x[i*3+1] = nodes[i][1];
     x[i*3+2] = nodes[i][2];

     m[i*3+0] = 0.001;
     m[i*3+1] = 0.001;
     m[i*3+2] = 0.001;

  }

}

Vector BendingNeedleModel::calcF(int i, double k)
{
  std::vector<Vector>& x = nodes;
  Vector f = Vector(
    needle_Fx(x[i-1],x[i],x[i+1], k),
    needle_Fy(x[i-1],x[i],x[i+1], k),
    needle_Fz(x[i-1],x[i],x[i+1], k)
  );
  return  f;
}

Vector BendingNeedleModel::calcFNext(int i, double k)
{
  std::vector<Vector>& x = nodes;
  Vector f = Vector(
    needle_Fx_next(x[i],x[i+1],x[i+2], k),
    needle_Fy_next(x[i],x[i+1],x[i+2], k),
    needle_Fz_next(x[i],x[i+1],x[i+2], k)
  );
  return  f;
}

Vector BendingNeedleModel::calcFPrev(int i, double k)
{
  std::vector<Vector>& x = nodes;
  Vector f = Vector(
    needle_Fx_next(x[i],x[i-1],x[i-2], k),
    needle_Fy_next(x[i],x[i-1],x[i-2], k),
    needle_Fz_next(x[i],x[i-1],x[i-2], k)
  );
  return  f;
}

Vector BendingNeedleModel::calcSpring(Vector a, Vector b, double k)
{
  Vector f = Vector(
    spring_FSx(a,b,k),
    spring_FSy(a,b,k),
    spring_FSz(a,b,k)
  );
  return  f;
}




void BendingNeedleModel::cg()
{
  cml::vectord& x = ap;
  cml::vectord r = b - A * x;

  cml::vectord p = r;

  double rsold=cml::dot(r,r);
  double rsnew;

  double eps = 0.000001;
  int i = 0;
  do {
    cml::vectord Ap = A * p;
    std::cout<<"1:"<<p<<std::endl;
    std::cout<<"2:"<<Ap<<std::endl;
    double alpha = rsold/cml::dot(p,Ap);
    x = x + alpha * p;
    r = r-alpha*Ap;
    rsnew = cml::dot(r,r);
    p=r+p*(rsnew/rsold);
    rsold=rsnew; 
    i++;

  } while (sqrt(rsnew) > eps && i < 1000);
  std::cout<<"cg needed "<<i<<" iterations"<<std::endl;
}

void BendingNeedleModel::updateJacobianForce()
{
  dF_dx.zero();
  int n = numNodes;
  for(int i = 0; i < numNodes; i++)
  {
    std::vector<Vector>& x = nodes;

    // pF_pxi
    if(i>=2) needle_JacobianCCMatrixAdd( i*3, i*3, dF_dx, x[i-2], x[i-1], x[i], kNeedle);
    if(i>=1 && i < n-1) needle_JacobianBBMatrixAdd( i*3, i*3, dF_dx, x[i-1], x[i], x[i+1], kNeedle);
    if(i>=0 && i< n-2) needle_JacobianAAMatrixAdd( i*3, i*3, dF_dx, x[i], x[i+1], x[i+2], kNeedle);

    // pF_pxi+1
    if(i>=1 && i < n-1) needle_JacobianBCMatrixAdd( i*3, i*3+3, dF_dx, x[i-1], x[i], x[i+1], kNeedle);
    if(i>=0 && i < n-2) needle_JacobianABMatrixAdd( i*3, i*3+3, dF_dx, x[i], x[i+1], x[i+2], kNeedle);

    // pF_pxi+2
    if(i>=0 && i < n-2) needle_JacobianACMatrixAdd( i*3, i*3+6, dF_dx, x[i], x[i+1], x[i+2], kNeedle);

    // pF_pxi-1
    if(i>=1 && i < n-1) needle_JacobianBAMatrixAdd( i*3, i*3-3, dF_dx, x[i-1], x[i], x[i+1], kNeedle);
    if(i>=2 && i < n) needle_JacobianCBMatrixAdd( i*3, i*3-3, dF_dx, x[i-2], x[i-1], x[i], kNeedle);

    // pF_pxi+2
    if(i>=2 && i < n) needle_JacobianCCMatrixAdd( i*3, i*3-6, dF_dx, x[i-2], x[i-1], x[i], kNeedle);

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

  //dF_dx *= 1.0;
  //dF_dx.transpose();
  if(debugOut) std::cout<<"dF_dx: "<<std::endl<<dF_dx<<std::endl;

}

void BendingNeedleModel::updateJacobianVelocity()
{
  //dF_dv = alphaC * M + betaC * dF_dx;
  // OK dF_dv = -50.5 * M + 0.005 * dF_dx;
  dF_dv.zero();

  double alphaC = 0.0;
  double betaC = 0.0;
  if(debugOut) std::cout<<"dF_dv: "<<std::endl<<dF_dv<<std::endl;

}

void BendingNeedleModel::updateSystemMatrix_A()
{
  // fill A: 
  //    = M - dF_dx * dt * dt * beta - dF_dv* dt * gamma;
  for(int j = 0; j < A.getSystemMatrix().getSize(); j++)
  {
    for(int i = 0; i < A.getSystemMatrix().getBandwidth(); i++)
    {
      double a = -(dF_dx(i,j)* dt * dt * 0.25) - (dF_dv(i,j) * dt * 0.5);
      A.getSystemMatrix()._at(i,j) = a;
    }
    A.getSystemMatrix()._at(j,j) += m[j];
  }
  if(debugOut) std::cout<<"A: "<<std::endl<<A.getSystemMatrix()<<std::endl;

}

void BendingNeedleModel::updateResultVector_b()
{
  // fill b: F + dt * dF_dx * ( v + dt * ( 0.5 - beta) * a ) + dF_dv * dt * (1-gamma) * a
  cml::vectord tmp = (( v + dt * ( 0.25 ) * ao ) * dt);
  b  = dF_dx * tmp;
  tmp = (ao * dt * 0.5);
  b += dF_dv * tmp;

  int n = numNodes;

  for(int i = 0; i < n; i++)
  {
    Vector f(0,0,0);

    if( i > 0 && i < n-1 ) f += calcF(i, kNeedle);
    if( i > 1 ) f += calcFPrev(i, kNeedle);
    if( i < n-2 ) f += calcFNext(i, kNeedle);

    f[1] += -9.81 * m[i*3+1];

    /*if(i==n-1) {

      Vector FP = calcSpring(x[i], sX, kSpring);
      if(true) std::cout<<"offset: "<<sX<<std::endl;
      f += FP;
    }*/

    b[i*3 + 0] += f[0]; 
    b[i*3 + 1] += f[1]; 
    b[i*3 + 2] += f[2];

  }
  if(debugOut) std::cout<<"b: "<<std::endl<<b<<std::endl;
}

double BendingNeedleModel::updateStep()
{
  // length enforcement
  bool useEnforcement = false;

  double error = 0;
  int n = numNodes;
  
  if(useEnforcement) 
  {
    ap[0] = 0;
    ap[1] = 0;
    ap[2] = 0;
    ap[3] = 0;
    ap[4] = 0;
    ap[5] = 0;
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
    double ux = dt * v[ix] + dt*dt*(0.25*ao[ix] + 0.25*ap[ix]);
    double uy = dt * v[iy] + dt*dt*(0.25*ao[iy] + 0.25*ap[iy]);
    double uz = dt * v[iz] + dt*dt*(0.25*ao[iz] + 0.25*ap[iz]);

    x[ix] += ux;
    x[iy] += uy;
    x[iz] += uz;

    error += ux*ux + uy*uy + uz*uz;

    // update std::vector v (giving v_plus)
    v[ix] += dt*(0.5*ao[ix] + 0.5*ap[ix]);
    v[iy] += dt*(0.5*ao[iy] + 0.5*ap[iy]);
    v[iz] += dt*(0.5*ao[iz] + 0.5*ap[iz]);

    nodes[i][0] = x[ix];
    nodes[i][1] = x[iy];
    nodes[i][2] = x[iz];
  }
  // todo readd errorcalc  r += pow(a[i][0]-ap[i*3 + 0], 2) + pow(a[i][1]-ap[i*3 + 1], 2) + pow(a[i][2]-ap[i*3 + 2], 2);

  // length enforcement
  /* todo readd
  if(useEnforcement) 
  {
  for(int i = 1; i < n; i++)
  {
    Vector N = x[i]-x[i-1];
    N.normalize();
    Vector p = cml::dot(N,v[i]) * N;
    v[i] -= p;

    double l = 1.0;
    x[i] -= N * ((x[i]-x[i-1]).length() - l);
    
  }
  }
  */
  /*for(int i = 0; i < n; i++)
  {
    ap[i] *= 0.99;
    v[i] *= 0.99;
  }*/

  return error;
}

void BendingNeedleModel::simulateImplicitChentanez( double _dt )
{
  dt = _dt;
  totaltime += dt;
  int n = numNodes;

  /*
      double offset = (sin(totaltime-3.141/2.0)+1.0)*2.5;
      Vector sX = Vector(10,5,0);
      kSpring = offset / 10.0;
      */

  //cml::vectord F(n*3+n); // sum of forces
  //cml::vectord b(n*3+n); // resulting vector for LSE

  //double m = 0.001;

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

  double error = updateStep();

  if(debugOut) std::cout<<"x: "<<std::endl<<x<<std::endl;
  if(debugOut) std::cout<<"v: "<<std::endl<<x<<std::endl;
  if(debugOut) std::cout<<"ap: "<<std::endl<<ap<<std::endl;


  std::cout<<"time: "<<totaltime<<" Error: "<<error<<std::endl; 
}

const std::vector<Vector>& BendingNeedleModel::getX() const
{
  return nodes;
}
