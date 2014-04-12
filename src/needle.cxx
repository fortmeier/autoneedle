
#include <iostream>

#include "needle.h"
#include "generatedCode.h"


SparseDiagonalMatrix::SparseDiagonalMatrix( int m, int b ) :
  size(m),
  bandwidth(b),
  bandwidth_2(b/2)
{

  values = new double[m*bandwidth];
  zero();

}

SparseDiagonalMatrix::~SparseDiagonalMatrix()
{
	delete[] values;
}

void SparseDiagonalMatrix::zero()
{
  for( int i = 0; i < bandwidth; i++ )
  {
    for( int j = 0; j < size; j++ )
    {
      _at(i,j) = 0;
    }
  }
}

std::ostream& operator<< ( std::ostream &out, SparseDiagonalMatrix &matrix )
{
  for( int j = 0; j < matrix.size; j++ )
  {
    out << "[ ";
    for( int k = matrix.bandwidth_2; k < j && k < matrix.size - matrix.bandwidth_2-1; k++) 
    {
      out << ". ";
    }

    for( int i = 0; i < matrix.bandwidth; i++ )
    {
      out << matrix._at(i,j) << " ";
    }

    for( int k = matrix.bandwidth; k < matrix.size - j -1 && k < matrix.size - matrix.bandwidth; k++) 
    {
      out << ". ";
    }

    out << " ]\n";
  }
  return out;
}

double& SparseDiagonalMatrix::_at( int i, int j )
{
  return values[ i + bandwidth*j ];
}

double& SparseDiagonalMatrix::operator() ( int i, int j )
{
  if( j <= bandwidth_2) return _at(i,j);
  else if( j >= size - bandwidth_2) return _at(i-size+bandwidth,j);
  return _at(i-j+bandwidth_2,j);
}

cml::vectord SparseDiagonalMatrix::operator* (cml::vectord x)
{
  cml::vectord r(x.size());

  for(int j = 0; j < bandwidth_2+1; j++ )
  {
    r[j] = 0;
    for( int i = 0; i < bandwidth; i++ )
    {
      r[j] += (*this)(i,j) * x[i];
    }

  }

  for(int j = bandwidth_2+1; j < size - bandwidth_2-1; j++ )
  {
    r[j] = 0;
    for( int i = 0; i < bandwidth; i++ )
    {
      r[j] += _at(i,j) * x[i+j-bandwidth_2];
    }

  }

  for(int j = size - bandwidth_2-1; j < size; j++ )
  {
    r[j] = 0;
    for( int i = 0; i < bandwidth; i++ )
    {
      r[j] += _at(i,j) * x[i+size-bandwidth];
    }

  }

  return r;
}

void testMatrix() {
  SparseDiagonalMatrix A( 10, 5 );
  A(0,0) = 1;
  A(1,1) = 2;
  A(0,2) = 3;
  A(3,3) = 4;
  A(4,4) = 5;
  A(5,5) = 6;
  A(6,6) = 7;
  A(9,7) = 8;
  A(8,8) = 9;
  A(9,9) = 3;

  cml::vectord r(10);
  r[0] = 11;
  r[1] = 22;
  r[2] = 33;
  r[3] = 44;
  r[4] = 55;
  r[5] = 66;
  r[6] = 77;
  r[7] = 88;
  r[8] = 99;
  r[9] = 11;
  std::cout<<A<<std::endl;
  std::cout<<r<<std::endl;
  std::cout<<A*r<<std::endl;
  exit(0);

}

BendingNeedleModel::BendingNeedleModel() :
  numNodes( 10 ),
  dF_dx_new( numNodes*3, 5*3 ),
  x( numNodes ),
  f( numNodes ),
  v( numNodes ),
  a( numNodes ),
  totaltime( 0 ),
  debugOut( true ),
  kSpring(1.0),
  kNeedle(1000.0)
{

  //testMatrix();
  for(int i = 0; i < numNodes; i++)
  {
     x[i] = Vector(i, 0.0f,0);
     f[i] = Vector(0,0,0);
     v[i] = Vector(0,0,0);
     a[i] = Vector(0,0,0);
  }

}

Vector BendingNeedleModel::calcF(int i, double k)
{
  Vector f = Vector(
    needle_Fx(x[i-1],x[i],x[i+1], k),
    needle_Fy(x[i-1],x[i],x[i+1], k),
    needle_Fz(x[i-1],x[i],x[i+1], k)
  );
  return  f;
}

Vector BendingNeedleModel::calcFNext(int i, double k)
{
  Vector f = Vector(
    needle_Fx_next(x[i],x[i+1],x[i+2], k),
    needle_Fy_next(x[i],x[i+1],x[i+2], k),
    needle_Fz_next(x[i],x[i+1],x[i+2], k)
  );
  return  f;
}

Vector BendingNeedleModel::calcFPrev(int i, double k)
{
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



void BendingNeedleModel::cg(const cml::matrixd_r& A, const cml::vectord& b, cml::vectord& x )
{
  cml::vectord r = b - A * x;

  cml::vectord p = r;

  double rsold=cml::dot(r,r);
  double rsnew;

  double eps = 0.000001;
  int i = 0;
  do {
    cml::vectord Ap = A * p;
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

void BendingNeedleModel::simulateImplicitChentanez( double dt )
{
  totaltime += dt;
  int n = numNodes;

  cml::vectord xo(numNodes*3+numNodes); // positions from last step
  cml::vectord vo(numNodes*3+numNodes); // velocities from last step
  cml::vectord ao(numNodes*3+numNodes); // accelerations from last step

  xo.zero();
  vo.zero();
  ao.zero();

      double offset = (sin(totaltime-3.141/2.0)+1.0)*2.5;
      Vector sX = Vector(10,5,0);
      kSpring = offset / 10.0;

  for(int i = 0; i < n; i++)
  {
    xo[i*3+0]=x[i][0];
    xo[i*3+1]=x[i][1];
    xo[i*3+2]=x[i][2];
    vo[i*3+0]=v[i][0];
    vo[i*3+1]=v[i][1];
    vo[i*3+2]=v[i][2];
    ao[i*3+0]=a[i][0];
    ao[i*3+1]=a[i][1];
    ao[i*3+2]=a[i][2];
  }

  cml::vectord F(n*3+n); // sum of forces
  //cml::vectord xp(n*3); // new positions
  //cml::vectord vp(n*3); // new velocities
  cml::vectord ap(n*3+n); // new accelerations, thats what we solve the LSE for

  cml::matrixd_r dF_dx(n*3+n, n*3+n); // Jacobian of force with respect to position
  cml::matrixd_r dF_dv(n*3+n, n*3+n); // Jacobian of force with respect to velocity
  cml::matrixd_r M(n*3+n, n*3+n); // Jacobian of force with respect to velocity


  cml::matrixd_r A(n*3+n, n*3+n); // resulting matrix for LSE
  cml::vectord b(n*3+n); // resulting vector for LSE

  double m = 0.001;

  double beta = 0.25;
  double gamma = 0.5;

  // mass matrix
  M.zero();
  for(int i = 0; i<n; i++)
  {
    M(i*3+0,i*3+0) = m;
    M(i*3+1,i*3+1) = m;
    M(i*3+2,i*3+2) = m;
    M(n*3 + i,n*3 + i) = 0;

  } 
  M(3*(n-1)+0,3*(n-1)+0) *=0.5;
  M(3*(n-1)+1,3*(n-1)+1) *=0.5;
  M(3*(n-1)+2,3*(n-1)+2) *=0.5;
  if(debugOut) std::cout<<"M"<<M<<std::endl;

  // fill the Jacobian dF_dx
  dF_dx.zero();
  dF_dx_new.zero();
  for(int i = 0; i < n; i++)
  {
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

    // testcopy
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

    /*if(i>0 && i < n-1) {


      needle_JacobianMatrixAdd( i*3, i*3, dF_dx, x[i-1],x[i],x[i+1]);
      needle_JacobianprevMatrixAdd( i*3, i*3-3, dF_dx, x[i+1], x[i], x[i-1]);
      needle_JacobianprevMatrixAdd( i*3, i*3+3, dF_dx, x[i-1], x[i], x[i+1]);
    }*/
  } 

  spring_JacobianMatrixAdd( (n-1)*3, (n-1)*3, dF_dx, x[n-1], sX, kSpring ); 

  //dF_dx *= 1.0;
  //dF_dx.transpose();
  if(debugOut) std::cout<<"dF_dx: "<<std::endl<<dF_dx<<std::endl;
  if(debugOut) std::cout<<"dF_dx_new: "<<std::endl<<dF_dx_new<<std::endl;

  // fill the Jacobian dF_dv
  /*
  for(int i = 0; i < n; i++)
  {
    if(i>1 && i < n-1) {
      needle_dfMatrixSetter( i*3, i*3, dF_dx, x[i-1],x[i],x[i+1]);
    }
    if(i>1 && i < n) {
      needle_dfprevMatrixSetter( i*3, i*3-3, dF_dx, x[i], x[i-1], x[i-2]);
    }
    if(i>1 && i < n-2) {
      needle_dfprevMatrixSetter( i*3, i*3+3, dF_dx, x[i], x[i+1], x[i+2]);
 
    }
    if(i==n-1) {
      needle_df2MatrixSetter( i*3-3, i*3, dF_dx, x[i-1], x[i] );
      needle_df2MatrixSetter( i*3, i*3, dF_dx, x[i], x[i-1] );
    } 
  }*/


  // OK dF_dv = -50.5 * M + 0.005 * dF_dx;
  dF_dv = -50.0 * M + 0.005 * dF_dx;
  //dF_dv.zero();

  if(debugOut) std::cout<<"dF_dv: "<<std::endl<<dF_dv<<std::endl;

  // fill A: A = ( M - dF_dx * dt**2 * beta - dF_dv* dt * gamma)
  A = M - dF_dx * dt * dt * beta - dF_dv* dt * gamma;





  // fill b: F + dt * dF_dx * ( v + dt * ( 0.5 - beta) * a ) + dF_dv * dt * (1-gamma) * a
  // but first, F
  F.zero();
  for(int i = 0; i < n; i++)
  {
    Vector f(0,0,0);

    if( i > 0 && i < n-1 ) f += calcF(i, kNeedle);
    if( i > 1 ) f += calcFPrev(i, kNeedle);
    if( i < n-2 ) f += calcFNext(i, kNeedle);

    f[1] += -9.81 * M(i*3+1,i*3+1);

    if(i==n-1) {

      Vector FP = calcSpring(x[i], sX, kSpring);
      if(true) std::cout<<"offset: "<<sX<<std::endl;
      f += FP;
    }

    F[i*3 + 0] = f[0]; 
    F[i*3 + 1] = f[1]; 
    F[i*3 + 2] = f[2]; 
  }

  if(debugOut) std::cout<<"F: "<<F<<std::endl;

  cml::vectord b1 = dt * dF_dx * ( vo + dt * ( 0.5 - beta) * ao );
  cml::vectord b2 = dF_dv * dt * (1.0-gamma) * ao;

  if(debugOut) std::cout<<"b1: "<<b1<<std::endl;
  if(debugOut) std::cout<<"b2: "<<b2<<std::endl;

  b = F + b1 + b2;

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
    b[n*3+i] = 0; //-4*cml::dot(N, v[i])*dt - dt*dt*cml::dot(N,a[i]);

    //if(i==0 || i==1) b[n*3+i] += 10.1*sin(totaltime);

  }

  if(debugOut) std::cout<<"A: "<<std::endl<<A<<std::endl<<std::endl;


  if(debugOut) std::cout<<"b: "<<b<<std::endl;

  if(true)
  {
    cg(A,b,ap);
  } else {
    ap = inverse(A) * b;
  }

  if(debugOut) std::cout<<"A^-1 * b = ap = "<<ap<<std::endl;


  // length enforcement
  bool useEnforcement = true;
  
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
    Vector N = x[i]-x[i-1];
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

  // update std::vector x (giving x_plus)
  for(int i = 0; i < n; i++)
  {
    x[i][0] += dt * v[i][0] + dt*dt*((0.5-beta)*a[i][0] + beta*ap[i*3 + 0]);
    x[i][1] += dt * v[i][1] + dt*dt*((0.5-beta)*a[i][1] + beta*ap[i*3 + 1]);
    x[i][2] += dt * v[i][2] + dt*dt*((0.5-beta)*a[i][2] + beta*ap[i*3 + 2]);
  }
  // update std::vector v (giving v_plus)
  for(int i = 0; i < n; i++)
  {
    v[i][0] += dt*((1.0-gamma)*a[i][0] + gamma*ap[i*3 + 0]);
    v[i][1] += dt*((1.0-gamma)*a[i][1] + gamma*ap[i*3 + 1]);
    v[i][2] += dt*((1.0-gamma)*a[i][2] + gamma*ap[i*3 + 2]);
  }
  // update std::vector a (giving a_plus)

  double r = 0;
  for(int i = 0; i < n; i++)
  {
    r += pow(a[i][0]-ap[i*3 + 0], 2) + pow(a[i][1]-ap[i*3 + 1], 2) + pow(a[i][2]-ap[i*3 + 2], 2);
    a[i][0] = ap[i*3 + 0];
    a[i][1] = ap[i*3 + 1];
    a[i][2] = ap[i*3 + 2];
  }

  // length enforcement
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
  for(int i = 0; i < n; i++)
  {
    a[i] *= 0.99;
    v[i] *= 0.99;
  }
  std::cout<<"time: "<<totaltime<<" Error: "<<r<<std::endl; 
}

const std::vector<Vector>& BendingNeedleModel::getX() const
{
  return x;
}
