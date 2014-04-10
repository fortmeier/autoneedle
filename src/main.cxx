#include <iostream>

#include "generatedCode.h"
#include "mathheader.h"
#include "Rendering.h"


#define numNodes 10
  
std::vector<Vector> x(numNodes);
std::vector<Vector> f(numNodes);
std::vector<Vector> v(numNodes);
std::vector<Vector> a(numNodes);



Rendering* r;

Vector calcF(int i, bool curr = true, bool next = true, bool prev = true)
{
      Vector fb = Vector(
         needle_df_dx(x[i-1],x[i],x[i+1]),
         needle_df_dy(x[i-1],x[i],x[i+1]),
         needle_df_dz(x[i-1],x[i],x[i+1])
      );
      Vector fa = Vector(
         needle_df_dxnext(x[i-2],x[i-1],x[i]),
         needle_df_dynext(x[i-2],x[i-1],x[i]),
         needle_df_dznext(x[i-2],x[i-1],x[i])
      );
      Vector fc = Vector(
         needle_df_dxprev(x[i],x[i+1],x[i+2]),
         needle_df_dyprev(x[i],x[i+1],x[i+2]),
         needle_df_dzprev(x[i],x[i+1],x[i+2])
      );
  Vector f(0,0,0);
  if(curr) f += fb;
  if(next) f += fa;
  if(prev) f += fc;
  return  f;
}

void simulateExplicit( double dt )
{
  for(int i = 2; i < x.size(); i++)
  {
    double m = 0.001;



    if(i != x.size()-1 ) 
    {
      f[i] = calcF(i);
    } else {
      f[i] = calcF(i, false, true, false );
    }
    f[i] *= -1.0;
  
    // double g = 9.81;
    double g = 9.81;

    f[i] += Vector(0,-1,0) * (m * g);
 
    // damping
    v[i] *= 0.99;
 
    // explicit update
    a[i] = f[i] / m * 1.0;

    v[i] = v[i] + a[i] * dt;

    x[i] = x[i] + f[i]*0.1;//v[i] * dt;// + a[i] * dt * dt;
  }
}

void cg(const cml::matrixd_r& A, const cml::vectord& b, cml::vectord& x )
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

  } while (sqrt(rsnew) > eps );
  std::cout<<"cg needed "<<i<<" iterations"<<std::endl;
}

void simulateImplicit( double dt )
{
  int n = numNodes;

  cml::vectord b(n*3);
  cml::vectord x_old(n*3);
  cml::vectord f_old(n*3);

  double m = 10.0;

  // fill A
  cml::matrixd_r A(n*3, n*3);
  cml::identity_transform( A );
  for(int i = 0; i < n; i++)
  {
    if(i>1 && i < n-1) {
      /*A(i*3 + 0,i*3 + 0) = needle_df_dx_dx(x[i-1],x[i],x[i+1]); 
      A(i*3 + 1,i*3 + 1) = needle_df_dy_dy(x[i-1],x[i],x[i+1]); 
      A(i*3 + 2,i*3 + 2) = needle_df_dz_dz(x[i-1],x[i],x[i+1]); */
      needle_dfMatrixSetter( i*3, i*3,A, x[i-1],x[i],x[i+1]);
    }
    if(i>1 && i < n) {
      /*A(i*3 + 0,i*3 + 0 - 3) = needle_df_dx_dxprev(x[i],x[i-1],x[i-2]); 
      A(i*3 + 1,i*3 + 1 - 3) = needle_df_dy_dyprev(x[i],x[i-1],x[i-2]); 
      A(i*3 + 2,i*3 + 2 - 3) = needle_df_dz_dzprev(x[i],x[i-1],x[i-2]); */
      needle_dfprevMatrixSetter( i*3, i*3-3,A, x[i], x[i-1], x[i-2]);
    }
    if(i>1 && i < n-2) {
      /*A(i*3 + 0,i*3 + 0 + 3) = needle_df_dx_dxprev(x[i],x[i+1],x[i+2]); 
      A(i*3 + 1,i*3 + 1 + 3) = needle_df_dy_dyprev(x[i],x[i+1],x[i+2]); 
      A(i*3 + 2,i*3 + 2 + 3) = needle_df_dz_dzprev(x[i],x[i+1],x[i+2]); */
      needle_dfprevMatrixSetter( i*3, i*3+3,A, x[i], x[i+1], x[i+2]);
 
    }
    if(i==n-1) {
      /*
      A(i*3 + 0 - 3, i*3 + 0) = needle_df2_dx_dx(x[i-1],x[i]); 
      A(i*3 + 1 - 3, i*3 + 1) = needle_df2_dy_dy(x[i-1],x[i]); 
      A(i*3 + 2 - 3, i*3 + 2) = needle_df2_dz_dz(x[i-1],x[i]);*/
      needle_df2MatrixSetter( i*3-3, i*3,A, x[i-1], x[i] );

      /*A(i*3 + 0,i*3 + 0 - 3) = needle_df2_dx_dx(x[i-1],x[i]); 
      A(i*3 + 1,i*3 + 1 - 3) = needle_df2_dy_dy(x[i-1],x[i]); 
      A(i*3 + 2,i*3 + 2 - 3) = needle_df2_dz_dz(x[i-1],x[i]);*/
      /*A(i*3 + 0,i*3 + 0) = needle_df2_dx_dx(x[i],x[i-1]); 
      A(i*3 + 1,i*3 + 1) = needle_df2_dy_dy(x[i],x[i-1]); 
      A(i*3 + 2,i*3 + 2) = needle_df2_dz_dz(x[i],x[i-1]);  */
      needle_df2MatrixSetter( i*3, i*3,A, x[i], x[i-1] );

    } 
//    else if(i==0) {
      //A(i*3 + 0,i*3 + 0) = -needle_df2_dx_dx(x[1],x[0]); 
      //A(i*3 + 1,i*3 + 1) = -needle_df2_dy_dy(x[1],x[0]); 
      //A(i*3 + 2,i*3 + 2) = -needle_df2_dz_dz(x[1],x[0]); 
//    }
    // add mass as damping term
    if(i>0 && i<n-1) {
      A(i*3 + 0,i*3 + 0) += m; 
      A(i*3 + 1,i*3 + 1) += m; 
      A(i*3 + 2,i*3 + 2) += m; 
    } else {
      A(i*3 + 0,i*3 + 0) += m;1 * 0.5; 
      A(i*3 + 1,i*3 + 1) += m;1 * 0.5; 
      A(i*3 + 2,i*3 + 2) += m;1 * 0.5; 
    }

  } 

  std::cout<<A<<std::endl;


  // fill b
  for(int i = 0; i < n; i++)
  {
    b[i*3 + 0] = 0; 
    b[i*3 + 1] = -9.81 * 0.001; 
    b[i*3 + 2] = 0; 

    //if(i==n-1) b[i*3+1]*=0.5;

    x_old[i*3 + 0] = x[i][0]; 
    x_old[i*3 + 1] = x[i][1]; 
    x_old[i*3 + 2] = x[i][2]; 

    Vector f;
    if( i > 0 && i < n-1 ) f = +calcF(i, true, true, true) * 1.0;
    if( i == n-1)
    {
      f = Vector (
        needle_df2_dx(x[n-2],x[n-1]), 
        needle_df2_dy(x[n-2],x[n-1]), 
        needle_df2_dz(x[n-2],x[n-1]) 
      );       
      f += calcF(i, false, true, false);
    }
    
    f_old[i*3 + 0] = f[0]; 
    f_old[i*3 + 1] = f[1]; 
    f_old[i*3 + 2] = f[2]; 
  }

  for(int j = 0; j < 3 * 2; j++) {
    f_old[j] = 0;
    b[j] = 0;
  }
  for(int j = 0; j < 3 * 0; j++) {
    f_old[n*3-1-j] = 0;
    b[n*3-1-j] = 0;
  }

  cml::vectord ax = A * x_old;

  b = b + ax - f_old;



  // calc r
  cml::vectord r = b * 0;
  if(false)
  {
    cg(A,b,r);
  } else {
    r = inverse(A) * b;
  }

  std::cout<<"f_old: "<<f_old<<std::endl;
  std::cout<<"b: "<<b<<std::endl;
  std::cout<<"A * x_old: "<<ax<<std::endl;
  std::cout<<"r: "<<r<<std::endl;

  // extract x
  for(int i = 0; i < n; i++)
  {
    x[i][0] = r[i*3 + 0];
    x[i][1] = r[i*3 + 1];
    x[i][2] = r[i*3 + 2];
  }
}

/*
  compute with method of [1]:
  1. Chentanez N, Alterovitz R, Ritchie D, Cho L, Hauser KK, Goldberg K, et al. 
  Interactive simulation of surgical needle insertion and steering. 
  ACM Trans Graph;28(3):1-10. Available from: http://portal.acm.org/citation.cfm?id=1531394

  result should be in the form A ap = b with

  A = ( M - dF_dx * dt**2 * beta - dF_dv* dt * gamma)
  b = F + dF_dx * dt * v + dF_dx * dt**2 * ( 0.5 - beta) * a + dF_dv * dt * (1-gamma) * a
   ?= F + dF_dx * ( dt * v + dt**2 * ( 0.5 - beta) * a ) + dF_dv * dt * (1-gamma) * a
   ?= F + dt * dF_dx * ( v + dt * ( 0.5 - beta) * a ) + dF_dv * dt * (1-gamma) * a

*/




void simulateImplicitChentanez( double dt )
{
  int n = numNodes;

  cml::vectord xo(numNodes*3+numNodes); // positions from last step
  cml::vectord vo(numNodes*3+numNodes); // velocities from last step
  cml::vectord ao(numNodes*3+numNodes); // accelerations from last step

  xo.zero();
  vo.zero();
  ao.zero();

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
  //M(3*(n-1)+0,3*(n-1)+0) *=0.5;
  //M(3*(n-1)+1,3*(n-1)+1) *=0.5;
  //M(3*(n-1)+2,3*(n-1)+2) *=0.5;
  std::cout<<"M"<<M<<std::endl;

  // fill the Jacobian dF_dx
  dF_dx.zero();
  for(int i = 0; i < n; i++)
  {
    if(i>1 && i < n-1) {
      needle_dfMatrixSetter( i*3, i*3, dF_dx, x[i-1],x[i],x[i+1]);
      needle_dfprevMatrixSetter( i*3, i*3-3, dF_dx, x[i+1], x[i], x[i-1]);
      needle_dfprevMatrixSetter( i*3, i*3+3, dF_dx, x[i-1], x[i], x[i+1]);
    }
    if(i==n-1) {
      //needle_dfprevMatrixSetter( i*3, i*3-3, dF_dx, x[i], x[i-1], x[i-2]);
      //needle_dfMatrixSetter( i*3, i*3, dF_dx, x[i], x[i-1], x[i-2]);
      //needle_dfMatrixSetter( i*3+3, i*3, dF_dx, x[i-1],x[i],x[i+1]);
      //needle_df2MatrixSetter( i*3+3, i*3, dF_dx, x[i], x[i+1] );
      //needle_df2MatrixSetter( i*3, i*3, dF_dx, x[i], x[i-1] );
    }
  } 

  dF_dx *= -1.0;
  std::cout<<"dF_dx: "<<std::endl<<dF_dx<<std::endl;

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


  dF_dv = 0.0 * M - 0.005 * dF_dx;
  dF_dv.zero();

  std::cout<<"dF_dv: "<<std::endl<<dF_dv<<std::endl;

  // fill A: A = ( M - dF_dx * dt**2 * beta - dF_dv* dt * gamma)
  A = M - dF_dx * dt * dt * beta - dF_dv* dt * gamma;





  // fill b: F + dt * dF_dx * ( v + dt * ( 0.5 - beta) * a ) + dF_dv * dt * (1-gamma) * a
  // but first, F
  F.zero();
  for(int i = 0; i < n; i++)
  {
    Vector f(0,0,0);

    if( i > 1 && i < n-1 ) f = +calcF(i, true, false, false) * -0.4;
    //if( i > 1 && i < n-1 ) f = +calcF(i, false, true, false) * -1.0;
    //if( i > 1 && i < n-2 ) f = +calcF(i, false, false, true) * -1.0;

    if( i > 1 ) f[1] += -9.81 * M(i*3+1,i*3+1);

    F[i*3 + 0] = f[0]; 
    F[i*3 + 1] = f[1]; 
    F[i*3 + 2] = f[2]; 
  }

  std::cout<<"F: "<<F<<std::endl;

  cml::vectord b1 = dt * dF_dx * ( vo + dt * ( 0.5 - beta) * ao );
  cml::vectord b2 = dF_dv * dt * (1.0-gamma) * ao;

  std::cout<<"b1: "<<b1<<std::endl;
  std::cout<<"b2: "<<b2<<std::endl;

  b = F + b1 + b2;

  // add lagrange multipliers
  for(int i = 0; i < n; i++) 
  {
    Vector N(1,0,0);

    if(i>1) N = (x[i]-x[i-1]).normalize();
    A(n*3+i, i*3+0) = N[0];
    A(n*3+i, i*3+1) = N[1];
    A(n*3+i, i*3+2) = N[2];
    A(i*3+0, n*3+i) = N[0];
    A(i*3+1, n*3+i) = N[1];
    A(i*3+2, n*3+i) = N[2];
    A(n*3+i, n*3+i) = 1;
    b[n*3+i] = 0;

  }

  std::cout<<"A: "<<std::endl<<A<<std::endl<<std::endl;


  std::cout<<"b: "<<b<<std::endl;

  if(false)
  {
    cg(A,b,ap);
  } else {
    ap = inverse(A) * b;
  }

  std::cout<<"A^-1 * b = ap = "<<ap<<std::endl;


  // length enforcement
  bool useEnforcement = true;
  
  if(useEnforcement) 
  {
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
  for(int i = 0; i < n; i++)
  {
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
  for(int i = 1; i < n; i++)
  {
    a[i] *= 0.95;
    v[i] *= 0.99;
  } 
}


bool ex = false;
void simulate()
{

  double dt = 0.01;

  static int c = 0;

  for(int i = 0; i < 1 && c < 30000; i++)
  {
    c++;
    //simulateExplicit(dt);
    //simulateImplicit(dt);
    simulateImplicitChentanez(dt);

  for(int i = 0; i < numNodes; i++)
  {
    std::cout<<i<<":"<<std::endl;
    std::cout<<"x: "<<x[i]<<std::endl;
    std::cout<<"f: "<<f[i]<<std::endl;
    std::cout<<"v: "<<v[i]<<std::endl;
    std::cout<<"a: "<<a[i]<<std::endl;
  }

  }

  r->update(x);
  if(ex) exit(0);
 
}

int main(int argi, char** argv)
{

  if(argi>1) ex = true;

  r = new Rendering();

  for(int i = 0; i < numNodes; i++)
  {
     x[i] = Vector(i, 0.0f,0);
     f[i] = Vector(0,0,0);
     v[i] = Vector(0,0,0);
     a[i] = Vector(0,0,0);
  }

  x[5][1] = 0.0;

  r->setup();
  r->setCallback( simulate );
  r->update(x);
  r->run();  
  

  return 0;
}
