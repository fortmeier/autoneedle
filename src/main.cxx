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

void simulateImplicit( double dt )
{
  int n = (numNodes-3);

  cml::vectord b(n*3);
  cml::vectord x_old(n*3);
  cml::vectord f_old(n*3);


  // fill A
  cml::matrixd_r A(n*3, n*3);
  cml::identity_transform( A );
  for(int i = 0; i < n; i++)
  {
    /*if(i<n-1)
    {
      A(i*3 + 0,i*3 + 0 + 3) = -needle_df_dx_dxprev(x[i+0+2],x[i+1+2],x[i+2+2]); 
      A(i*3 + 1,i*3 + 1 + 3) = -needle_df_dy_dyprev(x[i+0+2],x[i+1+2],x[i+2+2]); 
      A(i*3 + 2,i*3 + 2 + 3) = -needle_df_dz_dzprev(x[i+0+2],x[i+1+2],x[i+2+2]); 
    }
    if(i>0)
    {
      A(i*3 + 0,i*3 + 0 - 3) = -needle_df_dx_dxnext(x[i-2+2],x[i-1+2],x[i+0+2]); 
      A(i*3 + 1,i*3 + 1 - 3) = -needle_df_dy_dynext(x[i-2+2],x[i-1+2],x[i+0+2]); 
      A(i*3 + 2,i*3 + 2 - 3) = -needle_df_dz_dznext(x[i-2+2],x[i-1+2],x[i+0+2]); 
    }*/
    A(i*3 + 0,i*3 + 0) = -needle_df_dx_dx(x[i-1+2],x[i+2],x[i+1+2]); 
    A(i*3 + 1,i*3 + 1) = -needle_df_dy_dy(x[i-1+2],x[i+2],x[i+1+2]); 
    A(i*3 + 2,i*3 + 2) = -needle_df_dz_dz(x[i-1+2],x[i+2],x[i+1+2]); 
  }

  std::cout<<A<<std::endl;


  // fill b
  for(int i = 0; i < n; i++)
  {
    b[i*3 + 0] = 0; 
    b[i*3 + 1] = 0.0981; 
    b[i*3 + 2] = 0; 

    x_old[i*3 + 0] = x[i+2][0]; 
    x_old[i*3 + 1] = x[i+2][1]; 
    x_old[i*3 + 2] = x[i+2][2]; 

    Vector f = -calcF(i+2, true, false, false);
    f_old[i*3 + 0] = f[0]; 
    f_old[i*3 + 1] = f[1]; 
    f_old[i*3 + 2] = f[2]; 
  }

  b = b + A * x_old - f_old;


  // calc r
  cml::vectord r = inverse(A) * b;
  std::cout<<"f_old: "<<f_old<<std::endl;
  std::cout<<"b: "<<b<<std::endl;
  std::cout<<"r: "<<r<<std::endl;

  // extract x
  for(int i = 0; i < n; i++)
  {
    x[i+2][0] = r[i*3 + 0];
    x[i+2][1] = r[i*3 + 1];
    x[i+2][2] = r[i*3 + 2];
  }
}


void simulate()
{
  std::cout<<"updatesimulation"<<std::endl;

  double dt = 0.001;

  for(int i = 0; i < 1; i++)
  {
    //simulateExplicit(dt);
    simulateImplicit(dt);

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
//  exit(0);
  
}

int main()
{

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
