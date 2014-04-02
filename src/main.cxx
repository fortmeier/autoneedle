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

void simulateExplicit( double dt )
{
  for(int i = 2; i < x.size(); i++)
  {
    double m = 0.001;

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

    if(i != x.size()-1 ) 
    {
      f[i] = fa + fb + fc;
    } else {
      Vector f2b = Vector(
         needle_df2_dx(x[i-1],x[i]),
         needle_df2_dy(x[i-1],x[i]),
         needle_df2_dz(x[i-1],x[i])
      );
      f[i] = fa;// + f2b;
      std::cout<<f2b<<std::endl;
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

    x[i] = x[i] + v[i] * dt;// + a[i] * dt * dt;
  }
}

void simulate()
{
  std::cout<<"updatesimulation"<<std::endl;

  double dt = 0.001;

  for(int i = 0; i < 10; i++)
  {
    simulateExplicit(dt);
  }

  std::cout<<"f: "<<f[5]<<std::endl;
  std::cout<<"v: "<<v[5]<<std::endl;
  std::cout<<"a: "<<a[5]<<std::endl;
  r->update(x);
  
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

  r->setup();
  r->setCallback( simulate );
  r->update(x);
  r->run();  
  

  return 0;
}
