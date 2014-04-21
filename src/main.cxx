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

#include "needle.h"
#include "Rendering.h"


#define numNodes 10


//BendingNeedleModel needle( 9.0, 10, 10.0 );
//BendingNeedleModel needle( 9.0, 10, 10.0 );
BendingNeedleModel needle( 150.0, 31, 10000.0 );



Rendering* r;



bool ex = false;

double yoffset = 0.05;
double xrot = 0;
double error = 0;

void simulate()
{

  double dt = 0.01;

  static int c = 0;

  if( error > 100 || error != error ) return;

  for(int i = 0; i < 1 && c < 30000; i++)
  {
    c++;
    //simulateExplicit(dt);
    //simulateImplicit(dt);
    error = needle.simulateImplicitDynamic(dt);
    std::cout<<c<<" : "<<yoffset<<" : "<<error<<std::endl;
    switch(0) 
    {
      case 0:
        yoffset -= c < 100 ? 1.1 : 0;
        needle.setBasePosition( Vector(0,yoffset,0));
      break;
      case 10:
        xrot += c < 100 ? 0.001 : -0.001;
        needle.setBaseDirection( Vector(1,xrot,0));
      break;
      case 20:
        if(c%99 == 0 && c < 1600)
        {
          yoffset += yoffset;
          needle.setBasePosition( Vector(0,yoffset,0));
        }
      break;

    }

  // for(int i = 0; i < numNodes; i++)
  // {
  // if(debugOut) {
  //   std::cout<<i<<":"<<std::endl;
  //   std::cout<<"x: "<<x[i]<<std::endl;
  //   std::cout<<"f: "<<f[i]<<std::endl;
  //   std::cout<<"v: "<<v[i]<<std::endl;
  //   std::cout<<"a: "<<a[i]<<std::endl;
  // }
  // }

  }

  std::cout<<"needle length: "<<needle.getTotalLength()<<std::endl;

  r->update(needle.getX());
  if(ex) exit(0);
 
}

int main(int argi, char** argv)
{

  if(argi>1) ex = true;

  r = new Rendering();

  //needle.addLagrangeModifier(0, Vector(0,1,0));
  //needle.addLagrangeModifier(1, Vector(0,1,0));
  int last = needle.getX().size()-1;
  double s = needle.getSegmentLength();
  //needle.setSpring( last, Vector( last * s, 0, 0 ), 0.1 );
  //needle.setSpring( last-1, Vector( (last-1) * s, 0, 0 ), 0.1 );
  //needle.setSpring( last-2, Vector( (last-2) * s, 0, 0 ), 0.1 );

  needle.setSpring( 0, Vector( 0, 0, 0 ), 200.002 );
  needle.setSpring( 1, Vector( s, 0, 0 ), 200 );

  r->setup();
  r->setCallback( simulate );
  r->update(needle.getX());
  r->run();  
  

  return 0;
}
