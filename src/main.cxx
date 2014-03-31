#include <iostream>
#include "generatedCode.h"
#include "mathheader.h"

int main()
{
  Vector a(-1,0,0);
  Vector b(0,1,0);
  Vector c(1,0,0);

  double x = needle_jacobian_df_dx(a,b,c);
  double y = needle_jacobian_df_dy(a,b,c);
  double z = needle_jacobian_df_dz(a,b,c);

  Vector J(x,y,z);

  std::cout<<"jacobian: "<<J<<std::endl;
  return 0;
}
