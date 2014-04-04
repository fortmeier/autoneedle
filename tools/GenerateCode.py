from sympy import *
from sympy.physics.mechanics import *
from sympy.printing import print_ccode
from sympy.utilities.codegen import codegen

from sympy.simplify.cse_main import *

import numpy
from sympy.utilities.iterables import *

import CodeToC

a = symbols('a[0] a[1] a[2]')
b = symbols('b[0] b[1] b[2]')
c = symbols('c[0] c[1] c[2]')

u = ReferenceFrame('u')

u1=(u.x*a[0] + u.y*a[1] + u.z*a[2])
u2=(u.x*b[0] + u.y*b[1] + u.z*b[2])
u3=(u.x*c[0] + u.y*c[1] + u.z*c[2])

#s1=(u1-u2).normalize()
#s2=(u2-u3).normalize()
#v=cross(s1, s2)
#f=(dot(v,v))

s1=(u1-u2)
s2=(u2-u3)

s1n=s1/s1.magnitude()
s2n=s2/s2.magnitude()

alpha = acos(dot(s1n,s2n))
alpha =  1.0 - dot(s1n,s2n)
beta = (((u1+u3)*0.5) - u2).magnitude() ** 2
l1 = (s1.magnitude()-1.0)**2
l2 = (s2.magnitude()-1.0)**2 
f = alpha + l1 + l2

f2 = l1

code = ""

# 1st derivative of energy function
df2_dx=diff(f2, b[0])
df2_dy=diff(f2, b[1])
df2_dz=diff(f2, b[2])

df_dx=diff(f, b[0])
df_dy=diff(f, b[1])
df_dz=diff(f, b[2])

df_dxprev=diff(f, a[0])
df_dyprev=diff(f, a[1])
df_dzprev=diff(f, a[2])

df_dxnext=diff(f, c[0])
df_dynext=diff(f, c[1])
df_dznext=diff(f, c[2])

code += CodeToC.sympyToCMulti( [("df_dx", df_dx), ("df_dy", df_dy), ("df_dz", df_dz)], ["a", "b", "c"], prefix = "needle_" )

code += CodeToC.sympyToCMulti( [("df_dxprev", df_dxprev), ("df_dyprev", df_dyprev), ("df_dzprev", df_dzprev)], ["a", "b", "c"], prefix = "needle_" )
code += CodeToC.sympyToCMulti( [("df_dxnext", df_dxnext), ("df_dynext", df_dynext), ("df_dznext", df_dznext)], ["a", "b", "c"], prefix = "needle_" )

code += CodeToC.sympyToCMulti( [("df2_dx", df2_dx), ("df2_dy", df2_dy), ("df2_dz", df2_dz)], ["a", "b"], prefix = "needle_" )


# 2nd derivative of energy function
df2_dx_dx=diff(f2, b[0], b[0])
df2_dy_dy=diff(f2, b[1], b[1])
df2_dz_dz=diff(f2, b[2], b[2])

df_dx_dx=diff(f, b[0], b[0])
df_dy_dy=diff(f, b[1], b[1])
df_dz_dz=diff(f, b[2], b[2])

df_dx_dxprev=diff(f, a[0], a[0])
df_dy_dyprev=diff(f, a[1], a[1])
df_dz_dzprev=diff(f, a[2], a[2])

df_dx_dxnext=diff(f, c[0], c[0])
df_dy_dynext=diff(f, c[1], c[1])
df_dz_dznext=diff(f, c[2], c[2])

#df_dxx=diff(f, b[0], b[0])
#df_dyy=diff(f, b[1], b[1])
#df_dzz=diff(f, b[2], b[2])

#df_dxxx=diff(f, b[0], b[0], b[0])
#df_dyyy=diff(f, b[1], b[1], b[1])
#df_dzzz=diff(f, b[2], b[2], b[2])


code += CodeToC.sympyToCMulti( [("df_dx_dx", df_dx_dx), ("df_dy_dy", df_dy_dy), ("df_dz_dz", df_dz_dz)], ["a", "b", "c"], prefix = "needle_" )

code += CodeToC.sympyToCMulti( [("df_dx_dxprev", df_dx_dxprev), ("df_dy_dyprev", df_dy_dyprev), ("df_dz_dzprev", df_dz_dzprev)], ["a", "b", "c"], prefix = "needle_" )
code += CodeToC.sympyToCMulti( [("df_dx_dxnext", df_dx_dxnext), ("df_dy_dynext", df_dy_dynext), ("df_dz_dznext", df_dz_dznext)], ["a", "b", "c"], prefix = "needle_" )

code += CodeToC.sympyToCMulti( [("df2_dx_dx", df2_dx_dx), ("df2_dy_dy", df2_dy_dy), ("df2_dz_dz", df2_dz_dz)], ["a", "b"], prefix = "needle_" )
#code += CodeToC.sympyToCMulti( [("df_dxx", df_dxx), ("df_dyy", df_dyy), ("df_dzz", df_dzz)], ["a", "b", "c"], prefix = "needle_" )

#code += CodeToC.sympyToCMulti( [("df_dxxx", df_dxxx), ("df_dyyy", df_dyyy), ("df_dzzz", df_dzzz)], ["a", "b", "c"], prefix = "needle_" )

f = open('gen_src/generatedCode.h', 'w')
f.write("#include \"mathheader.h\"\n\n")
f.write(code)

f.close()



