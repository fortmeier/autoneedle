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
beta = (((u1+u3)*0.5) - u2).magnitude() ** 2
l1 = (s1.magnitude()-1.0)**2 
l2 = (s2.magnitude()-1.0)**2 
f = alpha

df_dx=diff(f, b[0])
df_dy=diff(f, b[1])
df_dz=diff(f, b[2])

code = CodeToC.sympyToCMulti( [("df_dx", df_dx), ("df_dy", df_dy), ("df_dz", df_dz)], ["a", "b", "c"], prefix = "needle_jacobian_" )

print code

f = open('gen_src/generatedCode.h', 'w')
f.write("#include \"mathheader.h\"\n\n")
f.write(code)

f.close()



