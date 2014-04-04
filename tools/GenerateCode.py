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


def secondDeriv( func, funcname, var, vectors, postfix =""):
  code=""
  syms = ['x','y','z']
  funcMat = [['','',''],['','',''],['','','']];
  for i in range(0,3):
    for j in range(0,3):
       name='needle_'+funcname+'_d'+syms[i]+'_d'+syms[j]+postfix
       print "Generating "+name
       deriv = diff(func, var[i], var[j])
       code += CodeToC.sympyToC( name, deriv, vectors )
       funcMat[i][j] = name

  code += CodeToC.sympyMatrixSetter( 'needle_'+funcname+postfix, funcMat, vectors ) 

  return code


code += secondDeriv( f2, "df2", b, ["a", "b"] )
code += secondDeriv( f, "df", b, ["a", "b", "c"] )
code += secondDeriv( f, "df", a, ["a", "b", "c"], "prev" )
code += secondDeriv( f, "df", c, ["a", "b", "c"], "next" )


f = open('gen_src/generatedCode.h', 'w')
f.write("#include \"mathheader.h\"\n\n")
f.write(code)

f.close()



