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

k, k1 = symbols('k k1');

u = ReferenceFrame('u')

u1=(u.x*a[0] + u.y*a[1] + u.z*a[2])
u2=(u.x*b[0] + u.y*b[1] + u.z*b[2])
u3=(u.x*c[0] + u.y*c[1] + u.z*c[2])

m1=(u.y*1)
m2=(u.z*1)

tt= u3/2 + u1/2 - u2
w1=dot(tt,m1)
w2=dot(tt,m2)

#s1=(u1-u2).normalize()
#s2=(u2-u3).normalize()
#v=cross(s1, s2)
#f=(dot(v,v))

s1=(u1-u2)
s2=(u2-u3)

s1n=s1/s1.magnitude()
s2n=s2/s2.magnitude()

#k1 = 1000
k2 = 0

e0 = u2-u1
e1 = u3-u2

e0d = u.x
e1d = u.x


kbi = (2 * cross(e0, e1)) / (e0d.magnitude()*e1d.magnitude() + e0.dot(e1))

#alpha = acos(dot(s1n,s2n))
alpha =  1.0 - dot(s1n,s2n)
beta = (((u1+u3)*0.5) - u2).magnitude() ** 2
l1 = (s1.magnitude()-1.0)**2*k2
l2 = (s2.magnitude()-1.0)**2*k2 
ld = e0d.magnitude()+e1d.magnitude();
E = dot(kbi,kbi)*k1/ld + l1 + l2

Espring = (s1.magnitude())**2*k

code = ""

# 1st derivative of energy function
Fx=-diff(E, b[0])
Fy=-diff(E, b[1])
Fz=-diff(E, b[2])

Fx_next=-diff(E, c[0])
Fy_next=-diff(E, c[1])
Fz_next=-diff(E, c[2])

code += CodeToC.sympyToCMulti( [("Fx_next", Fx_next), ("Fy_next", Fy_next), ("Fz_next", Fz_next)], ["a", "b", "c"], ["k1"], prefix = "needle_" )

code += CodeToC.sympyToCMulti( [("Fx", Fx), ("Fy", Fy), ("Fz", Fz)], ["a", "b", "c"], ["k1"], prefix = "needle_" )



def secondDeriv( func, funcname, var1, var2, vectors, scalars, postfix =""):
  code=""
  syms = ['x','y','z']
  funcMat = [['','',''],['','',''],['','','']];
  for i in range(0,3):
    for j in range(0,3):
       name=funcname+'_d'+syms[i]+'_d'+syms[j]+postfix
       print "Generating "+name
       deriv = diff(func, var1[i], var2[j])
       code += CodeToC.sympyToC( name, deriv, vectors, scalars )
       funcMat[i][j] = name

  code += CodeToC.sympyMatrixAdder( funcname+postfix, funcMat, vectors, scalars ) 

  return code


code += secondDeriv( -E, "needle_Jacobian", b, b, ["a", "b", "c"], ["k1"], "BB")
code += secondDeriv( -E, "needle_Jacobian", b, a, ["a", "b", "c"], ["k1"], "BA" )
code += secondDeriv( -E, "needle_Jacobian", b, c, ["a", "b", "c"], ["k1"], "BC" )

code += secondDeriv( -E, "needle_Jacobian", a, b, ["a", "b", "c"], ["k1"], "AB")
code += secondDeriv( -E, "needle_Jacobian", a, a, ["a", "b", "c"], ["k1"], "AA" )
code += secondDeriv( -E, "needle_Jacobian", a, c, ["a", "b", "c"], ["k1"], "AC" )

code += secondDeriv( -E, "needle_Jacobian", c, b, ["a", "b", "c"], ["k1"], "CB")
code += secondDeriv( -E, "needle_Jacobian", c, a, ["a", "b", "c"], ["k1"], "CA" )
code += secondDeriv( -E, "needle_Jacobian", c, c, ["a", "b", "c"], ["k1"], "CC" )


FSx=-diff(Espring, a[0])
FSy=-diff(Espring, a[1])
FSz=-diff(Espring, a[2])

code += CodeToC.sympyToCMulti( [("FSx", FSx), ("FSy", FSy), ("FSz", FSz)], ["a", "b"], ["k"], prefix = "spring_" )

code += secondDeriv( -Espring, "spring_Jacobian", a, a, ["a", "b"], ["k"], "")

f = open('gen_src/generatedCode.h', 'w')
f.write("#include \"mathheader.h\"\n\n")
f.write(code)

f.close()



