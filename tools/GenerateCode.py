from sympy import *
from sympy.physics.mechanics import *
from sympy.printing import print_ccode
from sympy.utilities.codegen import codegen

from sympy.simplify.cse_main import *

import numpy
from sympy.utilities.iterables import * 

x1, x2, x3 = symbols('x1 x2 x3')
y1, y2, y3 = symbols('y1 y2 y3')
z1, z2, z3 = symbols('z1 z2 z3')

u = ReferenceFrame('u')

u1=(u.x*x1 + u.y*y1 + u.z*z1)
u2=(u.x*x2 + u.y*y2 + u.z*z2)
u3=(u.x*x3 + u.y*y3 + u.z*z3)

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

df_dx2=diff(f, x2)
df_dy2=diff(f, y2)
df_dz2=diff(f, z2)




x=0.0
y=0.5
z=0.0

print "df_dx:"
print(df_dx2)
print str(df_dx2)

def sympyToC( symname, symfunc ):
	tmpsyms = numbered_symbols("tmp")
	symbols, simple = cse(symfunc, symbols=tmpsyms)
	symbolslist = map(lambda x:str(x), list(symfunc.atoms(Symbol)) )
	symbolslist.sort()
	varstring=",".join( " double "+x for x in symbolslist )

	c_code = "double "+str(symname)+"("+varstring+" )\n"
	c_code +=  "{\n"
	for s in symbols:
		#print s
		c_code +=  "  double " +ccode(s[0]) + " = " + ccode(s[1]) + ";\n"
	c_code +=  "  double r = " + ccode(simple[0])+";\n"
	c_code +=  "  return r;\n"
	c_code += "}\n"
	return c_code

print " "
code = sympyToC( "df_dx", df_dx2 )
print code

f = open('gen_src/generatedCode.h', 'w')
f.write("#include <math.h>\n\n")
f.write(code)
f.close()



