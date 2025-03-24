from pseudo_division import PseudoDivision as Pd
import sympy as sp

t = sp.symbols("t")
f5 = t**2 - 6*t + 5
f1 = "u1*x1 - u1*u3"
f2 = "u3*x2 - (u2-u1)*x1 + 7s"
f3 = "(u3*x2 - u2*x1 - u1*u3)*x3 + u1*u3*x1"
f4 = "u3*x4-u2*x3"

h1 = "u3*x2 - (u2-u1)*x1"
h2 = "u1*x1 - x2*u3"

g = "2*u2*x4 + 2*u3*x3 - u3*u3 - u2*u2"

D = Pd(g)
D.add_var(f1, f2, f5)
print(D.vars)
