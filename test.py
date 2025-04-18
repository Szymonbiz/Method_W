from pseudo_division import PseudoDivision
from polynomial import Poly
from chain import Chain

f1 = "u1*x1 - u1*u3"
f2 = "u3*x2 - (u2-u1)*x1"
f3 = "(u3*x2 - u2*x1 - u1*u3)*x3 + u1*u3*x1"
f4 = "u3*x4-u2*x3"

g = "2*u2*x4 + 2*u3*x3 - u3*u3 - u2*u2"

F = [f1, f2, f3, f4]
G = [g]
pd = PseudoDivision(*F)
pd1 = PseudoDivision()

for gi in G:
    R = gi
    for i in range(pd.find_index(pd.vars) - 1, -1, -1):
        print(pd.find_index(pd.vars))
        print(i)
        R = pd.divide_str(R, F[i], f"x{i + 1}", True)["r"]
    print(R)

F_p = [Poly(i) for i in F]
S = Chain(*F_p)
print(S)
print(Poly(g).vars)
print(pd1.prem(Poly(g), S))
