import sympy as sp
import re


x, y, z, u1, u2, u3, x1, x2, x3, x4 = sp.symbols('x y z u1 u2 u3 x1 x2 x3 x4')
# g_ = x*y**3 - x**2*y +1
# f_ = (x+1)*y -x
# g_ = 2*x*y**2 + x*y - x**2*y + y**3 -2 + 10*x**2 + 3*x**3*y
# f_ = x*y**2
# g_ = 2*x**3*y**2 - x**2*y + 10*x - 2 + 6*x**2*y**2 + 6*y + 4*x**4
# f_ = 5*x**2*y - 2*x
# g_str = "2*x**3*y**2 - x**2*y + 10*x - 2 + 6*x**2*y**2 + 6*y + 4*x**4"
# f_str = "5*x**2*y - 2*x"
# g_ = x**2*y*z + 3*x**3*y + y**2*z**2 + 3*x**4*z + 7*x**2*z
# f_ = 2*x**2*z**2 +y

f1 = u1*x1 - u1*u3
f2 = u3*x2 - (u2-u1)*x1
f3 = (u3*x2 - u2*x1 - u1*u3)*x3 + u1*u3*x1
f4 = u3*x4-u2*x3

h1 = u3*x2 - (u2-u1)*x1
h2 = u1*x1 - x2*u3

g = 2*u2*x4 + 2*u3*x3 - u3*u3 - u2*u2
print(type(g))


def pseudo_division(g: sp.core, f: sp.core, variable=x, verbose=False):

    f_v = sp.Poly(f, variable)
    g_v = sp.Poly(g, variable)

    bm = f_v.LC()
    m = f_v.degree()

    r = g_v.copy()
    q = sp.S(0)
    r_tem = r.copy()
    q_tem = q.copy()

    t = 0
    while not (r == 0) and r.degree() >= m:
        r = bm * r_tem - r_tem.LC() * f * variable ** (r_tem.degree() - m)
        q = bm * q_tem + r_tem.LC() * variable ** (r_tem.degree() - m)
        q = sp.Poly(q, variable)
        r = sp.Poly(r, variable)
        r_tem = r.copy()
        q_tem = q.copy()
        t += 1

    product = g_v.copy()
    combination = q * f_v + r

    s = 0
    while not (product == combination):
        product = product * bm
        s += 1

    if verbose:
        print("successful division" if q * f_v + r == bm ** s * g_v else "error occurred")
        print(f"Iteration done: {t}")
        print(f"q*f + r = ({bm})^{s}*g" if s != 1 else f"q*f + r = ({bm})*g")
    dict1 = {"q": q, "r": r, "s": s}

    return dict1


def pseudo_division_str(g: str, f: str, variable="x", verbose=False):

    variable = sp.symbols(variable)
    variables = sorted(set(re.findall(r'[a-zA-Z]+', g)))
    variables2 = sorted(set(re.findall(r'[a-zA-Z]+', f)))
    variables_unit = variables + [var for var in variables2 if var not in variables]

    symbol_dict = {var: sp.symbols(var) for var in variables_unit}

    f = sp.sympify(f, locals=symbol_dict)
    g = sp.sympify(g, locals=symbol_dict)

    f_v = sp.Poly(f, variable)
    g_v = sp.Poly(g, variable)

    bm = f_v.LC()
    m = f_v.degree()

    r = g_v.copy()
    q = sp.S(0)
    r_tem = r.copy()
    q_tem = q.copy()

    t = 0
    while not (r == 0) and r.degree() >= m:
        r = bm * r_tem - r_tem.LC() * f * variable ** (r_tem.degree() - m)
        q = bm * q_tem + r_tem.LC() * variable ** (r_tem.degree() - m)
        q = sp.Poly(q, variable)
        r = sp.Poly(r, variable)
        r_tem = r.copy()
        q_tem = q.copy()
        t += 1

    product = g_v.copy()
    combination = q * f_v + r

    s = 0
    while not (product == combination):
        product = product * bm
        s += 1

    if verbose:
        print("successful division" if q * f_v + r == bm ** s * g_v else "error occurred")
        print(f"Iteration done: {t}")
        print(f"q*f + r = ({bm})^{s}*g" if s != 1 else f"q*f + r = ({bm})*g")
    dict1 = {"q": q, "r": r, "s": s}

    return dict1


# print(pseudo_division(g_, f_, y))
# print(pseudo_division_str(g_str, f_str, "y"))
print(g)
print(pseudo_division(g, f4, x4))
R3 = pseudo_division(g, f4, x4)["r"]

print(pseudo_division(R3,f3,x3))
R2 = pseudo_division(R3,f3,x3)["r"]
print(pseudo_division(R2,f2,x2))
R1 = pseudo_division(R2,f2,x2)["r"]
R0 = pseudo_division(R1,f1,x1)["r"]
print(R0)

print(pseudo_division(h1, h2, x1))
