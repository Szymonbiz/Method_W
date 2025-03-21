import sympy as sp


x, y, z = sp.symbols('x y z')

# g_ = x*y**3 - x**2*y +1
# f_ = (x+1)*y -x
# g_ = 2*x*y**2 + x*y - x**2*y + y**3 -2 + 10*x**2 + 3*x**3*y
# f_ = x*y**2
g_ = 2*x**3*y**2 - x**2*y + 10*x - 2 + 6*x**2*y**2 + 6*y + 4*x**4
f_ = 5*x**2*y - 2*x
# g_ = x**2*y*z + 3*x**3*y + y**2*z**2 + 3*x**4*z + 7*x**2*z
# f_ = 2*x**2*z**2 +y


def pseudo_division(g: sp.core, f: sp.core, variable=x):

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

    print(f"Iteration done: {t}")
    # print(f"q: {q.as_expr()}")
    # print(f"r: {r.as_expr()}")

    product = g_v.copy()
    combination = q * f_v + r

    s = 0
    while not (product == combination):
        product = product * bm
        s += 1

    # print(q * f_v + r == bm ** s * g_v)
    print(f"q*f + r = ({bm})^{s}*g" if s != 1 else f"q*f + r = ({bm})*g")
    dict1 = {"q": q, "r": r, "s": s}

    return dict1


print(pseudo_division(g_, f_, y))
