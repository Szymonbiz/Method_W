import sympy as sp


x, y, z = sp.symbols('x y z')

# g = x*y**3 - x**2*y +1
# f = (x+1)*y -x
# g = 2*x*y**2 + x*y - x**2*y + y**3 -2 + 10*x**2 + 3*x**3*y
# f = x*y**2
# g = 2*x**3*y**2 - x**2*y + 10*x -2 + 6*x**2*y**2 + 6*y + 4*x**4
# f = 5*x**2*y - 2*x
g = x**2*y*z + 3*x**3*y + y**2*z**2 + 3*x**4*z + 7*x**2*z
f = 2*x**2*z**2 +y

variable = x

f_y = sp.Poly(f, variable)
g_y = sp.Poly(g, variable)


bm = f_y.LC()
m = f_y.degree()


r_in = g_y.copy()
q_in = sp.S(0)
r_tem = r_in.copy()
q_tem = q_in.copy()


t = 0
while not( r_in == 0) and r_in.degree() >= m:
    r_in = bm*r_tem- r_tem.LC()*f*variable**(r_tem.degree()-m)
    q_in = bm*q_tem + r_tem.LC()*variable**(r_tem.degree()-m)
    r_tem = r_in.copy()
    q_tem = q_in.copy()
    t += 1

print(f"Iteration done: {t}")
print(f"q: {q_in.as_expr()}")
print(f"r: {r_in.as_expr()}")

product = g_y.copy()
combination = q_in * f_y + r_in

s = 0
while not(product == combination):
    product = product*bm
    s += 1

print(q_in*f_y+r_in == bm**s*g_y)
print(f"q*f + r = ({bm})^{s}*g" if s != 1 else f"q*f + r = ({bm})*g")