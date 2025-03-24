import sympy as sp
import re

x, y, z, u1, u2, u3, x1, x2, x3, x4 = sp.symbols('x y z u1 u2 u3 x1 x2 x3 x4')


class PseudoDivision:

    def __init__(self, var=sp.S(0)):
        variables: list
        if type(var) == str:
            variables = sorted(set(re.findall(r'[a-zA-Z]\d*', var)))
        else:
            variables = var.free_symbols
            variables = [str(i) for i in variables]
        self.vars: list = variables

    def add_var(self, *string):

        for str_ in string:
            if type(str_) == str:
                variables1 = sorted(set(re.findall(r'[a-zA-Z]\d*', str_)))
            else:
                variables1 = str_.free_symbols
                variables1 = [str(i) for i in variables1]
            variables_unit = variables1 + [var for var in variables1 if var not in variables1]
            self.vars.extend([var for var in variables_unit if var not in self.vars])

    def divide(self, g: sp.core, f: sp.core, variable=x, verbose=False):

        PseudoDivision.add_var(g)
        PseudoDivision.add_var(f)

        for var in self.vars:
            sp.symbols(var)
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
        dict1 = {"q": q, "r": r, "s": s, "bm": bm}

        return dict1

    def divide_str(self, g: str, f: str, variable="x", verbose=False):

        variable = sp.symbols(variable)
        PseudoDivision.add_var(self, g)
        PseudoDivision.add_var(self, f)

        symbol_dict = {var: sp.symbols(var) for var in self.vars}

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
        dict1 = {"q": q, "r": r, "s": s, "bm": bm}

        return dict1
