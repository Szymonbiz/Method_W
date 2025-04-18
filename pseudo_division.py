import sympy as sp
import re
from polynomial import Poly
from chain import Chain


class PseudoDivision:

    """
    The primary tool in Wu’s Method is a variation on the division algorithm for multivariable polynomials called pseudo
    division. Let f and g belong to ring k[x1,...,xn,y], where k is a field, with
    g=a_p*y^p+...+a_0 and f=b_m*y^m+...+b_0, where the a_i, b_j are polynomials in the x1, . . . , xn. Then we have the
    following result.

        Let f, g be as above and assume that m ≤ p and that f != 0. Then,

            (i) There is an equation
                                                (b_m)^s*g=q*f+r

        where q, r ∈ k[x1, . . . , xn, y], s ≥ 0, and r is either the zero polynomial or its degree in y is less than m.

            (ii) r is in the ideal〈f, g〉in the ring k[x1, . . . , xn, y].
    """

    def __init__(self, *var):
        variables_unit = []
        for str_ in var:
            if type(str_) == str:
                variables1 = sorted(set(re.findall(r'[a-zA-Z]\d*', str_)))
            else:
                variables1 = str_.free_symbols
                variables1 = [str(i) for i in variables1]
            variables_unit = variables_unit + [var for var in variables1 if var not in variables_unit]

        self.vars: list = variables_unit

    def add_var(self, *string):

        for str_ in string:
            if type(str_) == str:
                variables1 = sorted(set(re.findall(r'[a-zA-Z]\d*', str_)))
            else:
                variables1 = str_.free_symbols
                variables1 = [str(i) for i in variables1]
            variables_unit = variables1 + [var for var in variables1 if var not in variables1]
            self.vars.extend([var for var in variables_unit if var not in self.vars])

    @staticmethod
    def find_index(F):
        max_ = 0
        for xi in F:
            if xi[0] == "x" and xi[-1] != "x":
                max_ = int(xi[1]) if int(xi[1]) > max_ else max_
        return max_

    """
    input: f, g
    Output: r, q
    
    r :=g, q := 0
    While r != 0 and deg(r,y) ≥ m:
        r := b_m*r - LC(r,y)*f*y^(deg(r,y) -m)
        q := b_m*q - LC(r,y)*y^(deg(r,y) -m)
        
    where notations deg(f, y) and LC(f, y) to denote the degree of f in the variable y and the leading coefficient of f 
    as a polynomial in y.
    """

    def divide_str(self, g: str, f: str, variable="x", verbose=False) -> dict:

        # Creating ring with given variables.
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

        # Start of the algorithm.
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

        # Finding power s.
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

    def prem(self, g: Poly, S: Chain):
        matrix = [f.vars for f in S.chain]
        domain = [var for vector in matrix for var in vector]
        self.add_var(*domain)

        i = PseudoDivision.find_index(self.vars)
        R = g.expr
        while i >= 0:
            R = PseudoDivision.divide_str(self, R, S.chain[i-1].expr, f"x{i + 1}", verbose=False)["r"]
            i -= 1
        return R

    def divide(self, g: sp.core, f: sp.core, variable, verbose=False) -> dict:

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
