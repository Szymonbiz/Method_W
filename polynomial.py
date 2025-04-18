import sympy as sp
import re


class Poly:

    def __init__(self, expression: str):

        self.vars: list = sorted(set(re.findall(r'[a-zA-Z]\d*', expression)))
        self.expr: str = expression
        self.Class = Poly.get_Class(self)
        self.LV = Poly.get_LV(self)
        self.LCx = Poly.get_LCx(self)
        self.LD = Poly.get_LD(self)

    def get_Class(self) -> int:

        """
        Definition: We call that the class of a polynomial f , denoted Class(f),
        the smallest integer c such that f ∈ k[x1, . . . , xc].
        """
        x_vars = re.findall(r'x(\d+)', self.expr)  # -> ["1", "2", ...]

        if not x_vars:
            return 0

        indices = [int(i) for i in x_vars]
        return max(indices)

    def get_LV(self):
        return f'x{Poly.get_Class(self)}'

    def get_LCx(self):

        """
        Definition: We call leading coefficient of f, denoted LV(f), as a polynomial in xc, where c = Class(f).
        """
        var = sp.Symbol(Poly.get_LV(self))
        f_poly = sp.sympify(self.expr, locals={Poly.get_LV(self): var})
        f_v = sp.Poly(f_poly, var)
        return f_v.LC()

    def get_deg_in_v(self, variable: str) -> int:
        var = sp.Symbol(variable)
        f_poly = sp.sympify(self.expr, locals={variable: var})
        f_v = sp.Poly(f_poly, var)
        return f_v.degree()

    def get_LD(self) -> int:

        """
        Definition: We call leading degree of f, denoted LD(f), as a polynomial in xc, where c = Class(f).
        """
        f_v = sp.Poly(self, Poly.get_LV(self))
        return f_v.degree()

    def is_reduced_to(self, other) -> bool:

        """
        Definition: A polynomial g is reduced with respect to f if deg(g, xc) < deg(f, xc), where Class(f) = c > 0.
        """
        if other.Class == 0:
            raise ValueError("Expected a polynomial with class > 0 for comparison.")
        if not isinstance(other, Poly):
            raise TypeError(f"Cannot compare object of type {type(other).__name__}. Expected type: 'Poly'.")
        return self.get_deg_in_v(other.LV) < other.get_deg_in_v(other.LV)

    def __lt__(self, other):
        """
        Definition: Given f, g ∈ k[x1, . . . , xn] we say that f < g (g is higher, or of
        higher rank) if either
            (i) class(f ) < class(g), or
            (ii) class(f ) = class(g) and ld(f ) < ld(g).
        """
        if not isinstance(other, Poly):
            raise TypeError(f"Cannot compare object of type {type(other).__name__}. Expected type: 'Poly'.")
        return (self.Class < other.Class) or \
            (self.Class == other.Class and self.LD < other.LD)

    def __eq__(self, other):
        """
        :param self:
        :param f, other:
        :return: f = other
        """
        if not isinstance(other, Poly):
            raise TypeError(f"Cannot compare object of type {type(other).__name__}. Expected type: 'Poly'.")
        return self.Class == other.Class and self.LD == other.LD

    def __gt__(self, other):
        return other < self

    def __le__(self, other):
        return self < other or self == other

    def __ge__(self, other):
        return self > other or self == other

    def __ne__(self, other):
        return not self == other
