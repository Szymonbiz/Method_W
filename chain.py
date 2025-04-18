import sympy as sp
from polynomial import Poly


class Chain:

    def __init__(self, *C):

        if not all(isinstance(f, Poly) for f in C):
            raise TypeError(f"Expected elements of type 'Poly'. Got: {[type(e).__name__ for e in C]}")

        self._chain = list(C)

    def __str__(self):
        inner = "\033[1m,\033[0m".join([f"\x1B[3m {str(f)} \x1B[0m" for f in self.chain])
        return "\033[1m(\033[0m{}\033[1m)\033[0m".format(inner)

    @property
    def chain(self):
        return self._chain

    @property
    def len(self):
        return len(self._chain)

    def is_seq_quasi(self):
        return (len(self.chain) == 1 and self.chain[0] != sp.S(0)) or (self.len > 1 and all(
            [0 < self.chain[i].Class < self.chain[i+1].Class for i in range(self.len - 1)]))

    def is_seq_ascending(self):
        return self.is_seq_quasi and \
            all([self.chain[j].is_reduced_to(self.chain[i]) for i in range(self.len) for j in range(i + 1, self.len)])

    def __lt__(self, other):

        """
        Definition: Let C = f1, . . . , fr and C1 = g1, . . . , gm be ascending chains. We
        say that C < C1 if either,
            (i) ∃s ≤ min(r, m) such that fi, gi are of the same rank for i < s and fs < gs, or
            (ii) m < r and fi and gi are of the same rank for i ≤ m.
        """

        # Requirements
        if not isinstance(other, Chain):
            raise TypeError(f"Cannot compare object of type {type(other).__name__}. Expected type: 'Chain'.")
        if not (self.is_seq_ascending and other.is_seq_ascending):
            raise TypeError("Expected ascending chains of polynomials.")

        # Case (i)
        if any([self.chain[s] < other.chain[s] and
                all([self.chain[i] == other.chain[i] for i in range(s - 1)])
                for s in range(min(self.len, other.len))]):
            bool2 = True
        else:
            bool2 = False

        # Case (ii)
        if other.len < self.len and all([self.chain[i] == other.chain[i] for i in range(other.len)]):
            bool1 = True
        else:
            bool1 = False

        return bool1 or bool2


# S = Chain(*[Poly("x1**2"), Poly("y**3")])
# print(S.chain)
# print(S.len)
# print(S)
# print("\n\n")
# print("[x1**2, y**3]")
