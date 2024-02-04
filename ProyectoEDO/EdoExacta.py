import sympy as sp

class EdoExacta:
    def __init__(self, M, N, x, y):
        self.M = M
        self.N = N
        self.x = x
        self.y = y

    def is_exact(self):
        # Compute the partial derivatives
        dM_dy = sp.diff(self.M, self.y)
        dN_dx = sp.diff(self.N, self.x)

        # Check if the partial derivatives are equal
        return dM_dy == dN_dx

    def solve(self):
        # Check if the equation is exact
        if not self.is_exact():
            return None

        # Integrate M with respect to x
        int_M = sp.integrate(self.M, self.x)

        # Integrate N with respect to y
        int_N = sp.integrate(self.N, self.y)

        # Check if the integrals are equal
        if int_M.subs(self.y, 0) != int_N.subs(self.x, 0):
            return None

        # Solve for the general solution
        x_sol = sp.integrate(int_M - self.N, self.x)
        y_sol = sp.integrate(self.M - int_M, self.y)

        return x_sol, y_sol