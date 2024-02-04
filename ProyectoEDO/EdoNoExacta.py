import sympy as sp

class EdoNoExacta:
    def __init__(self, M, N, x, y):
        self.M = M
        self.N = N
        self.x = x
        self.y = y

    def is_not_exact(self, tol=1e-5):
        # Compute the partial derivatives
        dM_dy = sp.diff(self.M, self.y)
        dN_dx = sp.diff(self.N, self.x)

        # Check if the partial derivatives are not equal up to a constant factor
        return sp.Abs(dM_dy - dN_dx) > tol

    def solve(self, method="integrating_factor"):
        # Check if the equation is not exact
        if not self.is_not_exact():
            return None

        factor = self.integrating_factor()

        M_factor = factor * self.M
        N_factor = factor * self.N

        # Integrate M with respect to x
        int_M = sp.integrate(M_factor, self.x)

        # Integrate N with respect to y
        int_N = sp.integrate(N_factor, self.y)

        # Solve for the general solution
        x_sol = sp.integrate(int_N - int_M, self.x)
        y_sol = sp.integrate(factor * (self.M - int_M), self.y)

        return x_sol, y_sol

    def integrating_factor(self, method="multiplicative"):
        # Compute the partial derivatives
        dM_dy = sp.diff(self.M, self.y)
        dN_dx = sp.diff(self.N, self.x)

        # Compute the multiplicative integrating factor
        if method == "multiplicative":
            return sp.exp(sp.integrate(dM_dy - dN_dx, self.x))

        # Compute the additive integrating factor
        if method == "additive":
            return sp.integrate(dM_dy - dN_dx, self.x)

        raise ValueError("Invalid method for integrating factor")