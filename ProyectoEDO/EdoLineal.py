import sympy as sp

class EdoLineal:
    def __init__(self, M, N, x, y):
        self.M = M
        self.N = N
        self.x = x
        self.y = y

    def is_linear(self, y):
        # Check if the EDO is linear
        return self.M.is_linear(y) and self.N.is_linear(y)

    def convert_to_linear_form(self, y):
        # Convert the EDO to linear form
        if self.is_linear(y):
            return self.M, self.N

        # Compute the coefficients of the linear form
        a = sp.diff(self.M, y)
        b = sp.diff(self.N, y)
        c = self.M.coeff(y.diff())
        d = self.N

        # Convert the EDO to linear form
        M_linear = (c*d - b*a) / c
        N_linear = d - M_linear * a

        return M_linear, N_linear

    def solve(self):
        # Convert the EDO to linear form
        M_linear, N_linear = self.convert_to_linear_form(self.y)

        # Integrate the linear EDO
        integral_N = sp.integrate(N_linear, self.x)
        integral_M = sp.integrate(M_linear, self.x)

        # Solve for y
        y_sol = sp.exp(-integral_M) * (integral_N + sp.Symbol('C'))

        return y_sol