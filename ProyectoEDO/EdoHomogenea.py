import sympy as sp

class EdoHomogenea:
    def __init__(self, M, N, x, y):
        self.M = M
        self.N = N
        self.x = x
        self.y = y

    def is_homogeneous(self):
        # Check if M and N have the same degree
        M_args = sp.degree(self.M, (self.x, self.y))
        N_args = sp.degree(self.N, (self.x, self.y))

        if M_args == N_args:
            return True
        else:
            return False

    def degree(self):
        # Calculate the degree of the ODE
        if not self.is_homogeneous():
            return None

        M_args = sp.degree(self.M, (self.x, self.y))
        N_args = sp.degree(self.N, (self.x, self.y))

        return M_args

    def change_to_homogeneous_form(self):
        # Change the ODE to homogeneous form
        M = self.M
        N = self.N

        if not self.is_homogeneous():
            return None

        x = self.x
        y = self.y

        homogeneous_M = M / (x**sp.degree(M, x) * y**sp.degree(M, y))
        homogeneous_N = N / (x**sp.degree(N, x) * y**sp.degree(N, y))

        return homogeneous_M, homogeneous_N
    
    def solve(self):
        # Change the ODE to homogeneous form
        homogeneous_M, homogeneous_N = self.change_to_homogeneous_form()

        # Interchange variables
        if sp.degree(homogeneous_M, self.x) > sp.degree(homogeneous_M, self.y):
            x, y = self.y, self.x
            M, N = homogeneous_N, homogeneous_M
        else:
            x, y = self.x, self.y
            M, N = homogeneous_M, homogeneous_N

        # Change variables to U and dU
        U = y / x
        dU = (y.diff() - U * x.diff()) / x

        # Change the ODE to an ODE in U and dU
        U_diff_eq = sp.Eq(M * dU + U * M.diff() - N * dU - U * N.diff(), 0)

        # Separate the ODE
        U_sol = sp.dsolve(U_diff_eq, U)

        # Convert U back to x or y
        if y == self.y:
            x_sol = sp.solve(U_sol.rhs - U, x)[0]
            y_sol = sp.solve(U_sol, y)[0]
        else:
            x_sol = sp.solve(U_sol.rhs - U, y)[0]
            y_sol = sp.solve(U_sol, x)[0]

        return x_sol, y_sol