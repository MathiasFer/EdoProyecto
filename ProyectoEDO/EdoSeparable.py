import sympy as sp

class EdoSeparable:
    def __init__(self, M, N):
        self.M = M
        self.N = N

    def is_separable(self):
        # Comprobar si la ecuación ya está separada
        if self.M.free_symbols.intersection(self.N.free_symbols) == set():
            return True

        # Despejar y factorizar para separarla
        try:
            y = sp.Symbol('y')
            x = sp.Symbol('x')
            if self.M.diff(x) == self.N.diff(y):
                # Despejar y
                y_solved = sp.solve(self.M - self.N, y)
                if len(y_solved) == 1:
                    y_solved = y_solved[0]
                    self.M = self.M.subs(y, y_solved)
                    self.N = self.N.subs(y, y_solved)
                    if self.M.free_symbols.intersection(self.N.free_symbols) == set():
                        return True
                else:
                    raise Exception("Error while solving for y")
            else:
                raise Exception("The equation is not separable")
        except Exception as e:
            print(f"Error while checking separability: {e}")
            return False

    def solve(self, method="integral"):
        y = sp.Symbol('y')
        x = sp.Symbol('x')
        if not self.is_separable():
            print("The equation is not separable")
            return None

        # Integrar
        if method == "integral":
            integral_M = sp.integrate(self.M, 'x')
            integral_N = sp.integrate(self.N, 'y')

            # Despejar y
            y = sp.solve(integral_M - integral_N, 'y')
            if len(y) == 1:
                y = y[0]
                return y
            else:
                return None
        elif method == "dsolve":
            y = sp.Function('y')
            eq = sp.Eq(self.M*sp.Derivative(y(x), x), self.N*sp.Derivative(y(x), x))
            try:
                y_sol = sp.dsolve(eq, y(x))
                return y_sol
            except sp.exc.SympifyError as e:
                print(f"Error while solving using dsolve: {e}")
                return None
        else:
            print(f"Method {method} is not supported")
            return None

        
    