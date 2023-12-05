import numpy as np
from sympy import symbols, Eq, solve

A1 = np.array([2586564.018, 274040.686])
B1 = np.array([2561432.582, 271464.628])
C1 = np.array([2672140.277, 278091.433])

A2 = np.array([2586357.074, 274869.216])
B2 = np.array([2561226.156, 272292.830])
C2 = np.array([2671934.535, 278920.676])

A, B, C, D = symbols('A B C D')

eq_AB_1 = Eq((A1[0]+B1[0])*A - (A1[1]+B1[1])*B - (A2[0]+B2[0]))
eq_AC_1 = Eq((A1[0]+C1[0])*A - (A1[1]+C1[1])*B - (A2[0]+C2[0]))
eq_BC_1 = Eq((B1[0]+C1[0])*A - (B1[1]+C1[1])*B - (B2[0]+C2[0]))

sol_AB = solve((eq_AB_1, eq_AC_1), (A, B))
print(sol_AB)

sol_AC = solve((eq_AB_1, eq_BC_1), (A, B))
print(sol_AC)

sol_BC = solve((eq_AC_1, eq_BC_1), (A, B))
print(sol_BC)
