import sympy as sp
import os

# Define the variables

os.system('cls') if os.name == 'nt' else os.system('clear')

s, r, Lb, dt = sp.symbols('s, r, L, t')

x = sp.symbols('x')

ka = dt / (Lb * r)
k1 = ka * (2*Lb + dt) / (2 * (Lb - dt))

curv_x = (ka - k1) * x / (Lb / 2) + k1

ang_s = sp.integrate(curv_x, (x,0,s))

ang_simple = sp.simplify(ang_s)
print(sp.latex(ang_simple))


hor_tbi = sp.sin(ang_simple)
ver_tbi = sp.cos(ang_simple)

horpos = sp.integrate(hor_tbi,(s,0,Lb))

print(sp.latex(sp.simplify(horpos)))


exit()
x, k1, ka, r, Lb, dt = sp.symbols('x \\kappa_1 \\kappa_avg r L t')

kb = 2 * (ka-k1) * x / Lb + k1

kt = kb / (1 - r * kb)

kt_subs = sp.Subs(kt, (ka), (dt/(Lb*r)))
print("kt------")
print(sp.latex(sp.simplify(kt)))
# print(sp.latex(sp.simplify(kt_subs)))
# Perform the definite integration with symbolic limits
integral = sp.simplify(sp.integrate(sp.simplify(kt), (x, 0, Lb-dt)))

# Output the result in LaTeX format
int_equation = integral - dt/r

latex_output = sp.latex(sp.simplify(integral))

print("integral------")
print(latex_output)


print("int_equation = 0 --- == -")
print(sp.latex(sp.simplify(int_equation)))

print("int_equation, subs kavg as t/Lr")
print(sp.latex(sp.simplify(sp.Subs(int_equation, (ka), (dt/(Lb*r))))))

eqn_simplever = sp.exp(2*r*(k1-ka)) - 1 - \
    (2*r*(Lb-dt)*(k1-ka) / Lb*(1-r*k1))

print("eqn_simplever = 0 --- == -")
print(sp.latex(eqn_simplever))

eqn_solve = sp.solve(eqn_simplever, k1, dict = True)

print("eqn_solve")
sp.pprint(eqn_solve)
exit()