import sympy as sp

# Define the variables

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