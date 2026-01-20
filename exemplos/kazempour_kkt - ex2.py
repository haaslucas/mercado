# exemplo retirado de Lecture 1 = 1h30m
# https://www.youtube.com/watch?v=1h3b0j2k4a8&list=PL6X9d7c5e2f8f7c4b&index=1

from pyomo.environ import *

m = ConcreteModel()

# Variáveis primais
m.x1, m.x2, m.x3, m.x4 = Var(), Var(), Var(), Var()

# Variáveis duais
m.mu1, m.mu2, m.mu3, m.mu4, m.mu5, m.mu6 = Var(), Var(), Var(), Var(), Var(), Var()

# Função objetivo
m.obj = Objective(expr = 18*m.x1 + 8*m.x2 + 12*m.x3 + 16*m.x4, sense=minimize)

# Restrições primais
m.restr1 = Constraint(expr = 2/3*m.x1 + 2*m.x2 + m.x3 >= 1) # mu1
m.restr2 = Constraint(expr = m.x1 + m.x2 + m.x4 >= 1) # mu2
m.restr3 = Constraint(expr = m.x1 >= 0) # mu3
m.restr4 = Constraint(expr = -m.x2 >= 0) # mu4
m.restr5 = Constraint(expr = m.x3 >= 0)# mu5
m.restr6 = Constraint(expr = m.x4 >= 0) # mu6

# Restrições duais: 2a condição KKT -> variáveis duais associadas às inequações primais são não-negativas
m.dual_restr1 = Constraint(expr = m.mu1 >= 0)
m.dual_restr2 = Constraint(expr = m.mu2 >= 0) 
m.dual_restr3 = Constraint(expr = m.mu3 >= 0)
m.dual_restr4 = Constraint(expr = m.mu4 >= 0)
m.dual_restr5 = Constraint(expr = m.mu5 >= 0)
m.dual_restr6 = Constraint(expr = m.mu6 >= 0)

# LaGrangiano
'''
Lagr(x,mu) = f(x) + mu*g(x)

Lagr(x,mu) = 18*m.x1 + 8*m.x2 + 12*m.x3 + 16*m.x4
           - m.mu1*(2/3*m.x1 + 2*m.x2 + m.x3 - 1)
           - m.mu2*(m.x1 + m.x2 + m.x4 - 1)
           - m.mu3*m.x1
           - m.mu4*(-m.x2)
           - m.mu5*m.mu3
           = m.mu6*m.x4
'''

# 1a condição KKT: gradiente da Lagrangiana em relação a x == 0
m.kkt1_x1 = Constraint(expr = 18 - 2/3*m.mu1 - m.mu2 - m.mu3 == 0) # x1
m.kkt1_x2 = Constraint(expr = 8 - 2*m.mu1 - m.mu2 + m.mu4 == 0) # x2
m.kkt1_x3 = Constraint(expr = 12 - m.mu1 - m.mu5 == 0) # x3
m.kkt1_x4 = Constraint(expr = 16 - m.mu2 - m.mu6 == 0) # x4

solver = SolverFactory('GUROBI')
resultados = solver.solve(m, tee=True, symbolic_solver_labels=True)

#printa o resultado
print('x1 =', m.x1())
print('x2 =', m.x2())
print('x3 =', m.x3())
print('x4 =', m.x4())
print('mu1 =', m.mu1())
print('mu2 =', m.mu2())
print('mu3 =', m.mu3())
print('mu4 =', m.mu4())
print('mu5 =', m.mu5())
print('mu6 =', m.mu6())
print('Objetivo =', m.obj())
assert m.obj() == 1.0, "Objective value does not match expected value."
a=1