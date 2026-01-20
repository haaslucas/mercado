import pyomo.kernel as pmo

h = 0.75  # custo de manter um item por um ano
c = 500.0  # custo de processamento de um pedido
d = 10000.0  # demanda anual

m = pmo.block()

# definir variáveis para restrições cônicas
m.u = pmo.variable(lb=0)
m.v = pmo.variable(lb=0)
m.t = pmo.variable(lb=0)

# relações para restrições cônicas em relação às variáveis de decisão
m.u_eq = pmo.constraint(m.u == 2)
m.q = pmo.conic.quadratic(m.t, [m.u, m.v])

# objetivo linear
m.eoq = pmo.objective(((h + 2 * c * d) * m.t + (h - 2 * c * d) * m.v) / 4)

# resolver com 'mosek_direct' ou 'gurobi_direct'
SOLVER = pmo.SolverFactory("mosek_direct")
SOLVER.solve(m)
a=1