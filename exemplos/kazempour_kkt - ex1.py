# exemplo retirado de Slide12 da Lecture 1
# https://www.youtube.com/watch?v=1h3b0j2k4a8&list=PL6X9d7c5e2f8f7c4b&index=1

from pyomo.environ import *
import numpy as np

m = ConcreteModel()


# Conjuntos
m.G = Set(initialize=['G1', 'G2'])  # Geradores
m.N = Set(initialize=['N1', 'N2', 'N3'])  # Nós
m.D = Set(initialize=['D1', 'D2'])  # Demandas
m.MapG = Set(initialize=[('G1', 'N1'), ('G2', 'N2')], within = m.G*m.N)  # Geradores por nó
m.MapD = Set(initialize=[('D1', 'N2'), ('D2', 'N3')], within = m.D*m.N)  # Demandas por nó
m.MapN = Set(initialize=[('N1','N2'),('N1','N3'),('N2','N3'),('N2','N1'),
('N3','N1'),('N3','N2')], within=m.N*m.N)  # Topologia da rede 

# Parâmetros
m.PGMAX  = Param(m.G, initialize = {'G1':100, 'G2':150})
m.PGCOST = Param(m.G, initialize = {'G1':12, 'G2':20})
m.PDMAX  = Param(m.D, initialize = {'D1':100, 'D2':50})
m.PDCOST = Param(m.D, initialize = {'D1':40, 'D2':35})
m.FMAX   = Param(m.MapN, initialize = {('N1','N2'):100, ('N2','N1'):100, ('N1','N3'):100, ('N3','N1'):100,
('N2','N3'):100, ('N3','N2'):100})
m.SUSCEP = Param(m.MapN, initialize = {('N1','N2'):500, ('N2','N1'):500, ('N1','N3'):500, ('N3','N1'):500,
('N2','N3'):500, ('N3','N2'):500})

# Variáveis
m.Pg = Var(m.G)  # Geração
m.Pd = Var(m.D)  # Demanda
m.Flow = Var(m.MapN)  
m.Theta = Var(m.N)  

# Restrições
m.pdMin    = Constraint(m.D, rule = lambda m,d: m.Pd[d] >= 0)
m.pdMax    = Constraint(m.D, rule = lambda m,d: m.Pd[d] <= m.PDMAX[d])
m.pgMin    = Constraint(m.G, rule = lambda m,g: m.Pg[g] >= 0)
m.pgMax    = Constraint(m.G, rule = lambda m,g: m.Pg[g] <= m.PGMAX[g])
m.thetaMin = Constraint(m.N, rule = lambda m,n: m.Theta[n] >= -np.pi/2)
m.thetaMax = Constraint(m.N, rule = lambda m,n: m.Theta[n] <= np.pi/2)

m.fluxo     = Constraint(m.MapN, rule = lambda m,i,j: m.Flow[i,j] == m.SUSCEP[i,j] *  (m.Theta[i] - m.Theta[j] )) 

m.balanco = Constraint(
    m.N,
    rule=lambda m, n:
        -sum(m.Pg[g]        for g in m.G if (g, n) in m.MapG) #
        +sum(m.Pd[d]        for d in m.D if (d, n) in m.MapD)
        +sum(m.Flow[n,  m2]     for m2 in m.N if (n, m2) in m.MapN)
        -sum(m.Flow[m1, n]      for m1 in m.N if (m1, n) in m.MapN)
        == 0
)
'''-sum(m.Pg[g] for (g,nn) in m.GN if nn == n)
+sum(m.Pd[d] for (d,nn) in m.DN if nn == n)
+sum(m.Flow[n,m2] for (n2,m2) in m.TOPO if n2 == n)
-sum(m.Flow[m1,n] for (m1,n2) in m.TOPO if n2 == n)
== 0'''



m.obj = Objective(
    expr = sum(m.PDCOST[d]*m.Pd[d] for d in m.D)
    -sum(m.PGCOST[g]*m.Pg[g] for g in m.G),
    sense=maximize
)

solver = SolverFactory('GUROBI')
resultados = solver.solve(m, tee=True, symbolic_solver_labels=True)
# Imprime o resultado
print('Geração:')
for g in m.G:
    print(f'{g}: {m.Pg[g].value}')

print('Demanda:')
for d in m.D:
    print(f'{d}: {m.Pd[d].value}')
    
print('Fluxos:')
for (i, j) in m.MapN:
    print(f'Fluxo de {i} para {j}: {m.Flow[i, j].value}')   
    
print('Ângulos de fase:')
for n in m.N:
    print(f'Ângulo em {n}: {m.Theta[n].value}')
    
print('Bem-estar social:', m.obj())


