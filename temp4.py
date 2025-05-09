import pyomo.environ as pyo
import pyomo.opt as opt

class DSO_Model:
    def __init__(self):
        # Modelo Bilevel 
        self.model = pyo.ConcreteModel()
        
        # Conjuntos
        self.model.SCENARIOS = pyo.Set()  # s ∈ S
        self.model.PERIODS = pyo.Set()    # t ∈ T
        self.model.GENERATORS = pyo.Set() # g ∈ G
        self.model.NODES = pyo.Set()      # i ∈ N
        
        # Parâmetros
        self.model.scenario_weights = pyo.Param(self.model.SCENARIOS, 
                                                self.model.PERIODS)  # p_{s,t}
        self.model.generation_costs = pyo.Param(self.model.GENERATORS) # C^D_{g,b}
        self.model.carbon_price = pyo.Param(self.model.SCENARIOS, 
                                            self.model.PERIODS)  # ψ_{t,s}
        
        # Variáveis de Primeiro Nível (DSO)
        self.model.power_injection = pyo.Var(self.model.GENERATORS, 
                                             self.model.PERIODS, 
                                             self.model.SCENARIOS)
        self.model.carbon_allowances = pyo.Var(self.model.SCENARIOS, 
                                               self.model.PERIODS)
        
        # Objetivo de Primeiro Nível (Minimização de Custos)
        def objective_rule(model):
            return sum(
                model.scenario_weights[s,t] * (
                    sum(model.generation_costs[g] * 
                        model.power_injection[g,t,s] 
                        for g in model.GENERATORS) +
                    model.carbon_price[s,t] * model.carbon_allowances[s,t]
                )
                for s in model.SCENARIOS 
                for t in model.PERIODS
            )
        
        self.model.objective = pyo.Objective(rule=objective_rule, 
                                             sense=pyo.minimize)
        
        # Adicionar restrições conforme equações do artigo
        # (2)-(23) Upper-level constraints
        # (25)-(36) ISO's primal constraints
        
    def solve(self):
        # Configurar solver
        solver = pyo.SolverFactory('gurobi')
        results = solver.solve(self.model)
        return results

# Instanciar e resolver
dso_model = DSO_Model()
results = dso_model.solve()
