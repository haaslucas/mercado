from contextlib import redirect_stdout
from amplpy import AMPL
import pathlib

# (1) Onde estão Novo.mod, input.dat, etc…
basedir = pathlib.Path(".")
ampl = AMPL()
ampl.setOption('solver', 'ipopt')

ampl.eval('include execute_linear.run;')             # chama o .run

# Função genérica para extrair nomes de objetos (Variable, Parameter, Constraint)
def get_names(obj_list):
    names = []
    for item in obj_list:
        if hasattr(item, 'name'):
            names.append(item.name())
        elif isinstance(item, tuple):
            for part in item:
                if hasattr(part, 'name'):
                    names.append(part.name())
                    break
    return names

# Captura nomes de variáveis, parâmetros e restrições
var_names   = get_names(ampl.getVariables())
param_names = get_names(ampl.getParameters())
con_names   = get_names(ampl.getConstraints())
set_names   = get_names(ampl.getSets())

# 1) Arquivo com todas as variáveis expandidas
with open('Variables.txt', 'w', encoding='utf-8') as f, redirect_stdout(f):
    for name in var_names:
        ampl.eval(f"expand {name};")

# 2) Arquivo com todos os parâmetros expandidos
with open('Parameters.txt', 'w', encoding='utf-8') as f, redirect_stdout(f):
    for name in param_names:
        try:
            ampl.eval(f"expand {name};")
        except Exception as e:
            print(f"Erro ao expandir {name}: {e}")
            continue

# 3) Arquivo com todas as restrições expandidas
with open('Constraints.txt', 'w', encoding='utf-8') as f, redirect_stdout(f):
    for name in con_names:
        ampl.eval(f"expand {name};")

# 4) Arquivo com todos os conjuntos expandidos
with open('Sets.txt', 'w', encoding='utf-8') as f, redirect_stdout(f):
    for name in set_names:
        ampl.eval(f"display {name};")
a=1