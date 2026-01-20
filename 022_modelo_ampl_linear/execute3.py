from contextlib import redirect_stdout
from amplpy import AMPL
import pathlib
import pandas as pd
import re
import sys
# (1) Onde estão Novo.mod, input.dat, etc…
basedir = pathlib.Path(".")
ampl = AMPL()
ampl.setOption('solver', 'gurobi')

ampl.eval('include execute_linear.run;')             # chama o .run

# 2. Utilidades --------------------------------------------------------
def safe_sheet_name(name: str, existing: set) -> str:
    """Nome de aba válido para Excel (≤31 caracteres, único e sem  \ / * ? : [ ] )."""
    name = re.sub(r'[\\/*?:\[\]]', '_', name)[:31]
    base = name
    k = 1
    while name in existing or not name:
        suf = f"_{k}"
        name = (base[:31 - len(suf)] + suf) if base else f"sheet{k}"
        k += 1
    existing.add(name)
    return name

def add_df(writer, df, sheet_name, used_sheets):
    """Grava df se ele não estiver vazio; devolve True se gravou algo."""
    if df.empty and len(df.columns) == 0:
        return False
    df.to_excel(writer, sheet_name=sheet_name, index=True)
    used_sheets.add(sheet_name)
    return True

# 1) Arquivo com todas as variáveis expandidas
with open('Variables.txt', 'w', encoding='utf-8') as f, redirect_stdout(f):
    for name in var_names:
        ampl.eval(f"display {name};")

# 4. ExcelWriter (confirme se o motor está instalado) ------------------
try:
    import xlsxwriter                     # assegura presença
    engine = "xlsxwriter"
except ImportError:
    print("xlsxwriter não instalado; usando openpyxl.", file=sys.stderr)
    engine = "openpyxl"

sheet_names, wrote_something = set(), False

# 4) Arquivo com todos os conjuntos expandidos
with open('Sets.txt', 'w', encoding='utf-8') as f, redirect_stdout(f):
    for name in set_names:
        ampl.eval(f"display {name};")
a=1
