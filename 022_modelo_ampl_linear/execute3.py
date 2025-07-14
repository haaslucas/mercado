from contextlib import redirect_stdout
from amplpy import AMPL
import pathlib
import pandas as pd
import re
import sys
# (1) Onde estão Novo.mod, input.dat, etc…
basedir = pathlib.Path(".")
ampl = AMPL()
#ampl.setOption('solver', 'ipopt')

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

# 3. Caminho de saída --------------------------------------------------
out_file = pathlib.Path("Resultados_AMPL.xlsx")

# 4. ExcelWriter (confirme se o motor está instalado) ------------------
try:
    import xlsxwriter                     # assegura presença
    engine = "xlsxwriter"
except ImportError:
    print("xlsxwriter não instalado; usando openpyxl.", file=sys.stderr)
    engine = "openpyxl"

sheet_names, wrote_something = set(), False

with pd.ExcelWriter(out_file, engine=engine, mode="w") as writer:

    # ----------- VARIÁVEIS -----------
    for name, var_obj in ampl.getVariables():
        df = var_obj.getValues().toPandas()
        if df.empty:
            df = pd.DataFrame({"value": [var_obj.value()]})
        wrote_something |= add_df(writer, df,
                                  safe_sheet_name(name, sheet_names),
                                  sheet_names)

    # ----------- PARÂMETROS ----------
    for name, par_obj in ampl.getParameters():
        df = par_obj.getValues().toPandas()
        if df.empty:
            df = pd.DataFrame({"value": [par_obj.value()]})
        wrote_something |= add_df(writer, df,
                                  safe_sheet_name(name, sheet_names),
                                  sheet_names)

    # ----------- RESTRIÇÕES ----------
    for name, con_obj in ampl.getConstraints():
        df = con_obj.getValues().toPandas()        # body, dual, slack
        wrote_something |= add_df(writer, df,
                                  safe_sheet_name(f"con_{name}", sheet_names),
                                  sheet_names)

    # ----------- OBJETIVOS -----------
    for objfun in ampl.getObjectives():
        df = pd.DataFrame({
            "value": [objfun.value()],
            "sense": [objfun.sense()]              # "minimize" ou "maximize"
        })
        wrote_something |= add_df(writer, df,
                                  safe_sheet_name(f"obj_{objfun.name()}", sheet_names),
                                  sheet_names)

    # Caso nada tenha sido gravado até aqui, cria aba dummy -------------
    if not wrote_something:
        pd.DataFrame({"msg": ["Sem dados para exportar"]}) \
          .to_excel(writer, sheet_name="Empty", index=False)

print(f"Arquivo '{out_file}' gerado com {len(sheet_names)} abas visíveis.")