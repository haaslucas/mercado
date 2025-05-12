from amplpy import AMPL
ampl = AMPL()
ampl.eval("""
reset;
model Novo.mod;
data input.dat;
""")

# ------------------------------------------------------------------
# 1) PELA API PYTHON
# ------------------------------------------------------------------
param_objs  = ampl.getParameters()            # lista de objetos Parameter
param_names = [p.name for p in param_objs]    # só os nomes

print(f"\nForam encontrados {len(param_names)} parâmetros:")
for n in param_names:
    p = ampl.getParameter(n)
    ndim = p.numIndices()                     # 0 = escalar, 1 = 1-D, …
    print(f" – {n:<20s}  (dim={ndim})")

# Se quiser salvar em arquivo:
# with open("lista_parametros.txt","w") as f:
#     f.write("\n".join(param_names))

# ------------------------------------------------------------------
# 2) VIA COMANDO AMPL  (mesma lista, estilo “console AMPL”)
# ------------------------------------------------------------------
ampl.eval("print _PARMS;")                    # ou 'show *;' se quiser tudo
print("\nSaída do AMPL ↓")
print(ampl.getOutput())                       # texto produzido pelo print

#sets
ampl.getSet('Trans_Nodes').get().to_pandas()   

#parametros
ampl.getValue('SE_Capacity')   
ampl.getParameter('load2').to_pandas()