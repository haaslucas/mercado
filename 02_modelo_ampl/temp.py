from amplpy import AMPL

SOLVER = "highs"

ampl = AMPL()
ampl.read("bim_dual.mod")
ampl.solve(solver=SOLVER)
assert ampl.solve_result == "solved", ampl.solve_result

print(
    f'y = ({ampl.var["y1"].value():.1f}, {ampl.var["y2"].value():.1f}, {ampl.var["y3"].value():.1f}, {ampl.var["y4"].value():.1f})'
)
print(f'optimal value = {ampl.obj["obj"].value():.2f}')