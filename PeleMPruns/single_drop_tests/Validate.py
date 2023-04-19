

import os
from CaseInfo import *
from ExtractData import *
import matplotlib.pyplot as plt

case = Tonini_4_33()

[refdvals, reftvals, refyvals] = ExtractRefVals(case)
# Set end time based on reference values
time = refdvals[-1, 0] / case.xconv
case.set_end_time(time)

os.system("mkdir -p {}".format(case.name))
os.system("rm -r {}/plt*".format(case.name))
params = CreateInputParams(case)

executable = "None"
for f in os.listdir("./"):
    if (f.startswith("Pele") and f.endswith(".ex")):
        executable = f
if (not os.path.exists(executable)):
    error = "Pele executable not found"
    raise ValueError(error)

input_file = "gen-input"
run_cmd = "mpiexec -np 4 ./"
os.system("{}{} {} {}".format(run_cmd, executable, input_file, params))

outfile = case.name + "/pele_vals.csv"
pele_vals = ExtractData(case, outfile)

numplots = 2
if (refyvals is not None):
    numplots += 1
plotnames = ["Diameter", "Temperature", "Mass fraction"]
fig, axs = plt.subplots(numplots)
for i in range(numplots):
    axs[i].set_title(plotnames[i])
    axs[i].plot(pele_vals[:,0], pele_vals[:,i+1], label="Pele", color='black')
    if (i == 0):
        refarr = refdvals
    elif (i == 1):
        refarr = reftvals
    elif (i == 2):
        refarr = refyvals
    if (refarr is not None):
        axs[i].scatter(refarr[:,0], refarr[:,1], label="Ref", color='red', s = 2.)
plt.legend()
plt.show()
plt.savefig(case.name + "/results.png")
