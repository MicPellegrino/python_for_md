import numpy as np
import matplotlib.pyplot as plt

def read_xvg_file(file_name) :
    x = []
    y = []
    with open(file_name) as f:
        for line in f:
            cols = line.split()
            if len(cols) == 2 and cols[0][0] != '@':
                x.append(float(cols[0]))
                y.append(float(cols[1]))
    return x, y

file_names = [  'FreeEnergyCorrugated/H4/FE1/dhdl.xvg', 
                'FreeEnergyCorrugated/H4/FE2/dhdl.xvg', 
                'FreeEnergyCorrugated/H4/FE3/dhdl.xvg'  ]

t, _ = read_xvg_file( file_names[0] )

lam = np.array(t) / max(t)
delta_lam = lam[1]-lam[0]

dhdl = np.zeros( (len(file_names), len(t)), dtype=float )

for k in range(len(file_names)) :
    _, y = read_xvg_file( file_names[k] )
    dhdl[k,:] = np.array(y)

dhdl_avg = np.mean(dhdl, axis=0)

delta_h = np.cumsum( delta_lam*dhdl_avg )

dH = delta_h[-1]
print("dH = "+str(dH)+" kJ/mol")

plt.plot(lam, dhdl_avg, 'b-', linewidth=0.25)
plt.plot(lam, delta_h, 'r-', linewidth=1.25)
plt.plot(lam, np.zeros(len(lam)), 'k--')
plt.show()
