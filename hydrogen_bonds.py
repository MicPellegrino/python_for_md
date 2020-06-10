import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Function to read .xvg data
def read_xvg_file (filename) :
    vars = {}
    n_vars = 0
    s = 0
    with open(filename) as f:
        for line in f:
            cols = line.split()
            if s==1:
                for i in range(n_vars):
                    vars[i].append(cols[i])
            elif cols[0][0]!='#' and cols[0][0]!='@':
                s = 1
                n_vars = len(cols)
                for i in range(n_vars):
                    vars[i] = [cols[i]]
    vars_dim = {}
    for i in range(n_vars):
        vars_dim[i] = np.asarray(vars[i]).astype(np.float)
    return vars_dim

# Main folder name
folder_name = '/home/michele/python_for_md/Droplet20nmExp/SubstrateSpreading/Q5'
file_names = [
    folder_name+'/Flat/hbnum.xvg',
    folder_name+'/Height2/hbnum.xvg',
    folder_name+'/Height1/hbnum.xvg',
    folder_name+'/Wave1/hbnum.xvg',
    folder_name+'/Wave2/hbnum.xvg',
    folder_name+'/Wave3/hbnum.xvg',
    folder_name+'/Wave4/hbnum.xvg',
    folder_name+'/Wave5/hbnum.xvg'
]

a = [0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75]
a_max = max(a)

time = read_xvg_file(file_names[0])[0]

h_bonds = dict()
for i in range(len(a)) :
    h_bonds[a[i]] = read_xvg_file(file_names[i])[2]

for a_val in a :
    if a_val==0.00 :
        plt.plot(time, h_bonds[a_val], label='a='+str(a_val), linewidth=3.00, c=cm.hot(a_val/(a_max+1)))
    else :
        plt.plot(time, h_bonds[a_val], label='a='+str(a_val), linewidth=1.50, c=cm.hot(a_val/(a_max+1)))
plt.legend(fontsize=20.0)
plt.xticks(fontsize=20.0)
plt.yticks(fontsize=20.0)
plt.xlabel('t [ps]', fontsize=20.0)
plt.ylabel('# H-bonds', fontsize=20.0)
plt.xlim([0,max(time)])
plt.title('Comparison relative displacement', fontsize=20.0)
plt.show()
