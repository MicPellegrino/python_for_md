# Reading input .xvg file

import sys

filename = sys.argv[1]
print('xvg_reader: reading ' + filename)
ext = filename.split('.')[-1]
if ext!="xvg":
    print("Python .xvg reader: ERROR invalid input file extension!")
    quit()

# Storing the input

vars = {}
n_vars = 0
s = 0

title = ''
xlabel = ''
ylabel = ''

with open(filename) as f:
    for line in f:
        cols = line.split()
        if s==1:
            for i in range(n_vars):
                vars[i].append(cols[i])
        elif cols[0][0]=='@':
            if cols[1] == 'title':
                title = ' '.join(cols[2:]).replace('"','')
            if cols[1] == 'xaxis':
                xlabel = ' '.join(cols[3:]).replace('"','')
            if cols[1] == 'yaxis':
                ylabel = ' '.join(cols[3:]).replace('"','')
        elif cols[0][0]!='#':
            s = 1
            n_vars = len(cols)
            for i in range(n_vars):
                vars[i] = [cols[i]]

# Matplotlib test

if n_vars == 0:
    print('Python .xvg reader: ERROR zero variables to plot!')
    quit()
elif n_vars == 1:
    print('Python .xvg reader: WARNING only one variable to plot!')
else:
    print('Python .xvg reader: '+str(n_vars-1)+' variable(s) will be plotted against the first one.')

import matplotlib.pyplot as plt
import numpy as np


vars_dim = {}
for i in range(n_vars):
    vars_dim[i] = np.asarray(vars[i]).astype(np.float)

# Rescaling: by default, the first variable (i.e. time) is not rescaled,
# while other variables are rescaled on their mean

vars_ref = {}
vars_ndm = {}

# Rescaling only if there are more than 2 dependent variables to plot
if n_vars > 2:
    for i in range(n_vars-1):
        vars_ref[i+1] = np.mean(vars_dim[i+1])
        vars_ndm[i+1] = vars_dim[i+1]/vars_ref[i+1]
        vars_ndm[0] = vars_dim[0]
else:
    vars_ndm = vars_dim

# Plotting

if n_vars > 2:
    for i in range(n_vars-1):
        plt.plot(vars_ndm[0], vars_ndm[i+1], label='ref = '+str(vars_ref[i+1]))
else:
    plt.plot(vars_ndm[0], vars_ndm[1])

plt.legend()
plt.title(title)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.show()
