import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy as sc
import scipy.special

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
# folder_name = '/home/michele/python_for_md/Droplet20nmExp/SubstrateSpreading/Q2'
# folder_name = '/home/michele/python_for_md/Droplet20nmExp/SubstrateRetract'
folder_name_adv = '/home/michele/python_for_md/Droplet20nmExp/SubstrateRetract'
folder_name_rec = '/home/michele/python_for_md/Droplet20nmExp/SubstrateSitting'

file_names_adv = [
    folder_name_adv+'/Flat/hbnum.xvg',
    folder_name_adv+'/Height2/hbnum.xvg',
    folder_name_adv+'/Height1/hbnum.xvg',
    folder_name_adv+'/Wave1/hbnum.xvg',
    folder_name_adv+'/Wave2/hbnum.xvg',
    folder_name_adv+'/Wave3/hbnum.xvg',
    folder_name_adv+'/Wave4/hbnum.xvg',
    folder_name_adv+'/Wave5/hbnum.xvg'
]
file_names_rec = [
    folder_name_rec+'/Flat/hbnum.xvg',
    folder_name_rec+'/Height2/hbnum.xvg',
    folder_name_rec+'/Height1/hbnum.xvg',
    folder_name_rec+'/Wave1/hbnum.xvg',
    folder_name_rec+'/Wave2/hbnum.xvg',
    folder_name_rec+'/Wave3/hbnum.xvg',
    folder_name_rec+'/Wave4/hbnum.xvg',
    folder_name_rec+'/Wave5/hbnum.xvg'
]


file_names_q5 = [
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
# a = [0.00, 0.75, 1.00, 1.25, 1.50, 1.75]
a_max = max(a)

r = lambda a : (2.0/np.pi) * np.sqrt(a+1.0) * sc.special.ellipe(a/(a+1.0))

time_adv = read_xvg_file(file_names_adv[1])[0]
time_rec = read_xvg_file(file_names_rec[1])[0]

time_q5 = read_xvg_file(file_names_q5[1])[0]

h_bonds_adv = dict()
h_bonds_rec = dict()
h_bonds_q5 = dict()
for i in range(len(a)) :
    h_bonds_adv[a[i]] = read_xvg_file(file_names_adv[i])[1]
    h_bonds_rec[a[i]] = read_xvg_file(file_names_rec[i])[1]
    h_bonds_q5[a[i]] = read_xvg_file(file_names_q5[i])[1]

# A little adjustment
h_bonds_adv[a[0]] = np.concatenate( (h_bonds_adv[a[0]], \
    np.NaN*np.ones( len(h_bonds_adv[a[1]]) - len(h_bonds_adv[a[0]]) ) ) )
# h_bonds_q5[a[0]] = np.concatenate( (h_bonds_q5[a[0]], \
#     np.NaN*np.ones( len(h_bonds_q5[a[1]]) - len(h_bonds_q5[a[0]]) ) ) )

for a_val in a :
    if a_val==0.00 :
        plt.plot(time_q5, h_bonds_q5[a_val], label='a='+str(a_val), linewidth=3.00, c=cm.hot(a_val/(a_max+1)))
    else :
        plt.plot(time_q5, h_bonds_q5[a_val], label='a='+str(a_val), linewidth=1.50, c=cm.hot(a_val/(a_max+1)))
plt.legend(fontsize=20.0)
plt.xticks(fontsize=20.0)
plt.yticks(fontsize=20.0)
plt.xlabel('t [ps]', fontsize=20.0)
plt.ylabel('# H-bonds', fontsize=20.0)
plt.xlim([0,max(time_adv)])
plt.title('Comparison relative displacement', fontsize=20.0)
plt.show()

n_init = 150

log_time = np.log(time_q5[n_init:-1])
log_nhb = dict()
for a_val in a :
    log_nhb[a_val] = np.log(h_bonds_q5[a_val][n_init:-1])

# Find power-law exponent
alpha = []
for a_val in a :
    alpha.append( np.polyfit(log_time, log_nhb[a_val], deg=1)[0] )

alpha = np.array(alpha)
alpha_avg = np.mean(alpha)
print(alpha_avg)

for a_val in a :
    if a_val==0.00 :
        plt.loglog(time_q5[n_init:-1], h_bonds_q5[a_val][n_init:-1], label='a='+str(a_val), linewidth=3.00, c=cm.hot(a_val/(a_max+1)))
    else :
        plt.loglog(time_q5[n_init:-1], h_bonds_q5[a_val][n_init:-1], label='a='+str(a_val), linewidth=1.50, c=cm.hot(a_val/(a_max+1)))
plt.legend(fontsize=20.0)
plt.xticks(fontsize=20.0)
plt.yticks(fontsize=20.0)
plt.xlabel('t [ps]', fontsize=20.0)
plt.ylabel('# H-bonds', fontsize=20.0)
plt.xlim([time_q5[n_init],max(time_q5)])
plt.title('Comparison relative displacement', fontsize=20.0)
plt.show()

N_avg = 20
hb_den_adv = np.zeros(len(a))
hb_den_rec = np.zeros(len(a))

for i in range(len(a)) :
    hb_den_adv[i] = (np.sum(h_bonds_adv[a[i]][-1-N_avg:-1])/N_avg) / r(a[i])
    hb_den_rec[i] = (np.sum(h_bonds_rec[a[i]][-1-N_avg:-1])/N_avg) / r(a[i])
    # hb_den_adv[i] = (np.sum(h_bonds_adv[a[i]][-1-N_avg:-1])/N_avg)
    # hb_den_rec[i] = (np.sum(h_bonds_rec[a[i]][-1-N_avg:-1])/N_avg)

plt.plot(a, hb_den_adv, 'b^--', linewidth=1.25, markersize=17.5, label='advancing')
plt.plot(a, hb_den_rec, 'rv--', linewidth=1.25, markersize=17.5, label='receding')
plt.legend(fontsize=20.0)
plt.xticks(fontsize=20.0)
plt.yticks(fontsize=20.0)
plt.xlabel('a [nondim.]', fontsize=20.0)
plt.ylabel('#HB / r [nondim.]', fontsize=20.0)
plt.title('Comparison final #HB', fontsize=20.0)
plt.show()
