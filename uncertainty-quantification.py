import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import numpy.random as rng
import os.path

from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

from matplotlib.ticker import MaxNLocator

fs_labels = 60
fs_ticks = 50

n_rep = 60
n_col = 3

list_file_names = []

# root = "Glycerol100p/ViscosityReplicas/"
# root = "Glycerol100p/ViscosityReplicas2ps/"
# y_upper = 35000
# t_min = 15000

# root = "Glycerol80p/ViscosityReplicas/BackupReplicas/"
# root = "Glycerol80p/ViscosityReplicas2ps/"
# y_upper = 400
# t_min = 1000

# root = "Glycerol60p/ViscosityReplicas/BackupReplicas/"
# root = "Glycerol60p/ViscosityReplicas2ps/"
# y_upper = 75
# t_min = 200

# root = "Glycerol40p/ViscosityReplicas/"
# root = "Glycerol40p/ViscosityReplicas2ps/"
# y_upper = 35
# t_min = 10

# root = "Glycerol20p/ViscosityReplicas/"
# root = "Glycerol20p/ViscosityReplicas2ps/"
# y_upper = 25
# t_min = 10

# root = "PureWater/ViscosityReplicas/"
# root = "/home/michele/python_for_md/PureWater/NVT-GK/"
root = "PureWater/ViscosityReplicas2ps/"
y_upper = 15
t_min = 10

# root = "/home/michele/python_for_md/HexaneOPLSAA/Viscosity/"
# y_upper = 30
# t_min = 100

# root = "/home/michele/python_for_md/PentanolmOPLSAAp/Viscosity/"
# root = "/home/michele/python_for_md/PentanolOPLSAA/Viscosity/"
# y_upper = 200
# t_min = 1000

for n in range(n_rep) :

    fn = root+"Rep"+str(n)+"/eviscoi.xvg"
    list_file_names.append(fn)

print("Analyzing Eintein integrals...")

times = []
integrals = []
n_min = 1e9
t_max = 0.0

# User-set parameters: tolerance on the max statistical error and number of initial steps to throw away
# toll = 0.05
toll = 0.075
# toll = 0.10

# Number of bootstrapping steps
kboot = 20

f, (ax1, ax2) = plt.subplots(1, 2)


"""
    Extracting the Einstein integrals from the
"""

n_corrupt = 0

for filename in list_file_names :

    print("Reading "+filename)

    vars = {}
    n_vars = 0
    s = 0

    title = ''
    xlabel = ''
    ylabel = '' 

    if os.path.isfile(filename) :

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

        vars_dim = {}
        for i in range(n_vars):
            vars_dim[i] = np.asarray(vars[i]).astype(np.float)

        n = len(vars_dim[0])
        t = np.array(vars_dim[0])
        
        if t[-1] > t_max :
            t_max = t[-1]

        for nnc in range(n_col) :
            I = np.array(vars_dim[1+nnc])    
            # ax1.plot((1e-3)*t, I, color="lightgray", linewidth=2)
            times.append(t)
            integrals.append(I)

        if n<n_min :
            n_min=n
    
    else :
        print("/!\ File "+str(filename)+" does not exist! Probaby the corresponding energy.edr file is corrupted")
        n_corrupt += 1

n_rep = n_rep-n_corrupt

# assert 4*n_rep == len(integrals), "! ERROR: number of replicas inconsistent with number of realizations of the observable I(t) !"

"""
    Comuting the average of the observables
"""

dt = times[0][1]-times[0][0]

idx_min = int(t_min/dt)

t_avg = times[0][0:n_min]
I_avg = np.zeros(n_min)
I2_avg = np.zeros(n_min)
for I in integrals :
    I_avg += I[0:n_min]
    I2_avg += I[0:n_min]*I[0:n_min]
I_avg /= len(integrals)
I2_avg /= len(integrals)

# Prepare plotting ...
ax1.plot((1e-3)*t_avg, I_avg, 'k-', linewidth=4.5, label='mean')

ax1.set_ylabel(r"$I^2$ [$10^6\times$Pa$\cdot$ns$^2$]", fontsize=fs_labels)
ax1.set_xlabel("$t$ [ns]", fontsize=fs_labels)
ax1.tick_params('both', labelsize=fs_ticks)
# ax1.ticklabel_format(axis='x', scilimits=[0,0])
# ax1.xaxis.get_offset_text().set_fontsize(fs_ticks)
ax1.set_ylim([0,y_upper])
ax1.set_xlim([0,(1e-3)*t_max])


"""
    Computing the normalized standard error and the max time given the tolerance
    (with a bit of data cleaning in case of nans)
"""

standard_deviation = np.sqrt(I2_avg-I_avg*I_avg)
standard_error = standard_deviation/np.sqrt(len(integrals))
ax1.plot((1e-3)*t_avg, I_avg+standard_error, 'k--', linewidth=3, label=r'$\pm$ standard error')
ax1.plot((1e-3)*t_avg, I_avg-standard_error, 'k--', linewidth=3)
ax1.fill_between((1e-3)*t_avg, I_avg+standard_deviation, I_avg-standard_deviation, color='lavender', label='standard deviation')

normalized_standard_error = standard_error/I_avg
norm_se = normalized_standard_error[~np.isnan(normalized_standard_error)]
t_avg_se = t_avg[~np.isnan(normalized_standard_error)]
idx = np.argmin(np.abs(norm_se-toll))
print("Max time for allowed tolerance: "+str(t_avg_se[idx])+"ns")
idx_max = np.argmin(np.abs(t_avg-t_avg_se[idx]))


"""
    Performing the linear fit of the Einstein integral in the time frame that balance
    systematic error and statistical accuracy
"""

# ein_fun = lambda t, p : p*t
# popt, pcov = opt.curve_fit(ein_fun, t_avg[idx_min:idx_max], I_avg[idx_min:idx_max])
popt, pcov = np.polyfit(t_avg[idx_min:idx_max], I_avg[idx_min:idx_max], deg=1, cov=True, w=1/norm_se[idx_min:idx_max])
# I_fit = ein_fun(t_avg[idx_min:idx_max], *popt)
I_fit = np.polyval(popt, t_avg[idx_min:idx_max])
eta_hat = (1e3)*popt[0]

"""
    Performing bootstrapping to quantify the uncertainty on viscosity
"""

print("Performing bootstrapping...")
n_take = int(0.8*n_col*n_rep)
samp_range = np.arange(n_col*n_rep)
eta_boot = []
for k in range(kboot) :
    np.random.shuffle(samp_range)
    I_rep_samp = []
    for ib in samp_range[0:n_take] :
        I_rep_samp.append(integrals[ib])
    I_avg_boot = np.zeros(n_min)
    for I in I_rep_samp :
        I_avg_boot += I[0:n_min]
    I_avg_boot /= n_take
    popt, pcov = np.polyfit(t_avg[idx_min:idx_max], I_avg_boot[idx_min:idx_max], deg=1, cov=True, w=1/norm_se[idx_min:idx_max])
    eta_boot.append((1e3)*popt[0])
eta_boot = np.array(eta_boot)
err = np.std(eta_boot)

print("--------------------------------------------------------")
print("Viscosity estimate: "+str(eta_hat)+" +/- "+str(err)+" cp")
print("--------------------------------------------------------")

# Prepare plotting ...
# ax1.plot((1e-3)*t_avg[idx_min:idx_max], I_fit, 'r--', linewidth=3.5, label='fit')
ax1.legend(fontsize=fs_labels, loc='upper left')
# ax1.xaxis.set_major_locator(MaxNLocator(prune='lower'))
ax1.set_xlim([0,(1e-3)*t_avg[-1]])

# Prepare plotting ...
ax2.plot(t_avg_se, norm_se, 'b-', linewidth=3.5)
ax2.plot([0, t_avg_se[idx]], [norm_se[idx], norm_se[idx]], 'b--', linewidth=2.0)
ax2.plot([t_avg_se[idx], t_avg_se[idx]], [norm_se[idx], 0], 'b--', linewidth=2.0)
ax2.set_ylabel(r"$\sigma^*/\sqrt{N}$", fontsize=fs_labels)
ax2.set_xlabel("t [ps]", fontsize=fs_labels)
ax2.tick_params('both', labelsize=fs_ticks)
ax2.ticklabel_format(axis='x', scilimits=[0,0])
ax2.xaxis.get_offset_text().set_fontsize(fs_ticks)
ax2.set_ylim([0,0.2])
ax2.set_xlim([0,t_max])
ax2.xaxis.set_major_locator(MaxNLocator(prune='lower'))

plt.subplots_adjust(left=0.075, bottom=None, right=0.975, top=None, wspace=0.3, hspace=None)
plt.show()

# plt.loglog(t_avg, I_avg, 'ko', markersize=10, linewidth=3.0, label='<I>(t)')
plt.loglog(t_avg, I_avg, 'k-', markersize=10, linewidth=3.0)
# plt.loglog(t_avg, I_avg[-1]*t_avg/t_avg[-1], 'b--', linewidth=3.0)
plt.loglog(t_avg, (1e-1)*t_avg, 'b--', linewidth=4.0)
plt.xlabel('t [ps]', fontsize=fs_labels)
plt.ylabel('<I> [(kg*ps)/(m*s)]', fontsize=fs_labels)
plt.tick_params('both', labelsize=fs_ticks)
plt.show()
