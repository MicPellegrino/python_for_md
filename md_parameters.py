"""
    Script to generate .mdp files
    See GROMACS documentation for further information
"""

from math import sqrt
from decimal import Decimal

# Phisical parametes [IS units]
mu = 8.77e-4        # viscosity [Pa*s]
rho = 986.0         # density [kg*m^-3]
D = 2.29e-9         # diffusion coefficient [m^2*s^-1]
gamma = 5.78e-2     # surface tension [Pa*m]
epsilon = 0.75e-9   # interface width [m]
R = 50.0e-9         # initial droplet radius [m]
T_ref = 300.0       # reference temperature [K]
p_ref = 1.0e5       # reference pressure [Pa]

###############################################################################
print("PHYSICAL PARAMETERS")
print("viscosity [Pa*s] = %E" % mu)
print("density [kg*m^-3] = %E" % rho)
print("diffusion coefficient [m^2*s^-1] = %E" % D)
print("surface tension [Pa*m] = %E" % gamma)
print("interface width [m] = %E" % epsilon)
print("initial droplet radius [m] = %E" % R)
print("reference temperature [K] = %E" % T_ref)
print("reference pressure [Pa] = %E" % p_ref)
###############################################################################

# ADD NON-DIMENSIONAL NUMBERS (IT MAY BE USEFUL)!

# Time scales [IS UNITS]
tau_mu = mu*R/gamma                 # viscous time scale
tau_rho = sqrt(rho*(R**3)/gamma)    # inertial time scale
tau_D = (R**2)/D                    # diffusion time scale

###############################################################################
print("TIME SCALES")
print("viscous t.s. [s] = %E" % tau_mu)
print("inertial t.s. [s] = %E" % tau_rho)
print("diffusive t.s. [s] = %E" % tau_D)
###############################################################################

# Ratio between total simulation time and viscous time
M_tau = 20
# Ratio between viscous time and time step
alpha_tau = 1e6

# Times [PICOSECONDS]
# t_tem = -1 : no temperature coupling
t_ini = 0.0
t_max = 1e12*M_tau*tau_mu
# NB: the time step should be small enough to capture hydrogen atoms vibrations;
#       typical values are on the order of 5fs; for water simulation better to stay
#       around 4fs or below
# dt = 1e12*tau_mu/alpha_tau
dt = 0.004
nsteps = M_tau*alpha_tau
t_tem = 0.1
# t_prs = 0.1
t_prs = -1.0

###############################################################################
print("MD TIME PARAMETERS (|!| picoseconds |!|)")
print("initial time [ps] = %f" % t_ini)
print("maximum time horizon [ps] = %f" % t_max)
print("time step [ps] = %f" % dt)
print("number of timesteps = %d" % nsteps)
print("temperature coupling time [ps] = %f" % t_tem)
print("pressure coupling time [ps] = %f" % t_prs)
###############################################################################

# Output time binning
nsample = 100           # number of sampling bins
t_sample = 100*dt       # sampling time window
nstxout = nsample       # number of steps that elapse between writing coordinates to output trajectory file, the last coordinates are always written
nstvout = nsample       # number of steps that elapse between writing velocities to output trajectory, the last velocities are always written
nstlog  = nsample       # number of steps that elapse between writing energies to the log file, the last energies are always written
nstenergy = nsample     # number of steps that else between writing energies to energy file, the last energies are always written

###############################################################################
print("sampling window [ps] = %f" % t_sample)
print("sample size = %d" % nsample)
###############################################################################

# Boundary conditions

pbc = "xyz"                     # xyz = p.b.c. on all edge
                                # no = use no p.b.c; see GROMACS documentation
                                # xy = use p.b.c. in x and y directions only; see GROMACS documentation
nwall = 0                       # 1 = wall at z=0; 2 = walls at z=0 and z=box_height
wall_type = "12-6"              # type of VdW interaction; see GROMACS documentation
wall_atomtype = "opls_740"      # need to be defined in a separate .itp file
wall_r_linpot = 0.0             # distance below which the wall potential is continued

if nwall > 0 :
    assert(pbc == "xy", "Walls are allowed if p.b.c. are set on x and y only")

###############################################################################
print("BOUNDARY CONDITIONS (|!| nanometers |!|)")
print("periodic boundary conditions: "+pbc)
print("number of walls = %d" % nwall)
print("type of wall: "+wall_type)
print("type of wall atoms: "+wall_atomtype)
print("below-wall potential distance [nm] = %f" % wall_r_linpot)
###############################################################################

# Neighbour searching
nstlist = 10                # frequency to update the neighbor list
ns_type = "grid"            # grid = make a grid in the box and only check atoms in neighboring grid cells when constructing a new neighbor list every nstlist steps
                            # simple = check every atom in the box when constructing a new neighbor list every nstlist steps
rlist = 1.0                 # cut-off distance for the short-range neighbor list [nm]
cutoff_scheme = "verlet"    # verlet = generate a pair list with buffering
                            # group = generate a pair list for groups of atoms

if pbc == "xy" :
    assert(ns_type == "grid", "Only neighbours search of type grid allowed with xy p.b.c.")

###############################################################################
print("NEIGHBOUR SEARCH")
print("search type: "+ns_type)
print("cutoff scheme: "+cutoff_scheme)
print("update frequency = %d" % nstlist)
print("cutoff distance [nm] = %f" % rlist)
###############################################################################

# Electrostatics and other non-bonded interactions
coulombtype = "pme"         # see GROMACS documentation
fourierspacing = 0.15       # spacing in FFT grid [nm]
pme_order = 6               # interpolation order for PME (4 equals cubic interpolation)
ewald_rtol = 1e-5           # relative strength of the Ewald-shifted direct potential at rcoulomb (?)
# ewald_geometry = "3dc"    # see documentations
ewald_geometry = "3d"
rcoulomb = 1.0              # distance for the Coulomb cut-off [nm]
vdw_type = "Cut-off"        # see GROMACS documentation
rvdw = 1.0                  # distance for the LJ or Buckingham cut-off [nm]
DispCorr = "EnerPres"       # no = don’t apply any correction
                            # EnerPres = apply long range dispersion corrections for Energy and Pressure
                            # Ener = apply long range dispersion corrections for Energy only

###############################################################################
print("ELECTROSTATICS AND VAN DER WAALS")
print("coulomb force computation: "+coulombtype)
print("spacing for fft [nm] = %f" % fourierspacing )
print("interpolation order for PME = %d" % pme_order)
print("pme relatiove strenght = %E" % ewald_rtol)
print("pme geometry: "+ewald_geometry)
print("cutoff distance coulomb [nm] = %f" % rcoulomb)
print("vdw force computation: "+vdw_type)
print("cutoff distnce vdw [nm] = %f" % rvdw)
print("long-range correction: "+DispCorr)
###############################################################################

# Termostat
tcoupl = "v-rescale"    # see GROMACS documentation
tc_grps = "System"      # groups to couple separately to temperature bath

###############################################################################
print("THERMOSTAT")
print("temperature coupling type: "+tcoupl)
print("coupled groups: "+tc_grps)
###############################################################################

# Barostat
# pcoupl = "berendsen"      # see GROMACS documentation
pcoupl = "no"             # see GROMACS documentation
pcoupltype = "isotropic"    # see GROMACS documentation
compressibility = 5e-5      # compressibility [bar^-1]
refcoord_scaling = "all"    # see GROMACS documentation

###############################################################################
print("BAROSTAT (|!| bar |!|)")
print("pressure coupling type: "+pcoupl)
print("coupling application type: "+pcoupltype)
print("compressibility [bar^-1] = %f" % compressibility)
print("reference coord. rescaling: "+refcoord_scaling)
###############################################################################

# Constraints
constraints = "all-bonds"   # see GROMACS documentation

###############################################################################
print("CONSTRAINTS")
print("constraint interactions: "+constraints)
###############################################################################

# Velocity generation
gen_vel = "no"          # no = do not generate velocities; the velocities are set to zero when there are no velocities in the input structure file
                        # yes = generate velocities in grompp according to a Maxwell distribution at temperature gen-temp [K], with random seed gen-seed
gen_seed = 1231234      # used to initialize random generator for random velocities; when gen-seed is set to -1, a pseudo random seed is used

###############################################################################
print("VELOCITY GENERATION")
print("generate velocities: "+gen_vel)
print("random seed = %d" % gen_seed)
###############################################################################

# GENERATE .mdp FILE

file_name = "mdpar.mdp"

###############################################################################
print("MD FILE NAME: "+file_name)
###############################################################################

md_param = open(file_name, 'w')

md_param.write("integrator\t\t = md\n")

md_param.write("dt\t\t\t = %f\n" % dt)
md_param.write("nsteps\t\t\t = %d\n" % nsteps)

md_param.write("nstxout\t\t\t = %d\n" % nstxout)
md_param.write("nstvout\t\t\t = %d\n" % nstvout)
md_param.write("nstlog\t\t\t = %d\n" % nstlog)
md_param.write("nstenergy\t\t = %d\n" % nstenergy)

md_param.write("pbc\t\t\t = "+pbc+"\n")
md_param.write("nwall\t\t\t = %d \n" % nwall)
md_param.write("wall-type\t\t = "+wall_type+" \n")
if nwall == 1 :
    md_param.write("wall-atomtype\t\t = "+wall_atomtype+" \n")
elif nwall == 2 :
    md_param.write("wall-atomtype\t\t = "+wall_atomtype+" "+wall_atomtype+" \n")
md_param.write("wall-r-linpot\t\t = %f\n" % wall_r_linpot)

md_param.write("nstlist\t\t\t = %d\n" % nstlist)
md_param.write("ns-type\t\t\t = "+ns_type+"\n")
md_param.write("rlist\t\t\t = %f\n" % rlist)
md_param.write("cutoff-scheme\t\t = "+cutoff_scheme+"\n")

md_param.write("coulombtype\t\t = "+coulombtype+"\n")
md_param.write("fourierspacing\t\t = %f\n" % fourierspacing)
md_param.write("pme-order\t\t = %d\n" % pme_order)
md_param.write("ewald-rtol\t\t = %e\n" % ewald_rtol)
md_param.write("rcoulomb\t\t = %f\n" % rcoulomb)
md_param.write("ewald-geometry\t\t = "+ewald_geometry+"\n")

md_param.write("vdw-type\t\t = "+vdw_type+"\n")
md_param.write("rvdw\t\t\t = %f\n" % rvdw)
md_param.write("DispCorr\t\t = "+DispCorr+"\n")

md_param.write("tcoupl\t\t\t = "+tcoupl+"\n")
md_param.write("tc-grps\t\t\t = "+tc_grps+"\n")
md_param.write("tau-t\t\t\t = %f\n" % t_tem )
md_param.write("ref-t\t\t\t = %d\n" % T_ref )

md_param.write("pcoupl\t\t\t = "+pcoupl+"\n")
md_param.write("pcoupltype\t\t = "+pcoupltype+"\n")
md_param.write("tau-p\t\t\t = %f\n" % t_prs )
md_param.write("compressibility\t\t = %e\n" % compressibility)
md_param.write("ref-p\t\t\t = %f\n" % p_ref)
md_param.write("refcoord-scaling\t = "+refcoord_scaling+"\n")

md_param.write("constraints\t\t = "+constraints+"\n")

md_param.write("gen-vel\t\t\t = "+gen_vel+"\n")
md_param.write("gen-seed\t\t = %d\n" % gen_seed)

md_param.close()
