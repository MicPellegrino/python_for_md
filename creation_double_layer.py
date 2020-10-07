import mdconf as md
import numpy as np

dist = 100.0    # [Å]
ni = int( 300/4.50 )
nj = int( 50/(4.50*0.5*np.sqrt(3.0)) )
nk = 1
output = '/home/michele/python_for_md/ChannelFlow/silica.pdb'

md.quad_double_layer(dist, ni, nj, nk, output)
