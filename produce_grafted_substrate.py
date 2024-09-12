import numpy as np
import numpy.random as rng
import scipy as sc
import mdconf_oop as md

def grafting_chain(n_chain=6, x0=0, y0=0, z0=0, file_name='graft-chain.gro') :

    n_atoms = 3*(n_chain-2)+4
    Lx = 1.0
    Ly = 1.0
    Lz = 1.0

    out_file = open(file_name, 'w')
    out_file.write('GRAFTING ALKANE CHAIN\n')
    out_file.write(str(n_atoms)+'\n')

    n = 1

    for ic in range(1,n_chain-1) :
        
        atomname='C'
        atomnum=n
        x=x0
        y=y0
        z=z0+ic*md.d_cc
        atom = md.Atom(atomname,atomnum,x,y,z,0,0,0)
        out_file.write(atom.gro_line(ic,'CH2'))
        n += 1

        atomname='H1'
        atomnum=n
        x=x0+md.d_ch
        y=y0
        z=z0+ic*md.d_cc
        atom = md.Atom(atomname,atomnum,x,y,z,0,0,0)
        out_file.write(atom.gro_line(ic,'CH2'))
        n += 1

        atomname='H2'
        atomnum=n
        x=x0-md.d_ch
        y=y0
        z=z0+ic*md.d_cc
        atom = md.Atom(atomname,atomnum,x,y,z,0,0,0)
        out_file.write(atom.gro_line(ic,'CH2'))
        n += 1

    atomname='C'
    atomnum=n
    x=x0
    y=y0
    z=z0+(n_chain-1)*md.d_cc
    atom = md.Atom(atomname,atomnum,x,y,z,0,0,0)
    out_file.write(atom.gro_line((n_chain-1),'CH3'))
    n += 1

    atomname='H1'
    atomnum=n
    x=x0+md.d_ch
    y=y0
    z=z0+(n_chain-1)*md.d_cc
    atom = md.Atom(atomname,atomnum,x,y,z,0,0,0)
    out_file.write(atom.gro_line((n_chain-1),'CH3'))
    n += 1

    atomname='H2'
    atomnum=n
    x=x0-md.d_ch
    y=y0
    z=z0+(n_chain-1)*md.d_cc
    atom = md.Atom(atomname,atomnum,x,y,z,0,0,0)
    out_file.write(atom.gro_line((n_chain-1),'CH3'))
    n += 1

    atomname='H3'
    atomnum=n
    x=x0
    y=y0
    z=z0+(n_chain-1)*md.d_cc+md.d_ch
    atom = md.Atom(atomname,atomnum,x,y,z,0,0,0)
    out_file.write(atom.gro_line((n_chain-1),'CH3'))
    n += 1

    out_file.write("%10.5f%10.5f%10.5f\n" % (Lx,Ly,Lz))

    out_file.close()