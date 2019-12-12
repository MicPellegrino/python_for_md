from math import log
from math import pi

import numpy as np

class weierstrass_mandelbrot :

    def __init__ ( self, G, D, gamma, n_0, n_inf, x_0, h_0 = 0, y_0 = 0 ) :
        self.scale = G
        self.fractal_dim = D
        self.mode_base = gamma
        self.low_cut = n_0
        self.high_cut = n_inf
        self.frequencies = dict()
        for n in range(n_0, n_inf) :
            self.frequencies[n] = 2*pi*(gamma**n)
        self.offset_x = x_0
        self.offset_y = y_0
        self.offset_h = h_0
        self.eval_x = np.empty(0, dtype=float)
        self.eval_y = np.empty(0, dtype=float)
        self.eval_h = np.empty(0, dtype=float)
        self.eval_structure = np.empty(0, dtype=float)

    def print( self ) :
        print("Parameters for Weierstrass-Mandelbrot function:")
        print("scale = %f" % self.scale)
        print("fractal dimension = %f" % self.fractal_dim)
        print("modes base = %f" % self.mode_base)
        print("low cutoff number = %d" % self.low_cut)
        print("high cutoff number = %d" % self.high_cut)
        print("offsets : x_0 = %f; y_0 = %f; h_0 = %f" % (self.offset_x, self.offset_y, self.offset_h))
        print("frequencies : ")
        print(self.frequencies)

    def evaluate_1D_vec ( self, x ) :
        self.eval_h = np.zeros(len(x), dtype=float)
        for n in self.frequencies :
            self.eval_h += (self.mode_base**(n*(self.fractal_dim-2))) * np.cos( (x + self.offset_x) * (self.frequencies[n]) )
        self.eval_h *= (self.scale)**(self.fractal_dim-1)
        self.eval_h += self.offset_h
        self.eval_x = x
        return self.eval_h

    def compute_structure_function( self ) :
        self.eval_structure = np.zeros(len(self.eval_h), dtype=float)
        for k in range(0,len(self.eval_h)) :
            self.eval_structure[k] = np.mean( (self.eval_h[0:-1-k]-self.eval_h[k:-1])**2 )
        return self.eval_structure

    def evaluate_2D_vec ( self, x, y ) :
        self.eval_h = np.zeros( shape=(len(x),len(y)) )
        for n in self.frequencies :
            self.eval_h += (self.mode_base**(n*(self.fractal_dim-3))) * np.outer( np.cos( (x + self.offset_x) * (self.frequencies[n]) ), np.cos( (y + self.offset_y) * (self.frequencies[n]) ) )
        # ANSATZ?
        self.eval_h *= (self.scale)**(self.fractal_dim-2)
        self.eval_x = x
        self.eval_y = y

    def evaluate_1D ( self, x ) :
        h = 0
        for n in self.frequencies :
            h += (self.mode_base**(n*(self.fractal_dim-2))) * np.cos( (x + self.offset_x) * (self.frequencies[n]) )
        h *= (self.scale)**(self.fractal_dim-1)
        h += self.offset_h
        return h
