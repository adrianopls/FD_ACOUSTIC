# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 16:54:16 2021

@author: Adriano
"""

import numpy as np


def rotatedStaggeredGrid(input_vec, frames, **kwargs):

    
    dx = kwargs['dx']
    dy = kwargs['dy']
    #dt = kwargs['dt']

    f0 = kwargs['f0']
    amp = kwargs['amp']


    # Discretization
    c1 = 20   # Number of grid points per dominant wavelength
    c2 = 0.5  # CFL-Number
    
    nx = input_vec.shape[0]  # Number of grid points in X
    ny = input_vec.shape[1]  # Number of grid points in Y
    T = 1.0     # Total propagation time
    
    # Source Signal
    #f0 = 500000       # Center frequency Ricker-wavelet
    #q0 = 30           # Maximum amplitude Ricker-Wavelet


    # Source Position
    sou_x = 128   # Source position (in grid points) in X
    sou_y = 128   # Source position (in grid points) in Y

    # No vetor grid, pontos com valor 1 estao associados a Fase 1 (rocha,
    # cor preta). Os pontos com valor 0 estao associados a Fase 2 (poros,
    # de cor branca).        
    rho1 = 3.0
    vp1 = 5000.0
    #
    rho2 = 2.2
    vp2 = 3000.0    

    # Velocity and density 
    vp_grid = np.where(input_vec > 0.0, vp1, vp2)
    rho_grid = np.where(input_vec > 0.0, rho1, rho2)
    
    
    ## Preparation
    
    # Init wavefields
    vx = np.zeros((ny, nx))
    vy = np.zeros((ny, nx))
    wavefield = np.zeros((frames, ny, nx))
        
 
    # Calculate first Lame-Paramter
    #self.lame_lambda = self._rho_grid * self._vp_grid * self._vp_grid
    
    
    lambda1 = rho1 * vp1 * vp1
    lambda2 = rho2 * vp2 * vp2
    
    lame_lambda = np.where(input_vec > 0.0, lambda1, lambda2)
    
    
    cmin = min(vp_grid.flatten())   # Lowest P-wave velocity
    cmax = max(vp_grid.flatten())   # Highest P-wave velocity
    
    fmax = 2 * f0                   # Maximum frequency
    



    
    #dx = 15.0 #0.000001 * 74 * 4
    #dy = dx                         # Spatial discretization (in m)
    dt = dx / (cmax) * c2           # Temporal discretization (in s)
    
    
    CFL_number = (cmax * dt) / dx        
    
    
    print ()
    print("CFL_number:")
    print(CFL_number)
    
    
    lampda_min = cmin / fmax        # Smallest wavelength
    
    
    
    
    dx = cmin/(fmax*c1)             # Spatial discretization (in m)
    dy = dx                         # Spatial discretization (in m)

    dr = np.sqrt((dx**2.0 + dy**2.0))
    
    
    # ## Create space and time vector
    # x = np.arange(0, dx * nx, dx) # Space vector in X
    # y = np.arange(0, dy * ny, dy) # Space vector in Y
    
    t = np.arange(0, T, dt)            # Time vector
    nt = np.size(t)                    # Number of time steps
    
            
    print()
    print('t.shape:')
    print(t.shape, T, dt, t[0], t[-1])
    print()    
    
    print(dx, dy, dt, dr)
        
    # ## Source signal - Ricker-wavelet
    tau = np.pi * f0 * (t - 1.5 / f0)
    
    print('tau')
    print(tau)
    
    wavelet = amp * (1.0 - 2.0 * tau**2.0) * np.exp(-tau**2)
            
    print()
    print('wavelet.shape:')
    print(wavelet.shape)
    print()
    
    
    
    # ## Source signal - Ricker-wavelet

#    tau = np.pi*f0*(t-1.5/f0)
#    q = q0*(1.0-2.0*tau**2.0)*np.exp(-tau**2)
    
    
    
    
    
    
    
    
    
    

    c1 = 1.0 / dr
    #
    c2_x = dr / (2*dy)
    c2_y = dr / (2*dy)
    # #
    # c1_y = c1 / dy
    # c2_y = c2 / dy
    #
    
    # The true calculation starts here...    
    #
    for it in range(frames):
        
        if (it % 10 == 0):
        
            print('Calculating SG [' + str(it+1) + '/' + str(frames) +']')
        
        if it<2:
            continue

        wavefield[it, :, :] = wavefield[it-1, :, :]
        
        
        
        # Update velocity
        for kx in range(5, nx-4):    
            for ky in range(5, ny-4):
                            
                # Calcula as derivadas do Stress, p_xr e p_yr 
                p_xr = (c1 * (wavefield[it, ky, kx+1] - wavefield[it, ky+1, kx]))   # eq. 8                     
                p_yr = (c1 * (wavefield[it, ky+1, kx+1] - wavefield[it, ky, kx]))   # eq. 9    
                       
                # Calcula as derivadas do Stress, p_x e p_y 
                p_x = c2_x * (p_yr + p_xr)    # eq. 11
                p_y = c2_y * (p_yr - p_xr)    # eq. 10
                
                # Extrapola velocidade
                vx[ky, kx] = vx[ky, kx] - dt / rho_grid[ky, kx] * p_x
                vy[ky, kx] = vy[ky, kx] - dt / rho_grid[ky, kx] * p_y



        # Inject source wavelet
        wavefield[it, sou_y, sou_x] = wavefield[it, sou_y, sou_x] + wavelet[it]

        # Update wavefield
        for kx in range(5, nx-4):     
            for ky in range(5, ny-4):
                
                
                vx_xr = (c1 * (vx[ky, kx+1] - vx[ ky+1, kx]))   # eq. 8   
                vy_yr = (c1 * (vx[ky+1, kx+1] - vx[ky, kx]))    # eq. 9  
                
                # Calcula as derivadas do Stress, p_x e p_y 
                vx_x = c2_x * (vy_yr + vx_xr)    # eq. 11
                vy_y = c2_y * (vy_yr - vx_xr)    # eq. 10               
                

                try:
                    termo = lame_lambda[ky, kx] * dt #* (vx_x + vy_y)
                    termo = (vx_x + vy_y) * termo
                except:
                    print()
                    print(it)
                    print(lame_lambda[ky, kx])
                    print(vx_x + vy_y)
                    print(dt)
                    print()
                    
                wavefield[it, ky, kx] = wavefield[it, ky, kx] - termo   # lame_lambda[ky, kx] * dt * (vx_x + vy_y)
 
    
    print("Wavefield Min-Max: ", np.min(wavefield), " - ",  np.max(wavefield))   
    print(dx, dy, dt)
    return wavefield, dx, dy, dt






