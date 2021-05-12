# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 16:54:16 2021

@author: Adriano
"""

import numpy as np


def staggeredGrid(input_vec, frames, **kwargs):

    
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
    sou_x = kwargs['sou_x']    # Source position (in grid points) in X
    sou_y = kwargs['sou_y']    # Source position (in grid points) in Y


    # No vetor grid, pontos com valor 1 estao associados a Fase 1 (rocha,
    # cor preta). Os pontos com valor 0 estao associados a Fase 2 (poros,
    # de cor branca).        
    rho1 = 3.0
    vp1 = 4000.0
    #
    rho2 = 2.2
    vp2 = 2500.0

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
    
    
    
    dt_dx = dt / dx  
        
    #CFL_number = (cmax * dt) / dx        
    CFL_number = cmax * dt_dx
    
    
    print ()
    print("CFL_number: ")
    print(CFL_number)
    
    
    lampda_min = cmin / fmax        # Smallest wavelength
    
   
    
 
    
    
    
    #dx = cmin/(fmax*c1)             # Spatial discretization (in m)
    #dy = dx                         # Spatial discretization (in m)

    
    
    
    # ## Create space and time vector
    x = np.arange(0, dx * nx, dx) # Space vector in X
    y = np.arange(0, dy * ny, dy) # Space vector in Y
    
    t = np.arange(0, T, dt)            # Time vector
    nt = np.size(t)                    # Number of time steps
    
            
    print()
    print('t.shape:')
    print(t.shape, T, dt, t[0], t[-1])
    print()    
    
    print(dx, dy, dt)
        
    # ## Source signal - Ricker-wavelet
    tau = np.pi * f0 * (t - 1.5 / f0)
    
    
    print('tau')
    print(tau.min(), tau.max())
    
    wavelet = amp * (1.0 - 2.0 * tau**2.0) * np.exp(-tau**2)
            
    print()
    print('wavelet.shape:')
    print(wavelet.shape)
    print()
    
    
    
    # ## Source signal - Ricker-wavelet

#    tau = np.pi*f0*(t-1.5/f0)
#    q = q0*(1.0-2.0*tau**2.0)*np.exp(-tau**2)
    
    
    
    
    
    
    
    
    
    
    # Calculation of coefficients c1 and c2. 
    # These are used to obtain a second-order accurate time, fourth-order 
    # accurate space. (Levander, 1988)
    c1 = 9.0 / 8.0
    c2 = 1.0 / 24.0
    #
    # c1_dx = c1 / dx
    # c2_dx = c2 / dx
    # #
    # c1_dy = c1 / dy
    # c2_dy = c2 / dy
    #
    
    
    # The true calculation starts here...    
    #
    for it in range(frames):
        
        # if (it % 10 == 0):
        
        #print('Calculating SG [' + str(it+1) + '/' + str(frames) +']')
        
        if it<2:
            continue

        wavefield[it, :, :] = wavefield[it-1, :, :]
        
        
        
        # Update velocity
        for kx in range(5, nx-4):    
            for ky in range(5, ny-4):
                
                # Stress derivatives, p_dx(+) e p_dy(+) 
                
                p_x = c1 * (wavefield[it, ky, kx+1] - wavefield[it, ky, kx]) - \
                      c2 * (wavefield[it, ky, kx+2] - wavefield[it, ky, kx-1])  # Eq. A-2 Lavender, 1988
                
                p_y = c1 * (wavefield[it, ky+1, kx] - wavefield[it, ky, kx]) - \
                      c2 * (wavefield[it, ky+2, kx] - wavefield[it, ky-1, kx])


                # Velocity extrapolation using Euler Method
                
#                print("vy[ky, kx]: ", ky, kx, vy[ky, kx] )

                try:
                    vx[ky, kx] -=  (dt_dx / rho_grid[ky, kx]) * p_x 
                    vy[ky, kx] -=  (dt_dx / rho_grid[ky, kx]) * p_y   
                except:
                    print("rho_grid[ky, kx]: ", ky, kx, rho_grid[ky, kx] )
                    print("vx[ky, kx]: ", ky, kx, vx[ky, kx] )
                    print("vy[ky, kx]: ", ky, kx, vy[ky, kx] )
                    
                    pass
                #vx[ky, kx] = vx[ky, kx] - dt / rho_grid[ky, kx] * p_x
                #vy[ky, kx] = vy[ky, kx] - dt / rho_grid[ky, kx] * p_y



        # Inject source wavelet
        if it < np.size(wavelet+2):
            
            wavefield[it, sou_y, sou_x] += wavelet[it]

        # Verificar a possibilidade abaixo - Curso W4V8 
        #vx[sou_y, sou_x] = vx[sou_y, sou_x]  + dt * wavelet[it] / (dt * rho_grid[sou_y, sou_x])
        #vy[sou_y, sou_x] = vy[sou_y, sou_x]  + dt * wavelet[it] / (dt * rho_grid[sou_y, sou_x])
        #



        # Update wavefield
        for kx in range(5, nx-4):     
            for ky in range(5, ny-4):
                
                # Velocity derivatives, vx_dx(-) e vy_dy(-)  # Qian et al 2013
                vx_x = c1 * (vx[ky, kx] - vx[ky, kx-1]) - c2 * (vx[ky, kx+1] - vx[ky, kx-2])  # Dx-, 
                vy_y = c1 * (vy[ky, kx] - vy[ky-1, kx]) - c2 * (vy[ky+1, kx] - vy[ky-2, kx])  # Dy-, 
                #
                wavefield[it, ky, kx] -=  (dt_dx * lame_lambda[ky, kx]) * (vx_x + vy_y)
                
                
                # try:
                #     termo = lame_lambda[ky, kx] * dt #* (vx_x + vy_y)
                #     termo = (vx_dx + vy_dy) * termo
                # except:
                #     print()
                #     print(it)
                #     print(lame_lambda[ky, kx])
                #     print(vx_dx + vy_dy)
                #     print(dt)
                #     print()
                    
                    
                #wavefield[it, ky, kx] = wavefield[it, ky, kx] - termo   # lame_lambda[ky, kx] * dt * (vx_x + vy_y)
 
    
    print("Wavefield Min-Max: ", np.min(wavefield), " - ",  np.max(wavefield))   
    print(dx, dy, dt)
    return wavefield, dx, dy, dt






