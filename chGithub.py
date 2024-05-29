"""
File Name: chGithub.py
Author: edgar avalos
Last Modified: 2024-05-27
Version: 1.0

Description:
This script solves cahn-hilliard coupled equations in 2D.

Dependencies:
- numpy
- matplotlib
- pyfftw

Permission Request:
If you would like to use or modify this code for your own projects, please contact
the author to request permission.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import integrate
import pyfftw
import multiprocessing
from pyfftw.interfaces.scipy_fftpack import fftn, ifftn

# Configure Matplotlib to use LaTeX for text rendering
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'


# Constants
N_STEPS = 2400
NS=200
DT = 0.002
N = 80
DX = 0.016
L = N * DX
NOISE = 0.4
C0 = 0.0
U0=-0.4
V0=0.0
EPSILON2 = 0.0025
SIGMA=20.0
B1=0.0 
B2=1.0 
TV=10.0
SEED = 12346 


def initialize_uv():
    rng = np.random.default_rng(SEED)
    u = np.empty((N_STEPS, N, N), dtype=np.float32)
    v = np.empty((N_STEPS, N, N), dtype=np.float32)
    u[0] = U0 + NOISE * rng.uniform(-1, 1, u[0].shape)
    v[0] = V0 + NOISE * rng.uniform(-1, 1, v[0].shape)

    # TotEnergy1= np.empty((N_STEPS), dtype=np.float32)
    # TotEnergy2= np.empty((N_STEPS), dtype=np.float32)
    # TotEnergy3= np.empty((N_STEPS), dtype=np.float32)
    # totu= np.empty((N_STEPS), dtype=np.float32)
    # totw= np.empty((N_STEPS), dtype=np.float32)

    return u, v

def update_uv(u,v, u_hat, v_hat, dfdu_hat, dfdv_hat, uv_hat, v_hat2):
    kx = ky = np.fft.fftfreq(N, d=DX) * 2 * np.pi
    #kx[0] = 10**(-12)
    kx[0] = kx[0] + 10**(-3)
    K = np.array(np.meshgrid(kx , ky ,indexing ='ij'), dtype=np.float32)
    K2 = np.sum(K*K, axis=0, dtype=np.float32)
    x = np.linspace(0, L, N, endpoint=False)
    y=x
    #X, Y = np.meshgrid(x, y, indexing='ij')
    Kx, Ky = np.meshgrid(kx, ky, indexing='ij')
    

    u_hat[:] = fftn(u[0])
    v_hat[:] = fftn(v[0])

    TotEnergy1 = []
    TotEnergy2 = []
    TotEnergy3 = []
    totu = []
    totw = []
    surfEnergy = None
    for i in range(1, N_STEPS):
        dfdu_hat[:] = fftn(u[i-1]**3 - u[i-1]) # for u
        dfdv_hat[:] = fftn(v[i-1]**3 - v[i-1]) # for v
        v_hat2[:]=fftn(v[i-1]**2)
        uv_hat[:] = fftn(np.multiply(u[i-1], v[i-1]))
        u_hat[:] = (u_hat - DT * (K2 * dfdu_hat)  - DT*(B1*K2 * v_hat - 0.5 * B2 * K2 * v_hat2)) /   \
                  (1 + DT * EPSILON2 * K2**2)
        v_hat[:] = (v_hat - (DT/TV) * (K2 * dfdv_hat)  - (DT/TV)*(B1*K2 * u_hat - B2 * K2 * uv_hat+ SIGMA * v_hat)) /   \
                  (1 + (DT/TV) * EPSILON2 * K2**2)
        # inverse fourier transform
        u[i] = ifftn(u_hat).real
        v[i] = ifftn(v_hat).real
        
        if i % NS == 0: 
            # ************************************************************************************************
            # energy calculation, E1, E2, E3
            # E1: energy of double well potential
            energy1 = (u[i]**2-1)**2 /4 + (v[i]**2-1)**2 /4+ \
                    B1 * u[i] *v[i] - 0.5 * B2 * u[i] * v[i]**2
            # Perform the double integral
            TotEnergy1.append(integrate.trapz(integrate.trapz(energy1, x, dx=DX), y, dx=DX))
            
            # E2: energy of the gradient
            #  term corresponding to u: abs(grad(u))^2
            du_dx_hat = -1j * Kx * u_hat
            du_dy_hat = -1j * Ky * u_hat

            du_dx = ifftn(du_dx_hat)  
            du_dy = ifftn(du_dy_hat) 

            # Compute the gradient magnitude squared of 
            grad_mag_squared_u = abs(du_dx)**2 + abs(du_dy)**2

            # term corresponding to v: abs(grad(v))^2
            dv_dx_hat = -1j * Kx * v_hat
            dv_dy_hat = -1j * Ky * v_hat
            
            dv_dx = ifftn(dv_dx_hat)
            dv_dy = ifftn(dv_dy_hat)

            # Compute the gradient magnitude squared of
            grad_mag_squared_v = abs(dv_dx)**2 + abs(dv_dy)**2

            energy2 = 0.5*EPSILON2*grad_mag_squared_u + 0.5* EPSILON2*grad_mag_squared_v

            # Compute the integral using the trapezoidal rule
            TotEnergy2.append(integrate.trapz(integrate.trapz(energy2, x, dx=DX), y, dx=DX))

            # E3:
            #energy3 = 0.5*SIGMA_term
            frctinlLaplcinfourier = 1.0 / np.sqrt((Kx)**2 + (Ky)**2) # bcos is  the norm

            inversetrans=ifftn(frctinlLaplcinfourier * v_hat)
            energy3=0.5*SIGMA*abs(inversetrans)**2

            # Compute the integral using the trapezoidal rule
            TotEnergy3.append(integrate.trapz(integrate.trapz(energy3, x, dx=DX), y, dx=DX))
            
            # integrand to plot energy density
            surfEnergy=energy1+energy2+energy3
            # total mass:

            # mass for uP
            totu.append(integrate.trapz(integrate.trapz(u[i], x, dx=DX), y, dx=DX))

            # mass for vP
            totw.append(integrate.trapz(integrate.trapz(v[i], x, dx=DX), y, dx=DX))

            # ************************************************************************************************
    TotEnergy1 = np.array(TotEnergy1)
    TotEnergy2 = np.array(TotEnergy2)
    TotEnergy3 = np.array(TotEnergy3)
    totu = np.array(totu)
    totw = np.array(totw)
    return u,v, TotEnergy1, TotEnergy2, TotEnergy3,  totu, totw, surfEnergy


# 3D plot (two components)
def plot_uv(u, v):
    x = np.arange(0, N, 1)
    y = np.arange(0, N, 1)
    X, Y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(12, 6))

    # Plot u
    ax1 = fig.add_subplot(121, projection='3d')  # 121 means 1 row, 2 columns, first plot
    im1 = ax1.plot_surface(X, Y, u[-1], cmap=cm.coolwarm, alpha=0.9)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('order parameter')
    cbar1 = plt.colorbar(im1, ax=ax1, orientation='vertical')
    cbar1.set_label('order parameter u')

    # Plot v
    ax2 = fig.add_subplot(122, projection='3d')  # 122 means 1 row, 2 columns, second plot
    im2 = ax2.plot_surface(X, Y, v[-1], cmap=cm.viridis, alpha=0.9)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_zlabel('order parameter')
    cbar2 = plt.colorbar(im2, ax=ax2, orientation='vertical')
    cbar2.set_label('order parameter v')

    plt.tight_layout()
    plt.show()  

    # plot integrand of energy functional:
    # 3D plot
def plot_surfEnergy(surfEnergy):
    x = np.arange(0, N, 1)
    y = np.arange(0, N, 1)
    X, Y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(12, 6))

    # Plot energy integrand
    ax1 = fig.add_subplot(111, projection='3d')  # 121 means 1 row, 2 columns, first plot
    im1 = ax1.plot_surface(X, Y, surfEnergy, cmap=cm.plasma, alpha=0.9)
    ax1.set_xlabel('x', fontsize=14)
    ax1.set_ylabel('y', fontsize=14)
    ax1.set_zlabel(r'$f[u[r], v[r], \nabla[u], \nabla[v]]$', fontsize=14)
    cbar1 = plt.colorbar(im1, ax=ax1, orientation='vertical')
    cbar1.set_label('free energy density')


def plot_energy(TotEnergy1, TotEnergy2, TotEnergy3):
    fig = plt.figure(figsize=(12, 6))
    #  start plotting from the second element (index 1)
    plt.plot(TotEnergy1[0:], label='Energy 1')
    plt.plot(TotEnergy2[0:], label='Energy 2')
    plt.plot(TotEnergy3[0:], label='Energy 3')
    plt.plot(TotEnergy1[0:] + TotEnergy2[0:] + TotEnergy3[0:], label='Total Energy')
    plt.xlabel('Time step')
    plt.ylabel('Energy')
    plt.ylim(0, 1)
    plt.title('Total Energy')
    plt.legend()
    plt.show()   

    # fucntion to plot mass
def plot_mass(totu, totw):
    fig = plt.figure(figsize=(12, 6))
    #  start plotting from the second element (index 1)
    plt.plot(totu[0:], label='mass u')
    plt.plot(totw[0:], label='mass v')
    plt.xlabel('Time step')
    plt.ylabel('mass')
    #plt.ylim(0, 1)
    plt.title('Total Mass')
    plt.legend()
    plt.show()  



def main():
    u,v =initialize_uv()
    u_hat = np.empty((N, N), dtype=np.complex64)
    v_hat = np.empty((N, N), dtype=np.complex64)
    dfdu_hat = np.empty((N, N), dtype=np.complex64)
    dfdv_hat = np.empty((N, N), dtype=np.complex64) 
    uv_hat = np.empty((N, N), dtype=np.complex64)
    v_hat2 = np.empty((N, N), dtype=np.complex64)
    u,v,TotEnergy1, TotEnergy2, TotEnergy3, totu, totw, surfEnergy= update_uv(u,v, u_hat, v_hat, dfdu_hat, dfdv_hat, uv_hat, v_hat2)
    plot_uv(u,v)
    plot_surfEnergy(surfEnergy)
    plot_energy(TotEnergy1,TotEnergy2, TotEnergy3)
    plot_mass(totu, totw)


if __name__ == "__main__":
    main()