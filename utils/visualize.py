#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.animation as manimation
import matplotlib.pyplot as plt
import sys

header_t = np.dtype([('nx', 'i4'), ('ny', 'i4'), ('num_iter', 'i4')])

def main(output='mpi.out', animation_file='diffusion_mpi.gif'):
    # Load data
    header = np.fromfile(output, dtype=header_t, count=1)
    nx = header['nx'][0]
    ny = header['ny'][0]
    num_iter = header['num_iter'][0] + 1
    print(f"nx={nx}, ny={ny}, num_iter={num_iter}")
    h = np.fromfile(output, offset=header_t.itemsize, dtype='f8').reshape(num_iter, nx, ny)

    # Create coordinate grid
    length = 1.0
    width = 1.0
    x = np.linspace(0, length, nx)
    y = np.linspace(0, width, ny)
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    # Find global maximum for scaling
    hmax = np.max(h)
    hmin = np.min(h)
    
    # Create figure with proper spacing
    fig = plt.figure(figsize=(15, 10))
    plt.subplots_adjust(left=0.1, right=0.85, bottom=0.1, top=0.9, wspace=0.4, hspace=0.4)
    
    # Create subplots
    ax1 = fig.add_subplot(221, projection='3d')
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    
    # Create colormap and normalization
    norm = plt.Normalize(hmin, hmax)
    
    # Create a single colorbar axis
    cax = plt.axes([0.88, 0.1, 0.02, 0.8])
    
    def update_anim(i):
        # Clear previous plots
        ax1.cla()
        ax2.cla()
        ax3.cla()
        cax.cla()
        
        Z = h[i]
        
        # 3D surface plot
        surf = ax1.plot_surface(X, Y, Z, cmap='viridis', norm=norm)
        ax1.set_title(f'Concentration field at step {i}')
        ax1.set_zlim(0, hmax)
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_zlabel('Concentration')
        
        # Contour plot
        levels = np.linspace(hmin, hmax, 20)
        cont = ax2.contourf(X, Y, Z, levels=levels, cmap='viridis', norm=norm)
        plt.colorbar(cont, cax=cax, label='Concentration')
        ax2.set_title(f'Concentration contours at step {i}')
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        
        # Cross-section plot
        ax3.plot(x, Z[:, ny//2], 'b-', label='y = 0.5')
        ax3.plot(y, Z[nx//2, :], 'r--', label='x = 0.5')
        ax3.set_title(f'Cross-sections at step {i}')
        ax3.set_xlabel('Position')
        ax3.set_ylabel('Concentration')
        ax3.set_ylim(0, hmax)
        ax3.grid(True)
        ax3.legend()
        
        # Add overall title
        plt.suptitle(f'Diffusion Simulation - Step {i}\nMax: {np.max(Z):.6f}, Min: {np.min(Z):.6f}', 
                    y=0.95)
    
    # Create animation
    animation_fig = manimation.FuncAnimation(
        fig, update_anim,
        frames=num_iter,
        interval=100
    )
    
    # Save animation
    animation_fig.save(animation_file, writer='pillow')
    print(f"Created animation: {animation_file}")

if __name__ == "__main__":
    main(*sys.argv[1:])