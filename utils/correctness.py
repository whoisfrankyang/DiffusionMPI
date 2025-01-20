import numpy as np
import os
import re
import sys

import matplotlib.animation as manimation
import matplotlib.pyplot as plt

header_t = np.dtype([('nx', 'i4'), ('ny', 'i4'), ('num_iter', 'i4')])


def main(file1='serial.out', file2='mpi.out', out='correctness.gif'):
    print(f"Comparing {file1} and {file2}...")
    header1 = np.fromfile(file1, dtype=header_t, count=1)
    header2 = np.fromfile(file2, dtype=header_t, count=1)
    if (np.array_equal(header1, header2) == False):
        print("Headers are not equal, so why compare these two runs? They can never be the same.")
        return
    
    nx = header1['nx'][0]
    ny = header1['ny'][0]
    num_iter = header1['num_iter'][0] + 1

    dx = 1 / (nx - 1.0);
    dy = 1 / (ny - 1.0);

    C1 = np.fromfile(file1, offset=header_t.itemsize, dtype='f8').reshape((num_iter, nx, ny))
    C2 = np.fromfile(file2, offset=header_t.itemsize, dtype='f8').reshape((num_iter, nx, ny))

    x = np.arange(num_iter)
    y = np.max(np.abs(C1[:] - C2[:]), axis=(1, 2))

    max_error = np.max(y)

    print(f"Max Error: {max_error:.6e}")


    fig = plt.figure(figsize=(15, 15))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    ax1.plot(x, y)

    hx = (-1/2 + dx/2.0 + np.arange(nx)*dx)[:, np.newaxis]
    hy = (-1/2 + dy/2.0 + np.arange(ny)*dy)[np.newaxis, :]

    X, Y = np.meshgrid(hx, hy)

    def update_anim(i):
        Z = np.abs(C1[i] - C2[i])

        ax2.cla()
        ax2.contourf(X, Y, Z, cmap=plt.cm.RdBu)
        ax2.set_title(f"Error for t = {i}")

    animation_fig = manimation.FuncAnimation(fig, update_anim, frames=num_iter, interval=100)
    animation_fig.save(out)


if __name__ == "__main__":
    main(*sys.argv[1:])