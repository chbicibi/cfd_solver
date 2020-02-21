import argparse
import glob
import os
import re
import shutil
import subprocess
from contextlib import contextmanager

import numpy as np
from matplotlib import pyplot as plt

import myutils as ut
import mk_grid as mg

'''NOTE: u, v, p, f, w'''

i = 0
n = [2, 4][i]
nm = ['images_p', 'image_w'][i]
force = False


def plot(dest):
    flag = False
    with ut.chdir(dest):
        ut.mkdir(nm)
        files = ut.iglobm('result/*.npy')
        for file in files:
            basename = ut.basename(file, '.*')
            image = f'{nm}/{basename}.png'
            if os.path.isfile(image) and not force:
                continue
            print(basename)

            data = np.load(file)
            # print(data.shape)
            val = data[:, :, n]
            if n == 2:
                a = val[val.shape[0]//2, 0]
                val -= a

            fig = plt.figure()
            fig.subplots_adjust(left=0.08, right=1, bottom=0, top=1)
            ax = fig.add_subplot(111)

            im = ax.imshow(val, vmin=-1, vmax=1)
            cax = fig.colorbar(im)

            ax_pos = ax.get_position()
            cax_pos0 = cax.ax.get_position()
            cax_pos1 = [cax_pos0.x0, ax_pos.y0, cax_pos0.x1 - cax_pos0.x0, ax_pos.y1 - ax_pos.y0]
            cax.ax.set_position(cax_pos1)

            # plt.show()
            fig.savefig(image)
            plt.close('all')
            flag = True
    return flag



def mk_v(dest):
    with ut.chdir(os.path.join(dest, nm)):
        subprocess.call(['ffmpeg', '-framerate', '30', '-y',
                        '-i', 'out_%04d.plt.png', '-vcodec', 'libx264',
                        '-pix_fmt', 'yuv420p', '-r', '30', 'out.mp4'])


def main():
    for d in os.listdir('.'):
        if os.path.isdir(d) and os.path.isfile(d+'/in2d.txt'):
            if plot(d):
                mk_v(d)


if __name__ == '__main__':
    main()

# set path=C:\Tools\ffmpeg-4.2.2-win64-static\bin;%path%
# ffmpeg -framerate 30 -i out_%04d.plt.png -vcodec libx264 -pix_fmt yuv420p -r 30 out.mp4
