import argparse
import csv
import ctypes
import os
import re
import time
from dataclasses import dataclass
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as plc
from tqdm import tqdm
import myutils as ut


def register_functions():
    libname = 'cfd.dll'
    loader_path = '.'
    cdll = np.ctypeslib.load_library(libname, loader_path)

    f_read_inputfile = cdll.f_read_inputfile
    f_read_inputfile.argtypes = None
    f_read_inputfile.restype = None

    f_initialize = cdll.f_initialize
    f_initialize.argtypes = None
    f_initialize.restype = None

    f_advance = cdll.f_advance
    f_advance.argtypes = [
        *(np.ctypeslib.ndpointer(dtype=np.float64) for i in range(6)),
        *(np.ctypeslib.ndpointer(dtype=np.int32) for i in range(2))]
    f_advance.restype = None

    f_calc_velociry = cdll.f_calc_velociry
    f_calc_velociry.argtypes = [
        *(np.ctypeslib.ndpointer(dtype=np.float64) for i in range(6)),
        *(np.ctypeslib.ndpointer(dtype=np.int32) for i in range(2))]
    f_calc_velociry.restype = None

    f_calc_pressure = cdll.f_calc_pressure
    f_calc_pressure.argtypes = [
        *(np.ctypeslib.ndpointer(dtype=np.float64) for i in range(3)),
        *(np.ctypeslib.ndpointer(dtype=np.int32) for i in range(4)),
        np.ctypeslib.ndpointer(dtype=np.uint8)]
    f_calc_pressure.restype = None

    f_bind_velocity = cdll.f_bind_velocity
    f_bind_velocity.argtypes = [
        *(np.ctypeslib.ndpointer(dtype=np.float64) for i in range(3)),
        *(np.ctypeslib.ndpointer(dtype=np.int32) for i in range(2))]
    f_bind_velocity.restype = None

    f_testsub = cdll.testsub
    f_testsub.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.uint8),
        np.ctypeslib.ndpointer(dtype=np.int32)
    ]
    f_testsub.restype = None

    globals().update({k: v for k, v in cdll.__dict__.items()
                     if k.startswith('f_')})


def plot_w(val, file=None, vr=0.1):
    fig, ax = plt.subplots(figsize=(6, 4))
    fig.subplots_adjust(left=0.08, right=1, bottom=0, top=1)

    colors = [(0, '#ff0000'), (0.5, '#000000'), (1, '#00ff00')]
    cmap = plc.LinearSegmentedColormap.from_list('custom_cmap', colors)
    im = ax.imshow(val, cmap=cmap, vmin=-vr, vmax=vr)
    cax = fig.colorbar(im)

    ax_pos = ax.get_position()
    pos0 = cax.ax.get_position()
    pos1 = [pos0.x0, ax_pos.y0, pos0.x1 - pos0.x0, ax_pos.y1 - ax_pos.y0]
    cax.ax.set_position(pos1)

    if file:
        fig.savefig(file, bbox_inches='tight', pad_inches=0.1)
    else:
        plt.show()
    plt.close('all')


def plot_nitr(itrs, file):
    num = 1000
    ave = [np.mean(itrs[max(i-num+1, 0):i+1]) for i in range(len(itrs))]

    plt.plot(itrs, lw=1)
    plt.plot(ave, lw=1)
    # plt.xlim(0, 20000)
    # plt.ylim(0, ymax)
    plt.yscale('log')
    plt.xlabel('cycle')
    plt.ylabel('itr')
    # plt.show()
    # plt.legend()
    plt.savefig(file, bbox_inches='tight', pad_inches=0.1, dpi=600)


def vorticity(u, v):
    w = (u[:-2, :-1] + u[:-2, 1:] - u[2:, :-1] - u[2:, 1:] +
         v[:-1, 2:] + v[1:, 2:] - v[1:, :-2] - v[:-1, :-2]) / 4
    return -w


def dump_data(file, u, v, p):
    u_ = ((u[1:-1, :-1] + u[1:-1, 1:]) / 2).astype(np.float32)
    v_ = ((v[:-1, 1:-1] + v[1:, 1:-1]) / 2).astype(np.float32)
    p_ = p[1:-1, 1:-1].astype(np.float32)
    w_ = vorticity(u, v).astype(np.float32)
    np.save(file, np.stack([u_, v_, p_, w_]))


def load_value(file):
    with open(file, 'rb') as f:
        nx, ny, nb = (int.from_bytes(f.read(4), 'little') for i in range(3))
        u = np.frombuffer(f.read(nx*ny*nb), dtype=f'f{nb}').reshape((ny, nx))

        nx, ny, nb = (int.from_bytes(f.read(4), 'little') for i in range(3))
        v = np.frombuffer(f.read(nx*ny*nb), dtype=f'f{nb}').reshape((ny, nx))

        nx, ny, nb = (int.from_bytes(f.read(4), 'little') for i in range(3))
        p = np.frombuffer(f.read(nx*ny*nb), dtype=f'f{nb}').reshape((ny, nx))
    return u, v, p


def read_inputfile(file='in2d.txt'):
    with open(file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    nitr_, cycle_ = map(int, re.findall(r'\d+', lines[15])[2:4])
    nx_, ny_ = map(int, re.findall(r'\d+', lines[24])[:2])
    save_ = int(re.findall(r'\d+', lines[30])[0])
    dest_ = re.findall(r'\S+', lines[30])[1].strip()

    @dataclass
    class Config:
        nitr: int = nitr_
        cycle: int = cycle_
        nx: int = nx_
        ny: int = ny_
        save: int = save_
        dest: str = dest_
    return Config()


def main(opts=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--cycle', '-n', type=int, default=None,
                        help='number of cycles')
    parser.add_argument('--out', '-o', type=str, default=None,
                        help='output directory name')
    parser.add_argument('--resume', '-r', type=str, default=None,
                        help='filename for resume calculation')
    args = parser.parse_args(args=opts)

    conf = read_inputfile('in2d.txt')
    if args.cycle:
        cycle = args.cycle
    else:
        cycle = conf.cycle
    if args.out:
        out = args.out
    else:
        out = conf.dest
    itr_max = conf.nitr
    save_interval = conf.save

    nx = np.array(conf.nx, dtype=np.int32)
    ny = np.array(conf.ny, dtype=np.int32)
    u = np.zeros((2, ny+2, nx+1), dtype=np.float64) # x方向速度
    v = np.zeros((2, ny+1, nx+2), dtype=np.float64) # y方向速度
    p = np.zeros((ny+2, nx+2), dtype=np.float64) # 圧力
    t = np.zeros((ny+2, nx+2), dtype=np.float64) # 温度
    f = np.ones((ny+2, nx+2), dtype=np.float64) # 流れ場情報（0=>物体上, 1=>流体）
    m = np.zeros((120,), dtype=np.uint8) # メッセージ格納用配列
    flg = np.array(0, dtype=np.int32) # 圧力計算収束確認用
    itr_hist = []
    save_count = 0

    register_functions()
    f_read_inputfile()
    f_initialize()

    with open('grid.csv') as fp:
        f[1:-1, 1:-1] = np.array(list(csv.reader(fp)), dtype=np.float64)

    if args.resume:
        u[0], v[0], p[:] = load_value(args.resume)

    for file in ut.iglobm('image/*.png'):
        os.remove(file)

    with ut.stopwatch('calc'):
        with ut.chdir(out):
            with tqdm(total=cycle, mininterval=1) as bar:
                for i in range(1, cycle+1):
                    f_calc_velociry(u[0], v[0], p, t, u[1], v[1], nx, ny)
                    f_bind_velocity(u[1], v[1], f, nx, ny)

                    for j in range(1, itr_max+1):
                        itr = np.array(j, dtype=np.int32)
                        m.fill(ord(' '))
                        f_calc_pressure(u[1], v[1], p, itr, flg, nx, ny, m)
                        f_bind_velocity(u[1], v[1], f, nx, ny)

                        if j % 100 == 0:
                            msg = ''.join(map(chr, m)).rstrip()
                            bar.write(f'cycle={i} {msg}')

                        if flg == 0:
                            if i % 100 == 0:
                                msg = ''.join(map(chr, m)).rstrip()
                                bar.write(f'cycle={i} {msg}')
                            break

                    itr_hist.append(j)
                    u[0] = u[1]
                    v[0] = v[1]

                    if i % save_interval == 0:
                        k = save_count
                        with ut.chdir('result'):
                            dump_data(f'out_{k:05d}.npy', u[0], v[0], p)

                        with ut.chdir('image'):
                            plot_w(vorticity(u[0], v[0]), f'out_{k:05d}.png')
                            plot_nitr(itr_hist, 'nitr.png')
                        save_count += 1

                    bar.update()


if __name__ == '__main__':
    main()
