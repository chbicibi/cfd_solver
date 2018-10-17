#! /usr/bin/env python3

import argparse
import csv
import math
import re
import numpy as np


################################################################################

def rotate(x, y, cx, cy, a):
    x_ = (x - cx) * math.cos(a) - (y - cy) * math.sin(a) + cx
    y_ = (x - cx) * math.sin(a) + (y - cy) * math.cos(a) + cy
    return x_, y_


################################################################################

def grid_plain(shape):
    grid = np.ones(shape, dtype=np.uint8)
    return grid


def grid_plate(shape):
    ny, nx = shape
    grid = np.ones(shape, dtype=np.uint8)
    grid[ny//2-2:ny//2+2, nx//8:nx//8+ny] = 0
    return grid


def grid_plate_a(shape, alpha):
    ny, nx = shape
    grid = np.ones(shape, dtype=np.uint8)

    left = nx // 8
    length = nx // 4
    thickness = 4

    a_rad = math.radians(alpha)
    cx = left + length // 2
    cy = ny // 2
    base_x = left, length
    base_y = cy - thickness // 2, thickness

    for i in range(1000): # x
        for j in range(1000): # y
            x = base_x[0] + base_x[1] * i / 999
            y = base_y[0] + base_y[1] * j / 999
            r = rotate(x, y, cx, cy, a_rad)
            x_, y_ = map(int, r)
            grid[y_, x_] = 0
    return grid


################################################################################

def get_info(input_file):
    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    nx, ny = [int(s) for s in re.findall(r'\d+', lines[24])]
    dest = re.findall(r'\S+', lines[30])[1].strip()
    info = {
        'nx': nx,
        'ny': ny,
        'dest': dest
    }
    return info


def dump_grid(grid, file):
    with open(file, 'w', encoding='utf-8') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerows(grid)


################################################################################

def make_grid(args, info):
    a = args.a
    print('alpha=', a, 'deg')
    shape = [info[k] for k in ('ny', 'nx')]
    grid = grid_plate_a(shape, a)
    return grid


################################################################################

def get_args():
    parser = argparse.ArgumentParser()
    # parser.add_argument('live_id', help='live id', nargs='*')
    # parser.add_argument('-f', help='get live ids from file', action='store_true')
    # parser.add_argument('-w', help='watch live page', action='store_true')
    # parser.add_argument('-c', help='getting comments', action='store_true')
    # parser.add_argument('-dl', help='call main', action='store_true')
    # parser.add_argument('-ck', help='check filename', action='store_true')
    # parser.add_argument('-op', help='open live page', action='store_true')
    # parser.add_argument('-debug', help='enable debug mode', action='store_true')
    # parser.add_argument('-concat', help='concatenate video files', action='store_true')

    parser.add_argument('-a', default=0, type=int, help='alpha')
    parser.add_argument('-test', action='store_true', help='test mode')
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    info = get_info('in2d.txt')
    grid = make_grid(args, info)
    dump_grid(grid, 'grid.csv')


if __name__ == '__main__':
    main()