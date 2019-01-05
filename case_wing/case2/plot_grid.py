#! /usr/bin/env python3

import argparse
import csv
import cv2
import numpy as np


def plot_shape(infile, outfile=None):
    print(infile)
    with open(infile, 'r') as f:
        reader = csv.reader(f)
        a = np.array(list(reader), dtype=np.uint8)
    print(a.shape)
    image = 255 * (1 - a)
    if outfile:
        cv2.imwrite(outfile, image)
    return image


def main():
    plot_shape('grid.csv', 'grid.png')


if __name__ == '__main__':
    main()
