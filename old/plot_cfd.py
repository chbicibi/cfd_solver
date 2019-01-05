# -*- coding: utf-8 -*-
import glob
import os
import re
import sys
import shutil
import cv2
import numpy as np
from contextlib import contextmanager

@contextmanager
def chdir(path):
  prev_path = os.getcwd()
  os.makedirs(path, exist_ok=True)
  os.chdir(path)
  yield
  os.chdir(prev_path)

def clamp255(c):
  return min(255, max(0, int(255 * c)))

def color_w(x, y, u, v, p, f, w):
  b = 0
  g = clamp255(-10 * w)
  r = clamp255(10 * w)
  return b, g, r

def color_p_log(x, y, u, v, p, f, w):
  plog = np.log10(abs(p)) * 5 if p > 0 else 0
  b = clamp255(plog) if p < 0 else 0
  g = 0
  r = clamp255(plog) if p > 0 else 0
  return b, g, r

def color_p(x, y, u, v, p, f, w):
  b = clamp255(-0.1 * p)
  g = 0
  r = clamp255(0.1 * p)
  return b, g, r

def color_u(x, y, u, v, p, f, w):
  b = clamp255(1 - abs(u))
  g = clamp255(u - 1)
  r = clamp255(-u)
  return b, g, r

def color_uv(x, y, u, v, p, f, w):
  b = 0
  g = clamp255(np.abs(u))
  r = clamp255(np.abs(v))
  return b, g, r

def plot(infile, fn, outfile=None):
  print("plot: " + infile + " -> " + outfile)
  with open(infile, "r") as file:
    # return [float(l.split()[5]) for l in file.readlines()[2:]]
    file.readline()
    line   = file.readline()
    nx, ny = (int(s) for s in re.findall(r"\d+", line))
    image  = np.empty((ny, nx, 3), np.uint8)
    for j in range(ny):
      for i in range(nx):
        line = file.readline()
        x, y, u, v, p, f, w = (float(s) for s in line.split())
        b, g, r = fn(x, y, u, v, p, f, w) if f > 0 else (255, 255, 255)
        image[j, i, 0] = b
        image[j, i, 1] = g
        image[j, i, 2] = r
  if outfile:
    cv2.imwrite(outfile, image)
  return image

def plot_shape(infile, outfile=None):
  print(infile)
  with open(infile, "r") as file:
    lines = file.readlines();
    nx, ny = (len(lines[0].split(",")), len(lines))
    image  = np.empty((ny, nx, 1), np.uint8)
    for j, l in enumerate(lines):
      for i, s in enumerate(l.split(",")):
        image[j, i, 0] = 255 * (1 - int(s))
  if outfile:
    cv2.imwrite(outfile, image)
  return image

def animation(infile, outfile):
  fps = 50
  os.system("..\\ffmpeg.exe -y -r %d -i %s -vcodec libx264 -pix_fmt yuv420p -r %d %s" % (fps, infile, fps, outfile))

def get_dir():
  dirs = [x for x in os.listdir('.') if os.path.isdir(x)]
  print("ディレクトリ番号を指定してください:")
  for i, d in enumerate(dirs):
    print(f"[{i}] {d}")
  no = int(sys.stdin.readline())
  if not no in range(len(dirs)):
    print("該当なし")
    return
  return dirs[no]

def main1(argv):
  outdir = ""
  color_fn = None
  table = {
    "uv": color_uv,
    "p": color_p,
    "u": color_u,
    "w": color_w,
    "plog": color_p_log
  }
  for a in argv:
    m = re.match(r"-m(.+)", a)
    if m:
      color_fn = table[m.group(1)]
      outdir = m.group(1) + "/"
      break
  for a in argv:
    m = re.match(r"-d(.+)", a)
    if m:
      outdir = m.group(1) + "/"
      break
  animfile = outdir + "flow.mp4"
  animkey = outdir + "out_%04d.png"
  dirname = get_dir()
  if not color_fn: # デフォルト
    color_fn = color_w
  if not (dirname and color_fn):
    return
  with chdir(dirname):
    os.makedirs(outdir, exist_ok=True)
    if not "-a" in argv: # 動画作成のみ
      for file in glob.glob("*.plt")[::-1]:
        outfile = outdir + file.replace("plt", "png")
        if not os.path.isfile(outfile) or "-f" in argv: # 全画像生成
          try:
            plot(file, color_fn, outfile)
          except ValueError:
            print("Error")
            pass
    animation(animkey, animfile)
    # sub1()

def main(argv):
  if "-s" in argv:
    plot_shape("grid.csv", "grid.png")
  else:
    main1(argv)

main(sys.argv[1:])
