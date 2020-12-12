import math
import os
import re
import subprocess
import time
from zipfile import ZipFile, ZIP_DEFLATED
import numpy as np
import myutils as ut


def comp_main():
    def read(files):
        for file in files:
            name = re.search(r'out_\d+', file)[0]
            print(name, end=' \r')
            yield name, np.load(file)

    for d in ut.fsort(os.listdir('.')):
        loc = f'{d}/result'
        if not os.path.isdir(loc):
            continue

        try:
            with ut.chdir(loc):
                # if any(ut.iglobm('*.npz')):
                #     continue

                if any(ut.iglobm('*.npy')):
                    print(loc)

                for i in range(10):
                    out = f'out_{i:03d}.npz'
                    if os.path.isfile(out):
                        continue

                    l = i * 100
                    u = (i + 1) * 100 if i < 9 else 1001
                    files = [f'out_{j:04d}.plt.npy' for j in range(l, u)]

                    if not all(map(os.path.isfile, files)):
                        continue

                    print(f'{l}=>{u}', ' '*20, end='\r')
                    h = dict(read(files))

                    try:
                        np.savez_compressed(out, **h)

                    except:
                        if os.path.isfile(out):
                            os.remove(out)
                        raise

                    ck = np.load(out)
                    for nm in ck.files:
                        f = f'{nm}.plt.npy'
                        if os.path.isfile(f):
                            os.remove(f)

        except OSError:
            print('OSError')
            time.sleep(10)
            continue


def zip_main():
    for d in ut.fsort(os.listdir('.')):
        loc = f'{d}/image_w'
        if not os.path.isdir(loc):
            continue

        with ut.chdir(loc):
            file = f'image_w.zip'
            if os.path.isfile(file):
                continue

            print(loc)
            # subprocess.run('zip image_w *.png', shell=True)

            try:
                def f_():
                    for png in ut.fsort(ut.iglobm('*.png')):
                        print(png, os.getsize(png), end='\r')
                        with open(png, 'rb'):
                            yield png

                pngs = list(f_())

                with ZipFile(file, 'w', compression=ZIP_DEFLATED) as z:
                    for png in pngs:
                        z.write(png)
            except:
                if os.path.isfile(file):
                    os.remove(file)
                raise

            assert os.path.isfile(file)

            if os.path.isfile(file):
                with ZipFile(file) as z:
                    for png in z.namelist():
                        if os.path.isfile(png):
                            print(png)
                            os.remove(png)


def main():
    comp_main()


if __name__ == '__main__':
    main()
