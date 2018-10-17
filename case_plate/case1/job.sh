#!/bin/bash
#PJM -L "rscunit=ito-a"
#PJM -L "rscgrp=ito-ss"
#PJM -L "vnode=1"
#PJM -L "vnode-core=36"
#PJM -L "elapse=10:00:00"
#PJM --mail-list "yunibo0614g@gmail.com"
#PJM -m b
#PJM -m e
#PJM -j
#PJM -X

module load python/3.6.2
module load cuda/8.0

# ulimit -v 33554432
# ulimit -u 512
# ulimit -a > ulimit.log
export PYTHONUSERBASE=~/yamada/local
export PYTHONPATH=~/yamada/main/lib:$PYTHONPATH
export OMP_NUM_THREADS=36

python3 cfd_main.py -run
