#!/bin/bash
export CFD_HOME=$HOME/yamada/main/cfd

ln -s -f $CFD_HOME/cfd_main.py cfd_main.py
ln -s -f $CFD_HOME/plot_grid.py plot_grid.py
ln -s -f $CFD_HOME/mk_grid.py mk_grid.py
ln -s -f $CFD_HOME/job.sh job.sh
# ln -s -f ../src/ bin
cd bin
ln -s -f $CFD_HOME/src/a.out a.out
cd ../
