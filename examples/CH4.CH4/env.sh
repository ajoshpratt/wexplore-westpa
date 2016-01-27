#! /bin/bash -l

module() { eval `/usr/local/Modules/3.2.10/bin/modulecmd bash $*`; }
module use /home/judas/modules
module load westpa/current
module load anaconda
module load gromacs

export WEST_PYTHON=$(which python2.7)
export WM_WORK_MANAGER=serial

