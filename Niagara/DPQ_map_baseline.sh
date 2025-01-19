#!/bin/bash

########SBATCH -A ctb-steinman
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=5:59:00
########SBATCH --time=1:00:00
##########SBATCH -p debug
#SBATCH --job-name DPQ_Baseline
#SBATCH --output=DPQ_%j.txt

export MPLCONFIGDIR=/scratch/s/steinman/ahaleyyy/.config/mpl
export PYVISTA_USERDATA_PATH=/scratch/s/steinman/ahaleyyy/.local/share/pyvista
export XDG_RUNTIME_DIR=/scratch/s/steinman/ahaleyyy/.local/temp
export TEMPDIR=$SCRATCH/.local/temp
export TMPDIR=$SCRATCH/.local/temp
export PYVISTA_OFF_SCREEN=true
export PYVISTA_USE_PANEL=true

module load CCEnv StdEnv/2020 gcc/9.3.0 python/3.9.6 #petsc/3.20.0
source $HOME/.virtualenvs/pv_updated/bin/activate 

#python slices_heat.py
#mpirun -np 40 python DPQ.py A notref noskip
#mpirun -np 40 python DPQ.py B notref noskip
mpirun -np 40 python DPQ.py C notref noskip
python mapping.py