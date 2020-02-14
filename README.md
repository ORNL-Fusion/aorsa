This directory contains the stable release of AORSA. 

src/ source for building AORSA

src/FFTPACK , src/CQL3DSETUP : source for additional AORSA dependencies

src/JAGERHP : special version of AORSA used for lower hybrid and testing odd order derivative

src/SAVE* , src/NEW , src/BACKUP : other versions of AORSA

# Cori
## Build
source env.cori
module unload darshan
module load cray-netcdf
module load dfftpack
mkdir -p obj/cori
make -f makefile_aorsa_v.24_simple_load.cori
## Run
source env.cori
ulimit -s unlimited
cp -r examples $SCRATCH/
cd $SCRATCH/examples
cd DIIID_SPONG_ICE
sbatch cori.batchscript
## Run interative
salloc -N 1 -C haswell -q interactive -t 01:00:00
srun -n 1 /path/to/xaorsa2d.cori 
