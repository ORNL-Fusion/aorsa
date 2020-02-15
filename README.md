# Dependencies

* netcdf
* scalapack
* mpi
* pgplot

# Build

## Cori
```
git clone https://github.com/ORNL-Fusion/aorsa.git
cd aorsa
source env.cori
module unload darshan
module load cray-netcdf
module load dfftpack
mkdir -p obj
make -f makefile_aorsa_v.24_simple_load.cori
```

## fusiont6
```
git clone https://github.com/ORNL-Fusion/aorsa.git
cd aorsa
mkdir -p obj
make -f makefile.fusiont6
```

## osx-mojave
```
brew install open-mpi
brew install scalapack
brew install pgplot
git clone https://github.com/ORNL-Fusion/aorsa.git
git checkout osx-mojave
mkdir obj
make -f makefile.osx-mojave

```

# Run

## Cori
```
source env.cori
ulimit -s unlimited
cp -r examples $SCRATCH/
cd $SCRATCH/examples
cd DIIID_SPONG_ICE
```
### batchscript
```
sbatch batchscript.cori
```
### interative
```
salloc -N 1 -C haswell -q interactive -t 01:00:00
srun -n 1 /path/to/xaorsa2d
```
## fusiont6
```
ulimit -s unlimited
cd examples/DIIID_SPONG_ICE
mpirun -n 1 ../../xaorsa2d
```

# Other Notes
This directory contains the stable release of AORSA. 

src/ source for building AORSA

src/FFTPACK , src/CQL3DSETUP : source for additional AORSA dependencies

src/JAGERHP : special version of AORSA used for lower hybrid and testing odd order derivative

src/SAVE* , src/NEW , src/BACKUP : other versions of AORSA


