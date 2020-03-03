[![Actions Status](https://github.com/ORNL-Fusion/aorsa/workflows/CI%20on%20Mac/badge.svg)](https://github.com/ORNL-Fusion/aorsa/actions)
[![Actions Status](https://github.com/ORNL-Fusion/aorsa/workflows/CI%20on%20Ubuntu/badge.svg)](https://github.com/ORNL-Fusion/aorsa/actions)

# Build & Run
## Dependencies

* netcdf / netcdff
* scalapack / blas / blacs
* mpi
* pgplot

## Cori
### Build
```
git clone https://github.com/ORNL-Fusion/aorsa.git
cd aorsa
source env.cori
make clean
make
```
### Run
```
source env.cori
ulimit -s unlimited
cp -r examples $SCRATCH/
cd $SCRATCH/examples
cd DIIID-helicon
```
#### batchscript
```
sbatch batchscript.cori
```
#### interative
```
salloc -N 1 -C haswell -q interactive -t 01:00:00
srun -n 1 /path/to/xaorsa2d
```

## Ubuntu 18.04
### Build
```
apt install gfortran libopenmpi-dev pgplot5 libnetcdff-dev libscalapack-openmpi-dev libpng-dev libopenblas-dev libx11-dev
```
### Run
```
cd examples/DIIID-helicon
mpirun -n 1 ../../xaorsa2d
```

## fusiont6
### Build
```
git clone https://github.com/ORNL-Fusion/aorsa.git
cd aorsa
make 
```
### Run
```
ulimit -s unlimited
cd examples/DIIID-helicon
mpirun -n 1 ../../xaorsa2d
```

## osx-mojave
### Build
```
brew install open-mpi scalapack pgplot
git clone https://github.com/ORNL-Fusion/aorsa.git
cd aorsa
make 
```

# Inputs
  - `ZTABLE.TXT` - numerical lookup table for Plasma Dispersion function 
  - `grfont.dat` - needed for graphics in both AORSA and CQL3D
  - `Ztable` - numerical lookup table for Plasma Dispersion function
  - `geqdsk` - EFIT G-EQDSK file for tokamak equilibrium, name is given in aorsa2d.in
  - `aorsa2d.in` - fortran namelist file for code parameters
    + processors: nprow,npcol - processor grid. Typically equal. nprow x npcol = nproc requested for batch job
    + profiles: iprofile=3 doubly parabolic, iprofile=5 numerical in STATE block
      - doubly parabolic options

# Outputs
  - `out15` - run time logging messages
  - `log_aorsa2d` - run time logging messages
  - vtk files (viewed with Visit or paraview or mayavi). 2D and 3D plots
    - `Bql_avg_2D.vtk`, Cql_avg_2D.vtk - Quasilinear diffusion coefficients, optional with flag
    - `E_kicks_2D.vtk` - for monte carlo kicks. optional with flag, 
    - `Eb_spectrum.vtk` - RF electric field Poynting_2D.vtk
    - `Bfield_2D.vtk` - magnetic equilibrium or RF?
    - `Efield_2D.vtk`  - electric field from RF 
    - `capd.vtk` - contour plots of dispersion relation D(k,x)
  - Postscript output (from pgplot)
    - `aorsa2d.ps` - main output. 1D and 2D field and power plots
    
      Colors for curves for species specific profiles are electrons: red, ions: cyan (majority), blue, green, magenta, orange, yellow 
      in order of species index from aorsa2d.in
      + slide p1+  - contour plots of the Ztable
      + slide p15+ - density and temperature profiles
      + slide p80+ - 2D contour plots of alpha, beta, and b components of rf electric fields and currents
      + slide p92+ - 1D plots of Eplus,Eminus,Eparallel. 2D contours of J.E and Wdot (power deposition)
      + slide p111+ - 1D and 2D plots of power spectra; useful for testing convergence
      + slide p123+ - 2D plots of absolute value of components of E
      + slide p130+ - 2D Btor and Bpol equilibria and dl/B
      + slide p133+ - power from quasilinear, J.E, and Wdot comparisons
    - `movie.ps` - frames of 2D electric field as a time harmonic movie to show phase velocity
    - `disper.ps` - dielectric properties of the wave. kperp for branches, location of resonances, etc.
    - `eqdsk_setup.ps` - show equilibrium properties from EFIT G-EQDSK file

# Other Notes
This directory contains the stable release of AORSA. 

src/ source for building AORSA

src/FFTPACK , src/CQL3DSETUP : source for additional AORSA dependencies

src/JAGERHP : special version of AORSA used for lower hybrid and testing odd order derivative

src/SAVE* , src/NEW , src/BACKUP : other versions of AORSA


