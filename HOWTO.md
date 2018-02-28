
Inputs:
  - ZTABLE.TXT - numerical lookup table for Plasma Dispersion function 
  - grfont.dat - needed for graphics in both AORSA and CQL3D
  - Ztable - numerical lookup table for Plasma Dispersion function
  - geqdsk - EFIT G-EQDSK file for tokamak equilibrium, name is given in aorsa2d.in
  - aorsa2d.in - fortran namelist file for code parameters
    + nprow,npcol - processor grid. Typically equal. nprow x npcol = nproc requested for batch job

Outputs:
  - out15 - run time logging messages
  - log_aorsa2d - run time logging messages
  - vtk files (viewed with Visit or paraview or mayavi). 2D and 3D plots
    - Bql_avg_2D.vtk, Cql_avg_2D.vtk - Quasilinear diffusion coefficients, optional with flag
    - E_kicks_2D.vtk - for monte carlo kicks. optional with flag, 
    - Eb_spectrum.vtk - RF electric field Poynting_2D.vtk
    - Bfield_2D.vtk - magnetic equilibrium or RF?
    - Efield_2D.vtk  - electric field from RF 
    - capd.vtk - contour plots of dispersion relation D(k,x)
  - Postscript output (from pgplot)
    - aorsa2d.ps - main output. 1D and 2D field and power plots
      + slide p80+ - 2D contour plots of alpha, beta, and b components of rf electric fields and currents
      + slide p92+ - 1D plots of Eplus,Eminus,Eparallel. 2D contours of J.E and Wdot (power deposition)
      + slide p111+ - 1D and 2D plots of power spectra; useful for testing convergence
      + slide p123+ - 2D plots of absolute value of components of E
      + slide p130+ - 2D Btor and Bpol equilibria and dl/B
      + slide p133+ - power from quasilinear, J.E, and Wdot comparisons
    - movie.ps - frames of 2D electric field as a time harmonic movie to show phase velocity
    - disper.ps - dielectric properties of the wave. kperp for branches, location of resonances, etc.
    - eqdsk_setup.ps - show equilibrium properties from EFIT G-EQDSK file
  
