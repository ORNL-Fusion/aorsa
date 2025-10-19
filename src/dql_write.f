
      subroutine Dql_write_nc(FILE_txt, FILE_nc)
      
C     This program writes the four quasilinear diffusion coefficients 
C     calculated by AORSA.  The companion program (dql_read.f) 
C     reads the netCDF data file created by this program.
      
      use netcdf
      implicit none
!     include 'netcdf.inc'

C     This is the name of the data file we will create.
      character(*):: FILE_nc
      character(*):: FILE_txt

      integer ncid

C     We are writing 3D data on a uperp, upara, rho  grid.  
c     We will need three netCDF dimensions
      integer NDIMS
      parameter (NDIMS = 3)            
      integer nuper, nupar, nnoderho
      
      character(*):: uperp_NAME, upara_NAME, rhon_NAME      
      parameter (uperp_NAME = 'uperp') 
      parameter (upara_NAME = 'upara')
      parameter (rhon_NAME = 'rho')               
      integer rhon_dimid, uperp_dimid, upara_dimid      

C     In addition to the uperp, upara, rhon dimensions, we will also
C     create uperp, upara, rhon netCDF variables which will hold the
C     actual uperp, upara, rhon. Since they hold data about the
C     coordinate system, the netCDF term for these is: "coordinate
C     variables."      

      real, dimension(:),  allocatable :: uperp
      real, dimension(:),  allocatable :: upara
      real, dimension(:),  allocatable :: rhon  

     
      integer rhon_varid, uperp_varid, upara_varid      

C     We will create four netCDF variables, 
      character*(*) bqlavg_NAME, cqlavg_NAME, eqlavg_name, fqlavg_name
      character*(*) xmi_NAME
      character*(*) vc_cgs_NAME
      character*(*) UminPara_NAME
      character*(*) UmaxPara_NAME      

                       
      parameter (bqlavg_NAME='B_ql')
      parameter (cqlavg_NAME='C_ql')
      parameter (eqlavg_NAME='E_ql')
      parameter (fqlavg_NAME='F_ql')
      parameter (xmi_NAME='mass')
      parameter (vc_cgs_NAME='vnorm') 
      parameter (UminPara_NAME='upara_min')
      parameter (UmaxPara_NAME='upara_max')            
                
      integer bqlavg_varid, cqlavg_varid, eqlavg_varid, fqlavg_varid 
      integer xmi_varid
      integer vc_cgs_varid
      integer UminPara_varid
      integer UmaxPara_varid      
      
      
      integer dimids(NDIMS)                 

C     It's good practice for each variable to carry a "units" attribute.
      character*(*) UNITS
      parameter (UNITS = 'units')
                  
      character*(*) bqlavg_UNITS, cqlavg_UNITS,eqlavg_UNITS,fqlavg_UNITS      
      character*(*) xmi_UNITS
      character*(*) vc_cgs_UNITS
      character*(*) UminPara_UNITS
      character*(*) UmaxPara_UNITS
                              
      character*(*) rhon_UNITS, uperp_UNITS, upara_UNITS                        
      parameter (uperp_UNITS = 'dimensionless')
      parameter (upara_UNITS = 'dimensionless') 
      parameter (rhon_UNITS = 'dimensionless')       
      parameter (bqlavg_UNITS = '(cm/sec)^5') 
      parameter (cqlavg_UNITS = '(cm/sec)^4') 
      parameter (eqlavg_UNITS = '(cm/sec)^4') 
      parameter (fqlavg_UNITS = '(cm/sec)^3')             
      parameter (xmi_UNITS = 'kg')
      parameter (vc_cgs_UNITS = 'cm/sec')
      parameter (UminPara_UNITS = 'dimensionless')
      parameter (UmaxPara_UNITS = 'dimensionless')      
                         

C     We will create some quasilinear data to write to FILE_nc.      

      double precision, dimension(:,:,:), allocatable :: bqlavg
      double precision, dimension(:,:,:), allocatable :: cqlavg 
      double precision, dimension(:,:,:), allocatable :: eqlavg
      double precision, dimension(:,:,:), allocatable :: fqlavg  
      
      double precision xmi 
      double precision vc_cgs  
      double precision UminPara 
      double precision UmaxPara                  

C     Loop indices.
      integer i_uperp, i_upara, n

C     Error handling.
      integer retval
                                 
C     Create test data by reading quasilinear diffusion coefficients 
c     from the text file FILE_txt

      open(unit=42, file = FILE_txt, status='unknown', form='formatted')
      read (42, 309) nuper
      read (42, 309) nupar
      read (42, 309) nnoderho
      
      allocate( uperp(nuper) )
      allocate( upara(nupar) )
      allocate( rhon(nnoderho) )       
      allocate( bqlavg(nuper, nupar, nnoderho) ) 
      allocate( cqlavg(nuper, nupar, nnoderho) ) 
      allocate( eqlavg(nuper, nupar, nnoderho) ) 
      allocate( fqlavg(nuper, nupar, nnoderho) )   

      read (42, 3310) vc_cgs
      read (42, 3310) UminPara, UmaxPara
      read (42, 3310) (rhon(n), n = 1, nnoderho)
      read (42, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
      read (42, 3310) (upara(i_upara), i_upara = 1, nupar)
      read (42, 3310) (((bqlavg(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar),
     &        n = 1, nnoderho)
      read (42, 3310) (((cqlavg(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), 
     &        n = 1, nnoderho)
      read (42, 3310) (((eqlavg(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), 
     &        n = 1, nnoderho)
      read (42, 3310) (((fqlavg(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), 
     &        n = 1, nnoderho)     
      read (42, 3310) xmi
      
      close (42)
                 
       
       write (6, *) 'nuper = ', nuper
       write (6, *) 'nupar = ', nupar
       write (6, *) 'nnoderho = ', nnoderho                          

       write (6, *) "bqlavg(midpt) = ",
     &       bqlavg(int(nuper/2), int(nupar/2), int(nnoderho/2))
       write(6,*) 'xmi =', xmi                 
       write(6,*) 'vc_cgs =', vc_cgs
       write(6,*) 'UminPara =', UminPara 
       write(6,*) 'UmaxPara =', UmaxPara               
       
             

C     Create the netcdf file. 
      retval = nf90_create(FILE_nc, nf90_clobber, ncid)
      if (retval .ne. nf90_noerr) call handle_err(retval)

C     Define the dimensions.       
      retval = nf90_def_dim(ncid, uperp_NAME, nuper, uperp_dimid)
      if (retval .ne. nf90_noerr) call handle_err(retval)
      retval = nf90_def_dim(ncid, upara_NAME, nupar, upara_dimid)
      if (retval .ne. nf90_noerr) call handle_err(retval)
      retval = nf90_def_dim(ncid, rhon_NAME, nnoderho, rhon_dimid)
      if (retval .ne. nf90_noerr) call handle_err(retval)
      
      write(6,*) "uperp_dimid = ", uperp_dimid
      write(6,*) "upara_dimid = ", upara_dimid
      write(6,*) "rhon_dimid = ", rhon_dimid            


C     Define the coordinate variables. They will hold the coordinate
C     information, that is, the latitudes and longitudes. A varid is
C     returned for each.                              
      retval = nf90_def_var(ncid, uperp_NAME, NF90_REAL, uperp_dimid, 
     &     uperp_varid)     
      if (retval .ne. nf90_noerr) call handle_err(retval)      
      retval = nf90_def_var(ncid, upara_NAME, NF90_REAL, upara_dimid, 
     &     upara_varid)     
      if (retval .ne. nf90_noerr) call handle_err(retval) 
      retval = nf90_def_var(ncid, rhon_NAME, NF90_REAL, rhon_dimid, 
     +     rhon_varid)          
      if (retval .ne. nf90_noerr) call handle_err(retval) 
      
      write(6,*) uperp_varid, upara_varid, rhon_varid
           

C     Assign units attributes to coordinate var data. This attaches a
C     text attribute to each of the coordinate variables, containing the
C     units.         
      retval = nf90_put_att(ncid, uperp_varid, UNITS, 
     .     uperp_UNITS)          
      if (retval .ne. nf90_noerr) call handle_err(retval) 
                 
      retval = nf90_put_att(ncid, upara_varid, UNITS, 
     .    upara_UNITS)
      if (retval .ne. nf90_noerr) call handle_err(retval) 
              
      retval = nf90_put_att(ncid, rhon_varid, UNITS, 
     .    rhon_UNITS)
      if (retval .ne. nf90_noerr) call handle_err(retval) 
      
      write(6, *)'uperp_units = ', uperp_units
      write(6, *)'upara_units = ', upara_units
      write(6, *)'rhon_units = ',  rhon_units
                                                   

C     Define the netCDF variables. The dimids array is used to pass the
C     dimids of the dimensions of the netCDF variables.
      dimids(1) = uperp_dimid
      dimids(2) = upara_dimid
      dimids(3) = rhon_dimid      

C     Define the netCDF variables for the D_quasi-linear data.            
      retval=nf90_def_var(ncid, bqlavg_NAME, nf90_double, dimids, 
     +     bqlavg_varid)
      if (retval .ne. nf90_noerr) call handle_err(retval)             
      retval=nf90_def_var(ncid, cqlavg_NAME, nf90_double, dimids, 
     +     cqlavg_varid)
      if (retval .ne. nf90_noerr) call handle_err(retval)            
      retval=nf90_def_var(ncid, eqlavg_NAME, nf90_double, dimids, 
     +     eqlavg_varid)
      if (retval .ne. nf90_noerr) call handle_err(retval)           
      retval=nf90_def_var(ncid, fqlavg_NAME, nf90_double, dimids, 
     +     fqlavg_varid)
      if (retval .ne. nf90_noerr) call handle_err(retval)              
      retval = nf90_def_var(ncid, xmi_Name, nf90_double, 
     +     xmi_varid)      
      if (retval .ne. nf90_noerr) call handle_err(retval)  
      retval = nf90_def_var(ncid, vc_cgs_Name, nf90_double,
     +     vc_cgs_varid)      
      if (retval .ne. nf90_noerr) call handle_err(retval)  
      
      retval = nf90_def_var(ncid, UminPara_Name, nf90_double,
     +     UminPara_varid)      
      if (retval .ne. nf90_noerr) call handle_err(retval)  
      retval = nf90_def_var(ncid, UmaxPara_Name, nf90_double,
     +     UmaxPara_varid)      
      if (retval .ne. nf90_noerr) call handle_err(retval)                        
                  
        
C     Assign units attributes to the netCDF variables.      
      
      retval = nf90_put_att(ncid, bqlavg_varid, UNITS, 
     .   bqlavg_UNITS)
      if (retval .ne. nf90_noerr) call handle_err(retval)         
      retval = nf90_put_att(ncid, cqlavg_varid, UNITS, 
     .   cqlavg_UNITS)
      if (retval .ne. nf90_noerr) call handle_err(retval)        
      retval = nf90_put_att(ncid, eqlavg_varid, UNITS, 
     .   eqlavg_UNITS)
      if (retval .ne. nf90_noerr) call handle_err(retval)      
      retval = nf90_put_att(ncid, fqlavg_varid, UNITS, 
     .   fqlavg_UNITS)
      if (retval .ne. nf90_noerr) call handle_err(retval) 
      
      retval = nf90_put_att(ncid, xmi_varid, UNITS, 
     .   xmi_UNITS)
      if (retval .ne. nf90_noerr) call handle_err(retval)   
      retval = nf90_put_att(ncid, vc_cgs_varid, UNITS, 
     .   vc_cgs_UNITS)
      if (retval .ne. nf90_noerr) call handle_err(retval)  
      
      retval = nf90_put_att(ncid, UminPara_varid, UNITS, 
     .   UminPara_UNITS)
      if (retval .ne. nf90_noerr) call handle_err(retval)   
      retval = nf90_put_att(ncid, UmaxPara_varid, UNITS, 
     .   UmaxPara_UNITS)
      if (retval .ne. nf90_noerr) call handle_err(retval)                   
      
      write(6, *)'bqlavg_units = ', bqlavg_units
      write(6, *)'cqlavg_units = ', cqlavg_units
      write(6, *)'eqlavg_units = ', eqlavg_units
      write(6, *)'fqlavg_units = ', fqlavg_units
      
      write(6, *)'xmi_units = ', xmi_units 
      write(6, *)'vc_cgs_units = ', vc_cgs_units 
      
      write(6, *)'UminPara_units = ', UminPara_units 
      write(6, *)'UmaxPara_units = ', UmaxPara_units                               
       
                             

C     End define mode.
      retval = nf90_enddef(ncid)
      if (retval .ne. nf90_noerr) call handle_err(retval)
            
C     Write the coordinate variable data. This will put uperp, upara,
c     and rhon's of our data grid into the netCDF file.      
      retval = nf90_put_var(ncid, uperp_varid, uperp)
      if (retval .ne. nf90_noerr) call handle_err(retval)
      retval = nf90_put_var(ncid, upara_varid, upara)
      if (retval .ne. nf90_noerr) call handle_err(retval) 
      retval = nf90_put_var(ncid, rhon_varid, rhon)
      if (retval .ne. nf90_noerr) call handle_err(retval) 
      
      write(6,*)uperp_varid, upara_varid, rhon_varid        
           
      write (6, 3310) (rhon(n), n = 1, nnoderho)
      write (6, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
      write (6, 3310) (upara(i_upara), i_upara = 1, nupar)  

C     Write the data. This will write our D_ql data. The arrays only hold one timestep worth
C     of data. We will just rewrite the same data for each timestep. In
C     a real application, the data would change between timesteps.
      
      retval = nf90_put_var(ncid, bqlavg_varid, bqlavg)
      if (retval .ne. nf90_noerr) call handle_err(retval)  
      retval = nf90_put_var(ncid, cqlavg_varid, cqlavg)
      if (retval .ne. nf90_noerr) call handle_err(retval)  
      retval = nf90_put_var(ncid, eqlavg_varid, eqlavg)
      if (retval .ne. nf90_noerr) call handle_err(retval)          
      retval = nf90_put_var(ncid, fqlavg_varid,  fqlavg)
      if (retval .ne. nf90_noerr) call handle_err(retval)       
      retval = nf90_put_var(ncid, xmi_varid, xmi)
      if (retval .ne. nf90_noerr) call handle_err(retval) 
      retval = nf90_put_var(ncid, vc_cgs_varid, vc_cgs)
      if (retval .ne. nf90_noerr) call handle_err(retval)  
      
      retval = nf90_put_var(ncid, UminPara_varid, UminPara)
      if (retval .ne. nf90_noerr) call handle_err(retval) 
      retval = nf90_put_var(ncid, UmaxPara_varid, UmaxPara)
      if (retval .ne. nf90_noerr) call handle_err(retval)                    
           

C     Close the file. 
      retval = nf90_close(ncid)
      if (retval .ne. nf90_noerr) call handle_err(retval)
   
C     If we got this far, everything worked as expected. Yipee!
      print *,'*** SUCCESS writing netcdf file, ', FILE_nc, '!'
      
      
      deallocate( uperp )
      deallocate( upara )
      deallocate( rhon ) 
      deallocate( bqlavg) 
      deallocate( cqlavg) 
      deallocate( eqlavg) 
      deallocate( fqlavg)         
      
 3310 format(1p,6e18.10)
  309 format(10i10)
  
      return      
      end

      subroutine handle_err(errcode)
      use netcdf
      implicit none
      !include 'netcdf.inc'
      integer errcode

      print *, 'Error: ', nf90_strerror(errcode)
      stop 2
      end
