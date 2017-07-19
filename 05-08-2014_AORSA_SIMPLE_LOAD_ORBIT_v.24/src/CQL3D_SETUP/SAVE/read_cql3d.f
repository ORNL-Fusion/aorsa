       module read_CQL3D

       USE CQL_kinds_m

       implicit none

       integer :: n_theta_max  ! max number of pitch angles
       integer :: n_u          ! Number of normalized velocities
       integer :: n_psi        ! Number of flux surfaces
       integer :: n_t          ! Number of cql3d time steps 

       integer, allocatable :: n_theta_(:) ! number of pitch angles for flux surface
       real(kind=real_kind), allocatable :: u(:) ! normalized velocity
       real(kind=real_kind), allocatable ::  theta(:,:) ! pitch angles theta(i_theta, i_psi)
       real(kind=real_kind), allocatable :: rho_a(:)    !normalized small radius
       real(kind=real_kind), allocatable :: f_CQL(:,:,:) ! f(i_theta, i_u, i_psi)
       real(kind=real_kind), allocatable :: f_cql_2d(:,:) ! f(i_u, i_theta)
       
       real(kind=real_kind), allocatable :: wperp_cql(:,:) ! wperp(i_psi, n_t) 
       real(kind=real_kind), allocatable :: wpar_cql(:,:)  ! wpar(i_psi, n_t)          


       real(kind=real_kind) :: vc



       CONTAINS

       subroutine netcdfr3d(netcdfnm)

!c-----Reads a
!c distribution and the mesh from a netCDF file created
!c     with variables in the primary CQL3D output netcdf file.

      implicit none

      INTEGER :: istat

!      implicit integer (i-n), real*8 (a-h,o-z)
!      save

       character*(*) netcdfnm ! input filename

!c --- include file for netCDF declarations
!c --- (obtained from NetCDF distribution)
       include 'netcdf.inc'

!c-----Need to set the dimensions through a parameter statement
!c       [dynamic dimensioning could get around this limitation].
!c       The user can use ncdump [associated with the netCDF distributn]
!c       to determine parameters and definitions of data in the
!c       netCDF file.
!c     parameters (iya,jxa,lrza,ngena) set with:
!c      include 'param.i'

!c-----output
!c      [Much more data is generally available in the input
!c       netcdf file from CQL3D.  ncdump can be used to view
!c       available data, dimensions, and short descriptions.
!c       Some additional data can be read below by removing
!c       line comments.  Other data can be obtained by
!c       emulation of the coding below.]

!c original data specifications.  Will make arrays allocatable below
!c      real*8 vnorm
!c      integer iy,jx,lrz
!c      integer iy_(lrza)
!c      real*8 dx(jxa)                   !0.5(x(j+1)-x(j-1))
!c      real*8 cint2(jxa)                !x**2dx
!c      real*8 x(jxa)
!c      real*8 dy(iya,lrza)              !0.5*(y(i+1,l)-y(i-1,l))
!c      real*8 cynt2(iya,lrza)           !2*pi*sin(y)*dy
!c      real*8 y(iya,lrza)
!c      real*8 rya(lrza)                  !normalized small radius
!c      real*8 f(iya,jxa,lrza)

       real*8 vnorm
       integer iy,jx,lrz, nt, nt_id
       integer, allocatable :: iy_(:)
       real*8, allocatable :: x(:)
       real*8, allocatable :: y(:,:)
       real*8, allocatable :: rya(:)                  !normalized small radius
       real*8, allocatable :: f(:,:,:)
       real*8, allocatable, dimension(:,:) :: wperp   !perp energy/particle
                                                      !tdim, rdim
       real*8, allocatable, dimension(:,:) :: wpar    !par energy/particle
                                                      !tdim, rdim
       

!c     y is pitch angle (radians).  It varies from one flux surface to
!c       the next, for example, due to increased resolution near the
!c       trapped-passing boundary.
!c     In general, the number of pitch angle points iy_() can vary
!c       with flux surface.  In this subroutine, it is presently
!c       assumed that the number of pitch angle points, iy,
!c       is constant as a function of flux
!c       surface.
!c     x(1:jx) is normalized momentum-per-mass from 0. to the
!c       maximum momentum-per-mass vnorm[cgs units].  The mesh
!c       does not vary with radius.
!c     rya(1:lrz) is the normalized radial flux surface grid.
!c       lrz is .le.lrza.
!c     The distribution function f(1:iy,1:jx,1:lrz) is normalized
!c       such that the integral[f x**2*dx*2*pi*sin(y)*dy] over
!c       momentum space gives species density per cm**3.
!c


!c --- some stuff for netCDF file ---
       character*128 name
       integer ncid,istatus
       integer xdim,ydim,rdim,kdim,vid
       integer ngen,ntotal
       integer ll,j,i
       integer start(3),count(3),start_y(2),count_y(2)

       data start/1,1,1/,start_y/1,1/

!.......................................................................
!     Open cql3d netcdf file
!.......................................................................

       write(*,*)'before ncopn netcdfnm=',netcdfnm
       ncid = ncopn(TRIM(netcdfnm),NCNOWRIT,istatus)
       write(*,*)'after ncopn ncid=',ncid,'istatus',istatus
       
!c.......................................................................
!c     read in dimension IDs and sizes

       write(*,*)'before ncdid xdim'
       xdim = ncdid(ncid,'xdim',istatus)
       write(*,*)'after ncdid xdim=',xdim,'istatus',istatus
       ydim = ncdid(ncid,'ydim',istatus)
       write(*,*)'after ncdid ydim=',ydim,'istatus',istatus
       rdim = ncdid(ncid,'rdim',istatus)
       write(*,*)'after ncdid rdim=',rdim,'istatus',istatus
       kdim=ncdid(ncid,'species_dim',istatus)
       write(*,*)'afterncid kim=',kdim,'istatus',istatus
       
       istatus = nf_inq_dimid(ncid,'tdim',nt_id)
       write(*,*)'proc_cql3d_op: after ncdid nt_id = ',nt_id,'istatus = ',istatus 
       

       istatus = nf_inq_dimlen(ncid, nt_id, nt)
       !call ncdinq(ncid, nt_id,'tdim', nt, istatus)
       write(*,*)'proc_cql3d_op: after ncdinq, # of t steps = ',nt, ' istatus=',istatus      

!c --- inquire about dimension sizes ---
!c     ncdinq(netCDF_id, dimension_id_from_ncdid, returned_dim_name,
!c     returned_dim_size)
!c     Note: for unlimited dimension, returned_dim_size=current maximum
!c     which is the same as the maximum record number

       call ncdinq(ncid,ydim,name,iy,istatus)
       call ncdinq(ncid,xdim,name,jx,istatus)
       call ncdinq(ncid,rdim,name,lrz,istatus)
       call ncdinq(ncid,kdim,name,ntotal,istatus)
       write(*,*)'iy,jx,lrz,ntotal',iy,jx,lrz,ntotal
!c      write(*,*)'netcdfr3d 2 iya,jxa,lrza',iya,jxa,lrza

!c********** Skip dimension checking since are allocating arrays DBB 8/27/03
!c-----check the dimensions

!c      if (lrz.ne.lrza) then
!c      if (lrz.gt.lrza) then
!c         write(*,*)'netcdfrw2 lrz.ne.lrza'
!c         write(*,*)'netcdfr3d lrz.gt.lrza'
!c         write(*,*)'lrz from netcdfnm.nc file =',lrz
!c         write(*,*)'lrza from param.i file =',lrza
!c         write(*,*)'Attention!!! It should be lrz=lrza'
!c         write(*,*)'Please change lrza in param.i and recomplile codes'
!c         stop
!c      endif

!c      write(*,*)'netcdfr3d 3 iya,jxa,lrza',iya,jxa,lrza
!c      write(*,*)'netcdfr3d 3 iy,jx,lrz',iy,jx,lrz

!c      if (jx.ne.jxa) then
!c       if (jx.gt.jxa) then
!c         write(*,*)'netcdfr3d jx.ne.jxa'
!c         write(*,*)'netcdfr3d jx.gt.jxa'
!c         write(*,*)'jx from netcdfnm.nc file =',jx
!c         write(*,*)'jxa from param.i file =',jxa
!c         write(*,*)'Please change jxa in param.i'
!c         stop
!c      endif

!c      if (iy.ne.iya) then
!c      if (iy.gt.iya) then
!c         write(*,*)'netcdfr3d iy.ne.iya'
!c         write(*,*)'netcdfr3d iy.gt.iya'
!c         write(*,*)'iy from netcdfnm.nc file =',iy
!c         write(*,*)'iya from param.i file =',iya
!c         write(*,*)'Please change iya in param.i'
!c         stop
!c      endif
!c ***************end dimension checking ***************************

!c ************* Allocate arrays

!c ************* Allocate space for Harvey arrays.

	ALLOCATE( iy_(lrz), stat=istat )
       IF (istat /= 0 ) THEN
       WRITE (*,'("read_CQL3D: allocate failed for iy_")')
       PAUSE
       END IF

	ALLOCATE( y(iy, lrz), stat=istat )
       IF (istat /= 0 ) THEN
       WRITE (*,'("read_CQL3D: allocate failed for y")')
       PAUSE
       END IF

	ALLOCATE( x(jx), stat=istat )
       IF (istat /= 0 ) THEN
       WRITE (*,'("read_CQL3D: allocate failed for x")')
       PAUSE
       END IF

	ALLOCATE( rya(lrz), stat=istat )
       IF (istat /= 0 ) THEN
       WRITE (*,'("read_CQL3D: allocate failed for rya")')
       PAUSE
       END IF

	ALLOCATE( f(iy, jx, lrz), stat=istat )
       IF (istat /= 0 ) THEN
       WRITE (*,'("read_CQL3D: allocate failed for f")')
       PAUSE
       END IF
       
	ALLOCATE( wperp(lrz, nt), stat=istat )
       IF (istat /= 0 ) THEN
       WRITE (*,'("read_CQL3D: allocate failed for wperp")')
       PAUSE
       END IF 
       
	ALLOCATE( wpar(lrz, nt), stat=istat )
       IF (istat /= 0 ) THEN
       WRITE (*,'("read_CQL3D: allocate failed for wpar")')
       PAUSE
       END IF               
       
       

       count(1)=iy
       count(2)=jx
       count(3)=lrz

!c-----normalized momentum x (momentum/mass/vnorm) variables

!c     vnorm - character velocity (momentum-per-mass)[cms/sec]
       vid = ncvid(ncid,'vnorm',istatus)
       call ncvgt(ncid,vid,1,1,vnorm,istatus)
       write(*,*)'after ncvgp vnorm=',vnorm

       vid = ncvid(ncid,'x',istatus)
       call ncvgt(ncid,vid,1,jx,x,istatus)
!c      write(*,*)'x',(x(j),j=1,jx)

!c      vid = ncvid(ncid,'dx',istatus)
!c      call ncvgt(ncid,vid,1,jx,dx,istatus)
!c      write(*,*)'dx',(dx(j),j=1,jx)

!c      vid = ncvid(ncid,'cint2',istatus)
!c      call ncvgt(ncid,vid,1,jx,cint2,istatus)
!c      write(*,*)'cint2',(cint2(j),j=1,jx)

!c-----pitch angle variavles y

       vid = ncvid(ncid,'iy_',istatus)
       call ncvgt(ncid,vid,1,lrz,iy_,istatus)
!c      write(*,*)'iy_',(iy_(ll),ll=1,lrz)

       count_y(1)=iy
       count_y(2)=lrz
       vid = ncvid(ncid,'y',istatus)
       call ncvgt(ncid,vid,start,count_y,y,istatus)
!c      do ll=1,lrza
!c         write(*,*)'ll=',ll,'iy_(ll)=',iy_(ll)
!c         write(*,*)'y',(y(i,ll),i=1,iy_(ll))
!c      enddo

!c      vid = ncvid(ncid,'dy',istatus)
!c      call ncvgt(ncid,vid,start,count_y,dy,istatus)
!c      do ll=1,lrza
!c         write(*,*)'ll=',ll,'iy_(ll)=',iy_(ll)
!c         write(*,*)'dy',(dy(i,ll),i=1,iy_(ll))
!c      enddo

!c      vid = ncvid(ncid,'dy',istatus)
!c      call ncvgt(ncid,vid,start,count_y,cynt2,istatus)
!c      do ll=1,lrza
!c         write(*,*)'ll=',ll,'iy_(ll)=',iy_(ll)
!c         write(*,*)'cynt2',(cynt2(i,ll),i=1,iy_(ll))
!c      enddo

!c-----normalized small radius
       vid = ncvid(ncid,'rya',istatus)
       call ncvgt(ncid,vid,1,lrz,rya,istatus)
!c      write(*,*)'rya',(rya(ll),ll=1,lrz)


!c-----distribution function f(i,j,ll) [vnorm**3/(cm**3*(cm/sec)**3)]
       vid=ncvid(ncid,'f',istatus)
       call ncvgt(ncid,vid,start,count,f,istatus)
       
!c      do ll=1,lrz
!c         write(*,*)' netcdfr3d ll=',ll
!c         do j=1,jx
!c            write(*,*)'ll=',ll,'j=',j,'f(i=1,...,iy_(ll))'
!c            write(*,*)(f(i,j,ll),i=1,iy_(ll))
!c         enddo
!c      enddo


      ! the energies wperp--wpar
      write(*,*)'shape of wperp ', shape(wperp)  
      istatus = nf_inq_varid(ncid, 'wperp', vid)
      istatus = nf_get_var_double(ncid, vid, wperp)       
      write(*,*)'proc_cql3d_op: after ncvgt, wperp = ', wperp(:, nt)

      write(*,*)'shape of wpar ', shape(wpar)    
      istatus = nf_inq_varid(ncid, 'wpar', vid)
      istatus = nf_get_var_double(ncid, vid, wpar)            
      write(*,*)'proc_cql3d_op: after ncvgt, wpar = ', wpar(:, nt)
      
       
!c      do ll=1,lrz
!c         write(*,*)' netcdfr3d ll=',ll
!c         do j=1,jx
!c            write(*,*)'ll=',ll,'j=',j,'f(i=1,...,iy_(ll))'
!c            write(*,*)(f(i,j,ll),i=1,iy_(ll))
!c         enddo
!c      enddo

!c-----Close netCDF file
       call ncclos(ncid,istatus)
       call check_err(istatus)

!c ************* Allocate space for module arrays.
!c **************If previously allocated, release space first.

	n_theta_max = iy
	n_u = jx
	n_psi = lrz
	n_t = nt

	IF ( ALLOCATED(n_theta_) ) DEALLOCATE (n_theta_, u, theta, rho_a, f_CQL, f_cql_2d, wperp_cql, wpar_cql)


	ALLOCATE( n_theta_(n_psi), stat=istat )
       IF (istat /= 0 ) THEN
       WRITE (*,'("read_CQL3D: allocate failed for n_theta_")')
       PAUSE
       END IF

	ALLOCATE( theta(n_theta_max, n_psi), stat=istat )
       IF (istat /= 0 ) THEN
       WRITE (*,'("read_CQL3D: allocate failed for theta")')
       PAUSE
       END IF

	ALLOCATE( u(n_u), stat=istat )
       IF (istat /= 0 ) THEN
       WRITE (*,'("read_CQL3D: allocate failed for u")')
       PAUSE
       END IF

	ALLOCATE( rho_a(n_psi), stat=istat )
       IF (istat /= 0 ) THEN
       WRITE (*,'("read_CQL3D: allocate failed for rho_a")')
       PAUSE
       END IF

	ALLOCATE( f_CQL(n_theta_max, n_u, n_psi), stat=istat )
      ALLOCATE( f_cql_2d(n_u, n_theta_max),     stat=istat )
       IF (istat /= 0 ) THEN
       WRITE (*,'("read_CQL3D: allocate failed for f_CQL")')
       PAUSE
       END IF
       
       

       ALLOCATE( wperp_cql(n_psi, n_t),    stat=istat )
       IF (istat /= 0 ) THEN
       WRITE (*,'("read_CQL3D: allocate failed for wperp_CQL")')
       PAUSE
       END IF  
       
       ALLOCATE( wpar_cql(n_psi, n_t),    stat=istat )
       IF (istat /= 0 ) THEN
       WRITE (*,'("read_CQL3D: allocate failed for wpar_CQL")')
       PAUSE
       END IF              
       

! ********** load Harvey arrays into module arrays and deallocate harvey arrays

       vc = vnorm
       n_theta_ = iy_
       theta = y
       u = x
       rho_a = rya
       f_CQL = f       

       wperp_cql = wperp
       wpar_cql = wpar

       DEALLOCATE (iy_, x, y, rya, f, wperp, wpar)

       return
       end subroutine netcdfr3d
!c
!c
       subroutine check_err(iret)
       integer iret
       include 'netcdf.inc'
       if (iret .ne. NF_NOERR) then
!c      print *, nf_strerror(iret)
          print *, 'netCDF error'
       stop 'check_err:'
       endif
       end subroutine check_err
!c

end module read_CQL3D
