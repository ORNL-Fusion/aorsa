       module read_CQL3D

       USE CQL_kinds_m

       implicit none

       integer :: n_theta_max  ! max number of pitch angles
       integer :: n_u          ! Number of normalized velocities
       integer :: n_psi        ! Number of flux surfaces
       integer :: n_t          ! Number of cql3d time steps

       integer, allocatable :: n_theta_(:) ! number of pitch angles for flux surface
       real(kind=real_kind), allocatable :: u(:) ! normalized velocity
       real(kind=real_kind), allocatable :: theta(:,:) ! pitch angles theta(i_theta, i_psi)
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

         use netcdf
         implicit none

         INTEGER :: istat

!      implicit integer (i-n), real*8 (a-h,o-z)
!      save

         character(*):: netcdfnm ! input filename

!c --- include file for netCDF declarations
!c --- (obtained from NetCDF distribution)

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
!c      real:: vnorm
!c      integer iy,jx,lrz
!c      integer iy_(lrza)
!c      real:: dx(jxa)                   !0.5(x(j+1)-x(j-1))
!c      real:: cint2(jxa)                !x**2dx
!c      real:: x(jxa)
!c      real:: dy(iya,lrza)              !0.5*(y(i+1,l)-y(i-1,l))
!c      real:: cynt2(iya,lrza)           !2*pi*sin(y)*dy
!c      real:: y(iya,lrza)
!c      real:: rya(lrza)                  !normalized small radius
!c      real:: f(iya,jxa,lrza)

         real:: vnorm
         integer:: iy,jx,lrz, nt, nt_id
         integer, allocatable :: iy_(:)
         real, allocatable :: x(:)
         real, allocatable :: y(:,:)
         real, allocatable :: rya(:)                  !normalized small radius
         real, allocatable :: f(:,:,:,:)
         real, allocatable, dimension(:,:,:) :: wperp, wpar   !perp, par energy/particle
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
         character(nf90_max_name):: name
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
         istatus=nf90_open(TRIM(netcdfnm),NF90_NOWRITE,ncid)
         write(*,*)'after ncopn ncid=',ncid,'istatus',istatus

!c.......................................................................
!c     read in dimension IDs and sizes

         write(*,*)'before ncdid xdim'

         istatus = nf90_inq_dimid(ncid,'xdim',xdim)
!       xdim = ncdid(ncid,'xdim',istatus)
         write(*,*)'after inq_dimid xdim=',xdim,'istatus',istatus

         istatus = nf90_inq_dimid(ncid,'ydim',ydim)
!       ydim = ncdid(ncid,'ydim',istatus)
         write(*,*)'after inq_dimid ydim=',ydim,'istatus',istatus

         istatus = nf90_inq_dimid(ncid,'rdim',rdim)
!       rdim = ncdid(ncid,'rdim',istatus)
         write(*,*)'after inq_dimid rdim=',rdim,'istatus',istatus

         istatus = nf90_inq_dimid(ncid,'species_dim',kdim)
!       kdim=ncdid(ncid,'species_dim',istatus)
         write(*,*)'after inq_dimid kdim=',kdim,'istatus',istatus

         istatus = nf90_inq_dimid(ncid,'tdim',nt_id)
         write(*,*)'proc_cql3d_op: after nf90_inq_dimid nt_id = ',nt_id,'istatus = ',istatus

         istatus = nf90_inq_dimid(ncid,'gen_species_dim',gkdim_id)
         write(*,*)'proc_cql3d_op: after nf90_inq_dimid ngen_id = ',gkdim_id,'istatus = ',istatus

         istatus = nf90_inquire_dimension(ncid, nt_id, len = nt)
       !call ncdinq(ncid, nt_id,'tdim', nt, istatus)
         write(*,*)'proc_cql3d_op: after inquire dimenstion, # of t steps = ',nt, ' istatus=',istatus

!c --- inquire about dimension sizes ---
!c     ncdinq(netCDF_id, dimension_id_from_ncdid, returned_dim_name,
!c     returned_dim_size)
!c     Note: for unlimited dimension, returned_dim_size=current maximum
!c     which is the same as the maximum record number

         istatus = nf90_inquire_dimension(ncid, ydim, len = iy)
         istatus = nf90_inquire_dimension(ncid, xdim, len = jx)
         istatus = nf90_inquire_dimension(ncid, rdim, len = lrz)
         istatus = nf90_inquire_dimension(ncid, gkdim, len = ngen)
         istatus = nf90_inquire_dimension(ncid, kdim, len = ntotal)

!c ************* Allocate space for Harvey arrays.

         ALLOCATE( iy_(lrz), stat=istat )
         IF (istat /= 0 ) THEN
            WRITE (*,'("read_CQL3D: allocate failed for iy_")')
         END IF

         ALLOCATE( y(iy, lrz), stat=istat )
         IF (istat /= 0 ) THEN
            WRITE (*,'("read_CQL3D: allocate failed for y")')
         END IF

         ALLOCATE( x(jx), stat=istat )
         IF (istat /= 0 ) THEN
            WRITE (*,'("read_CQL3D: allocate failed for x")')
         END IF

         ALLOCATE( rya(lrz), stat=istat )
         IF (istat /= 0 ) THEN
            WRITE (*,'("read_CQL3D: allocate failed for rya")')
         END IF

         ALLOCATE( f(iy, jx, lrz, ngen), stat=istat )
         IF (istat /= 0 ) THEN
            WRITE (*,'("read_CQL3D: allocate failed for f")')
         END IF

         ALLOCATE( wperp(lrz, ngen, nt), stat=istat )
         IF (istat /= 0 ) THEN
            WRITE (*,'("read_CQL3D: allocate failed for wperp")')
         END IF

         ALLOCATE( wpar(lrz, ngen, nt), stat=istat )
         IF (istat /= 0 ) THEN
            WRITE (*,'("read_CQL3D: allocate failed for wpar")')
         END IF


         count(1)=iy
         count(2)=jx
         count(3)=lrz

  310  format(1p,6e12.4)
  311  format(1p,6i6)
!c-----normalized momentum x (momentum/mass/vnorm) variables

!c     vnorm - character velocity (momentum-per-mass)[cms/sec]
       !       vid = ncvid(ncid,'vnorm',istatus)
       !       call ncvgt(ncid,vid,1,1,vnorm,istatus)
         istatus = nf90_inq_varid(ncid,'vnorm', vid)
         istatus = nf90_get_var(ncid, vid, vnorm)
         write(*,*)'after ncvgp vnorm=',vnorm

       !       vid = ncvid(ncid,'x',istatus)
       ! call ncvgt(ncid,vid,1,jx,x,istatus)
         istatus = nf90_inq_varid(ncid, 'x', vid)
         istatus = nf90_get_var(ncid, vid, x) !, start=1, count=jx)
         write(*,*) 'x status',istatus
         write(*,310)(x(j),j=1,jx)

!c      vid = ncvid(ncid,'dx',istatus)
!c      call ncvgt(ncid,vid,1,jx,dx,istatus)
!c      write(*,*)'dx',(dx(j),j=1,jx)

!c      vid = ncvid(ncid,'cint2',istatus)
!c      call ncvgt(ncid,vid,1,jx,cint2,istatus)
!c      write(*,*)'cint2',(cint2(j),j=1,jx)

!c-----pitch angle variables y
!       vid = ncvid(ncid,'iy_',istatus)
!       call ncvgt(ncid,vid,1,lrz,iy_,istatus)
         istatus = nf90_inq_varid(ncid, 'iy_', vid)
         istatus = nf90_get_var(ncid, vid, iy_)!, 1, lrz)
         write(*,*)'iy_',istatus
         write(*,311)(iy_(ll),ll=1,lrz)

         count_y(1)=iy
         count_y(2)=lrz
!       vid = ncvid(ncid,'y',istatus)
!       call ncvgt(ncid,vid,start,count_y,y,istatus)
         istatus = nf90_inq_varid(ncid, 'y', vid)
         istatus = nf90_get_var(ncid, vid, y)
         write(*,*)'y',istatus
         write(*,311)(y(1,ll),ll=1,lrz)
         
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
!       vid = ncvid(ncid,'rya',istatus)
!       call ncvgt(ncid,vid,1,lrz,rya,istatus)
         istatus = nf90_inq_varid(ncid, 'rya', vid)
         istatus = nf90_get_var(ncid, vid, rya)
         write(*,*)'rya',istatus, shape(rya)
         write(*,310) (rya(ll),ll=1,lrz)


!c-----distribution function f(i,j,ll) [vnorm**3/(cm**3*(cm/sec)**3)]
!       vid=ncvid(ncid,'f',istatus)
!       call ncvgt(ncid,vid,start,count,f,istatus)
         istatus = nf90_inq_varid(ncid, 'f', vid)
         istatus = nf90_get_var(ncid, vid, f)!, start, count)

!c      do ll=1,lrz
!c         write(*,*)' netcdfr3d ll=',ll
!c         do j=1,jx
!c            write(*,*)'ll=',ll,'j=',j,'f(i=1,...,iy_(ll))'
!c            write(*,*)(f(i,j,ll),i=1,iy_(ll))
!c         enddo
!c      enddo


      ! the energies wperp--wpar
         write(*,*)'shape of wperp ', shape(wperp)
         istatus = nf90_inq_varid(ncid, 'wperp', vid)
         istatus = nf90_get_var(ncid, vid, wperp)
         write(*,*)'proc_cql3d_op: after ncvgt, wperp = ', wperp(:,1,nt)

         write(*,*)'shape of wpar ', shape(wpar)
         istatus = nf90_inq_varid(ncid, 'wpar', vid)
         istatus = nf90_get_var(ncid, vid, wpar)
         write(*,*)'proc_cql3d_op: after ncvgt, wpar = ', wpar(:,1,nt)


!c      do ll=1,lrz
!c         write(*,*)' netcdfr3d ll=',ll
!c         do j=1,jx
!c            write(*,*)'ll=',ll,'j=',j,'f(i=1,...,iy_(ll))'
!c            write(*,*)(f(i,j,ll),i=1,iy_(ll))
!c         enddo
!c      enddo

!c-----Close netCDF file
      !       call ncclos(ncid,istatus)
         istatus = nf90_close(ncid)
         call check_err(istatus)

!c ************* Allocate space for module arrays.
!c **************If previously allocated, release space first.

         n_theta_max = iy
         n_u = jx
         n_psi = lrz
         n_t = nt

         IF ( ALLOCATED(n_theta_) ) THEN
            DEALLOCATE (n_theta_, u, theta, rho_a, f_CQL, f_cql_2d, &
       &     wperp_cql, wpar_cql)
         END IF

         ALLOCATE( n_theta_(n_psi), stat=istat )
         IF (istat /= 0 ) THEN
            WRITE (*,'("read_CQL3D: allocate failed for n_theta_")')
         END IF

         ALLOCATE( theta(n_theta_max, n_psi), stat=istat )
         IF (istat /= 0 ) THEN
            WRITE (*,'("read_CQL3D: allocate failed for theta")')
         END IF

         ALLOCATE( u(n_u), stat=istat )
         IF (istat /= 0 ) THEN
            WRITE (*,'("read_CQL3D: allocate failed for u")')
         END IF

         ALLOCATE( rho_a(n_psi), stat=istat )
         IF (istat /= 0 ) THEN
            WRITE (*,'("read_CQL3D: allocate failed for rho_a")')
         END IF

         ALLOCATE( f_CQL(n_theta_max, n_u, n_psi), stat=istat )
         ALLOCATE( f_cql_2d(n_u, n_theta_max),     stat=istat )
         IF (istat /= 0 ) THEN
            WRITE (*,'("read_CQL3D: allocate failed for f_CQL")')
         END IF


         ALLOCATE( wperp_cql(n_psi, n_t),    stat=istat )
         IF (istat /= 0 ) THEN
            WRITE (*,'("read_CQL3D: allocate failed for wperp_CQL")')
         END IF

         ALLOCATE( wpar_cql(n_psi, n_t),    stat=istat )
         IF (istat /= 0 ) THEN
            WRITE (*,'("read_CQL3D: allocate failed for wpar_CQL")')
         END IF


! ********** load Harvey arrays into module arrays and deallocate harvey arrays

         vc = vnorm
         n_theta_ = iy_
         theta = y
         u = x
         rho_a = rya
         f_CQL = f

         wperp_cql = wperp(:,1,:)
         wpar_cql = wpar(:,1,:)

         DEALLOCATE (iy_, x, y, rya, f, wperp, wpar)

         return
       end subroutine netcdfr3d


       subroutine check_err(iret)
         use netcdf
         integer iret
         if (iret .ne. NF90_NOERR) then
!c      print *, nf_strerror(iret)
            print *, 'netCDF error'
            stop 'check_err:'
         end if
       end subroutine check_err


end module read_CQL3D
