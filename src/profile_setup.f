C subroutine added by LAB to convert from flux to rz for numerical profiles.
c
c***************************************************************************
c
      subroutine flux_to_rz(nnodex, nnodey, profile_in, 
     &   profile_out, rho, nrho, rho_ij)
     
      implicit none

      integer,intent(in) :: nnodex, nnodey, nrho
      real,intent(in) :: profile_in(nrho)
      real,intent(out) :: profile_out(nnodex,nnodey)      
      real,intent(in) :: rho(nrho), rho_ij(nnodex,nnodey)
      integer :: i, j, k

      do i = 1, nnodex
         do j = 1, nnodey
            do k = 1, nrho

               if(rho_ij(i,j) .le. rho(1)) then
               !  need one sided to the left

c                 call flush(6)
                  profile_out(i,j) = profile_in(1) + 
     &            (profile_in(2) - profile_in(1)) /
     &            (rho(2) - rho(1)) * (rho_ij(i,j) - rho(1))
       
               elseif(rho_ij(i,j) .ge.rho(nrho)) then
               ! need one sidded to the right.
                   
c                 profile_out(i,j) = profile_in(nrho) +
c     &            (profile_in(nrho) - profile_in(nrho - 1)) /
c     &            (rho(nrho) - rho(nrho-1)) * (rho_ij(i,j) - rho(nrho))

                  profile_out(i,j) = profile_in(nrho)     
      
               elseif(rho_ij(i,j) - rho(k) .ge. 0.0 .and.
     &             rho_ij(i,j) -rho(k+1) .lt. 0.0 ) then
                  ! standard  caseprint*, 'last if'
                  profile_out(i,j) = profile_in(k)  +
     &            (profile_in(k+1) - profile_in(k))/
     &            (rho(k+1) - rho(k)) * (rho_ij(i,j) - rho(k))
               endif
               
               if (profile_out(i,j) .lt. 0.0) profile_out(i,j) = 1.0e-10
         
            enddo
         end do
      end do
         
      end subroutine flux_to_rz
      
c
c***************************************************************************
c
      subroutine fluxavg3(f, favg, rho, nxdim, nydim, nrhodim,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, r0, vol)

      implicit none

      integer nxdim, nydim, nrhodim, nnodex, nnodey, nnoderho
      integer n, i, j

      real f(nxdim, nydim), favg(nrhodim), rho(nxdim, nydim), r0,
     &   drho, dx, dy, fvol(nrhodim), vol(nrhodim), capr(nxdim)


      do n = 1, nnoderho
          fvol(n) = 0.0
          vol(n) = 0.0
          favg(n) = 0.0
      end do

      do i = 1, nnodex
         do j = 1, nnodey
            n = int(rho(i,j) / drho) + 1
            if(n .le. nnoderho)then
               fvol(n) = fvol(n) + dx * dy * capr(i) / r0 * f(i,j)
               vol(n) =  vol(n) + dx * dy * capr(i) / r0
            end if
         end do
      end do

      do n = 1, nnoderho
         favg(n) = 0.0
         if(vol(n) .ne. 0.0)favg(n) = fvol(n) / vol(n)
      end do

      do n = 2, nnoderho -1
         if(favg(n) .eq. 0.0)then
            favg(n) = (favg(n-1) + favg(n+1)) / 2.0
         end if
      end do

      if (favg(1) .eq. 0.0) favg(1) = favg(2)
      if (favg(nnoderho) .eq. 0.0) favg(nnoderho) = favg(nnoderho - 1)


  100 format (1i10, 1p,8e12.4)
  102 format (2i10)
      return
      end

c
c***************************************************************************
c

      subroutine plasma_state(ndata, ndim, rhod, xneavg, tekev, tikev,
     &     xn_beam, e_beam, zeffavg, nsmooth)

      implicit none

      integer n, ndata, ndim, ismooth, nsmooth

      real q

      real rhod(ndim), xneavg(ndim), tekev(ndim), tikev(ndim),
     &     omega(ndim), zeffavg(ndim), xn_beam(ndim), e_beam(ndim)

      q = 1.6e-19

*     ---------------------
*     read the profile data
*     ---------------------

      open(unit=70, file='n_electron.txt', status='unknown',
     &                                          form='formatted')
      open(unit=71, file='T_electron.txt', status='unknown',
     &                                          form='formatted')
      open(unit=72, file='T_ion.txt', status='unknown',
     &                                          form='formatted')
      open(unit=73, file='n_beam.txt', status='unknown',
     &                                          form='formatted')
      open(unit=74, file='beam_energy.txt', status='unknown',
     &                                          form='formatted')
      open(unit=75, file='Z_eff.txt', status='unknown',
     &                                          form='formatted')

      read(70, 1312) ndata
      write(6, 1312) ndata

      do n = 1, ndata
         read (70, 2314)rhod(n), xneavg(n)
         read (71, 2314)rhod(n), tekev(n)
         read (72, 2314)rhod(n), tikev(n)
         read (73, 2314)rhod(n), xn_beam(n)
         read (74, 2314)rhod(n), e_beam(n)
         read (75, 2314)rhod(n), zeffavg(n)

         xneavg(n) = xneavg(n) * 1.0e+19
         tekev(n)  = tekev(n)  * 1.0e+03 * q
         tikev(n)  = tikev(n)  * 1.0e+03 * q
         xn_beam(n)= xn_beam(n) * 1.0e+19
c         e_beam(n) = e_beam(n) * 1.0e+03 * q * 2./3.
         e_beam(n) = e_beam(n) * 1.0e+03 * q 


*     ---------------------------------------
*     For LDRD project, turn off beam n and T
*     ---------------------------------------
c         xn_beam(n) = 0.0
c         e_beam(n)  = 0.0

         write(6, 2313)rhod(n), xneavg(n), tekev(n), tikev(n),
     &      xn_beam(n), e_beam(n), zeffavg(n)
     
      end do
      
*     -----------
*     Smooth data
*     -----------
      
      ismooth = 1
      
      do n = 1, nsmooth
         call smooth1d(xneavg, ismooth, ndim)
         call smooth1d(tekev,  ismooth, ndim)
         call smooth1d(tikev,  ismooth, ndim)
         call smooth1d(xneavg, ismooth, ndim)
         call smooth1d(xn_beam,ismooth, ndim)
         call smooth1d(e_beam, ismooth, ndim)
         call smooth1d(zeffavg,ismooth, ndim)
      end do
      

 1312 format(1i10, 1p,8e12.4)
 2313 format(1p,8e12.4)
 2314 format(1p,e11.4, 1p,e11.4)

      return
      end

c
c***************************************************************************
c


c     1D smoothing algorithm

      subroutine smooth1d(x, ismooth, nr)

   
      real x(nr)
      real, dimension(:), allocatable :: y

   
   
      integer :: ismooth
      integer :: ir,nr,k,kmin,kmax
   

 
      allocate(y(nr))
      
c     Apply smoothing operation

   
      do ir = 1, nr
        y(ir) = 0.0

        kmin = max(1, ir-ismooth)  
        kmax = min(nr, ir+ismooth)
        
        do k = kmin, kmax
          y(ir) = y(ir) + x(k)
        enddo
        
        y(ir) = y(ir) / real((1+kmax - kmin))
      enddo

      x = y
      
      deallocate(y)
   
      end subroutine smooth1d








