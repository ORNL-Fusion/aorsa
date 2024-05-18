      module size_mod
      implicit none

      integer, parameter:: nmodesmax=16, mmodesmax=16
!JCW AORSA presently not happy if sizes are not 2x nmodex,
!     ---------------------------------------------------
!     Maximum number of modes allowed:
!     IMPORTANT NOTE: when changing these parameters, you 
!     must touch all files and recomile everything!!! 
!     ---------------------------------------------------

      end module size_mod

