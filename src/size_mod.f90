      module size_mod
      implicit none

      integer, parameter:: nmodesmax=64 , mmodesmax=256
!JCW AORSA presently not happy if sizes are not 2x nmodex,
!     ---------------------------------------------------
!     Maximum number of modes allowed:
!     IMPORTANT NOTE: when changing these parameters, you 
!     must touch all files and recompile everything!!! 
!     ---------------------------------------------------

      end module size_mod

