      module size_mod
      implicit none

      integer nmodesmax, mmodesmax
	
!     ---------------------------------------------------
!     Maximum number of modes allowed:
!     IMPORTANT NOTE: when changing these parameters, you 
!     must touch all files and recomile everything!!! 
!     ---------------------------------------------------
      parameter (nmodesmax = 256)
      parameter (mmodesmax = 256)

      end module size_mod

