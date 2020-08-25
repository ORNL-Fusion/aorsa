      module size_mod
      implicit none

      integer nmodesmax, mmodesmax
	
!     ---------------------------------------------------
!     Maximum number of modes allowed:
!     IMPORTANT NOTE: when changing these parameters, you 
!     must touch all files and recomile everything!!! 
!     ---------------------------------------------------
      parameter (nmodesmax = 600)
      parameter (mmodesmax = 600)

      end module size_mod

