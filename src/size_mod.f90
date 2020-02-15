      module size_mod
      implicit none

      integer nmodesmax, mmodesmax
	
!     ---------------------------------------------------
!     Maximum number of modes allowed:
!     IMPORTANT NOTE: when changing these parameters, you 
!     must touch all files and recomile everything!!! 
!     ---------------------------------------------------
      parameter (nmodesmax = 500)
      parameter (mmodesmax = 500)

      end module size_mod

