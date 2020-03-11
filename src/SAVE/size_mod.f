      module size_mod
      implicit none

      integer nmodesmax, mmodesmax
	
!     ---------------------------------------------------
!     Maximum number of modes allowed:
!     IMPORTANT NOTE: when changing these parameters, you 
!     must touch all files and recomile everything!!! 
!     ---------------------------------------------------
      parameter (nmodesmax = 800)
      parameter (mmodesmax = 400)

      end module size_mod

