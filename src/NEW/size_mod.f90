      module size_mod
      implicit none

      integer nmodesmax, mmodesmax
	
!     ---------------------------------------------------
!     Maximum number of modes allowed:
!     IMPORTANT NOTE: when changing these parameters, you 
!     must touch all files and recomile everything!!! 
!     ---------------------------------------------------
      parameter (nmodesmax = 450)
      parameter (mmodesmax = 450)

      end module size_mod

