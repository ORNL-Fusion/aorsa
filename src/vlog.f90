        subroutine vlog( vtarget, vsrc, n )
        implicit none
        integer:: n
        real:: vtarget(n),vsrc(n)
        integer:: i
        do i=1,n
          vtarget(i) = log( vsrc(i) )
        enddo
        return
        end subroutine vlog
