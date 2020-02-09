      subroutine ezplot7(title, titll, titlr, titlb, x1, 
     .   y1, y2, y3, y4, y5, y6, ye,
     .   nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax), y2(nrmax), y3(nrmax), y4(nrmax),
     .     y5(nrmax), y6(nrmax), ye(nrmax)
      real y1max, y2max, y3max, y4max, y5max, y6max, yemax
      real y1min, y2min, y3min, y4min, y5min, y6min, yemin
      real ymin, ymax

      character*32 title
      character*32 titll
      character*32 titlr
      character*32 titlb
      
      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan,
     1    ncolelec, ncolln2, ncollin, ncolbrd
     
      integer nblack,nred,nyellow, ngreen,naqua,npink,
     1   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     1   nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      call a1mnmx(y1, nrmax, nr, y1min, y1max)
      call a1mnmx(y2, nrmax, nr, y2min, y2max)
      call a1mnmx(y3, nrmax, nr, y3min, y3max)
      call a1mnmx(y4, nrmax, nr, y4min, y4max)
      call a1mnmx(y5, nrmax, nr, y5min, y5max)
      call a1mnmx(y6, nrmax, nr, y6min, y6max)
      call a1mnmx(ye, nrmax, nr, yemin, yemax)

      ymax = amax1(y1max, y2max, y3max, y4max, y5max, y6max, yemax)
      ymin = amin1(y1min, y2min, y3min, y4min, y5min, y6min, yemin)
c      ymin=0.0

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax=ymax*1.1
      
c      ymax = 2.8e+06
      
c Advance plotter to a new page, define coordinate range of graph and draw axes
c      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)      
      
c Label the axes (note use of \u and \d for raising exponent).

      call pglab(titlb, titll, title)
      
c Plot the line graph.

      call pgsci(nred) 
      call pgline(nr, x1, ye)  

      call pgsci(ncyan)      
      call pgline(nr, x1, y1)

      call pgsci(nblue)       
      call pgline(nr, x1, y2)
      
      call pgsci(ngreen) 
      call pgline(nr, x1, y3)
      
      call pgsci(nyellow) 
      call pgline(nr, x1, y4)
      
      call pgsci(nmagenta) 
      call pgline(nr, x1, y5)
      
      call pgsci(norange) 
      call pgline(nr, x1, y6)
      
      call pgsci(nblack) 
	    
  300 format (1p9e11.3)

      return
      end
