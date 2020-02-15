      SUBROUTINE SETUPBLACS( M, N, MB, NB, WANTNPROCS, MYID, NPROC,
     $                       MYROW, MYCOL, NPROW, NPCOL, CONTXT )
* intent(in) :: m,n,  mb,nb
* intent(inout) :: myid,nproc
* myrow,mycol     coordinate of processor
* nprow,npcol  processor grid is nprow by npcol
*     .. Scalar Arguments ..
      INTEGER            CONTXT, M, MB, MYCOL, MYID, MYROW, N, NB,
     $                   NPCOL, NPROC, NPROW, WANTNPROCS
*     ..
*     .. Local Scalars ..
      LOGICAL            ISROOT
      CHARACTER          ITOP
      INTEGER            COLSRC, ROWSRC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ASSERT, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,
     $                   BLACS_SETUP, IGEBR2D, IGEBS2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, MAX, MIN, MOD, SQRT
*     ..
*     .. Executable Statements ..
      CALL BLACS_PINFO( MYID, NPROC )
      IF( NPROC.LT.1 ) THEN
         NPROC = WANTNPROCS
         CALL BLACS_SETUP( MYID, NPROC )
      ENDIF
      CALL BLACS_GET( 0, 0, CONTXT )
      ISROOT = ( MYID.EQ.0 )
      ITOP = ' '
      ROWSRC = 0
      COLSRC = 0
* broadcast parameters for matrix
      CALL BLACS_GRIDINIT( CONTXT, 'row', 1, NPROC )
      IF( ISROOT ) THEN
         CALL IGEBS2D( CONTXT, 'all', ITOP, 1, 1, M, 1 )
         CALL IGEBS2D( CONTXT, 'all', ITOP, 1, 1, N, 1 )
         CALL IGEBS2D( CONTXT, 'all', ITOP, 1, 1, MB, 1 )
         CALL IGEBS2D( CONTXT, 'all', ITOP, 1, 1, NB, 1 )
      ELSE
         CALL IGEBR2D( CONTXT, 'all', ITOP, 1, 1, M, 1, ROWSRC, COLSRC )
         CALL IGEBR2D( CONTXT, 'all', ITOP, 1, 1, N, 1, ROWSRC, COLSRC )
         CALL IGEBR2D( CONTXT, 'all', ITOP, 1, 1, MB, 1, ROWSRC,
     $                 COLSRC )
         CALL IGEBR2D( CONTXT, 'all', ITOP, 1, 1, NB, 1, ROWSRC,
     $                 COLSRC )
      ENDIF
      CALL BLACS_GRIDEXIT( CONTXT )
* initialize blacs for rectangular mesh
      IF( NPROW*NPCOL.NE.NPROC ) THEN
* calculate a reasonable vaue for nprow,npcol
         NPROW = MIN( NPROC, MAX( 1, 1+INT( SQRT( 0.01d0+
     $           DBLE( NPROC ) ) ) ) )
   10    CONTINUE
         IF( ( MOD( NPROC, NPROW ).NE.0 ) .AND. ( NPROW.GE.1 ) ) THEN
            NPROW = NPROW - 1
            GOTO 10
         ENDIF
   20    CONTINUE
         CALL ASSERT( ( NPROW.GE.1 ) .AND. ( MOD( NPROC, NPROW ).EQ.0 ),
     $                '** setupblacs: invalid nprow ', NPROW )
         NPCOL = NPROC / NPROW
         CALL ASSERT( NPROW*NPCOL.EQ.NPROC, 'invalid nprow*npcol ',
     $                NPROW*NPCOL )
      ENDIF
      CALL BLACS_GET( 0, 0, CONTXT )
      CALL BLACS_GRIDINIT( CONTXT, 'Rowmajor', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( CONTXT, NPROW, NPCOL, MYROW, MYCOL )
      RETURN
      END
