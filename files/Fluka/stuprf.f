*$ CREATE STUPRF.FOR
*COPY STUPRF
*
*=== stuprf ===========================================================*
*
      SUBROUTINE STUPRF ( IJ, MREG, XX, YY, ZZ, NPSECN, NPPRMR )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1997-2005      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     SeT User PRoperties for Fluka particles:                         *
*                                                                      *
*     Created on  09 october 1997  by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on  14-jul-05    by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(EVTFLG)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(TRACKR)'

      INCLUDE '(GENSTK)'
*
*      LOUSE   (NPFLKA)  = LLOUSE
      louse(npflka) = jtrack
      sparek(1,npflka) = etrack
      DO 100 ISPR = 2, MKBMX1
         SPAREK (ISPR,NPFLKA) = SPAUSR (ISPR)
  100 CONTINUE
      DO 200 ISPR = 1, MKBMX2
         ISPARK (ISPR,NPFLKA) = ISPUSR (ISPR)
  200 CONTINUE


      IF ( LDECAY .OR. LINEVT ) THEN
         ISPARK (1,NPFLKA) = IJ
         ISPARK (2,NPFLKA) = ISPUSR(MKBMX2)
         ISPARK (3,NPFLKA) = KPART (NPSECN)
         ISPARK (4,NPFLKA) = MREG
         IF(LDECAY) THEN
            ISPARK (5,NPFLKA) = 10
         ELSE IF(LINEVT) THEN
            ISPARK (5,NPFLKA) = 20
         END IF
         SPAREK (1,NPFLKA)  = ETRACK
         SPAREK (2,NPFLKA)  = CXTRCK;
         SPAREK (3,NPFLKA)  = CYTRCK; 
         SPAREK (4,NPFLKA)  = CZTRCK;
         SPAREK (5,NPFLKA)  = XX
         SPAREK (6,NPFLKA)  = YY
         SPAREK (7,NPFLKA)  = ZZ
         SPAREK (8,NPFLKA)  = PLR(NPSECN) * CXR (NPSECN)
         SPAREK (9,NPFLKA)  = PLR(NPSECN) * CYR (NPSECN)
         SPAREK (10,NPFLKA) = PLR(NPSECN) * CZR (NPSECN)
      END IF
      
*  Increment the track number and put it into the last flag:
      IF ( NPSECN .GT. NPPRMR ) THEN
         IF ( NTRCKS .EQ. 2000000000 ) NTRCKS = -2000000000
         NTRCKS = NTRCKS + 1
         ISPARK (MKBMX2,NPFLKA) = NTRCKS
      END IF
      RETURN
*=== End of subroutine Stuprf =========================================*
      END

