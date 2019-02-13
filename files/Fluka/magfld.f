*$ CREATE MAGFLD.FOR
*COPY MAGFLD
*
*===magfld=============================================================*
*
      SUBROUTINE MAGFLD ( X, Y, Z, BTX, BTY, BTZ, B, NREG, IDISC )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'

*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1988-2009      by Alberto Fasso` & Alfredo Ferrari *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     Created  in 1988    by     Alberto Fasso`, CERN - TIS            *
*                                                                      *
*     Last change on 15-oct-09     by    Alfredo Ferrari               *
*                                                                      *
*     Input variables:                                                 *
*            x,y,z = current position                                  *
*            nreg  = current region                                    *
*     Output variables:                                                *
*            btx,bty,btz = cosines of the magn. field vector           *
*            B = magnetic field intensity (Tesla)                      *
*            idisc = set to 1 if the particle has to be discarded      *
*                                                                      *
*----------------------------------------------------------------------*
*
*     Last change on 5-oct-10     by    Advanced FLUKA course teacher  *
*
      INCLUDE '(RTDFCM)'
      INCLUDE '(LTCLCM)'
*
      DIMENSION IROTIN(10)
      LOGICAL LFIRST
      SAVE LFIRST
      SAVE  DIPFLD, GRADIENT, NDIPOLE, IROTIN
      SAVE  NDIPOL1, NDIPOL2, NDIPOL3, NDIPOL4 
      SAVE  NDIPOL5, NDIPOL6, NDIPOL7, NDIPOL8 
      SAVE  NDIPOL9, NDIPO10, NDIPO11, NDIPO12
      SAVE  NDIPO13, TORFLD2
*
      DATA LFIRST / .TRUE. /, IROTIN / 10*0 /
*
      IDISC = 0
      IF (LFIRST) THEN
*     gradient in tesla per meter for 1 GeV protons
         GRADIENT = 1.D+00 
*     field in tesla for 14 GeV protons
         DIPFLD   = 3.0D+00 
* ATLAS TOROID
*         TORFLD2  = 0.16D+00 
* ACM TOROID
*         TORFLD2  = 32.D+00 
         TORFLD2 = 0.25*100.D+00
*         DIPFLD   = 0.00
*         DIPFLD2  = 0.00
         CALL GEON2R("VacBox  ",NDIPOL1,IERR)
         CALL GEON2R("PixStn0 ",NDIPOL2,IERR)
         CALL GEON2R("Dum0    ",NDIPOL3,IERR)
         CALL GEON2R("PixStn1 ",NDIPOL4,IERR)
         CALL GEON2R("Dum1    ",NDIPOL5,IERR)
         CALL GEON2R("PixStn2 ",NDIPOL6,IERR)
         CALL GEON2R("Dum2    ",NDIPOL7,IERR)
         CALL GEON2R("PixStn3 ",NDIPOL8,IERR)
         CALL GEON2R("Dum3    ",NDIPOL9,IERR)
         CALL GEON2R("PixStn4 ",NDIPO10,IERR)
         CALL GEON2R("Dum4    ",NDIPO11,IERR)
*         CALL GEON2R("MS2     ",NDIPO12,IERR)
         CALL GEON2R("MSBox   ",NDIPO12,IERR)
         LFIRST  = .FALSE.
      END IF
*
      IF (NREG.EQ.NDIPOL1.OR.NREG.EQ.NDIPOL2.OR.
     &    NREG.EQ.NDIPOL3.OR.NREG.EQ.NDIPOL4.OR.
     &    NREG.EQ.NDIPOL5.OR.NREG.EQ.NDIPOL6.OR.
     &    NREG.EQ.NDIPOL7.OR.NREG.EQ.NDIPOL8.OR.
     &    NREG.EQ.NDIPOL9.OR.NREG.EQ.NDIPO10.OR.
     &    NREG.EQ.NDIPO11) THEN
         BTX = -ONEONE
         BTY = ZERZER
         BTZ = ZERZER
         B   = DIPFLD 
      ELSE IF (NREG.EQ.NDIPO12) THEN
         R = SQRT(X**2 + Y **2)
         BTX = -Y / R 
         BTY = X / R
         BTZ = ZERZER
* ATLAS TOROID
*         B   = TORFLD2  
* ACM TOROID
         B   = TORFLD2 / R 
      ELSE
         CALL FLABRT ( 'MAGFLD', 'NO MAGNETIC FIELD IN THIS REGION!' )
      END IF
      RETURN
*=== End of subroutine magfld =========================================*
      END
