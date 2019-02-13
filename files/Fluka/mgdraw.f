*$ CREATE MGDRAW.FOR
*COPY MGDRAW
*                                                                      *
*=== mgdraw ===========================================================*
*                                                                      *
      SUBROUTINE MGDRAW ( ICODE, MREG )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1990-2006      by        Alfredo Ferrari           *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     MaGnetic field trajectory DRAWing: actually this entry manages   *
*                                        all trajectory dumping for    *
*                                        drawing                       *
*                                                                      *
*     Created on   01 march 1990   by        Alfredo Ferrari           *
*                                              INFN - Milan            *
*     Last change  05-may-06       by        Alfredo Ferrari           *
*                                              INFN - Milan            *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(CASLIM)'
      INCLUDE '(COMPUT)'
      INCLUDE '(SOURCM)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(GENSTK)'
      INCLUDE '(MGDDCM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(QUEMGD)'
      INCLUDE '(SUMCOU)'
      INCLUDE '(TRACKR)'
*
      DIMENSION DTQUEN ( MXTRCK, MAXQMG )
*
      CHARACTER*20 FILNAM
      LOGICAL LFCOPE
      SAVE LFCOPE
      DATA LFCOPE / .FALSE. /

      CHARACTER*8 MRGNAM, NRGNAM
*
*----------------------------------------------------------------------*
*                                                                      *
*     Icode = 1: call from Kaskad                                      *
*     Icode = 2: call from Emfsco                                      *
*     Icode = 3: call from Kasneu                                      *
*     Icode = 4: call from Kashea                                      *
*     Icode = 5: call from Kasoph                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      IF ( .NOT. LFCOPE ) THEN
*         LFCOPE = .TRUE.
*         IF ( KOMPUT .EQ. 2 ) THEN
*            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
*         ELSE
*            FILNAM = CFDRAW
*         END IF
*         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'UNKNOWN')
**         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
**     &          'UNFORMATTED' )
*      END IF
**      IF(JTRACK .EQ. 10 .OR. JTRACK .EQ. 11 .OR.
**     &   JTRACK .EQ. 13 .OR. JTRACK .EQ. 14) THEN
*      IF(JTRACK .EQ. 14) THEN
*         WRITE (IODRAW,*)  NTRACK, MTRACK, JTRACK, SNGL (ETRACK),
*     &        SNGL (WTRACK), ATRACK, LTRACK
**     100  FORMAT(5E25.15)
*         WRITE (IODRAW,*) ( SNGL (XTRACK (I)), SNGL (YTRACK (I)),
*     &        SNGL (ZTRACK (I)), I = 0, NTRACK ),
*     &        ( SNGL (DTRACK (I)), I = 1, MTRACK ),
*     &        SNGL (CTRACK)
**     101  FORMAT(5E25.15)
*      END IF
**  +-------------------------------------------------------------------*
**  |  Quenching is activated
*      IF ( LQEMGD ) THEN
*         IF ( MTRACK .GT. 0 ) THEN
*            RULLL  = ZERZER
*            CALL QUENMG ( ICODE, MREG, RULLL, DTQUEN )
*            WRITE (IODRAW) ( ( SNGL (DTQUEN (I,JBK)), I = 1, MTRACK ),
*     &                         JBK = 1, NQEMGD )
*         END IF
*      END IF
**  |  End of quenching
**  +-------------------------------------------------------------------*
      RETURN
*
*======================================================================*
*                                                                      *
*     Boundary-(X)crossing DRAWing:                                    *
*                                                                      *
*     Icode = 1x: call from Kaskad                                     *
*             19: boundary crossing                                    *
*     Icode = 2x: call from Emfsco                                     *
*             29: boundary crossing                                    *
*     Icode = 3x: call from Kasneu                                     *
*             39: boundary crossing                                    *
*     Icode = 4x: call from Kashea                                     *
*             49: boundary crossing                                    *
*     Icode = 5x: call from Kasoph                                     *
*             59: boundary crossing                                    *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY BXDRAW ( ICODE, MREG, NEWREG, XSCO, YSCO, ZSCO )
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
*         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
*     &          'UNFORMATTED' )
         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'UNKNOWN')
      END IF

      CALL GEOR2N(MREG,MRGNAM,IERR1)
      CALL GEOR2N(NEWREG,NRGNAM,IERR2)


      IF(IERR1 .NE. 0 .OR. IERR2 .NE. 0) STOP 
      IF(MRGNAM .EQ. "VOID" .AND. NRGNAM .EQ. "TrigStn0" .OR.
     &   MRGNAM .EQ. "VOID" .AND. NRGNAM .EQ. "TrigStn1" .OR.
     &   MRGNAM .EQ. "VacBox" .AND. NRGNAM .EQ. "PixStn0" .OR.
     &   MRGNAM .EQ. "VacBox" .AND. NRGNAM .EQ. "PixStn1" .OR.
     &   MRGNAM .EQ. "VacBox" .AND. NRGNAM .EQ. "PixStn2" .OR.
     &   MRGNAM .EQ. "VacBox" .AND. NRGNAM .EQ. "PixStn3" .OR.
     &   MRGNAM .EQ. "VacBox" .AND. NRGNAM .EQ. "PixStn4" .OR.
     &   MRGNAM .EQ. "VOID" .AND. NRGNAM .EQ. "MS0" .OR.
     &   MRGNAM .EQ. "VOID" .AND. NRGNAM .EQ. "MS1" .OR.
     &   MRGNAM .EQ. "VOID" .AND. NRGNAM .EQ. "MS2" .OR.
*     &   MRGNAM .EQ. "VOID" .AND. NRGNAM .EQ. "MS3" .OR.
     &   MRGNAM .EQ. "VOID" .AND. NRGNAM .EQ. "MS3") THEN
        
         IF (JTRACK .EQ. 10 .OR. JTRACK .EQ. 11 .OR.
     &       JTRACK .EQ. 13 .OR. JTRACK .EQ. 14 .OR.
     &       JTRACK .EQ. 15 .OR. JTRACK .EQ. 16 .OR.
     &       JTRACK .EQ. 1  .OR. JTRACK .EQ. 2  .OR.
     &       JTRACK .EQ. 3  .OR. JTRACK .EQ. 4) THEN
*            WRITE(IODRAW,*) 'EVENT = ', NCASE
*            WRITE(IODRAW,*) NRGNAM
* LLOUSE is the particle mother, see http://www.fluka.org/fluka.php?id=faq&sub=5
* If zero use the same ID as the one of the particle being considered
         if (llouse .eq. 0) then
	     llouse=jtrack
	 endif
            WRITE (IODRAW,*) NRGNAM,JTRACK,llouse,
     &	    SNGL (ETRACK), SNGL (XSCO), SNGL (YSCO), SNGL (ZSCO), 
     &           SNGL (CXTRCK), SNGL (CYTRCK), SNGL (CZTRCK)
*     &           XTRACK(NTRACK), YTRACK(NTRACK), ZTRACK(NTRACK)
*            IF(JTRACK .NE. 3 .AND. JTRACK .NE. 4) THEN
*               WRITE(IODRAW,*) 'PARTICLE ORIGIN'
*               WRITE(IODRAW,*)(ISPUSR(I),I=1,5),(SNGL(SPAUSR(I)),i=1,10)
*            END IF
*     100  FORMAT(5E25.15)
*     WRITE (IODRAW,*) ( SNGL (XTRACK (I)), SNGL (YTRACK (I)),
*     &           SNGL (ZTRACK (I)), I = 0, NTRACK ),
*     &           ( SNGL (DTRACK (I)), I = 1, MTRACK ),
*     &           SNGL (CTRACK)
*     101  FORMAT(5E25.15)
         END IF
         
      END IF
      
      RETURN
*
*======================================================================*
*                                                                      *
*     Event End DRAWing:                                               *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY EEDRAW ( ICODE )
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
*         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
*     &          'UNFORMATTED' )
         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'UNKNOWN')
      END IF
      WRITE (IODRAW,*) '**************** end of event **************'
      RETURN
*
*======================================================================*
*                                                                      *
*     ENergy deposition DRAWing:                                       *
*                                                                      *
*     Icode = 1x: call from Kaskad                                     *
*             10: elastic interaction recoil                           *
*             11: inelastic interaction recoil                         *
*             12: stopping particle                                    *
*             13: pseudo-neutron deposition                            *
*             14: escape                                               *
*             15: time kill                                            *
*     Icode = 2x: call from Emfsco                                     *
*             20: local energy deposition (i.e. photoelectric)         *
*             21: below threshold, iarg=1                              *
*             22: below threshold, iarg=2                              *
*             23: escape                                               *
*             24: time kill                                            *
*     Icode = 3x: call from Kasneu                                     *
*             30: target recoil                                        *
*             31: below threshold                                      *
*             32: escape                                               *
*             33: time kill                                            *
*     Icode = 4x: call from Kashea                                     *
*             40: escape                                               *
*             41: time kill                                            *
*             42: delta ray stack overflow                             *
*     Icode = 5x: call from Kasoph                                     *
*             50: optical photon absorption                            *
*             51: escape                                               *
*             52: time kill                                            *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY ENDRAW ( ICODE, MREG, RULL, XSCO, YSCO, ZSCO )
*      IF ( .NOT. LFCOPE ) THEN
*         LFCOPE = .TRUE.
*         IF ( KOMPUT .EQ. 2 ) THEN
*            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
*         ELSE
*            FILNAM = CFDRAW
*         END IF
*         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
*     &          'UNFORMATTED' )
*      END IF
*      WRITE (IODRAW)  0, ICODE, JTRACK, SNGL (ETRACK), SNGL (WTRACK)
*      WRITE (IODRAW)  SNGL (XSCO), SNGL (YSCO), SNGL (ZSCO), SNGL (RULL)
*  +-------------------------------------------------------------------*
*  |  Quenching is activated : calculate quenching factor
*  |  and store quenched energy in DTQUEN(1, jbk)
*      IF ( LQEMGD ) THEN
*         RULLL = RULL
*         CALL QUENMG ( ICODE, MREG, RULLL, DTQUEN )
*         WRITE (IODRAW) ( SNGL (DTQUEN(1, JBK)), JBK = 1, NQEMGD )
*      END IF
*  |  end quenching
*  +-------------------------------------------------------------------*
      RETURN
*
*======================================================================*
*                                                                      *
*     SOurce particle DRAWing:                                         *
*                                                                      *
*======================================================================*
*
      ENTRY SODRAW
*      IF ( .NOT. LFCOPE ) THEN
*         LFCOPE = .TRUE.
*         IF ( KOMPUT .EQ. 2 ) THEN
*            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
*         ELSE
*            FILNAM = CFDRAW
*         END IF
*         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
*     &          'UNFORMATTED' )
*      END IF

      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
*         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
*     &          'UNFORMATTED' )
         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'UNKNOWN')
      END IF


*      WRITE (IODRAW) -NCASE, NPFLKA, NSTMAX, SNGL (TKESUM),
*     &                SNGL (WEIPRI)
*  +-------------------------------------------------------------------*
*  |  (Radioactive) isotope: it works only for 1 source particle on
*  |  the stack for the time being
*      IF ( ILOFLK (NPFLKA) .GE. 100000 .AND. LRADDC (NPFLKA) ) THEN
*         IARES  = MOD ( ILOFLK (NPFLKA), 100000  )  / 100
*         IZRES  = MOD ( ILOFLK (NPFLKA), 10000000 ) / 100000
*         IISRES = ILOFLK (NPFLKA) / 10000000
*         IONID  = ILOFLK (NPFLKA)
*         WRITE (IODRAW) ( IONID,SNGL(-TKEFLK(I)),
*     &                    SNGL (WTFLK(I)), SNGL (XFLK (I)),
*     &                    SNGL (YFLK (I)), SNGL (ZFLK (I)),
*     &                    SNGL (TXFLK(I)), SNGL (TYFLK(I)),
*     &                    SNGL (TZFLK(I)), I = 1, NPFLKA )
*  |
*  +-------------------------------------------------------------------*
*  |  Patch for heavy ions: it works only for 1 source particle on
*  |  the stack for the time being
*      ELSE IF ( ABS (ILOFLK (NPFLKA)) .GE. 10000 ) THEN
*         IONID = ILOFLK (NPFLKA)
*         CALL DCDION ( IONID )
*         WRITE (IODRAW) ( IONID,SNGL(TKEFLK(I)+AMNHEA(-IONID)),
*     &                    SNGL (WTFLK(I)), SNGL (XFLK (I)),
*     &                    SNGL (YFLK (I)), SNGL (ZFLK (I)),
*     &                    SNGL (TXFLK(I)), SNGL (TYFLK(I)),
*     &                    SNGL (TZFLK(I)), I = 1, NPFLKA )
*  |
*  +-------------------------------------------------------------------*
*  |  Patch for heavy ions: ???
*      ELSE IF ( ILOFLK (NPFLKA) .LT. -6 ) THEN
*         WRITE (IODRAW) ( IONID,SNGL(TKEFLK(I)+AMNHEA(-ILOFLK(NPFLKA))),
*     &                    SNGL (WTFLK(I)), SNGL (XFLK (I)),
*     &                    SNGL (YFLK (I)), SNGL (ZFLK (I)),
*     &                    SNGL (TXFLK(I)), SNGL (TYFLK(I)),
*     &                    SNGL (TZFLK(I)), I = 1, NPFLKA )
*  |
*  +-------------------------------------------------------------------*
*  |
*      ELSE
*      WRITE(IODRAW,*) 'PRIMARY PARTICLE KINEMATICS'
*      WRITE (IODRAW,*)   'Primary',NPFLKA, ILOFLK(NPFLKA),
*     &     SNGL(TKEFLK(NPFLKA)+AM(ILOFLK(NPFLKA))),
*     &     SNGL (XFLK (NPFLKA)),
*     &     SNGL (YFLK (NPFLKA)), SNGL (ZFLK (NPFLKA)),
*     &     SNGL (TXFLK(NPFLKA)), SNGL (TYFLK(NPFLKA)),
*     &                    SNGL (TZFLK(NPFLKA))

* Loop on the primaries, put ID of mother equal to ID of particle
      do i=1,npflka
      WRITE (IODRAW,*) 'Primary ', (ILOFLK(I)),(iloflk(i)),
     &     SNGL(TKEFLK(I)+AM(ILOFLK(I))),
     &     SNGL (XFLK (I)),
     &     SNGL (YFLK (I)), SNGL (ZFLK (I)),
     &     SNGL (TXFLK(I)), SNGL (TYFLK(I)),
     &     SNGL (TZFLK(I))
      enddo
      
*     END IF
*  |
*  +-------------------------------------------------------------------*
      RETURN
*
*======================================================================*
*                                                                      *
*     USer dependent DRAWing:                                          *
*                                                                      *
*     Icode = 10x: call from Kaskad                                    *
*             100: elastic   interaction secondaries                   *
*             101: inelastic interaction secondaries                   *
*             102: particle decay  secondaries                         *
*             103: delta ray  generation secondaries                   *
*             104: pair production secondaries                         *
*             105: bremsstrahlung  secondaries                         *
*             110: decay products                                      *
*     Icode = 20x: call from Emfsco                                    *
*             208: bremsstrahlung secondaries                          *
*             210: Moller secondaries                                  *
*             212: Bhabha secondaries                                  *
*             214: in-flight annihilation secondaries                  *
*             215: annihilation at rest   secondaries                  *
*             217: pair production        secondaries                  *
*             219: Compton scattering     secondaries                  *
*             221: photoelectric          secondaries                  *
*             225: Rayleigh scattering    secondaries                  *
*     Icode = 30x: call from Kasneu                                    *
*             300: interaction secondaries                             *
*     Icode = 40x: call from Kashea                                    *
*             400: delta ray  generation secondaries                   *
*  For all interactions secondaries are put on GENSTK common (kp=1,np) *
*  but for KASHEA delta ray generation where only the secondary elec-  *
*  tron is present and stacked on FLKSTK common for kp=npflka          *
*                                                                      *
*======================================================================*
*
      ENTRY USDRAW ( ICODE, MREG, XSCO, YSCO, ZSCO )
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
     &          'UNFORMATTED' )
      END IF
* No output by default:
      RETURN
*=== End of subrutine Mgdraw ==========================================*
      END

