*$ CREATE SOURCE.FOR
*COPY SOURCE
*
*=== source ===========================================================*
*
      SUBROUTINE SOURCE ( NOMORE )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1990-2010      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     New source for FLUKA9x-FLUKA20xy:                                *
*                                                                      *
*     Created on 07 January 1990   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on  17-Oct-10    by    Alfredo Ferrari               *
*                                                                      *
*  This is just an example of a possible user written source routine.  *
*  note that the beam card still has some meaning - in the scoring the *
*  maximum momentum used in deciding the binning is taken from the     *
*  beam momentum.  Other beam card parameters are obsolete.            *
*                                                                      *
*       Output variables:                                              *
*                                                                      *
*              Nomore = if > 0 the run will be terminated              *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(BEAMCM)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(IOIOCM)'
      INCLUDE '(LTCLCM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(SOURCM)'
      INCLUDE '(SUMCOU)'
*
      LOGICAL LFIRST
*
      SAVE LFIRST
      DATA LFIRST / .TRUE. /

      integer jpsi_flag
      save jpsi_flag
*======================================================================*
*                                                                      *
*                 BASIC VERSION                                        *
*                                                                      *
*======================================================================*
      NOMORE = 0
*  +-------------------------------------------------------------------*
*  |  First call initializations:
      IF ( LFIRST ) THEN
*  |  *** The following 3 cards are mandatory ***
         TKESUM = ZERZER
         LFIRST = .FALSE.
         LUSSRC = .TRUE.
*     |  *** User initialization ***
         jpsi_flag = int(whasou(1))
      END IF
*  |
*  +-------------------------------------------------------------------*
*  Push one source particle to the stack. Note that you could as well
*  push many but this way we reserve a maximum amount of space in the
*  stack for the secondaries to be generated
*  Npflka is the stack counter: of course any time source is called it
*     must be =0

      if(jpsi_flag.gt.0) then
         call signal_rootgen(nomore,jpsi_flag)
      else
         NPFLKA = NPFLKA + 1
*  Wt is the weight of the particle
         WTFLK  (NPFLKA) = ONEONE
         WEIPRI = WEIPRI + WTFLK (NPFLKA)
*  Particle type (1=proton.....). Ijbeam is the type set by the BEAM
*  card
*  +-------------------------------------------------------------------*
*  |  (Radioactive) isotope:
         IF ( IJBEAM .EQ. -2 .AND. LRDBEA ) THEN
            IARES  = IPROA
            IZRES  = IPROZ
            IISRES = IPROM
            CALL STISBM ( IARES, IZRES, IISRES )
            IJHION = IPROZ  * 1000 + IPROA
            IJHION = IJHION * 100 + KXHEAV
            IONID  = IJHION
            CALL DCDION ( IONID )
            CALL SETION ( IONID )
*  |
*  +-------------------------------------------------------------------*
*  |  Heavy ion:
         ELSE IF ( IJBEAM .EQ. -2 ) THEN
            IJHION = IPROZ  * 1000 + IPROA
            IJHION = IJHION * 100 + KXHEAV
            IONID  = IJHION
            CALL DCDION ( IONID )
            CALL SETION ( IONID )
            ILOFLK (NPFLKA) = IJHION
*  |  Flag this is prompt radiation
            LRADDC (NPFLKA) = .FALSE.
*  |  Group number for "low" energy neutrons, set to 0 anyway
            IGROUP (NPFLKA) = 0
*  |
*  +-------------------------------------------------------------------*
*  |  Normal hadron:
         ELSE
            IONID = IJBEAM
            ILOFLK (NPFLKA) = IJBEAM
*  Put by hand J/psi code
*         IONID = 443
*         ILOFLK (NPFLKA) = 443
*  |  Flag this is prompt radiation
            LRADDC (NPFLKA) = .FALSE.
*  |  Group number for "low" energy neutrons, set to 0 anyway
            IGROUP (NPFLKA) = 0
         END IF
*  |
*  +-------------------------------------------------------------------*
*  From this point .....
*  Particle generation (1 for primaries)
         LOFLK  (NPFLKA) = 1
*     User dependent flag:
         LOUSE  (NPFLKA) = 0
*     No channeling:
         LCHFLK (NPFLKA) = .FALSE.
         DCHFLK (NPFLKA) = ZERZER
*     User dependent spare variables:
         DO 100 ISPR = 1, MKBMX1
            SPAREK (ISPR,NPFLKA) = ZERZER
 100     CONTINUE
*     User dependent spare flags:
         DO 200 ISPR = 1, MKBMX2
            ISPARK (ISPR,NPFLKA) = 0
 200     CONTINUE
*  Save the track number of the stack particle:
         ISPARK (MKBMX2,NPFLKA) = NPFLKA
         NPARMA = NPARMA + 1
         NUMPAR (NPFLKA) = NPARMA
         NEVENT (NPFLKA) = 0
         DFNEAR (NPFLKA) = +ZERZER
*     ... to this point: don't change anything
*  Particle age (s)
         AGESTK (NPFLKA) = +ZERZER
         AKNSHR (NPFLKA) = -TWOTWO
* Rapidity of the particle
* Gaussian distribution with sigma=1
* 51   CONTINUE   
*      CALL FLNRRN(RGAUSS)
*      Y = 1.9 + RGAUSS
*      IF(Y > 3.7 .OR. Y < 1.8) GOTO 51
*      WRITE(99,*) "Y = "
*      WRITE(99,*) Y

* rapidity

         IF(ILOFLK (NPFLKA) .EQ. 13 .OR. ILOFLK (NPFLKA) .EQ. 14 .OR.
     &        ILOFLK (NPFLKA) .EQ. 15 .OR. ILOFLK (NPFLKA) .EQ. 16 .OR.
     &        ILOFLK (NPFLKA) .EQ.  1 .OR. ILOFLK (NPFLKA) .EQ. 8) THEN

*     Double Gaussian from NA49
 51         CONTINUE
*      R1 = -2. + 8 * FLRNDM(XDUMMY)
*     Not generating in full phase space!
            R1 = 1.5 + 3. * FLRNDM(XDUMMY)
*      R1 = 1.9. + 1.8 * FLRNDM(XDUMMY)
            IF(ILOFLK (NPFLKA) .EQ. 13 .OR. ILOFLK (NPFLKA) .EQ. 14)THEN
               Y0 = 0.666
               S0 = 0.872
               YMAX = 1.5
            END IF
*      IF(ILOFLK (NPFLKA) .EQ. 15 .OR. ILOFLK(NPFLKA) .EQ. 1) THEN
            IF(ILOFLK (NPFLKA) .EQ. 15) THEN
               Y0 = 0.694
               S0 = 0.725
               YMAX = 1.27
            END IF
            IF(ILOFLK (NPFLKA) .EQ. 1 .OR. ILOFLK (NPFLKA) .EQ. 8) THEN
               Y0 = 0.907
               S0 = 0.798
               YMAX = 1.1
            END IF
            IF(ILOFLK (NPFLKA) .EQ. 16) THEN
               Y0 = 0.569
               S0 = 0.635
               YMAX = 1.35
            END IF
            dNDY=EXP(-(R1-2.2-Y0)**2/(2.*S0**2)) + 
     &           EXP(-(R1-2.2+Y0)**2/(2.*S0**2))
            
            R2 = FLRNDM(XDUMMY)
            
            IF(R2 < DNDY/YMAX) THEN
               Y = R1
**      ELSE IF(R2 > PT_TEMP / MAX .OR. PT_TEMP > 3.) THEN
            ELSE
**         WRITE(*,*) "SONO QUI QUI"
               GOTO 51
            END IF
            
         END IF

         IF(ILOFLK (NPFLKA) .EQ. 10 .OR. ILOFLK (NPFLKA) .EQ. 11) THEN
* 52      CONTINUE
*         R1 = 1.5 + 3. * FLRNDM(XDUMMY)
*         S0 = 1.
*        dNDY=EXP(-(R1-1.9)**2/(2.*S0))
*        IF(R2 < dNDY) THEN
*           Y = R1
*        ELSE
*           GOTO 52
*        END IF
            Y = 1.5 + 3. * FLRNDM(XDUMMY)
         END IF

* Transverse momentum

         IF(ILOFLK (NPFLKA) .EQ. 13 .OR. ILOFLK (NPFLKA) .EQ. 14 .OR.
     &        ILOFLK (NPFLKA) .EQ. 15 .OR. ILOFLK (NPFLKA) .EQ. 16 .OR.
     &        ILOFLK (NPFLKA) .EQ.  1 .OR. ILOFLK (NPFLKA) .EQ. 8) THEN

* Sample transverse momentum of the particle
* NA49 parameters

            IF(ILOFLK (NPFLKA) .EQ. 13 .OR. ILOFLK (NPFLKA) .EQ. 14)THEN
               T = 0.170
               PTMAX = 0.049
            END IF
            IF(ILOFLK (NPFLKA) .EQ. 15 .OR. ILOFLK(NPFLKA) .EQ. 1 .OR.
     &           ILOFLK(NPFLKA)  .EQ. 8) THEN
               T = 0.232
               PTMAX = 0.074
            END IF
            IF(ILOFLK(NPFLKA) .EQ. 1 .OR.
     &           ILOFLK(NPFLKA)  .EQ. 8) THEN
               T = 0.257
               PTMAX = 0.083
            END IF
            IF(ILOFLK (NPFLKA) .EQ. 16) THEN
               T = 0.226
               PTMAX = 0.071
            END IF
            
 50         CONTINUE
            R1 = 3. * FLRNDM(XDUMMY)
            R2 = FLRNDM(XDUMMY)
            
            PT_TEMP = R1 * EXP(-SQRT(R1**2 + AM (IONID)**2)/T) / PTMAX
*      IF(R2 < R1) THEN
            IF(R2 < PT_TEMP) THEN
               Pt =R1
*     ELSE IF(R2 > PT_TEMP / MAX .OR. PT_TEMP > 3.) THEN
            ELSE
*         WRITE(*,*) "SONO QUI QUI"
               GOTO 50
            END IF
*     Pt = 0.5
*      WRITE(*,*) "SONO QUI"
*      WRITE(*,*) PT

*      WRITE(70,*) AM(IONID)

         END IF
         
         IF(ILOFLK (NPFLKA) .EQ. 10 .OR. ILOFLK (NPFLKA) .EQ. 11) THEN
            Pt = 3. * FLRNDM(XDUMMY)
         END IF
         
*     Azimuthal angle      
         
         Phi = FLRNDM(XDUMMY) * 2 * 3.14159
*     Components of momentum
         CALL SFECFE(SINT,COST)
         Px = Pt * COST
         Py = Pt * SINT
         Pz = SQRT(Pt**2 + AM (IONID)**2 ) * SINH(Y)
*     WRITE(*,*) Pz

*      Pz = SQRT(Pt**2 + AM (IONID)**2 ) * SINH(1.9)

*     Particle momentum
         PMOFLK (NPFLKA) = SQRT(Px**2 + Py**2 + Pz**2)
*     PMOFLK (NPFLKA) = 20.
         
*     Kinetic energy of the particle (GeV)
         TKEFLK (NPFLKA) =SQRT(PMOFLK (NPFLKA)**2+AM (IONID)**2)-
     &        AM (IONID)
*     Particle momentum
*     PMOFLK (NPFLKA) = PBEAM
*     PMOFLK (NPFLKA) = SQRT ( TKEFLK (NPFLKA) * ( TKEFLK (NPFLKA)
*    &                       + TWOTWO * AM (IONID) ) )
*  Cosines (tx,ty,tz)
         TNORM = PMOFLK (NPFLKA)
         
         TXFLK  (NPFLKA) = Px / TNORM
         TYFLK  (NPFLKA) = Py / TNORM
         TZFLK  (NPFLKA) = Pz / TNORM
         
*      TXFLK  (NPFLKA) = UBEAM
*      TYFLK  (NPFLKA) = VBEAM
*      TZFLK  (NPFLKA) = WBEAM
*      TNORM = SQRT(TXFLK(NPFLKA)**2+TYFLK(NPFLKA)**2+TZFLK(NPFLKA)**2)
*
*      TXFLK  (NPFLKA) = TXFLK  (NPFLKA) / TNORM
*      TYFLK  (NPFLKA) = TYFLK  (NPFLKA) / TNORM
*      TZFLK  (NPFLKA) = TZFLK  (NPFLKA) / TNORM

*     TZFLK  (NPFLKA) = SQRT ( ONEONE - TXFLK (NPFLKA)**2
*    &                       - TYFLK (NPFLKA)**2 )
*  Polarization cosines:
         TXPOL  (NPFLKA) = -TWOTWO
         TYPOL  (NPFLKA) = +ZERZER
         TZPOL  (NPFLKA) = +ZERZER
*  Particle coordinates
         XFLK   (NPFLKA) = XBEAM
         YFLK   (NPFLKA) = YBEAM
         ZFLK   (NPFLKA) = ZBEAM
*     Calculate the total kinetic energy of the primaries: don't change
         IF ( ILOFLK (NPFLKA) .EQ. -2 .OR. ILOFLK (NPFLKA) .GT. 100000 )
     &        THEN
            TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
         ELSE IF ( ILOFLK (NPFLKA) .NE. 0 ) THEN
            TKESUM = TKESUM + (TKEFLK (NPFLKA) + AMDISC(ILOFLK(NPFLKA)))
     &           * WTFLK (NPFLKA)
         ELSE
            TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
         END IF
         RADDLY (NPFLKA) = ZERZER
*  Here we ask for the region number of the hitting point.
*     NREG (NPFLKA) = ...
*  The following line makes the starting region search much more
*  robust if particles are starting very close to a boundary:
         CALL GEOCRS ( TXFLK (NPFLKA), TYFLK (NPFLKA), TZFLK (NPFLKA) )
         CALL GEOREG ( XFLK  (NPFLKA), YFLK  (NPFLKA), ZFLK  (NPFLKA),
     &        NRGFLK(NPFLKA), IDISC )
*  Do not change these cards:
         CALL GEOHSM ( NHSPNT (NPFLKA), 1, -11, MLATTC )
         NLATTC (NPFLKA) = MLATTC
         CMPATH (NPFLKA) = ZERZER
         CALL SOEVSV
      endif
      RETURN
*=== End of subroutine Source =========================================*
      END



c-------------------------------------------------------------------
      SUBROUTINE SIGNAL_ROOTGEN ( NOMORE, JPSI_FLAG)
*======================================================================*
*                                                                      *
*            Call GenMUONLMR for signal generation                     *
*                                                                      *
*======================================================================*

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
      INCLUDE '(BEAMCM)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(IOIOCM)'
      INCLUDE '(LTCLCM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(SOURCM)'
      INCLUDE '(SUMCOU)'
*
      LOGICAL LFIRST_JPSI
*
      SAVE LFIRST_JPSI
      DATA LFIRST_JPSI / .TRUE. /

      integer jpsi_flag
      integer NENE, ijmuon_plus, ijmuon_minus
      parameter (NENE=7)
      parameter (ijmuon_plus = 10, ijmuon_minus=11)
c      double precision amass_j, amass_muon, pmuon
c      parameter ( amass_j = 3.096, pmuon = 1.5444, amass_muon=0.105 )
      
c      double precision par1_pt(NENE), par2_pt(NENE), par3_pt(NENE)
c      double precision plab(NENE)
c      data par1_pt /11.03, 11.03, 3.86, 4.89, 3.73, 13.96, 5.34 /
c      data par2_pt /429.7, 429.7, 10.6, 17.3, 8.83, 201.1, 18.8 /
c      data par3_pt /3.27, 3.27, 2.92, 2.80, 2.84, 2.58, 2.70 / 
c      data plab /40., 50., 70., 90., 110., 130., 150. /

c      double precision par_pt(3)
c      double precision xmin_pt,xmax_pt,ymin_pt,ymax_pt
c      save par_pt, xmin_pt,xmax_pt,ymin_pt,ymax_pt

c      double precision par_y(4)
c      save par_y
c      double precision pl, el
c      double precision ylab, sqrts
c      double precision ay,by
c      double precision ymin,ymax

c      double precision pt_J, y_j, px_j, py_j, pz_j
c      double precision pmuon_plus(4), pmuon_minus(4), pjpsi(4)
c      double precision pmuon_plus_lab(4), pmuon_minus_lab(4)
c      double precision x_sam, y_sam, fun_pt, sint, cost
      
      double precision px1, py1, pz1, px2, py2, pz2

c      integer Nsam_max
c      parameter (nsam_max=100)
c
*  +-------------------------------------------------------------------*
      NOMORE = 0

c    pushing the 2 muons in the fluka stack
*
      call dimugenlmr(px1,py1,pz1,px2,py2,pz2)

      do imuon = 1,2
      
         NPFLKA = NPFLKA + 1
*     Wt is the weight of the particle
         WTFLK  (NPFLKA) = ONEONE
         WEIPRI = WEIPRI + WTFLK (NPFLKA)

*     Particle type

         if(imuon.eq.1) then
            IONID = ijmuon_plus
            ILOFLK (NPFLKA) = ijmuon_plus
         else
            IONID = ijmuon_minus
            ILOFLK (NPFLKA) = ijmuon_minus
         endif
         
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
*  |  Group number for "low" energy neutrons, set to 0 anyway
         IGROUP (NPFLKA) = 0
*     |

*  From this point .....
*  Particle generation (1 for primaries)
         LOFLK  (NPFLKA) = 1
*     User dependent flag:
         LOUSE  (NPFLKA) = 0
*     No channeling:
         LCHFLK (NPFLKA) = .FALSE.
         DCHFLK (NPFLKA) = ZERZER
*     User dependent spare variables:
         DO 100 ISPR = 1, MKBMX1
            SPAREK (ISPR,NPFLKA) = ZERZER
 100     CONTINUE
*     User dependent spare flags:
         DO 200 ISPR = 1, MKBMX2
            ISPARK (ISPR,NPFLKA) = 0
 200     CONTINUE
*  Save the track number of the stack particle:
         ISPARK (MKBMX2,NPFLKA) = NPFLKA
         NPARMA = NPARMA + 1
         NUMPAR (NPFLKA) = NPARMA
         NEVENT (NPFLKA) = 0
         DFNEAR (NPFLKA) = +ZERZER
*     ... to this point: don't change anything
         
*  Particle age (s)
         AGESTK (NPFLKA) = +ZERZER
         AKNSHR (NPFLKA) = -TWOTWO
c
c     muon kinematic: kinetic energies, momentum and cosines
c
         if(imuon.eq.1) then
            TKEFLK (NPFLKA) = sqrt(px1**2+py1**2+pz1**2
     &      +amass_muon**2) - amass_muon
            PMOFLK (NPFLKA) = sqrt(px1**2+py1**2+pz1**2)
            TXFLK  (NPFLKA) = px1/PMOFLK (NPFLKA)
            TYFLK  (NPFLKA) = py1/PMOFLK (NPFLKA)
         else
            TKEFLK (NPFLKA) = sqrt(px2**2+py2**2+pz2**2
     &      +amass_muon**2) - amass_muon
            PMOFLK (NPFLKA) = sqrt(px2**2+py2**2+pz2**2)
            TXFLK  (NPFLKA) = px2/PMOFLK (NPFLKA)
            TYFLK  (NPFLKA) = py2/PMOFLK (NPFLKA)
         endif
         TZFLK  (NPFLKA) = SQRT ( ONEONE - TXFLK (NPFLKA)**2
     &                       - TYFLK (NPFLKA)**2 )
         write (15,*) 'p, cosines', PMOFLK (NPFLKA), 
     &   TXFLK (NPFLKA), TYFLK (NPFLKA), TZFLK (NPFLKA)

*  Polarization cosines:
         TXPOL  (NPFLKA) = -TWOTWO
         TYPOL  (NPFLKA) = +ZERZER
         TZPOL  (NPFLKA) = +ZERZER
*  Particle coordinates

         XFLK   (NPFLKA) = XBEAM 
         YFLK   (NPFLKA) = YBEAM 
         ZFLK   (NPFLKA) = ZBEAM
         
*  Calculate the total kinetic energy of the primaries: don't change
         IF ( ILOFLK (NPFLKA) .EQ. -2 .OR. ILOFLK (NPFLKA) .GT. 100000 )
     &        THEN
            TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
         ELSE IF ( ILOFLK (NPFLKA) .NE. 0 ) THEN
            TKESUM = TKESUM + (TKEFLK(NPFLKA) + AMDISC(ILOFLK(NPFLKA)) )
     &           * WTFLK (NPFLKA)
         ELSE
            TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
         END IF
         RADDLY (NPFLKA) = ZERZER
*  Here we ask for the region number of the hitting point.
*     NREG (NPFLKA) = ...
*  The following line makes the starting region search much more
*  robust if particles are starting very close to a boundary:
         CALL GEOCRS ( TXFLK (NPFLKA), TYFLK (NPFLKA), TZFLK (NPFLKA) )
         CALL GEOREG ( XFLK  (NPFLKA), YFLK  (NPFLKA), ZFLK  (NPFLKA),
     &        NRGFLK(NPFLKA), IDISC )
*  Do not change these cards:
         CALL GEOHSM ( NHSPNT (NPFLKA), 1, -11, MLATTC )
         NLATTC (NPFLKA) = MLATTC
         CMPATH (NPFLKA) = ZERZER
         CALL SOEVSV
	 
	 write(15,*) 'imuon=',imuon, 'npflka=',npflka, 
     &   'particle code =',ILOFLK(NPFLKA),
     &   'p, px,py,pz cosines=', PMOFLK (NPFLKA), TXFLK (NPFLKA), 
     &                           TYFLK (NPFLKA), TZFLK (NPFLKA)

      end do
      
      RETURN

      END

