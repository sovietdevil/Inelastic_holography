********************************************************************
        FUNCTION FSCATT (G,UL,Z,SYMBOL,ACCVLT,ABSFLG,ACCFLG,DWFLG)

*       modified: 07.08.90
*
*       Author:
*       A.Weickenmeier
*       Technische Hochschule Darmstadt
*       Hochschulstr. 6
*       D-6100 Darmstadt
*       Germany
*       bitnet: XLTODA6L@DDATHD21.BITNET
*       Tel.: XX-49-6151-16-3381
*
*       CALCULATES THE COMPLEX SCATTERING AMPLITUDE
*
*
*       INPUT:
*       ======
*
*       G:      SCATTERING VECTOR 4*PI*G = S = SIN(THETA)/LAMBDA
*       UL:     THE RMS THERMAL DISPLACEMENT OF THE ATOM
*       Z:      ATOM
*       ACCVLT: ACCELERATION VOLTAGE IN KV
*                       NOT REQUIRED IF ACCFLG = .FALSE.
*       ABSFLG: IF TRUE IMAGINARY PART OF FSCATT WILL BE CALCULATED
*       ACCFLG: IF TRUE THEN F DEPENDS ON THE ACCELERATION VOLTAGE
*               IF FALSE REAL(F) MUST BE MULTIPLIED BY GAMMA,
*                        IMAG(F) BY GAMMA^2/K0
*       DWFLG:  IF TRUE REAL(F) WILL BE MULTIPLIED WITH THE
*                       DEBYE-WALLER FACTOR
*
*       OUTPUT:
*       =======
*
*       FSCATT: COMPLEX SCATTERING AMPLITUDE
*       SYMBOL: SYMBOL OF THE ELEMENT Z
*
*       UNITS:
*       ======
*
*       UNIT OF LENGTH IS ANGSTROEM
*
*       NOTE:
*       =====
*
*       THE SCATTERING AMPLITUDES ARE MULTIPLIED BY AN ADDITIONAL
*       FACTOR 4PI.
*
***************************************************************************

        PARAMETER (FOURPI=4.*3.1415927)

        LOGICAL         ABSFLG
        LOGICAL         ACCFLG
        LOGICAL         DWFLG
        INTEGER         Z
        REAL            K0
        REAL            A(2)
        REAL            B(6)
        CHARACTER*2     SYMBOL
        COMPLEX         FSCATT

C       CHECK INPUT
        IF (Z .LT. 2) THEN
           WRITE (6,1000) Z
           STOP
        END IF
        IF (UL .LT. 0.) THEN
           WRITE (6,1010) UL
           STOP
        END IF
        IF (G .LT. 0.) THEN
           WRITE (6,1020) G
           STOP
        END IF
        IF (ACCVLT .LT. 0.) THEN
           WRITE (6,1030) ACCVLT
           STOP
        END IF

C       GET FITTING COEFFICIENTS
        CALL GETWK (Z,SYMBOL,A,B)

        S      = G / FOURPI
        FREAL  = FOURPI * WEKO (A,B,S)

        IF (ABSFLG) THEN
           FIMA   = FIMAG (G,UL,A,B)
        ELSE
           FIMA   = 0.
        END IF
        IF (ACCFLG) THEN
C          CALCULATE WAVENUMBER AND GAMMA
           K0   = .5068 * SQRT(1022.*ACCVLT + ACCVLT*ACCVLT)
           GAMMA = (ACCVLT+511.) / 511.
           FREAL = FREAL * GAMMA
           FIMA  = FIMA  * GAMMA * GAMMA / K0
        END IF
        IF (DWFLG) THEN
C          CALCULATE DEBYE-WALLER FACTOR
           DEWA  = EXP(-.5*UL*UL*G*G)
           FREAL = FREAL * DEWA
        END IF

        FSCATT = CMPLX( FREAL,FIMA )

        RETURN
1000    FORMAT (2X, 'Z = ',I3,' THIS VALUE IS NOT SUPPLIED!')
1010    FORMAT (2X, 'UL = ',E13.4,' THIS VALUE IS NOT SUPPLIED!')
1020    FORMAT (2X, 'G = ',E13.4,' THIS VALUE IS NOT SUPPLIED!')
1030    FORMAT (2X, 'ACCVLT = ',E13.4,' THIS VALUE IS NOT SUPPLIED!')
        END

********************************************************************
        FUNCTION WEKO (A,B,S)

C       ELECTRON SCATTERING APMLITUDE F(S)

        REAL A (2)
        REAL B (6)

        WEKO=0.
        IF (S .GT. 0.) THEN
           S2 = 1./(S*S)
        END IF

        DO 10 I=1,6
           J = 1+(I-1)/3
           ARGU = B(I)*S*S
           IF (ARGU .LT. .1) THEN
              WEKO = WEKO + A(J)*B(I) * (1.-.5*ARGU)
           ELSE IF (ARGU .GT. 20.) THEN
              WEKO = WEKO + A(J)*S2
           ELSE
              WEKO = WEKO + A(J)*(1.-EXP(-ARGU))*S2
           END IF
10      CONTINUE

        RETURN
        END

C ***************************************************************
        FUNCTION FIMAG (G,UL,A,B)

        PARAMETER (FOURPI = 12.56636)
        PARAMETER (FP2    = FOURPI*FOURPI)

        INTEGER Z
        REAL A  (2)
        REAL B  (6)
        REAL A1 (2)
        REAL B1 (6)

        U2 = UL*UL

        DO 10 I=1,2
           A1(I) = A(I) * FP2
10      CONTINUE
        DO 20 I=1,6
           B1(I) = B(I) / FP2
20      CONTINUE

        FIMAG = 0.
        G2    = G*G
        DEWA  = EXP(-.5*U2*G2)

        DO 40 J=1,6
           JJ = 1+(J-1)/3
           DO 30 I=1,6
              II = 1+(I-1)/3
              FIMAG = FIMAG + A1(JJ)*A1(II)*(  DEWA * RI1(B1(I),B1(J),G)
     1                                       - RI2(B1(I),B1(J),G,UL)  )
30         CONTINUE
40      CONTINUE

        RETURN
        END

***************************************************************************

        FUNCTION RI1 (BI,BJ,G)

C       ERSTES INTEGRAL FUER DIE ABSORPTIONSPOTENTIALE

        PARAMETER (PI=3.1415927)
        PARAMETER (C =0.5772157)

        IF (G .EQ. 0.) THEN
           RI1 = BI * LOG( (BI+BJ)/BI ) + BJ * LOG( (BI+BJ)/BJ )
           RI1 = RI1 * PI
           RETURN
        END IF

        G2   = G*G
        BIG2 = BI*G2
        BJG2 = BJ*G2

        RI1 = 2.*C + LOG(BIG2) + LOG(BJG2) - 2.*EI( -BI*BJ*G2/(BI+BJ) )
        X1  = BIG2
        X2  = BIG2*BI/(BI+BJ)
        X3  = BIG2
        RI1 = RI1 + RIH1(X1,X2,X3)

        X1  = BJG2
        X2  = BJG2*BJ/(BI+BJ)
        X3  = BJG2
        RI1 = RI1 + RIH1(X1,X2,X3)

        RI1 = RI1 * PI / G2

        RETURN
        END

***************************************************************************

        FUNCTION RI2 (BI,BJ,G,U)

C       ZWEITES INTEGRAL FUER DIE ABSORPTIONSPOTENTIALE

        PARAMETER (PI=3.1415927)

        U2 = U*U

        IF (G .EQ. 0.) THEN
           RI2 = (BI+U2) * LOG( (BI+BJ+U2)/(BI+U2) )
           RI2 = RI2 + BJ * LOG( (BI+BJ+U2)/(BJ+U2) )
           RI2 = RI2 + U2 * LOG( U2/(BJ+U2) )
           RI2 = RI2 * PI
           RETURN
        END IF

        G2   = G*G
        BIUH = BI + .5*U2
        BJUH = BJ + .5*U2
        BIU  = BI + U2
        BJU  = BJ + U2

        RI2 = EI( -.5*U2*G2*BIUH/BIU ) + EI( -.5*U2*G2*BJUH/BJU )
        RI2 = RI2 - EI( -BIUH*BJUH*G2/(BIUH+BJUH) ) - EI( -.25*U2*G2 )
        RI2 = 2.*RI2
        X1  = .5*U2*G2
        X2  = .25*U2*G2
        X3  = .25*U2*U2*G2/BIU
        RI2 = RI2 + RIH1(X1,X2,X3)

        X1  = .5*U2*G2
        X2  = .25*U2*G2
        X3  = .25*U2*U2*G2/BJU
        RI2 = RI2 + RIH1(X1,X2,X3)

        X1  = BIUH*G2
        X2  = BIUH*BIUH*G2/(BIUH+BJUH)
        X3  = BIUH*BIUH*G2/BIU
        RI2 = RI2 + RIH1(X1,X2,X3)

        X1  = BJUH*G2
        X2  = BJUH*BJUH*G2/(BIUH+BJUH)
        X3  = BJUH*BJUH*G2/BJU
        RI2 = RI2 + RIH1(X1,X2,X3)

        RI2 = RI2 * PI / G2

        RETURN
        END

***************************************************************************

        FUNCTION RIH1 (X1,X2,X3)

C       WERTET DEN AUSDRUCK EXP(-X1) * ( EI(X2)-EI(X3) ) AUS

        IF (X2 .LE. 20.  .AND.  X3 .LE. 20.) THEN
           RIH1 = EXP(-X1) * ( EI(X2)-EI(X3) )
           RETURN
        END IF

        IF (X2 .GT. 20) THEN
           RIH1 = EXP(X2-X1)*RIH2(X2)/X2
        ELSE
           RIH1 = EXP(-X1)*EI(X2)
        END IF

        IF (X3 .GT. 20) THEN
           RIH1 = RIH1 - EXP(X3-X1)*RIH2(X3)/X3
        ELSE
           RIH1 = RIH1 - EXP(-X1)*EI(X3)
        END IF

        RETURN
        END

***************************************************************************

        FUNCTION RIH2 (X)

C       WERTET X*EXP(-X)*EI(X) AUS FUER GROSSE X
C       DURCH INTERPOLATION DER TABELLE ... AUS ABRAMOWITZ

        REAL F (0:20)

        DATA F / 1.000000,1.005051,1.010206,1.015472,1.020852,
     1           1.026355,1.031985,1.037751,1.043662,1.049726,
     2           1.055956,1.062364,1.068965,1.075780,1.082830,
     3           1.090140,1.097737,1.105647,1.113894,1.122497,
     4           1.131470 /

        X1 = 1./X
        I  = INT( 200.*X1 )
        I1 = I+1

        RIH2 = F(I) + 200.*( F(I1)-F(I) ) * ( X1-.5E-3*REAL(I) )

        RETURN
        END

C ***********************************************************
        FUNCTION EI (X)

C       EXPONENTIALINTEGRAL
C         GETESTET -60 < X < 60

        PARAMETER (A1=8.57332,A2=18.05901,A3=8.63476,A4=.26777)
        PARAMETER (B1=9.57332,B2=25.63295,B3=21.09965,B4=3.95849)

        IF (X .GT. 60.) THEN
           WRITE (6,*) '>>> EI FUER X= ',X,' NICHT GETESTET <<<'
           STOP
        END IF

        IF (X .LT. -60.) THEN
           EI = 0.
           RETURN
        END IF

        IF (X .LT. -1.) THEN
C          ABRAMOWITZ (5.1.56)
           XP = ABS(X)
           EI = -( A4+XP*(A3+XP*(A2+XP*(A1+XP))) ) /
     1           ( B4+XP*(B3+XP*(B2+XP*(B1+XP))) ) * EXP(-XP)/XP
           RETURN
        END IF

        EI   = .577216 + LOG( ABS(X) )

        I    = 1
        SI   = X
        SUMM = SI

10      CONTINUE
        RI   = REAL(I)
        RI1  = RI + 1.
        SI   = SI * X * RI/(RI1*RI1)
        SUMM = SUMM + SI

        IF (ABS(SI/X) .GT. 1.E-6) THEN
           I = I+1
           GOTO 10
        ELSE
           EI = EI + SUMM
           RETURN
        END IF

        RETURN
        END

***************************************************************************
        SUBROUTINE GETWK (Z,SYMBOL,A,B)

C       DATEN VON UND FUER DIE AUFRUFENDE ROUTINE
        INTEGER Z
        REAL    A       (2)
        REAL    B       (6)
        CHARACTER*2     SYMBOL

C       DATEN, DIE NUR INTERN BENOETIGT WERDEN
        REAL    V       (98)
        REAL    BB      (6,98)
        CHARACTER*2 SY  (98)

        SAVE V,BB,SY

        DATA SY /'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     1  'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     2  'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     3  'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     4  'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     5  'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     6  'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     7  'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     8  'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     9  'Pa','U ','Np','Pu','Am','Cm','Bk','Cf'/

        DATA V /0.0,0.5,0.5,0.3,0.5,0.5,0.5,0.5,0.5,0.5,
     1    0.5,0.5,0.4,0.5,0.5,0.5,0.5,0.5,0.2,0.3,
     2    0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
     3    0.5,0.5,0.5,0.5,0.5,0.5,0.2,0.3,0.5,0.5,
     4    0.5,0.5,0.5,0.4,0.5,0.5,0.5,0.3,0.4,0.6,
     5    0.6,0.6,0.4,0.4,0.1,0.1,0.3,0.3,0.2,0.2,
     6    0.2,0.2,0.1,0.2,0.1,0.2,0.1,0.2,0.1,0.1,
     7    0.1,0.1,0.4,0.2,0.5,0.4,0.5,0.5,0.4,0.4,
     8    0.4,0.3,0.4,0.4,0.4,0.4,0.1,0.2,0.2,0.3,
     9    0.2,0.2,0.2,0.2,0.2,0.3,0.2,0.3/
        DATA ((BB(I,Z),I=1,6),Z=1,10) /
     1      0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,
     2      2.54216,  8.74302, 12.69098,  0.43711,  5.29446, 28.25045,
     3      0.68454,  3.06497,  6.23974,126.17816,131.20160,131.76538,
     4      0.53996,  3.38752, 55.62340, 50.78098, 67.00502, 96.36635,
     5      0.33138,  2.97485, 34.01118, 35.98365, 36.68364, 60.80991,
     6      0.29458,  3.93381, 24.97836, 25.27916, 25.46696, 46.70328,
     7      0.23925,  4.93515, 18.11895, 15.69698, 15.81922, 40.24150,
     8      6.37582,  8.03744, 27.20649,  0.11157,  0.38686, 10.89944,
     9      0.21800,  6.76987,  7.05056,  6.67484, 12.38148, 28.08398,
     *      0.20055,  5.49814,  6.28052,  7.19211,  7.54763, 23.26388/

        DATA ((BB(I,Z),I=1,6),Z=11,20) /
     1      0.21902,  5.30022,  5.31938,  5.28281,  5.28546,128.18391,
     2      1.97633,  2.80902, 16.39184,  0.05494,  2.06121,121.70512,
     3      2.29692,  2.35822, 24.98576,  0.07462,  0.55953,128.50104,
     4      1.73656,  3.04329, 30.57191,  0.05070,  0.99181, 86.18340,
     5      0.17949,  2.63250,  2.67559, 34.57098, 36.77888, 54.06180,
     6      1.00609,  4.90414, 31.34909,  0.03699,  0.98700, 44.94354,
     7      0.18464,  1.47963,  5.20989, 24.79470, 32.06184, 39.09933,
     8      0.20060,  6.53262, 22.72092,  1.20022,  1.27398, 36.25907,
     9      0.44424,  3.36735, 19.63031,  0.01824, 23.51332,212.86819,
     *      0.18274,  2.06638, 16.99062, 11.57795, 13.97594,186.10446/
        DATA ((BB(I,Z),I=1,6),Z=21,30) /
     1      0.14245,  1.46588, 15.46955,  4.24287,  9.80399,121.46864,
     2      0.12782,  1.45591, 12.09738,  4.61747, 11.96791,105.00546,
     3      0.13126,  1.39923,  8.00762,  7.98129, 13.41408, 95.30811,
     4      0.12311,  2.38386,  9.92149,  1.64793, 11.00035, 68.45583,
     5      0.48173,  3.78306,  8.47337,  0.04690,  8.74544, 77.44405,
     6      0.44704,  6.89364,  6.90335,  0.05691,  3.02647, 70.86599,
     7      0.10705,  3.63573,  7.55825,  1.27986,  5.14045, 67.16051,
     8      0.11069,  1.61889,  6.00325,  5.97496,  6.06049, 59.41419,
     9      0.11293,  1.89077,  5.08503,  5.07335,  5.09928, 46.38955,
     *      0.10209,  1.73365,  4.78298,  4.80706,  5.64485, 51.21828/
        DATA ((BB(I,Z),I=1,6),Z=31,40) /
     1      0.10642,  1.53735,  5.13798,  4.74298,  4.99974, 61.42872,
     2      0.09583,  1.67715,  4.70275,  2.91198,  7.87009, 64.93623,
     3      0.09428,  2.21409,  3.95060,  1.52064, 15.81446, 52.41380,
     4      0.09252,  1.60168,  3.04917,  3.18476, 18.93890, 47.62742,
     5      0.09246,  1.77298,  3.48134,  1.88354, 22.68630, 40.69434,
     6      0.49321,  2.08254, 11.41282,  0.03333,  2.09673, 42.38068,
     7      0.15796,  1.71505,  9.39164,  1.67464, 23.58663,152.53635,
     8      0.36052,  2.12757, 12.45815,  0.01526,  2.10824,133.17088,
     9      0.09003,  1.41396,  2.05348, 10.25766, 10.74831, 90.63555,
     *      0.10094,  1.15419,  2.34669, 10.58145, 10.94962, 82.82259/
        DATA ((BB(I,Z),I=1,6),Z=41,50) /
     1      0.09243,  1.16977,  5.93969,  1.30554, 13.43475, 66.37486,
     2      0.43543,  1.24830,  7.45369,  0.03543,  9.91366, 61.72203,
     3      0.45943,  1.18155,  8.31728,  0.03226,  8.32296, 64.97874,
     4      0.08603,  1.39552, 11.69728,  1.39552,  3.45200, 55.55519,
     5      0.09214,  1.11341,  7.65767,  1.12566,  8.32517, 48.38017,
     6      0.09005,  1.12460,  9.69801,  1.08539,  5.70912, 33.48585,
     7      0.08938,  3.19060,  9.10000,  0.80898,  0.81439, 41.34453,
     8      0.28851,  1.61312,  8.99691,  0.01711,  9.46666, 58.13256,
     9      0.08948,  1.23258,  8.23129,  1.22390,  7.06201, 59.69622,
     *      0.07124,  0.85532,  6.40081,  1.33637,  6.38240, 50.92361/
        DATA ((BB(I,Z),I=1,6),Z=51,60) /
     1      0.35749,  1.32481,  6.51696,  0.03550,  6.51913, 50.80984,
     2      0.50089,  3.95301,  7.62830,  0.03005,  0.50737, 49.62628,
     3      0.08429,  1.12959,  8.86209,  1.12981,  9.13243, 56.01965,
     4      0.27796,  1.62147, 11.45200,  0.02032,  3.27497, 51.44078,
     5      0.12045,  1.53654,  9.81569, 41.21656, 42.62216,224.34816,
     6      0.12230,  1.44909,  9.50159, 49.40860, 74.94942,217.04485,
     7      0.08930,  1.26225,  8.09703,  1.20293, 17.65554,116.61481,
     8      0.08504,  1.28286, 11.22123,  1.32741,  4.61040,112.19678,
     9      0.09805,  1.52628,  8.58953,  1.23893, 22.49126,140.02856,
     *      0.09413,  1.26616,  5.98844, 17.78775, 18.14397,132.59305/
        DATA ((BB(I,Z),I=1,6),Z=61,70) /
     1      0.09447,  1.25111,  5.91205, 16.28675, 16.73089,127.90916,
     2      0.09061,  1.59281, 10.64077,  1.78861,  2.22148,124.56328,
     3      0.10485,  1.54396,  8.65223,  7.09290, 53.36537,183.69014,
     4      0.09338,  1.38681,  7.35883,  1.55122, 20.81916,111.03201,
     5      0.10190,  1.52368,  7.16923, 20.86269, 49.29465,166.09206,
     6      0.08402,  1.40890,  7.14042,  1.34848, 11.42203,108.01204,
     7      0.09441,  1.61807,  6.27142, 40.34946, 42.82722,130.59616,
     8      0.08211,  1.25106,  4.81241, 10.84493, 10.90164,100.07855,
     9      0.09662,  1.60236,  5.67480, 30.59014, 31.12732,138.69682,
     *      0.09493,  1.60220,  5.43916, 28.31076, 29.27660,138.08665/
        DATA ((BB(I,Z),I=1,6),Z=71,80) /
     1      0.09658,  1.56751,  5.32170, 34.18217, 35.25187,121.42893,
     2      0.09294,  1.55499,  5.25121, 37.51883, 38.88302,105.16978,
     3      0.06298,  0.81950,  2.89124,  5.54290,  5.98101, 54.42459,
     4      0.07902,  1.37096,  8.23364,  1.38300,  1.39219, 77.11813,
     5      0.05266,  0.90718,  4.43830,  0.94590,  4.37477, 43.97909,
     6      0.22700,  1.56975,  6.34451,  0.01564,  1.61769, 46.15815,
     7      0.05055,  0.86775,  5.09325,  0.88123,  3.56919, 39.77390,
     8      0.05253,  0.83773,  3.95899,  0.81515,  6.44217, 34.21146,
     9      0.54927,  1.72752,  6.71952,  0.02637,  0.07253, 35.45745,
     *      0.21941,  1.41611,  6.68241,  0.01472,  1.57578, 37.15826/
        DATA ((BB(I,Z),I=1,6),Z=81,90) /
     1      0.22459,  1.12822,  4.30289,  0.01485,  7.15607, 43.08737,
     2      0.06432,  1.19406,  7.39342,  1.14160,  1.28905, 51.13401,
     3      0.05380,  0.86719,  1.87540,  7.64796,  7.86794, 45.63897,
     4      0.50112,  1.63784,  6.78551,  0.02187,  0.08602, 46.72951,
     5      0.22321,  1.10827,  3.59116,  0.01011, 11.63732, 45.06839,
     6      0.21152,  1.14015,  3.41473,  0.01188, 13.41211, 43.11389,
     7      0.09435,  1.02649,  6.25480, 32.51444, 36.29119,149.11722,
     8      0.07300,  1.01825,  5.89629,  1.03089, 20.37389,115.34722,
     9      0.07515,  0.94941,  3.72527, 17.58346, 19.75388,109.12856,
     *      0.06385,  0.90194,  4.65715,  0.90253, 15.70771, 83.69695/
        DATA ((BB(I,Z),I=1,6),Z=91,98) /
     1      0.07557,  0.84920,  4.00991, 16.95003, 17.78767,100.20415,
     2      0.07142,  1.14907,  9.21231,  0.95923,  1.20275,104.32746,
     3      0.06918,  0.98102,  5.95437,  0.99086, 22.06437, 90.98156,
     4      0.07136,  0.95772,  6.13183,  0.97438, 15.67499, 89.86625,
     5      0.07301,  0.93267,  6.34836,  0.91032, 13.26179, 86.85986,
     6      0.05778,  0.72273,  3.01146,  9.21882,  9.53410, 65.86810,
     7      0.07088,  0.77587,  6.14295,  1.79036, 15.12379, 83.56983,
     8      0.06164,  0.81363,  6.56165,  0.83805,  4.18914, 61.41408/

        SYMBOL = SY(Z)

        A(1) = 0.02395*REAL(Z) / (3.*(1.+V(Z)))
        A(2) = V(Z) * A(1)

        DO 10 I=1,6
           B(I) = BB(I,Z)
10      CONTINUE

        RETURN
        END
