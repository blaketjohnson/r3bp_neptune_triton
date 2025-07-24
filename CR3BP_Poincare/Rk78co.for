
      SUBROUTINE RK78CO

C
C-----
C
C PURPOSE:
C          THE SUBROUTINE RK78CO SETS UP FEHLBERG COEFFI-
C          CIENTS FOR THE NUMERICAL INTEGRATION ROUTINE
C          RKF78.
C
C INPUT:
C          NONE.
C
C OUTPUT:
C          COMMON/COEF78/
C          A0,...  FEHLBERG COEFFICIENTS A0,...,A12, B10,...
C                  B1211,CH0,...,CH12, E0,...,E12.
C           B,...  STEPSIZE CONTROL FACTORS, AND NECESSARY
C                  CONSTANTS USEFUL TO CALCULATE THEM.
C
C SUBCALLS:
C          NONE.
C
C
      IMPLICIT REAL*8(A-H,O-Z)

      COMMON/COEF78/A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,
     *              A12,B10,B20,B21,B30,B31,B32,B40,B41,B42,
     1              B43,B50,B51,B52,B53,B54,B60,B61,B62,B63,
     2              B64,B65,B70,B71,B72,B73,B74,B75,B76,B80,
     3              B81,B82,B83,B84,B85,B86,B87,B90,B91,B92,
     4              B93,B94,B95,B96,B97,B98,B100,B101,B102,
     5              B103,B104,B105,B106,B107,B108,B109,B110,
     6              B111,B112,B113,B114,B115,B116,B117,B118,
     7              B119,B1110,B120,B121,B122,B123,B124,
     8              B125,B126,B127,B128,B129,B1210,B1211,
     9              CH0,CH1,CH2,CH3,CH4,CH5,CH6,CH7,CH8,CH9,
     A              CH10,CH11,CH12,E0,E1,E2,E3,E4,E5,E6,E7,
     B              E8,E9,E10,E11,E12,B,BLO,BUP,REMIN,DTINC,
     C              DTDEC,MAXREJ
C
C-----
C
      MAXREJ = 10
      DTINC = 20.D0
      DTDEC = 0.025D0
      REMIN = 3.0D-15
      B = 0.85D0

C
C     SET COEFFICIENTS
C

      A0 = 0.D0
      A1 = 2.D0/27.D0
      A2 = 1.D0/9.D0
      A3 = 1.D0/6.D0
      A4 = 5.D0/12.D0
      A5 = 1.D0/2.D0
      A6 = 5.D0/6.D0
      A7 = 1.D0/6.D0
      A8 = 2.D0/3.D0
      A9 = 1.D0/3.D0
      A10 = 1.D0
      A11 = 0.D0
      A12 = 1.D0

      B10 = 2.D0/27.D0
      B20 = 1.D0/36.D0
      B21 = 1.D0/12.D0
      B30 = 1.D0/24.D0
      B31 = 0.D0
      B32 = 1.D0/8.D0
      B40 = 5.D0/12.D0
      B41 = 0.D0
      B42 = -25.D0/16.D0
      B43 = 25.D0/16.D0
      B50 = 1.D0/20.D0
      B51 = 0.D0
      B52 = 0.D0
      B53 = 1.D0/4.D0
      B54 = 1.D0/5.D0
      B60 = -25.D0/108.D0
      B61 = 0.D0
      B62 = 0.D0
      B63 = 125.D0/108.D0
      B64 = -65.D0/27.D0
      B65 = 125.D0/54.D0
      B70 = 31.D0/300.D0
      B71 = 0.D0
      B72 = 0.D0
      B73 = 0.D0
      B74 = 61.D0/225.D0
      B75 = -2.D0/9.D0
      B76 = 13.D0/900.D0
      B80 = 2.D0
      B81 = 0.D0
      B82 = 0.D0
      B83 = -53.D0/6.D0
      B84 = 704.D0/45.D0
      B85 = -107.D0/9.D0
      B86 = 67.D0/90.D0
      B87 = 3.D0
      B90 = -91.D0/108.D0
      B91 = 0.D0
      B92 = 0.D0
      B93 = 23.D0/108.D0
      B94 = -976.D0/135.D0
      B95 = 311.D0/54.D0
      B96 = -19.D0/60.D0
      B97 = 17.D0/6.D0
      B98 = -1.D0/12.D0
      B100 = 2383.D0/4100.D0
      B101 = 0.D0
      B102 = 0.D0
      B103 = -341.D0/164.D0
      B104 = 4496.D0/1025.D0
      B105 = -301.D0/82.D0
      B106 = 2133.D0/4100.D0
      B107 = 45.D0/82.D0
      B108 = 45.D0/164.D0
      B109 = 18.D0/41.D0
      B110 = 3.D0/205.D0
      B111 = 0.D0
      B112 = 0.D0
      B113 = 0.D0
      B114 = 0.D0
      B115 = -6.D0/41.D0
      B116 = -3.D0/205.D0
      B117 = -3.D0/41.D0
      B118 = 3.D0/41.D0
      B119 = 6.D0/41.D0
      B1110 = 0.D0
      B120 = -1777.D0/4100.D0
      B121 = 0.D0
      B122 = 0.D00
      B123 = -341.D0/164.D0
      B124 = 4496.D0/1025.D0
      B125 = -289.D0/82.D0
      B126 = 2193.D0/4100.D0
      B127 = 51.D0/82.D0
      B128 = 33.D0/164.D0
      B129 = 12.D0/41.D0
      B1210 = 0.D0
      B1211 = 1.D0

      CH0 = 0.D0
      CH1 = 0.D0
      CH2 = 0.D0
      CH3 = 0.D0
      CH4 = 0.D0
      CH5 = 34.D0/105.D0
      CH6 = 9.D0/35.D0
      CH7 = 9.D0/35.D0
      CH8 = 9.D0/280.D0
      CH9 = 9.D0/280.D0
      CH10 = 0.D0
      CH11 = 41.D0/840.D0
      CH12 = 41.D0/840.D0

      E0 = 41.D0/840.D0
      E1 = 0.D0
      E2 = 0.D0
      E3 = 0.D0
      E4 = 0.D0
      E5 = 0.D0
      E6 = 0.D0
      E7 = 0.D0
      E8 = 0.D0
      E9 = 0.D0
      E10 = 41.D0/840.D0
      E11 = -41.D0/840.D0
      E12 = -41.D0/840.D0

C
C     SET STEP SIZE CONTROL FACTORS
C

      BUP = (B/DTDEC)**8
      BLO = (B/DTINC)**8
      RETURN
      END
