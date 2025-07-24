! Purpose:

! This program is used to generate Poincaré Surface of Section (or Poincaré maps)
! of systems in the Circular Restricted Three-Body Problem (CRTBP).

! Author: Diogo Merguizo Sanchez, 2015
! Addes to git-based source control in May, 2020
!
! Based on orginal routines written by Roger Broucke
! Some references:
! * Broucke, R. A., Numerical Integration of Periodic Orbits in
!   the Main Problem of Artificial Satellite Theory, Celestial
!   Mechanics & Dynamical Astronomy, Volume 58, Issue 2, pp.99-123
!   DOI: 10.1007/BF00695787
! * Broucke, R. A., Periodic Orbits in the Restricted Three-Body Problem
!   With Earth-Moon Masses, NASA Technical Report 32-1168, 1968.


! The folowing module reads the initial parameters and writes them on the screen
! for user conference
! This module also is used to define some constants and global variables
! Importante note: Since the initial conditions a enclose inside a module,
! changes in the initial conditions may require recompilation of the main
! program.

      MODULE input_parameters
        IMPLICIT NONE
        CHARACTER (100), PROTECTED :: folder
        REAL*8, PROTECTED :: C0, CF, DC, mu, TLIM, DT1, XI, XF, x_step, &
     &                       ZI, ZF, z_step, TOL, mu_star, d_p, r_max,  &
     &                       r1_min, r2_min
! Parameters
        INTEGER, PARAMETER :: NEQN = 6  ! Number of differential equations to be integrated
        REAL*8, PARAMETER :: PI = 4.d0*DATAN(1.d0)
        REAL*8, PARAMETER :: RAD = PI/180.d0
        REAL*8, PARAMETER :: DEG = 180.d0/PI
        INTEGER :: traj, traj_step
        PUBLIC :: read_input, write_input
      CONTAINS
        SUBROUTINE read_input()
          IMPLICIT NONE
          OPEN(20, FILE = 'input.in')
            READ(20,*)folder
            READ(20,*)C0, CF, DC
            READ(20,*)mu
            mu_star = 1.d0 - mu
            READ(20,*)d_p
            READ(20,*)r_max
            READ(20,*)r1_min
            READ(20,*)r2_min
            READ(20,*)TLIM, DT1
            READ(20,*)XI, XF, x_step
            READ(20,*)ZI, ZF, z_step
            READ(20,*)TOL
            READ(20,*)traj
            READ(20,*)traj_step
          CLOSE(20)
        END SUBROUTINE

        SUBROUTINE write_input()
          IMPLICIT NONE
          PRINT *
          CALL SYSTEM ('date')
          PRINT *
          PRINT *, '-----------------Initial Parameters----------------'
          PRINT *, 'Folder: ', TRIM(folder)
          PRINT '(3(1x, A, E))', 'C0 =', C0, ', CF =', CF, ', DC =', DC
          PRINT *, 'mu =', mu
          PRINT *, 'Distance between the primaries, km =', d_p
          PRINT *, 'Distance to consider as escape, km =', r_max
          PRINT *, 'Minimum distance from M1, km =', r1_min
          PRINT *, 'Minimum distance from M2, km =', r2_min
          PRINT '(2(1x, A, E))', 'TLIM =', TLIM, ', DT1 =', DT1
          PRINT '(3(1x, A, E))', 'XI =', XI, ', XF =', XF, ', x_step =',&
     &                           x_step
          PRINT '(3(1x, A, E))', 'ZI =', ZI, ', ZF =', ZF, ', z_step =',&
     &          z_step
          PRINT *, 'Save trajectory =', traj
          IF (traj /= 0) THEN
            PRINT *, 'Saving trajectory every', traj_step, 'points'
          END IF
          PRINT *, 'Double check all these values.'
          PRINT *, '---------------------------------------------------'
          PRINT *
        END SUBROUTINE
      END MODULE

      MODULE GLOBAL_CTRL
        IMPLICIT NONE
        REAL*8 :: iemax
      END MODULE

! Procedure modules

      MODULE critical
        USE input_parameters
        IMPLICIT NONE
        PUBLIC critical_events
      CONTAINS
        SUBROUTINE critical_events(X, cond)
          IMPLICIT NONE
          REAL*8, DIMENSION(6) :: X
          REAL*8 :: R1, R2, R
          INTEGER :: cond

          R1 = DSQRT((mu + X(1))**2 + X(2)**2 + X(3)**2)
          R2 = DSQRT((X(1) - mu_star)**2 + X(2)**2 + X(3)**2)
          R = DSQRT(X(1)**2 + X(2)**2 + X(3)**2)

          cond = 0
          IF(R1 <= r1_min/d_p)THEN
            ! PRINT *, 'COLLISION with M1, R1 =', R1, 'n.d.; x_0 =',  x(1)
            cond = 1
          ELSE IF(R2 <= r2_min/d_p)THEN
            ! PRINT *, 'COLLISION with M2, R2 =', R2, 'n.d.; x_0 =',  x(1)
            cond = 2
          ELSE IF(R >= r_max/d_p)THEN
            ! PRINT *, 'ESCAPE from the system, R =', R, 'n.d.; x_0 =',  x(1)
            cond = 9
          END IF

          RETURN
        END SUBROUTINE
      END MODULE

      MODULE system_rotations
        IMPLICIT NONE
        PUBLIC :: FIXtoROT, ROTtoFIX
      CONTAINS
        SUBROUTINE FIXtoROT(XF, T, XR)
          IMPLICIT NONE
          REAL*8, DIMENSION (6) :: XF, XR
          REAL*8 :: T, ST, CT
          CT = DCOS(T)
          ST = DSIN(T)

          XR(1) = XF(1)*CT + XF(2)*ST
          XR(2) = XF(2)*CT - XF(1)*ST
          XR(3) = XF(3)

          XR(4) = XF(4)*CT + XF(2)*CT - XF(1)*ST + XF(5)*ST
          XR(5) =-XF(1)*CT + XF(5)*CT - XF(4)*ST - XF(2)*ST
          XR(6)=XF(6)
          RETURN
        END SUBROUTINE

        SUBROUTINE ROTtoFIX(XR, T, XF)
          IMPLICIT NONE
          REAL*8, DIMENSION (6) :: XF, XR
          REAL*8 :: T, ST, CT
          CT = DCOS(T)
          ST = DSIN(T)
          XF(1) = XR(1)*CT - XR(2)*ST
          XF(2) = XR(2)*CT + XR(1)*ST
          XF(3) = XR(3)
          XF(4) =-XR(1)*ST - XR(5)*ST + XR(4)*CT - XR(2)*CT
          XF(5) = XR(4)*ST - XR(2)*ST + XR(1)*CT + XR(5)*CT
          XF(6) = XR(6)
          RETURN
        END
      END MODULE

      INCLUDE 'Rkf78-V20180418.for'
      INCLUDE 'Rk78co.for'

      PROGRAM POINCARE
!
! Variable declaration:
!        WORK (W) = real work matrix (W = 14*NEQN)
!        IFLAG    = (input) parameter of routine initialization
!        IFLAG    = (output) code for the type of return
!        I,J      = counters
!        ABSERR   = absolute error matrix
!        RELERR   = relative error matrix
!        DT       = initial stepsize
!        DT1      = variable stepsize
!        mu       = MSN/(MP + MSN)
!        mu_star  = 1 - mu
!        MSN      = mass of the secondary body
!        MP       = mass of the primary body
!        TI       = initial time
!        TF       = final time
!        TOUT     = current time
!        X(NEQN)  = matrix of variables
!        NEQN     = number of diff equations to be integrated

      USE input_parameters
      USE critical
      USE GLOBAL_CTRL
      USE system_rotations

      IMPLICIT NONE
      REAL*8, DIMENSION(6):: X, ABSERR, RELERR, XR, X_inertial
      REAL*8, DIMENSION(84):: WORK
      REAL*8:: DT, HMIN, H1, EPS, TI0, TI, X1, Z1, R1,                  &
     &         R2, ARG, VY, XANT, YANT, ZANT, VXANT, VYANT, VZANT, XM,  &
     &         YM, ZM, VXM, VYM, VZM, TOUT, C1, CJ, aux
      INTEGER:: I, IFLAG, J, ii, jj, nc, nx, nz, ij, ik, il, nt, im, ni,&
     &          ifl, cond, iaux
      CHARACTER(4):: sfl
      CHARACTER(100):: fname01, fname02, fmt01, fmt02, path
      LOGICAL :: fstatus

!  	  DERIVS = subrotina corespondente a F
!       F  = nome da subrotina dentro de RKF78
!            que fornece o valor das derivadas
!       RFK78  = subrotina de integraηγo do RUNGE-KUTTA
!
        EXTERNAL DERIVS

        CALL read_input
        CALL write_input

! The results will be stored in a folder whose name comes from input.in
! This part of the program creates this folther if the folder does not exist.
! Warning: this feature will only work for Intel Fortran (ifort)
        INQUIRE(DIRECTORY = TRIM(folder), EXIST = fstatus)
        IF(.NOT.fstatus)THEN
          CALL SYSTEM ('mkdir '//TRIM(folder))
        ENDIF
        path = TRIM(folder)//'/'

! 01/10/2015: loop to generate PS for several values of CJ
      ni = NINT((CF - C0)/DC)

      DO ij = 0, ni
        CJ = C0 + DBLE(ij)*DC

        PRINT '(1x,A,1x,D11.4)', 'Running C =', CJ
        PRINT *

        aux = 1000.d0*CJ
        ifl = NINT(aux)
        WRITE(sfl,'(I4.4)')ifl

        fname01 = TRIM(path)//'poincY'//sfl//'.dat'
        fname02 = TRIM(path)//'trajY'//sfl//'.dat'

        ! PRINT '(1X,A,1X,A)', 'File name =', TRIM(fname01)
        OPEN(20, FILE = TRIM(fname01), STATUS = 'new')
        IF(traj /= 0)OPEN(22, FILE = TRIM(fname02))

        fmt01 = '(8(1x,E21.14),1x,E23.16)'
        fmt01 = TRIM(fmt01)
        fmt02 = '(9(1x,E21.14),1x,I)'
        fmt02 = TRIM(fmt02)

        ! WRITE(20,fmt02)0.d0, 0.d0, 0.d0, 0.d0, 0

        DT = DT1/100.d0
        ! NEQN = 6
        HMIN = 0.d0
        H1 = 0.01d0
        EPS = TOL
        TI0 = 0.d0
        TI = TI0
        ! mu_star = 1.0D0 - mu
        ii = 0
        jj = 0
!      PRINT *, XI, XF, x_step
!
! Set the integrator stepsize.
        DO J=1, NEQN
          RELERR(J) = TOL
          ABSERR(J) = TOL
        END DO

        nx = NINT((XF - XI)/x_step)
        nz = NINT((ZF - ZI)/z_step)

        DO ik = 0, nx
          X1 = XI + DBLE(ik)*x_step

          DO il = 0, nz
            Z1 = ZI + DBLE(il)*z_step

            iflag = 1
            iemax = 0

            X(1) = X1
            X(2) = 0.D0
            X(3) = Z1
            R1 = DSQRT((mu + X(1))**2 + X(2)**2 + X(3)**2)
            R2 = DSQRT((X(1) - mu_star)**2 + X(2)**2 + X(3)**2)
            ARG = -CJ + (X(1)**2 + X(2)**2)                             &
     &            + 2.d0*((1.D0 - mu)/R1 + mu/R2)
            IF(ARG < 0.d0)THEN
              ! PRINT *, 'ARG < 0.d0', X1, Z1
              CYCLE
            ENDIF
            VY = DSQRT(ARG)
            X(4) = 0.d0
            X(5) = VY
            X(6) = 0.d0
            ! PRINT *, VY
            TI = TI0

            nt = NINT((TLIM - TI)/DT1)

            DO im = 0, nt
              TOUT = TI0 + DBLE(im)*DT1

              print *, TOUT, X
              pause

              XANT = X(1)
              YANT = X(2)
              ZANT = X(3)
              VXANT = X(4)
              VYANT = X(5)
              VZANT = X(6)

              CALL critical_events(X, cond)
              ! PRINT *, 'cond =', cond
              IF (cond /= 0)iemax = 2222

          CALL RKF78(DERIVS,NEQN,X,TI,TOUT,RELERR,ABSERR,IFLAG,WORK,DT)
              ! CALL JACOBI(X, mu, C1)
              ! IF(TOUT == TI0) THEN
                ! WRITE(20,fmt02)X1, Z1, X, CJ, cond
                ! write(3,fmt01)X1,Z1,XANT,YANT,ZANT,VXANT,VYANT,VZANT,CJ
                ! PRINT fmt01, XANT, YANT, ZANT, VXANT, VYANT, VZANT, C1
              ! ENDIF
              IF(jj == traj_step)THEN
                XR = X
                IF(traj == 1) THEN
                  CALL ROTtoFIX(XR, TOUT, X_inertial)
                  WRITE(22, fmt01)X1, Z1, X_inertial
                ELSE IF (traj == 2) THEN
                  WRITE(22, fmt01)X1, Z1, XR
                END IF
                ! PRINT fmt01, x(1), x(2), x(3), x(4), x(5), x(6), CJ
                jj = 0
              ENDIF
              IF (traj /= 0) jj = jj + 1
              ! PRINT *, jj
              ! PRINT '(4(1x,E))', tout, X(1), X(2), (YANT*X(2))
              IF (((YANT*X(2)) .LT. 0.d0) .AND. (X(2) > 0.d0)) THEN
!                IF(ii == 10)THEN
                XM = (X(1) + XANT)/2.d0
                YM = (X(2) + YANT)/2.d0
                ZM = (X(3) + ZANT)/2.d0
                VXM = (X(4) + VXANT)/2.d0
                VYM = (X(5) + VYANT)/2.d0
                VZM = (X(6) + VZANT)/2.d0

                WRITE(20,fmt02)X1,Z1,XM,YM,ZM,VXM,VYM,VZM,CJ,cond

                ! PRINT fmt01, XM, YM, ZM, VXM, VYM, VZM
!                ii = 0
!               ENDIF
!               ii = ii + 1
              ENDIF
            ENDDO ! DO for TOUT
          ENDDO ! DO for Z1
        ENDDO ! DO for X1
        CLOSE(20)
      ENDDO !DO for CJ

      CALL SYSTEM ('date')

      STOP
      END
!
      SUBROUTINE DERIVS (T, X, DX)

        USE input_parameters

        IMPLICIT NONE
        REAL*8, DIMENSION(6):: X, DX
        REAL*8:: R1, R2, T, PR1, PR2
!
        R1 = (X(1) + mu)**2 + X(2)**2 + X(3)**2
        PR1 = R1**1.5d0
        R2 = (X(1) - mu_star)**2 + X(2)**2 + X(3)**2
        PR2 = R2**1.5d0

        DX(1) = X(4)
        DX(2) = X(5)
        DX(3) = X(6)
        DX(4) = 2.d0*X(5) + X(1) - (mu_star*(X(1) + mu))/PR1            &
     &         -(mu*(X(1) - mu_star))/PR2
        DX(5) =-2.d0*X(4) + X(2) - (mu_star*X(2))/PR1 - (mu*X(2))/PR2
        DX(6) =-(mu_star*X(3))/PR1-(mu*X(3))/PR2

        RETURN
      END
!
      SUBROUTINE TWOB(U,XR,EN,CANG,T,CTOT,AINCL)
!
!	 CANG is the z-component of the angular momentum
!
        USE system_rotations

        IMPLICIT NONE
        REAL*8, DIMENSION(6):: XR, XF
        REAL*8:: R1, R2, U, G, EN, CANG, T, CTOT, AINCL, CX, CY, CZ

        R1=DSQRT((-XR(1) - U)**2 + XR(2)**2 + XR(3)**2)
        R2=DSQRT((-XR(1) - U + 1.0D0)**2 + XR(2)**2 + XR(3)**2)
        G=(1.0D0 - U)/R1 + U/R2
        EN =((XR(1) + XR(5))**2 + (XR(4) - XR(2))**2 + XR(6)**2)/2.0D0-G
        CANG = XR(1)**2 + XR(2)**2 + XR(1)*XR(5) - XR(2)*XR(4)

        CALL ROTtoFIX(XR, T, XF)
!
!	components of the angular momentum
!
        CX = (XF(2)*XF(6) - XF(3)*XF(5))
        CY = (XF(3)*XF(4) - XF(1)*XF(6))
  	    CZ = (XF(1)*XF(5) - XF(4)*XF(2))
        CTOT = DSQRT(CX**2 + CY**2 + CZ**2)
!
!			inclination
!
        AINCL = DACOS(CZ/CTOT)
        RETURN
      END

      SUBROUTINE JACOBI(XR, U, C1)

!	 This routine uses Broucke's system, like in
!	 Traveling Between the Lagrangian Points and the Moon
!
        IMPLICIT NONE
        REAL*8, DIMENSION(6) :: XR
        REAL*8 :: R1, R2, U, C, OME, C1
        R2 = DSQRT((-1.d0 + U + XR(1))**2 + XR(2)**2 + XR(3)**2)
        R1 = DSQRT((U + XR(1))**2 + XR(2)**2 + XR(3)**2)
        C = XR(5)**2 + XR(4)**2 + XR(6)**2
        OME = (XR(1)**2 + XR(2)**2)/2.d0 + (1.d0 - U)/R1 + U/R2
        C1 = C/2 - OME
        RETURN
      END
