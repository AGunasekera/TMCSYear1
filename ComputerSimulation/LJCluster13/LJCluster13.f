      PROGRAM LJ13CLUSTER
*
*
*
      IMPLICIT NONE
      INTEGER :: NATOMS, I, NSTEPS
      REAL*8 :: XMATRIX(3, 13), VMATRIX(3, 13), FMATRIX(3, 13), 
     &NEWFMATRIX(3, 13), DT
      REAL*8 :: POTENTIALENERGY, KINETICENERGY, 
     &TEMPERATURE, SITEMPERATURE

      NATOMS = 13
      NSTEPS = 1000
      DT = 0.01

      CALL INITIALISE(XMATRIX, VMATRIX, NATOMS)

      CALL GETFORCES(FMATRIX, XMATRIX, NATOMS)

      DO 1 I = 1, NSTEPS
          CALL V_VERLET_UPDATE_XMATRIX(XMATRIX, VMATRIX, 
     &FMATRIX, NATOMS, DT)
          CALL GETFORCES(NEWFMATRIX, XMATRIX, NATOMS)
          CALL V_VERLET_UPDATE_VMATRIX(VMATRIX, FMATRIX, 
     &NEWFMATRIX, NATOMS, DT)
          FMATRIX = NEWFMATRIX
          IF (MOD(I, 10) .EQ. 0) THEN
              CALL WRITECONFIG(XMATRIX, NATOMS)
              CALL TOTALPOTENTIALENERGY(XMATRIX, NATOMS, 
     &POTENTIALENERGY)
              CALL TOTALKINETICENERGY(VMATRIX, NATOMS, KINETICENERGY)
              CALL GETTEMP(KINETICENERGY, NATOMS, TEMPERATURE)
              CALL CONVERTTEMP(TEMPERATURE, SITEMPERATURE)
              WRITE (*, *) "Time step: ", I
              WRITE (*, *) "Potential energy: ", POTENTIALENERGY
              WRITE (*, *) "Kinetic energy: ", KINETICENERGY
              WRITE (*, *) "Total energy: ", POTENTIALENERGY 
     &+ KINETICENERGY
              WRITE (*, *) "Temperature (reduced units): ", TEMPERATURE
              WRITE (*, *) "Approximate eqivalent for Argon / K: ", 
     &SITEMPERATURE
              WRITE (*, *) ""
              OPEN (3, ACCESS = "APPEND", FILE = "thermo.dat")
              WRITE (3, *) I * DT, POTENTIALENERGY, KINETICENERGY, 
     &POTENTIALENERGY + KINETICENERGY, TEMPERATURE
              CLOSE (3)
              
          END IF
1     CONTINUE

      END PROGRAM

      SUBROUTINE SETMATRIXZERO(MATRIX, N)
      REAL*8 :: MATRIX(3, *)
      INTEGER :: N
*
* Sets the elements of a 3 x N matrix to zero
*
      INTEGER :: I, J

      DO 30 I = 1, 3
          DO 40 J = 1, N
              MATRIX(I, J) = 0.
40        CONTINUE
30    CONTINUE

      END SUBROUTINE

      SUBROUTINE GETDISP(XI, XJ, DISP)
      REAL*8 :: XI(3), XJ(3), DISP(3)
*
* Displacement between two vectors in 3D space
*
      INTEGER :: N

      DO 10 N = 1, 3
          DISP(N) = XJ(N) - XI(N)
10    CONTINUE

      END SUBROUTINE

      SUBROUTINE GETSQUAREMAGNITUDE(VEC, MAG2)
*
* Square magnitude of a given 3D vector
*
      REAL*8 :: VEC(3), MAG2
      INTEGER :: N

      MAG2 = 0

      DO 20 N = 1, 3
          MAG2 = MAG2 + VEC(N) ** 2
20    CONTINUE

      END SUBROUTINE

      SUBROUTINE GETFORCESCALAR(R2, LJ2FSCAL)
      REAL*8 :: R2, LJ2FSCAL
*
* Calculate 2-particle Lennard-Jones force divided by displacement vector
*
      REAL*8 :: RTOMIN6

      RTOMIN6 = R2 ** (-3)
      LJ2FSCAL = 24 * ((2 * (RTOMIN6 ** 2)) - RTOMIN6) / R2

      END SUBROUTINE

      SUBROUTINE GETLJFORCEIJ(XI, XJ, LJFIJ)
      REAL*8 :: XI(3), XJ(3), LJFIJ(3)
*
* Lennard-Jones force of particle j on particle i
*
      REAL*8 :: DISPIJ(3), R2, FSCAL
      INTEGER :: N

      CALL GETDISP(XI, XJ, DISPIJ)
      CALL GETSQUAREMAGNITUDE(DISPIJ, R2)
      CALL GETFORCESCALAR(R2, FSCAL)

      DO 90 N = 1, 3
          LJFIJ(N) = - DISPIJ(N) * FSCAL
90    CONTINUE

      END SUBROUTINE

      SUBROUTINE GET2ATOMPOTENTIAL(XI, XJ, LJPOTIJ)
      REAL*8 :: XI(3), XJ(3), LJPOTIJ
*
* Lennard-Jones pair potential for particles i and j
*
      REAL*8 :: DISPIJ(3), R2, RTOMIN6

      CALL GETDISP(XI, XJ, DISPIJ)
      CALL GETSQUAREMAGNITUDE(DISPIJ, R2)

      RTOMIN6 = R2 ** (-3)
      LJPOTIJ = 4 * (RTOMIN6 ** 2 - RTOMIN6)

      END SUBROUTINE

      SUBROUTINE V_VERLET_UPDATE_XMATRIX
     &(XMATRIX, VMATRIX, FMATRIX, NATOMS, DT)
      REAL*8 :: XMATRIX(3, *), VMATRIX(3, *), FMATRIX(3, *), DT
      INTEGER :: NATOMS
*
* Update positions by velocity Verlet algorithm
*
      INTEGER :: I, N

      DO 50 I = 1, NATOMS
          DO 60 N = 1, 3
              XMATRIX(N, I) = XMATRIX(N, I) + VMATRIX(N, I) 
     &* DT + (FMATRIX(N, I) / 2) * DT ** 2
60        CONTINUE
50    CONTINUE

      END SUBROUTINE

      SUBROUTINE V_VERLET_UPDATE_VMATRIX(VMATRIX, FMATRIX, 
     &NEWFMATRIX, NATOMS, DT)
      REAL*8 :: VMATRIX(3, *), FMATRIX(3, *), NEWFMATRIX(3, *), DT
      INTEGER :: NATOMS
*
* Update velocities by velocity Verlet algorithm
*
      INTEGER :: I, N

      DO 70 I = 1, NATOMS
          DO 80 N = 1, 3
              VMATRIX(N, I) = VMATRIX(N, I) + (FMATRIX(N, I) 
     &+ NEWFMATRIX(N, I)) * DT / 2
80        CONTINUE
70    CONTINUE

      END SUBROUTINE

      SUBROUTINE GETFORCES(FMATRIX, XMATRIX, NATOMS)
      REAL*8 :: XMATRIX(3, *), FMATRIX(3, *)
      INTEGER :: NATOMS
*
* Get matrix of forces for NATOMS atoms
*
      REAL*8 :: XI(3), XJ(3), FIJ(3)
      INTEGER :: I, J, N

      CALL SETMATRIXZERO(FMATRIX, NATOMS)

      DO 100 I = 1, NATOMS - 1
          XI = XMATRIX(:, I)
          DO 110 J = I + 1, NATOMS
              XJ = XMATRIX(:, J)
              CALL GETLJFORCEIJ(XI, XJ, FIJ)
              DO 120 N = 1, 3
                  FMATRIX(N, I) = FMATRIX(N, I) + FIJ(N)
                  FMATRIX(N, J) = FMATRIX(N, J) - FIJ(N)
120           CONTINUE
110       CONTINUE
100   CONTINUE

      END SUBROUTINE

      SUBROUTINE INITIALISE(XMATRIX, VMATRIX, NATOMS)
      REAL*8 :: XMATRIX(3, *), VMATRIX(3, *)
      INTEGER :: NATOMS
*
* Initialise positions and velocities of NATOMS atoms.
*
      INTEGER :: I

      OPEN(1, FILE = "initialconfig.txt")

      DO 130 I = 1, NATOMS
          READ (1, *) XMATRIX(1, I), XMATRIX(2, I), XMATRIX(3, I)
130   CONTINUE

      CLOSE(1)

      CALL SETMATRIXZERO(VMATRIX, NATOMS)

      END SUBROUTINE

      SUBROUTINE WRITECONFIG(XMATRIX, NATOMS)
      REAL*8 :: XMATRIX(3, *)
      INTEGER :: NATOMS
*
* Writes current configuration to .xyz file
*
      INTEGER :: I

      OPEN(1, ACCESS = "APPEND", FILE = "LJ13.xyz")
      WRITE (1, *) NATOMS
      WRITE (1, *) ""
      DO 140 I = 1, NATOMS
          WRITE (1, *) "Ar", XMATRIX(:,I)
140   CONTINUE
      CLOSE(1)

      END SUBROUTINE

      SUBROUTINE TOTALKINETICENERGY(VMATRIX, NATOMS, TOTKE)
      REAL*8 :: VMATRIX(3, *), TOTKE
      INTEGER :: NATOMS
*
* Get the total kinetic energy (in reduced units)
*
      REAL*8 :: VI(3), V2I
      INTEGER :: I

      TOTKE = 0

      DO 150 I = 1, NATOMS
          VI = VMATRIX(:, I)
          CALL GETSQUAREMAGNITUDE(VI, V2I)
          TOTKE = TOTKE + V2I / 2
150   CONTINUE

      END SUBROUTINE

      SUBROUTINE GETTEMP(TOTALKE, NATOMS, TEMP)
      REAL*8 :: TOTALKE, TEMP
      INTEGER :: NATOMS
*
* Get temperature (in reduced units)
* from average kinetic energy of particles
*
      TEMP = TOTALKE / DBLE(NATOMS) * (2. / 3.)

      END SUBROUTINE

      SUBROUTINE CONVERTTEMP(TEMP, SITEMP)
      REAL*8 :: TEMP, SITEMP
*
* Convert temperature from reduced units to SI, for Argon
*
      REAL*8 :: SIGMA, EPSILONBYKBK

      EPSILONBYKBK = 120
      SIGMA = 3.4 * 10 ** (-10)

      SITEMP = TEMP * EPSILONBYKBK

      END SUBROUTINE

      SUBROUTINE TOTALPOTENTIALENERGY(XMATRIX, NATOMS, TOTPE)
      REAL*8 :: XMATRIX(3, *), TOTPE
      INTEGER :: NATOMS
*
* Get total potential energy using Lennard-Jones pair potential
*
      REAL*8 :: XI(3), XJ(3), POTIJ
      INTEGER :: I, J

      TOTPE = 0

      DO 160 I = 1, NATOMS - 1
          XI = XMATRIX(:, I)
          DO 170 J = I + 1, NATOMS
              XJ = XMATRIX(:, J)
              CALL GET2ATOMPOTENTIAL(XI, XJ, POTIJ)
              TOTPE = TOTPE + POTIJ
170       CONTINUE
160   CONTINUE

      END SUBROUTINE
