      FUNCTION NORMALDIST(MEAN, VAR, PI)
      DOUBLE PRECISION :: MEAN, VAR, PI, NORMALDIST
*
* Generates pseudorandom normally distributed value
* using Box-Muller transform from uniform distribution
*
      DOUBLE PRECISION :: STDEV, U1, U2

      STDEV = SQRT(VAR)
      U1 = RAND(0)
      U2 = RAND(0)

      NORMALDIST = SQRT(-2 * LOG(U1)) * COS(2 * PI * U2)
      NORMALDIST = NORMALDIST * STDEV + MEAN

      RETURN
      END FUNCTION


      PROGRAM LJFLUID
*
*
*
      IMPLICIT NONE
      INTEGER :: NATOMS, I, NSTEPS, NFCCINLENSIM
      DOUBLE PRECISION :: XMATRIX(3, 500), VMATRIX(3, 500), FMATRIX(3
     &, 500), NEWFMATRIX(3, 500)
      DOUBLE PRECISION :: TEMP, NU, DENSITY, LATTICEPARAM, BOXLEN, DT
      DOUBLE PRECISION :: POTENTIALENERGY, KINETICENERGY, 
     &TEMPERATURE, SITEMPERATURE
      DOUBLE PRECISION, PARAMETER :: PI = 4. * ATAN(1.)

      TEMP = 2
      NU = 1
      NFCCINLENSIM = 5
      DENSITY = 1
      LATTICEPARAM = (4. / DENSITY) ** (1. / 3.)
      NATOMS = 4 * NFCCINLENSIM ** 3
      BOXLEN = NFCCINLENSIM * LATTICEPARAM
      NSTEPS = 20000
      DT = 0.001

      CALL INITIALISE(XMATRIX, VMATRIX, NFCCINLENSIM, LATTICEPARAM, 
     &TEMP, PI)
              CALL WRITECONFIG(XMATRIX, NATOMS, 0)
              CALL TOTALPOTENTIALENERGY(XMATRIX, NATOMS, BOXLEN, 
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

      CALL GETFORCES(FMATRIX, XMATRIX, BOXLEN, NATOMS)

      DO 1 I = 1, NSTEPS
          CALL V_VERLET_UPDATE_XMATRIX(XMATRIX, VMATRIX, 
     &FMATRIX, NATOMS, BOXLEN, DT)
          CALL GETFORCES(NEWFMATRIX, XMATRIX, BOXLEN, NATOMS)
          CALL V_VERLET_UPDATE_VMATRIX(VMATRIX, FMATRIX, 
     &NEWFMATRIX, NATOMS, DT)
          CALL ANDERSEN(VMATRIX, NATOMS, TEMP, NU, DT, PI)
          FMATRIX = NEWFMATRIX
          IF (MOD(I, 50) .EQ. 0) THEN
              CALL WRITECONFIG(XMATRIX, NATOMS, I)
              CALL TOTALPOTENTIALENERGY(XMATRIX, NATOMS, BOXLEN, 
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
      DOUBLE PRECISION :: MATRIX(3, *)
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

      SUBROUTINE GETMINDISP(XI, XJ, BOXLEN, DISP)
      DOUBLE PRECISION :: XI(3), XJ(3), BOXLEN, DISP(3)
*
* Displacement between two vectors in 3D space,
* applying the minimum image convention
* for periodic boundary conditions.
*
      INTEGER :: N

      DO 10 N = 1, 3
          DISP(N) = XJ(N) - XI(N)
          IF (DISP(N) .GT. 0.5 * BOXLEN) THEN
              DISP(N) = DISP(N) - BOXLEN
          ELSE IF (DISP(N) .LE. -0.5 * BOXLEN) THEN
              DISP(N) = DISP(N) + BOXLEN
          END IF
10    CONTINUE

      END SUBROUTINE

      SUBROUTINE GETSQUAREMAGNITUDE(VEC, MAG2)
*
* Square magnitude of a given 3D vector
*
      DOUBLE PRECISION :: VEC(3), MAG2
      INTEGER :: N

      MAG2 = 0

      DO 20 N = 1, 3
          MAG2 = MAG2 + VEC(N) ** 2
20    CONTINUE

      END SUBROUTINE

      SUBROUTINE GETFORCESCALAR(R2, LJ2FSCAL)
      DOUBLE PRECISION :: R2, LJ2FSCAL
*
* Calculate 2-particle Lennard-Jones force divided by displacement vector
*
      DOUBLE PRECISION :: RTOMIN6

      RTOMIN6 = R2 ** (-3)
      LJ2FSCAL = 24 * ((2 * (RTOMIN6 ** 2)) - RTOMIN6) / R2

      END SUBROUTINE

      SUBROUTINE GETLJFORCEIJ(XI, XJ, BOXLEN, LJFIJ)
      DOUBLE PRECISION :: XI(3), XJ(3), BOXLEN, LJFIJ(3)
*
* Lennard-Jones force of particle j on particle i
*
      DOUBLE PRECISION :: DISPIJ(3), R2, FSCAL
      INTEGER :: N

      CALL GETMINDISP(XI, XJ, BOXLEN, DISPIJ)
      CALL GETSQUAREMAGNITUDE(DISPIJ, R2)
      CALL GETFORCESCALAR(R2, FSCAL)

      DO 90 N = 1, 3
          LJFIJ(N) = - DISPIJ(N) * FSCAL
90    CONTINUE

      END SUBROUTINE

      SUBROUTINE GET2ATOMPOTENTIAL(XI, XJ, BOXLEN, LJPOTIJ)
      DOUBLE PRECISION :: XI(3), XJ(3), BOXLEN, LJPOTIJ
*
* Lennard-Jones pair potential for particles i and j
*
      DOUBLE PRECISION :: DISPIJ(3), R2, RTOMIN6

      CALL GETMINDISP(XI, XJ, BOXLEN, DISPIJ)
      CALL GETSQUAREMAGNITUDE(DISPIJ, R2)

      RTOMIN6 = R2 ** (-3)
      LJPOTIJ = 4 * (RTOMIN6 ** 2 - RTOMIN6)

      END SUBROUTINE

      SUBROUTINE V_VERLET_UPDATE_XMATRIX
     &(XMATRIX, VMATRIX, FMATRIX, NATOMS, BOXLEN, DT)
      DOUBLE PRECISION :: XMATRIX(3, *), VMATRIX(3, *), 
     &FMATRIX(3, *), BOXLEN, DT
      INTEGER :: NATOMS
*
* Update positions by velocity Verlet algorithm
* applying the minimum image convention
*
      INTEGER :: I, N

      DO 50 I = 1, NATOMS
          DO 60 N = 1, 3
              XMATRIX(N, I) = XMATRIX(N, I) + VMATRIX(N, I) 
     &* DT + (FMATRIX(N, I) / 2) * DT ** 2
              DO WHILE (XMATRIX(N, I) .GE. BOXLEN)
                  XMATRIX(N, I) = XMATRIX(N, I) - BOXLEN
              END DO
              DO WHILE (XMATRIX(N, I) .LT. 0)
                  XMATRIX(N, I) = XMATRIX(N, I) + BOXLEN
              END DO
60        CONTINUE
50    CONTINUE

      END SUBROUTINE

      SUBROUTINE V_VERLET_UPDATE_VMATRIX(VMATRIX, FMATRIX, 
     &NEWFMATRIX, NATOMS, DT)
      DOUBLE PRECISION :: VMATRIX(3, *), FMATRIX(3, *), 
     &NEWFMATRIX(3, *), DT
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

      SUBROUTINE GETFORCES(FMATRIX, XMATRIX, BOXLEN, NATOMS)
      DOUBLE PRECISION :: XMATRIX(3, *), FMATRIX(3, *), BOXLEN
      INTEGER :: NATOMS
*
* Get matrix of forces for NATOMS atoms
*
      DOUBLE PRECISION :: XI(3), XJ(3), FIJ(3)
      INTEGER :: I, J, N

      CALL SETMATRIXZERO(FMATRIX, NATOMS)

      DO 100 I = 1, NATOMS - 1
          XI = XMATRIX(:, I)
          DO 110 J = I + 1, NATOMS
              XJ = XMATRIX(:, J)
              CALL GETLJFORCEIJ(XI, XJ, BOXLEN, FIJ)
              DO 120 N = 1, 3
                  FMATRIX(N, I) = FMATRIX(N, I) + FIJ(N)
                  FMATRIX(N, J) = FMATRIX(N, J) - FIJ(N)
120           CONTINUE
110       CONTINUE
100   CONTINUE

      END SUBROUTINE

      SUBROUTINE INITIALISE(XMATRIX, VMATRIX, NBOX, 
     &LATTICEPARAM, TEMP, PI)
      DOUBLE PRECISION :: XMATRIX(3, *), VMATRIX(3, *), 
     &LATTICEPARAM, TEMP, PI
      INTEGER :: NBOX
*
* Initialise fcc lattice of 4 * NBOX^3 atoms.
*
      INTEGER :: I, J, K, N
      DOUBLE PRECISION :: NORMALDIST

      N = 1

      DO 130 I = 1, NBOX
          DO 131 J = 1, NBOX
              DO 132 K = 1, NBOX
                  XMATRIX(1, N) = (I - 1) * LATTICEPARAM
                  XMATRIX(2, N) = (J - 1) * LATTICEPARAM
                  XMATRIX(3, N) = (K - 1) * LATTICEPARAM
                  VMATRIX(1, N) = NORMALDIST(0.D0, TEMP, PI)
                  VMATRIX(2, N) = NORMALDIST(0.D0, TEMP, PI)
                  VMATRIX(3, N) = NORMALDIST(0.D0, TEMP, PI)
                  N = N + 1
                  XMATRIX(1, N) = (I - 0.5) * LATTICEPARAM
                  XMATRIX(2, N) = (J - 0.5) * LATTICEPARAM
                  XMATRIX(3, N) = (K - 1) * LATTICEPARAM
                  VMATRIX(1, N) = NORMALDIST(0.D0, TEMP, PI)
                  VMATRIX(2, N) = NORMALDIST(0.D0, TEMP, PI)
                  VMATRIX(3, N) = NORMALDIST(0.D0, TEMP, PI)
                  N = N + 1
                  XMATRIX(1, N) = (I - 0.5) * LATTICEPARAM
                  XMATRIX(2, N) = (J - 1) * LATTICEPARAM
                  XMATRIX(3, N) = (K - 0.5) * LATTICEPARAM
                  VMATRIX(1, N) = NORMALDIST(0.D0, TEMP, PI)
                  VMATRIX(2, N) = NORMALDIST(0.D0, TEMP, PI)
                  VMATRIX(3, N) = NORMALDIST(0.D0, TEMP, PI)
                  N = N + 1
                  XMATRIX(1, N) = (I - 1) * LATTICEPARAM
                  XMATRIX(2, N) = (J - 0.5) * LATTICEPARAM
                  XMATRIX(3, N) = (K - 0.5) * LATTICEPARAM
                  VMATRIX(1, N) = NORMALDIST(0.D0, TEMP, PI)
                  VMATRIX(2, N) = NORMALDIST(0.D0, TEMP, PI)
                  VMATRIX(3, N) = NORMALDIST(0.D0, TEMP, PI)
                  N = N + 1
132           CONTINUE
131       CONTINUE
130   CONTINUE

      WRITE(*,*) TEMP

      END SUBROUTINE

      SUBROUTINE WRITECONFIG(XMATRIX, NATOMS, STEP)
      DOUBLE PRECISION :: XMATRIX(3, *)
      INTEGER :: NATOMS, STEP
*
* Writes current configuration to .xyz file
*
      INTEGER :: I

      OPEN(1, ACCESS = "APPEND", FILE = "LJFluid.xyz")
      WRITE (1, *) NATOMS
      WRITE (1, *) "Timestep:", STEP
      DO 140 I = 1, NATOMS
          WRITE (1, *) "Ar", XMATRIX(:,I)
140   CONTINUE
      CLOSE(1)

      END SUBROUTINE

      SUBROUTINE TOTALKINETICENERGY(VMATRIX, NATOMS, TOTKE)
      DOUBLE PRECISION :: VMATRIX(3, *), TOTKE
      INTEGER :: NATOMS
*
* Get the total kinetic energy (in reduced units)
*
      DOUBLE PRECISION :: VI(3), V2I
      INTEGER :: I

      TOTKE = 0

      DO 150 I = 1, NATOMS
          VI = VMATRIX(:, I)
          CALL GETSQUAREMAGNITUDE(VI, V2I)
          TOTKE = TOTKE + V2I / 2
150   CONTINUE

      END SUBROUTINE

      SUBROUTINE GETTEMP(TOTALKE, NATOMS, TEMP)
      DOUBLE PRECISION :: TOTALKE, TEMP
      INTEGER :: NATOMS
*
* Get temperature (in reduced units)
* from average kinetic energy of particles
*
      TEMP = TOTALKE / DBLE(NATOMS) * (2. / 3.)

      END SUBROUTINE

      SUBROUTINE CONVERTTEMP(TEMP, SITEMP)
      DOUBLE PRECISION :: TEMP, SITEMP
*
* Convert temperature from reduced units to SI, for Argon
*
      DOUBLE PRECISION :: SIGMA, EPSILONBYKBK

      EPSILONBYKBK = 120
*     SIGMA = 3.4 * 10 ** (-10)

      SITEMP = TEMP * EPSILONBYKBK

      END SUBROUTINE

      SUBROUTINE TOTALPOTENTIALENERGY(XMATRIX, NATOMS, BOXLEN, TOTPE)
      DOUBLE PRECISION :: XMATRIX(3, *), BOXLEN, TOTPE
      INTEGER :: NATOMS
*
* Get total potential energy using Lennard-Jones pair potential
*
      DOUBLE PRECISION :: XI(3), XJ(3), POTIJ
      INTEGER :: I, J

      TOTPE = 0

      DO 160 I = 1, NATOMS - 1
          XI = XMATRIX(:, I)
          DO 170 J = I + 1, NATOMS
              XJ = XMATRIX(:, J)
              CALL GET2ATOMPOTENTIAL(XI, XJ, BOXLEN, POTIJ)
              TOTPE = TOTPE + POTIJ
170       CONTINUE
160   CONTINUE

      END SUBROUTINE

      SUBROUTINE ANDERSEN(VMATRIX, NATOMS, TEMP, NU, DT, PI)
      DOUBLE PRECISION :: VMATRIX(3, *), TEMP, NU, PI, DT
      INTEGER :: NATOMS
*
* Pick random particles and give them M-B velocities according to TEMP
*
      DOUBLE PRECISION :: NORMALDIST
      INTEGER :: I

      DO 180 I = 1, NATOMS
          IF (RAND(0) .LT. NU * DT) THEN
              VMATRIX(1, I) = NORMALDIST(0.D0, TEMP, PI)
              VMATRIX(2, I) = NORMALDIST(0.D0, TEMP, PI)
              VMATRIX(3, I) = NORMALDIST(0.D0, TEMP, PI)
          END IF
180   CONTINUE

      END SUBROUTINE

      SUBROUTINE RANDOMPARTICLE(NATOMS, I)
      INTEGER :: NATOMS, I
*
* Pick a random particle I from NATOMS
*

      I = INT(NATOMS * RAND(0)) + 1

      END SUBROUTINE
