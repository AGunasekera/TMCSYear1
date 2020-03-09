      PROGRAM LJFLUID
*
*
*
      IMPLICIT NONE
      INTEGER :: NATOMS, I, NSTEPS, NFCCINLENSIM
      DOUBLE PRECISION :: XMATRIX(3, 120), VMATRIX(3, 120), FMATRIX(3
     &, 120), NEWFMATRIX(3, 120), DENSITY, LATTICEPARAM, BOXLEN, DT
      DOUBLE PRECISION :: POTENTIALENERGY, KINETICENERGY, 
     &TEMPERATURE, SITEMPERATURE

      NFCCINLENSIM = 3
      DENSITY = 1.1
      LATTICEPARAM = (4 / DENSITY) ** (1. / 3.)
      NATOMS = 4 * NFCCINLENSIM ** 3
      BOXLEN = NFCCINLENSIM * LATTICEPARAM
      NSTEPS = 1000
      DT = 0.01

      CALL INITIALISE(XMATRIX, NFCCINLENSIM, LATTICEPARAM)
      CALL TOTALPOTENTIALENERGY(XMATRIX, POTENTIALENERGY)

      DO 1 I = 1, NSTEPS
          CALL MC_TRIAL(XMATRIX, NEWXMATRIX, NATOMS)
          CALL TOTALPOTENTIALENERGY(NEWXMATRIX, NEWPOTENTIALENERGY)
          CALL METROPOLIS(XMATRIX, NEWXMATRIX, POTENTIALENERGY, 
     &NEWPOTENTIALENERGY, NATOMS)
          CALL UPDATE_AVERAGES(XMATRIX)
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
     &(XMATRIX, VMATRIX, FMATRIX, NATOMS, DT)
      DOUBLE PRECISION :: XMATRIX(3, *), VMATRIX(3, *), 
     &FMATRIX(3, *), DT
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
     &LATTICEPARAM)
      DOUBLE PRECISION :: XMATRIX(3, *), VMATRIX(3, *), LATTICEPARAM
      INTEGER :: NBOX
*
* Initialise fcc lattice of 4 * NBOX^3 atoms.
*
      INTEGER :: I, J, K, N

      N = 1

      DO 130 I = 1, NBOX
          DO 131 J = 1, NBOX
              DO 132 K = 1, NBOX
                  XMATRIX(1, N) = (I - 1) * LATTICEPARAM
                  XMATRIX(2, N) = (J - 1) * LATTICEPARAM
                  XMATRIX(3, N) = (K - 1) * LATTICEPARAM
                  N = N + 1
                  XMATRIX(1, N) = (I - 0.5) * LATTICEPARAM
                  XMATRIX(2, N) = (J - 0.5) * LATTICEPARAM
                  XMATRIX(3, N) = (K - 1) * LATTICEPARAM
                  N = N + 1
                  XMATRIX(1, N) = (I - 0.5) * LATTICEPARAM
                  XMATRIX(2, N) = (J - 1) * LATTICEPARAM
                  XMATRIX(3, N) = (K - 0.5) * LATTICEPARAM
                  N = N + 1
                  XMATRIX(1, N) = (I - 1) * LATTICEPARAM
                  XMATRIX(2, N) = (J - 0.5) * LATTICEPARAM
                  XMATRIX(3, N) = (K - 0.5) * LATTICEPARAM
                  N = N + 1
132           CONTINUE
131       CONTINUE
130   CONTINUE

      CALL SETMATRIXZERO(VMATRIX, 4 * NBOX ** 3)

      END SUBROUTINE

      SUBROUTINE WRITECONFIG(XMATRIX, NATOMS)
      DOUBLE PRECISION :: XMATRIX(3, *)
      INTEGER :: NATOMS
*
* Writes current configuration to .xyz file
*
      INTEGER :: I

      OPEN(1, ACCESS = "APPEND", FILE = "LJFluid.xyz")
      WRITE (1, *) NATOMS
      WRITE (1, *) ""
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

      SUBROUTINE MC_TRIAL(XMATRIX, NEWXMATRIX, NATOMS, VAR)
      DOUBLE PRECISION :: XMATRIX(3, *), NEWXMATRIX(3, *), VAR
      INTEGER :: NATOMS
*
* Pick a random particle and give it a M-B velocity according to TEMP
*
      INTEGER :: I, J, N

      CALL RANDOMPARTICLE(NATOMS, N)

      DO 180 I = 1, NATOMS
          DO 181 J = 1, 3
              NEWXMATRIX(J, I) = XMATRIX(J, I)
181       CONTINUE
180   CONTINUE

      D0 190 J = 1, 3
          NEWXMATRIX(J, N) = NEWXMATRIX(J, N) + NORMALDIST(0, VAR)
190   CONTINUE

      END SUBROUTINE

      SUBROUTINE RANDOMPARTICLE(NATOMS, I)
      INTEGER :: NATOMS, I
*
* Pick a random particle I from NATOMS
*
      I = INT(NATOMS * RAND(0)) + 1

      END SUBROUTINE
