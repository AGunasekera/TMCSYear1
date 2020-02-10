      FUNCTION DISP(XI, XJ)
      REAL*8 XI(3), XJ(3), DISP(3)
*
* Displacement between two vectors in 3D space
*
      INTEGER N

      DO 10 N = 1, 3
          DISP(I) = XJ(N) - XI(N)
10    CONTINUE

      RETURN
      END FUNCTION

      FUNCTION MAG(VEC)
*
* Magnitude of a given 3D vector
*
      REAL*8 VEC(3), MAG
      INTEGER N

      MAG = 0

      DO 20 N = 1, 3
          MAG = MAG + VEC(I) ** 2
20    CONTINUE

      RETURN
      END FUNCTION

      FUNCTION V_VERLET_X(X, V, F, DT)
      REAL :: X(3), V(3), F(3), DT
      REAL :: V_VERLET_X(3)
*
* New position from velocity Verlet algorithm
*
      INTEGER N

      DO 30 N = 1, 3
          V_VERLET_X(N) = X(N) + V(N) * DT + (F(N) / 2) * DT ** 2
30    CONTINUE

      RETURN
      END FUNCTION

      FUNCTION V_VERLET_V(V, F, NEWF, DT)
      REAL :: V(3), F(3), NEWF(3), DT
      REAL :: V_VERLET_V(3)
*
* New velocity from velocity Verlet algorithm
*
      INTEGER N

      DO 40 N = 1, 3
          V_VERLET_V(N) = V(N) + (F(N) + NEWF(N)) * DT / 2
40    CONTINUE

      RETURN
      END FUNCTION

      FUNCTION LJ2FORCE(R)
      REAL*8 R, LJ2FORCE
*
* Calculate 2-particle Lennard-Jones force divided by displacement vector
*
      REAL*8 RTOMIN6

      RTOMIN6 = R ** (-6)
      LJ2FORCE = 24 * (RTOMIN6 - 2 * RTOMIN6 ** 2) / (R ** 2)

      RETURN
      END FUNCTION

      FUNCTION LJFIJ(XI, XJ)
      REAL*8 XI(3), XJ(3), LJFIJ(3)
*
* Lennard-Jones force of particle j on particle i
*
      REAL*8 DISPIJ(3), DISTIJ, FSCAL
      INTEGER N

      DISPIJ = DISP(XI, XJ)
      DISTIJ = MAG(DISPIJ)
      FSCAL = LJ2FORCE(DISTIJ)

      DO 30 N = 1, 3
          LJFIJ(N) = DISPIJ(N) * FSCAL
30    CONTINUE

      RETURN
      END FUNCTION

      SUBROUTINE GETINITIALCONFIG(CONFIG)
      REAL*8 CONFIG(3,*)

      PROGRAM LJ13CLUSTER
*
*
*
      IMPLICIT NONE
      REAL :: X0, V0, DT
      INTEGER :: I, NSTEPS
      REAL :: X, V, F, NEWF
      REAL :: V_VERLET_X, V_VERLET_V, FORCE

      NSTEPS = 500
      X0 = 1
      V0 = 0
      DT = 0.5

      I = 0
      X = X0 
      V = V0
      F = FORCE(X)
      OPEN(1, FILE = "1Dharmonic.out", STATUS = "NEW")
      WRITE (1,*) "#T             X           V"
      WRITE (1,*) I * DT, "    ", X, "    ", V

      DO 10 I = 1, NSTEPS
      X = V_VERLET_X(X, V, F, DT)
      NEWF = FORCE(X)
      V = V_VERLET_V(V, F, NEWF, DT)
      F = NEWF
      WRITE (1,*) I * DT, "    ", X, "    ", V
10    CONTINUE
      CLOSE(1)

      END PROGRAM
