      FUNCTION V_VERLET_X(X, V, F, DT)
      REAL :: X, V, F, DT
      REAL :: V_VERLET_X

      V_VERLET_X = X + V * DT + (F / 2) * DT ** 2
      RETURN
      END FUNCTION

      FUNCTION V_VERLET_V(V, F, NEWF, DT)
      REAL :: V, F, NEWF, DT
      REAL :: V_VERLET_V

      V_VERLET_V = V + (F + NEWF) * DT / 2
      RETURN
      END FUNCTION

      FUNCTION FORCE(X)
      REAL :: X
      REAL :: FORCE

      FORCE = -X
      END FUNCTION

      PROGRAM HARM1D
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
