      SUBROUTINE DLACN1( N, T, V, X, LDX, XOLD, LDXOLD, WRK,
     $                   H, IND, INDH, EST, KASE, ISEED, INFO )
*
*     .. Scalar Arguments ..
      INTEGER            INFO, KASE, LDXOLD, LDX, N, T
      DOUBLE PRECISION   EST
*     ..
*     .. Array Arguments ..
      INTEGER            IND( * ), INDH( * ), ISEED( 4 )
      DOUBLE PRECISION   H( * ),  V( * ), X( LDX, * ), WRK( * ),
     $                   XOLD( LDXOLD, * )
*     ..
*
*  Purpose
*  =======
*
*  DLACN1 estimates the 1-norm of a square, real matrix A.
*  Reverse communication is used for evaluating matrix-matrix products.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The order of the matrix.  N >= 1.
*
*  T      (input) INTEGER
*         The number of columns used at each step.
*
*  V      (output) DOUBLE PRECISION array, dimension (N).
*         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
*         (W is not returned).
*
*  X      (input/output) DOUBLE PRECISION array, dimension (N,T)
*         On an intermediate return, X should be overwritten by
*               A * X,   if KASE=1,
*               A' * X,  if KASE=2,
*         and DLACN1 must be re-called with all the other parameters
*         unchanged.
*
*  LDX    (input) INTEGER
*         The leading dimension of X.  LDX >= max(1,N).
*
*  XOLD   (workspace) DOUBLE PRECISION array, dimension (N,T)
*
*  LDXOLD (input) INTEGER
*         The leading dimension of XOLD.  LDXOLD >= max(1,N).
*
*  WRK    (workspace) DOUBLE PRECISION array, dimension (T)
*
*  H      (workspace) DOUBLE PRECISION array, dimension (N)
*
*  IND    (workspace) INTEGER array, dimension (N)
*
*  INDH   (workspace) INTEGER array, dimension (N)
*
*  EST    (output) DOUBLE PRECISION
*         An estimate (a lower bound) for norm(A).
*
*  KASE   (input/output) INTEGER
*         On the initial call to DLACN1, KASE should be 0.
*         On an intermediate return, KASE will be 1 or 2, indicating
*         whether X should be overwritten by A * X  or A' * X.
*         On the final return from DLACN1, KASE will again be 0.
*
*  ISEED  (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  INFO   (output) INTEGER
*         INFO describes how the iteration terminated:
*            INFO = 1: iteration limit reached.
*            INFO = 2: estimate not increased.
*            INFO = 3: repeated sign matrix.
*            INFO = 4: power method convergence test.
*            INFO = 5: repeated unit vectors.
*
*  ====================================================================
*
*     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IBEST, ITEMP, ITER, J, JUMP
      DOUBLE PRECISION   ESTOLD, TEMP
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DASUM, DDOT
      EXTERNAL           DASUM, DDOT, IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLACPY, DLAPST, DLARNV, DLARPC, DLASCL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, NINT, SIGN
*     ..
*     .. Save statement ..
      SAVE
*     ..
*     .. Executable Statements ..
*
*     
      IF( KASE .EQ. 0 ) THEN
*
         ESTOLD = ZERO
         ITER = 1
         ITEMP = 1
         INFO = 0
*
         DO 10 I = 1, N
            X( I, 1 ) = ONE
            IND( I ) = I
            INDH( I ) = 0
   10    CONTINUE    
*
         DO 30 J = 2, T
            CALL DLARNV( 2, ISEED, N, X( 1, J ) )
            DO 20 I = 1, N
               X( I, J ) = SIGN( ONE, X( I, J ) )
   20       CONTINUE        
   30    CONTINUE
*
         IF ( T .GT. 1 )
     $      CALL DLARPC( N, T, X, LDX, XOLD, LDXOLD, WRK, KASE, ISEED)
*        
         CALL DLASCL( 'G', 0, 0, DBLE(N), ONE, N, T, X, LDX, INFO )
*
         KASE = 1
         JUMP = 1
         RETURN
      END IF     
*
      GO TO ( 40, 100 ) JUMP
*
*     ................ ENTRY   (JUMP = 1)
*     FIRST HALF OF THE ITERATION: X HAS BEEN OVERWRITTEN BY A*X.
*
   40 CONTINUE
*
      IF ( ITER .EQ. 1  .AND.  N .EQ. 1 ) THEN
         V( 1 ) = X( 1, 1 )
         EST = ABS( V( 1 ) )
*        ... QUIT
         GO TO 210
      END IF

      EST = ZERO
      DO 50 J = 1, T
         TEMP = DASUM( N, X( 1, J ), 1 )
         IF ( TEMP .GT. EST ) THEN
             EST = TEMP
             ITEMP = J 
         END IF
   50 CONTINUE
*
      IF ( EST .GT. ESTOLD  .OR.  ITER .EQ. 2 ) THEN
         IBEST = IND( ITEMP )
      END IF
*
      IF ( EST .LE. ESTOLD  .AND.  ITER .GE. 2 ) THEN
         EST = ESTOLD
         INFO = 2
         GO TO 210
      END IF
*
      ESTOLD = EST
      CALL DCOPY( N, X( 1, ITEMP ), 1, V, 1 )
*
      IF ( ITER .GT. ITMAX ) THEN
         INFO = 1
         GO TO 210
      END IF
*
      DO 70 J = 1, T
         DO 60 I = 1, N
            X( I, J ) = SIGN( ONE, X( I, J ) )
   60    CONTINUE
   70 CONTINUE
*
      IF ( ITER .GT. 1 ) THEN
*
*        IF ALL COLUMNS of X PARALLEL TO XOLD, EXIT. 
*
         DO 80 J = 1, T
            CALL DGEMV( 'Transpose', N, T, ONE, XOLD, LDXOLD,
     $           X( 1, J ), 1, ZERO, WRK, 1)
            IF ( NINT(ABS(WRK(IDAMAX(T, WRK, 1)))) .LT. N ) GO TO 90
   80    CONTINUE
         INFO = 3
         GO TO 210
*     
   90    CONTINUE
*
         IF ( T. GT. 1 )
     $      CALL DLARPC( N, T, X, LDX, XOLD, LDXOLD, WRK, KASE, ISEED)
*
      ELSE 
*         
         IF ( T. GT. 1 )
     $      CALL DLARPC( N, T, X, LDX, XOLD, LDXOLD, WRK, 0, ISEED)
*
      END IF
*
      CALL DLACPY( 'Whole', N, T, X, LDX, XOLD, LDXOLD )
*
      KASE = 2
      JUMP = 2
      RETURN
*
*     ................ ENTRY   (JUMP = 2)
*     SECOND HALF OF THE ITERATION: 
*                      X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
*
  100 CONTINUE
*   
      DO 110 I = 1, N
         H( I ) = ABS( X ( I, IDAMAX( T, X( I, 1 ), N ) ) )
         IND( I ) = I
  110 CONTINUE
*
      IF ( ITER .GE. 2 .AND. H( IDAMAX(N, H, 1) ) .EQ. H(IBEST) ) THEN
         INFO = 4
         GO TO 210
      END IF
*
*     Sort so that h(i) >= h(j) for i < j
*
      CALL DLAPST( 'D', N, H, IND, ITEMP )
* 
      IF ( ITER .EQ. 1 ) THEN
         ITEMP = T
         GO TO 170
      END IF
*
*     IF IND(1:T) IS CONTAINED IN INDH, TERMINATE.
*
      IF ( T .GT. 1 ) THEN
         DO 130 J = 1, T
            DO 120 I = 1, (ITER-1)*T
               IF (I .GT. N .OR. IND( J ) .EQ. INDH( I )) GO TO 130
  120       CONTINUE
            GO TO 140
  130    CONTINUE
         INFO = 5
         GO TO 210
  140    CONTINUE   
*
*        REPLACE IND(1:T) BY THE FIRST T INDICES IN IND THAT
*        ARE NOT IN INDH. 
*
         ITEMP = 1
         DO 160 J = 1, N
            DO 150 I = 1, (ITER-1)*T
               IF ( I .GT. N .OR. IND( J ) .EQ. INDH( I ) ) GO TO 160
  150       CONTINUE
            IND( ITEMP ) = IND( J )
            IF ( ITEMP .EQ. T ) GO TO 170
            ITEMP = ITEMP + 1
  160    CONTINUE   
      END IF
*
      ITEMP = ITEMP - 1
*
  170 CONTINUE
*   
      IF ( (ITER-1)*T .GE. N ) THEN
         DO 180 J = 1, ITEMP
            INDH( (ITER-1)*T+J ) = IND( J )
  180    CONTINUE
      END IF
*
      DO 200 J = 1, T
         DO 190 I = 1, N
            X( I, J ) = ZERO
  190    CONTINUE
         X( IND( J ), J ) = ONE
  200 CONTINUE
*
      ITER = ITER + 1
*
      KASE = 1
      JUMP = 1
      RETURN
*
  210 CONTINUE
      KASE = 0
      RETURN
*
*     End of DLACN1
*
      END
