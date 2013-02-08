      SUBROUTINE ZLACN1( N, T, V, X, LDX, XOLD, LDXOLD,
     $                   H, IND, INDH, EST, KASE, ISEED, INFO )
*
*     .. Scalar Arguments ..
      INTEGER            INFO, KASE, LDXOLD, LDX, N, T
      DOUBLE PRECISION   EST
*     ..
*     .. Array Arguments ..

      INTEGER            IND( * ), INDH( * ), ISEED( 4 )
      DOUBLE PRECISION   H( * )
      COMPLEX*16         V( * ), X( LDX, * ), XOLD( LDXOLD, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLACN1 estimates the 1-norm of a square, real matrix A.
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
*  V      (output) COMPLEX*16 array, dimension (N).
*         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
*         (W is not returned).
*
*  X      (input/output) COMPLEX*16 array, dimension (N,T)
*         On an intermediate return, X should be overwritten by
*               A * X,   if KASE=1,
*               A' * X,  if KASE=2,
*         and ZLACN1 must be re-called with all the other parameters
*         unchanged.
*
*  LDX    (input) INTEGER
*         The leading dimension of X.  LDX >= max(1,N).
*
*  XOLD   (workspace) COMPLEX*16 array, dimension (N,T)
*
*  LDXOLD (input) INTEGER
*         The leading dimension of XOLD.  LDXOLD >= max(1,N).
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
*         On the initial call to ZLACN1, KASE should be 0.
*         On an intermediate return, KASE will be 1 or 2, indicating
*         whether X should be overwritten by A * X  or A' * X.
*         On the final return from ZLACN1, KASE will again be 0.
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
*            INFO = 3: power method convergence test.
*            INFO = 4: repeated unit vectors.
*
*  ====================================================================
*
*     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ),
     $                   CONE = ( 1.0D0, 0.0D0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IBEST, IDUMM, ITEMP, ITER, J, JUMP
      DOUBLE PRECISION   ABSXIJ, ESTOLD, SAFMIN, TEMP
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX, IZMAX1
      DOUBLE PRECISION   DLAMCH, DZSUM1
      EXTERNAL           DLAMCH, DZSUM1, IDAMAX, IZMAX1
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLACPY, DLAPST, ZLARNV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DIMAG
*     ..
*     .. Save statement ..
      SAVE
*     ..
*     .. Executable Statements ..
*
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      IF( KASE .EQ. 0 ) THEN
*
         ESTOLD = ZERO
         ITER = 1
         ITEMP = 1
         INFO = 0
*
         DO 10 I = 1, N
            X( I, 1 ) = DCMPLX( ONE )
            IND( I ) = I
            INDH( I ) = 0
   10    CONTINUE    
*
         DO 20 J = 2, T
            CALL ZLARNV( 5, ISEED, N, X( 1, J ) )
   20    CONTINUE         
*
         CALL ZLASCL( 'G', IDUMM, IDUMM, DBLE(N), ONE, N, T, X, LDX, 
     $                INFO )
         KASE = 1
         JUMP = 1
         RETURN
      END IF     
*
      GO TO ( 30, 70 ) JUMP
*
*     ................ ENTRY   (JUMP = 1)
*     FIRST HALF OF THE ITERATION: X HAS BEEN OVERWRITTEN BY A*X.
*
   30 CONTINUE
*
      IF ( ITER .EQ. 1  .AND.  N .EQ. 1 ) THEN
         V( 1 ) = X( 1, 1 )
         EST = ABS( V( 1 ) )
*        ... QUIT
         GO TO 180
      END IF
*     
      EST = ZERO
      DO 40 J = 1, T
         TEMP = DZSUM1( N, X( 1, J ), 1 )
         IF ( TEMP .GT. EST ) THEN
             EST = TEMP
             ITEMP = J 
         END IF
   40 CONTINUE
*
      IF ( EST .GT. ESTOLD  .OR.  ITER .EQ. 2 ) THEN
         IBEST = IND( ITEMP )
      END IF
*
      IF ( EST .LE. ESTOLD  .AND.  ITER .GE. 2 ) THEN
         EST = ESTOLD
         INFO = 2
         GO TO 180
      END IF
*
      ESTOLD = EST
      CALL ZCOPY( N, X( 1, ITEMP ), 1, V, 1 )
*
      IF ( ITER .GT. ITMAX ) THEN
         INFO = 1
         GO TO 180
      END IF
*
*     COMPUTING THE SIGN MATRIX 
*
      DO 60 J = 1, T
         DO 50 I = 1, N
            ABSXIJ = ABS( X( I, J ) ) 
            IF (ABSXIJ .GT. SAFMIN) THEN
               X( I, J ) = DCMPLX( DBLE( X( I, J ) ) / ABSXIJ,
     $                     DIMAG( X( I, J ) ) / ABSXIJ )
            ELSE
               X( I, J ) = CONE
            END IF  
   50    CONTINUE
   60 CONTINUE
*
      CALL ZLACPY( 'Whole', N, T, X, LDX, XOLD, LDXOLD )
*
      KASE = 2
      JUMP = 2
      RETURN
*
*     ................ ENTRY   (JUMP = 2)
*     SECOND HALF OF THE ITERATION: 
*                      X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
*
   70 CONTINUE
*
      DO 80 I = 1, N
         H( I ) = ABS( X ( I, IZMAX1( T, X( I, 1 ), N ) ) )
         IND( I ) = I
   80 CONTINUE
*
      IF ( ITER .GE. 2 .AND. H( IDAMAX(N, H, 1) ) .EQ. H(IBEST) ) THEN
         INFO = 3
         GO TO 180
      END IF
*
*     Sort so that h(i) >= h(j) for i < j
*
      CALL DLAPST( 'D', N, H, IND, ITEMP )
* 
      IF ( ITER .EQ. 1 ) THEN
         ITEMP = T
         GO TO 140
      END IF
*
*     IF IND(1:T) IS CONTAINED IN INDH, TERMINATE.
*
      IF ( T .GT. 1 ) THEN
         DO 100 J = 1, T
            DO 90 I = 1, (ITER-1)*T
               IF (I .GT. N .OR. IND( J ) .EQ. INDH( I )) GO TO 100
   90       CONTINUE
            GO TO 110
  100    CONTINUE
         INFO = 4
         GO TO 180
  110    CONTINUE   
*
*        REPLACE IND(1:T) BY THE FIRST T INDICES IN IND THAT
*        ARE NOT IN INDH. 
*
         ITEMP = 1
         DO 130 J = 1, N
            DO 120 I = 1, (ITER-1)*T
               IF ( I .GT. N .OR. IND( J ) .EQ. INDH( I ) ) GO TO 130
  120       CONTINUE
            IND( ITEMP ) = IND( J )
            IF ( ITEMP .EQ. T ) GO TO 140
            ITEMP = ITEMP + 1
  130    CONTINUE   
      END IF
*
      ITEMP = ITEMP - 1
*
  140 CONTINUE
*
      IF ( (ITER-1)*T .GE. N ) THEN
         DO 150 J = 1, ITEMP
            INDH( (ITER-1)*T+J ) = IND( J )
  150    CONTINUE
      END IF
*
      DO 170 J = 1, T
         DO 160 I = 1, N
            X( I, J ) = CZERO
  160    CONTINUE
         X( IND( J ), J ) = CONE
  170 CONTINUE
*
      ITER = ITER + 1
*
      KASE = 1
      JUMP = 1
      RETURN
*
  180 CONTINUE
      KASE = 0
      RETURN
*
*     End of ZLACN1
*
      END
