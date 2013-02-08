      SUBROUTINE DLARPC( N, T, X, LDX, XOLD, LDXOLD, WRK, KASE, ISEED)
*
*     .. Scalar Arguments ..
      INTEGER            N, T, LDX, LDXOLD, KASE, ISEED(4)
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   WRK( * ), X( LDX, * ), XOLD( LDXOLD, * )
*     ..
*
*  Purpose
*  =======
*
*  DLARPC looks for and replaces columns of X which are parallel to
*         columns of XOLD and itself.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The number of rows.  N >= 1.
*
*  T      (input) INTEGER
*         The number of columns used at each step.
*
*  X      (input/output) DOUBLE PRECISION array, dimension (N,T)
*         On return, X will have full rank.
*
*  LDX    (input) INTEGER
*         The leading dimension of X.  LDX >= max(1,N).
*
*  XOLD   (input/output) DOUBLE PRECISION array, dimension (N,T)
*         On return, XOLD will have full rank.
*
*  LDXOLD (input) INTEGER
*         The leading dimension of XOLD.  LDXOLD >= max(1,N).
*
*  WRK    (workspace) DOUBLE PRECISION array, dimension (T)
*
*  KASE   (input) INTEGER
*          Check parallel columns within X only when KASE = 0,
*          check both X and XOLD otherwise.
*
*  ISEED  (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, JSTART, PCOL
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DLARNV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, NINT, SIGN
*     ..
*     .. Executable Statements ..
*
      IF ( KASE .EQ. 0 ) THEN
         JSTART = 2
      ELSE 
         JSTART = 1
      END IF
*      
      DO 50 J = JSTART, T
*
         PCOL = 0
*
         IF ( KASE .EQ. 0 ) GO TO 30
   10    CALL DGEMV( 'Transpose', N, T, ONE, XOLD, LDXOLD,
     $                   X( 1, J ), 1, ZERO, WRK, 1)
         IF ( NINT ( ABS ( WRK ( IDAMAX ( T, WRK, 1 ) ) ) ) 
     $        .EQ. N ) THEN
            PCOL = PCOL + 1
            CALL DLARNV( 2, ISEED, N, X( 1, J ) )
            DO 20 I = 1, N
               X( I, J ) = SIGN( ONE, X( I, J ) )
   20       CONTINUE
            IF ( PCOL .GE. N/T ) GO TO 60
            GO TO 10
         END IF
*
         IF ( J .EQ. 1 ) GO TO 50
   30    CALL DGEMV( 'Transpose', N, J-1, ONE, X, LDX,
     $                   X( 1, J ), 1, ZERO, WRK, 1)
         IF ( NINT ( ABS ( WRK ( IDAMAX ( J-1, WRK, 1 ) ) ) ) 
     $        .EQ. N ) THEN
            PCOL = PCOL + 1
            CALL DLARNV( 2, ISEED, N, X( 1, J ) )
            DO 40 I = 1, N
               X( I, J ) = SIGN( ONE, X( I, J ) )
   40       CONTINUE
            IF ( PCOL .GE. N/T ) GO TO 60
            IF ( KASE .EQ. 0 ) THEN
               GO TO 30
            ELSE 
               GO TO 10
            END IF
         END IF
*
   50 CONTINUE
   60 CONTINUE
      RETURN
*
*     End of DLARPC
*
      END
