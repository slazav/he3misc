*UPTODATE X04AAETEXT
      SUBROUTINE X04AAE(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     *** NOTE ***
C     THIS ROUTINE ASSUMES THAT THE VALUE OF NERR1 IS SAVED
C     BETWEEN CALLS.  IN SOME IMPLEMENTATIONS IT MAY BE
C     NECESSARY TO STORE NERR1 IN A LABELLED COMMON
C     BLOCK /AX04AA/ TO ACHIEVE THIS.
C
C     .. SCALAR ARGUMENTS ..
      INTEGER I, NERR
C     ..
C     .. LOCAL SCALARS ..
      INTEGER NERR1
C     FOR ms/FORTRAN MUH
      NERR1=0
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END
