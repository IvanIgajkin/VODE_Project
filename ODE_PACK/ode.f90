SUBROUTINE DSOLVE(NEQ, RPAR, IPAR, Y, T, TSTEP)
    EXTERNAL FEX, JEX
    INTEGER IPAR, i, NRPD
    DOUBLE PRECISION ATOL, RPAR, RTOL, RWORK, T, TOUT, Y, TSTEP
    DIMENSION IPAR(1), Y(NEQ), ATOL(4), RWORK(90), IWORK(34), RPAR(IPAR(1) * (2 * NEQ + 1))
    TOUT = TSTEP
    NRPD = 4
    ITOL = 2
    RTOL = 1.D-4
    ATOL(1) = 1.D-8
    ATOL(2) = 1.D-14
    ATOL(3) = 1.D-6
    ATOL(4) = 1.D-6
    ITASK = 1
    ISTATE = 1
    IOPT = 0
    LRW = 90
    LIW = 34
    MF = 21
    DO IOUT = 1,120
        CALL DVODE(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE, IOPT,RWORK,LRW,IWORK,LIW, JEX, MF,RPAR,IPAR)
        print *, T,(Y(i), i = 1, NEQ)
        IF (ISTATE .LT. 0) GO TO 80
        TOUT = T+TSTEP
    ENDDO
    WRITE(6,60) IWORK(11),IWORK(12),IWORK(13),IWORK(19),IWORK(20),IWORK(21),IWORK(22)
    60  FORMAT(/' No. steps =',I4,'   No. f-s =',I4, &
            '   No. J-s =',I4,'   No. LU-s =',I4,&
            '  No. nonlinear iterations =',I4,&
            '  No. nonlinear convergence failures =',I4,&
            '  No. error test failures =',I4/)
    STOP
    80  WRITE(6,90)ISTATE
    90  FORMAT(///' Error halt: ISTATE =',I3)
    STOP
END

SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
     DOUBLE PRECISION RPAR, T, Y, YDOT
     INTEGER IPAR, M
     DOUBLE PRECISION K, v1, v2, W
     INTEGER cnum, n, i
     DIMENSION IPAR(1)
     DIMENSION Y(NEQ), YDOT(NEQ), RPAR((2 * NEQ + 1) * IPAR(1))
     DIMENSION v1(NEQ), v2(NEQ)
     M = 2 * NEQ + 1
     DO cnum = 1, NEQ
        YDOT(cnum) = 0
        DO n = 1, IPAR(1)
            v1 = RPAR((n - 1) * M + 1 : NEQ + (n - 1) * M)
            v2 = RPAR(NEQ + (n - 1) * M + 1 : n * M - 1)
            K = RPAR(M * n)
            W = K
            DO i = 1, NEQ
                IF (Y(i) > 0 .OR. v1(i) == 1) W = W * Y(i)**v1(i)
            ENDDO
            YDOT(cnum) = YDOT(cnum) + (v2(cnum + (n - 1) * M + NEQ + 1) - v1(cnum + (n - 1) * M)) * W
        ENDDO
     ENDDO
     RETURN
END SUBROUTINE

SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
    INTEGER M, IPAR, cnum, i, n
     DOUBLE PRECISION PD, RPAR, T, Y, K, v1, v2, W
     DIMENSION Y(NEQ), PD(NEQ, NRPD), IPAR(1), RPAR((2 * NEQ + 1) * IPAR(1)), v1(NEQ), v2(NEQ)
     M = 2 * NEQ + 1
     DO cnum = 1, NEQ
        DO i = 1, NRPD
            PD(cnum, i) = 0
        ENDDO
        DO n = 1, IPAR(1)
            v1 = RPAR((n - 1) * M + 1 : NEQ + (n - 1) * M)
            v2 = RPAR(NEQ + (n - 1) * M + 1 : n * M - 1)
            K = RPAR(M * (n - 1))
            W = K
            DO i = 1, NEQ
                IF (Y(i) > 0 .OR. v1(i) == 1) W = W * Y(i)**v1(i)
            ENDDO
            PD(cnum, 1) = PD(cnum, 1) + (v2(cnum + (n - 1) * M + NEQ + 1) - v1(cnum + (n - 1) * M)) * W
            PRINT *, (PD(i, 1), i = 1, NEQ)
        ENDDO
     ENDDO
     RETURN
END