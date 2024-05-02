!----------------------------------------------------------------------
!    UMAT FOR FINITE VISCOELASTIC THEORY - REESE & GOVINDJEE
!    WITH OGDEN HYPERELASTICITY MODEL
!    IMPLEMENTED AS PART OF MINI THESIS  
!    - ALAN J CORREA
!----------------------------------------------------------------------
!    PROPS(1) - MU
!    PROPS(2) - ALPHA
!    PROPS(3) - KELAS
!    PROPS(4) - MUVIS
!    PROPS(5) - ALPHAVIS
!    PROPS(6) - KVIS
!    PROPS(7) - ETADEV
!    PROPS(8) - ETAVOL
!----------------------------------------------------------------------
      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)
!
      INCLUDE 'ABA_PARAM.INC'
!
      CHARACTER*8 CMNAME
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2 PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DROT(3, 3),
     3 DFGRD0(3, 3), DFGRD1(3, 3)

!
! ---------------------------------------------------------------------
!     LOCAL ARRAYS
! ---------------------------------------------------------------------
!     BeOLD         - VISCOUS LEFT CAUCHY-GREEN DEFORMATION TENSOR (AT N-1)
!     BeOLDINV      - INVERSE(BeOLD) 
!     CiOLD         - VISCOUS RIGHT CAUCHY-GREEN DEFORMATION TENSOR (AT N-1)
!     CiOLDINV      - INVERSE(CiOLD)
!     BeTR          - NEQ TRIAL LEFT CAUCHY-GREEN TENSOR
!     BeTR_         - COPY OF BeTR FOR EIGEN VALUE\VECTOR DECOMPOSITION 
!     Be            - NEQ LEFT CAUCHY-GREEN TENSOR 
!     TAUNEQ        - NEQ KIRCHOFF STRESS TENSOR
!     PVBeTR        - EIGEN VALUES OF BeTR
!     PDBeTR        - EIGEN VECTORS OF BeTR
!     EPSeTR        - TRIAL NEQ PRINCIPAL(EIGEN) STRAIN
!     PVBe          - EIGEN VALUES OF Be 
!     EPSe          - NEQ LOGARITHMIC PRINCIPAL(EIGEN) STRAIN 
!     PVBeBAR       - EIGEN VALUES OF NEQ DEVIATORIC LEFT-CAUCHY GREEN TENSOR
!     DEVTAU        - EIGEN VALUES OF DEVIATORIC NEQ KIRCHOFF STRESS TENSOR
!     RESVEC        - RESIDUAL VECTOR FOR NEWTON RAPHSON'S ITERATION
!     DDEVTAUDEPSe  - D(DEVTAU) / D(EPSe)
!     KMAT          - K MATRIX FOR NEWTON RAPHSON'S ITERATION
!     KINV          - INVERSE OF K MATRIX
!     DELEPSe       - DELTA EPSe FOR NEWTON RAPHSON'S ITERATION 
!     PVTAU         - EIGEN VALUES OF TAUNEQ
!     PDTAU         - EIGEN VECTORS OF TAUNEQ
!     DPVTAUDEPSe   - D(PVTAU) / D(EPSe)
!     CALG          - NEQ ALGORITHMIC TANGENT MODULUS
!     SIGMANEQ      - NEQ CAUCHY STRESS TENSOR
!     L4NEQ         - NEQ INTERMEDIATE MATERIAL TANGENT STIFFNESS MODULUS
!     C4NEQ         - NEQ SPATIAL TANGENT STIFFNESS MODULUS
!     C4NEQJ        - NEQ SPATIAL TANGENT STIFFNESS MODULUS (FE)
!     CNEQ          - NEQ SPATIAL TANGENT STIFFNESS MODULUS (VOIGT)
!     BTOT          - EQ LEFT CAUCHY-GREEN TENSOR
!     BTOT_         - COPY OF BTOT FOR EIGEN VALUE/VECTOR DECOMPOSITION
!     PVBTOT        - EIGEN VALUES OF BTOT
!     PDBTOT        - EIGEN VECTORS OF BTOT
!     PVBBAR        - EIGEN VALUES OF EQ DEVIATORIC LEFT CAUCHY-GREEN TENSOR
!     PVTAUEQ       - EIGEN VALUES OF EQ KIRCHOFF STRESS TENSOR
!     SIGMAEQ       - EQ CAUCHY STRESS TENSOR
!     CAB           - CAB FOR EQ SPATIAL TANGENT STIFFNESS MODULUS
!     GAB           - GAB FOR EQ SPATIAL TANGENT STIFFNESS MODULUS
!     C4EQMAT       - COEFFICIENTS OF EQ SPATIAL TANGENT STIFFNESS MODULUS 
!     C4EQ          - EQ SPATIAL TANGENT STIFFNESS MODULUS
!     C4EQJ         - EQ SPATIAL TANGENT STIFFNESS MODULUS (FE)
!     CEQ           - EQ SPATIAL TANGENT STIFFNESS MODULUS (VOIGT)
!     STRESSTOT     - TOTAL CAUCHY STRESS
!     IDT2          - 2ND ORDER IDENTTITY TENSOR
! ---------------------------------------------------------------------
      DOUBLE PRECISION, DIMENSION(3, 3) :: BeOLD
      DOUBLE PRECISION, DIMENSION(3, 3) :: BeOLDINV
      DOUBLE PRECISION, DIMENSION(3, 3) :: CiOLD
      DOUBLE PRECISION, DIMENSION(3, 3) :: CiOLDINV
      DOUBLE PRECISION, DIMENSION(3, 3) :: BeTR
      DOUBLE PRECISION, DIMENSION(3, 3) :: BeTR_
      DOUBLE PRECISION, DIMENSION(6) :: BeTRV
      DOUBLE PRECISION, DIMENSION(3, 3) :: Be
      DOUBLE PRECISION, DIMENSION(3,3) :: TAUNEQ      
      DOUBLE PRECISION, DIMENSION(3) :: PVBeTR
      DOUBLE PRECISION, DIMENSION(3,3) :: PDBeTR
      DOUBLE PRECISION, DIMENSION(3) :: EPSeTR
      DOUBLE PRECISION, DIMENSION(3) :: PVBe
      DOUBLE PRECISION, DIMENSION(3) :: EPSe
      DOUBLE PRECISION, DIMENSION(3) :: PVBeBAR
      DOUBLE PRECISION, DIMENSION(3) :: DEVTAU
      DOUBLE PRECISION, DIMENSION(3) :: RESVEC
      DOUBLE PRECISION, DIMENSION(3,3) :: DDEVTAUDEPSe
      DOUBLE PRECISION, DIMENSION(3,3) :: KMAT
      DOUBLE PRECISION, DIMENSION(3,3) :: KINV
      DOUBLE PRECISION, DIMENSION(3) :: DELEPSe
      DOUBLE PRECISION, DIMENSION(3) :: PVTAU
      DOUBLE PRECISION, DIMENSION(3,3) :: PDTAU
      DOUBLE PRECISION, DIMENSION(3,3) :: DPVTAUDEPSe
      DOUBLE PRECISION, DIMENSION(3,3) :: CALG
      DOUBLE PRECISION, DIMENSION(3,3) :: SIGMANEQ
      DOUBLE PRECISION, DIMENSION(3,3,3,3) :: L4NEQ
      DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C4NEQ
      DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C4NEQJ
      DOUBLE PRECISION, DIMENSION(6,6) :: CNEQ
      DOUBLE PRECISION, DIMENSION(3,3) :: BTOT
      DOUBLE PRECISION, DIMENSION(3,3) :: BTOT_
      DOUBLE PRECISION, DIMENSION(3) :: PVBTOT
      DOUBLE PRECISION, DIMENSION(3,3) :: PDBTOT
      DOUBLE PRECISION, DIMENSION(3) :: PVBBAR
      DOUBLE PRECISION, DIMENSION(3) :: PVTAUEQ
      DOUBLE PRECISION, DIMENSION(3,3) :: SIGMAEQ
      DOUBLE PRECISION, DIMENSION(3,3) :: CAB
      DOUBLE PRECISION, DIMENSION(3,3) :: GAB
      DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C4EQMAT
      DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C4EQ
      DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C4EQJ
      DOUBLE PRECISION, DIMENSION(6,6) :: CEQ
      DOUBLE PRECISION, DIMENSION(3,3) :: STRESSTOT
      DOUBLE PRECISION, DIMENSION(3,3) :: IDT2
!
! ---------------------------------------------------------------------
!     LOCAL VARIABLES
! ---------------------------------------------------------------------
!     DET      - JACOBIAN OF TOTAL DEFORMATION GRADIENT   
!     Je       - JACOBIAN OF NEQ DEFORMATION GRADIENT
!     NORMRES  - NORM OF RESIDUAL VECTOR
!     RESTOL   - TOLERANCE FOR NEWTON RHAPSON CONVERGENCE
!     EPS      - TOLERANCE FOR EQUIVALENCE OF FLOAT VALUES   
!     ITER     - LOCAL ITERATION NUMBER
! ---------------------------------------------------------------------
      DOUBLE PRECISION MU
      DOUBLE PRECISION ALPHA
      DOUBLE PRECISION KELAS
      DOUBLE PRECISION MUVIS
      DOUBLE PRECISION ALPHAVIS
      DOUBLE PRECISION KVIS
      DOUBLE PRECISION ETADEV
      DOUBLE PRECISION ETAVOL
      DOUBLE PRECISION DET
      DOUBLE PRECISION Je
      DOUBLE PRECISION NORMRES
      INTEGER ITER, I, J, K, L, M, N
      LOGICAL OK_FLAG
      DOUBLE PRECISION :: RESTOL = 1.0D-12
      DOUBLE PRECISION :: EPS = 1.0D-12
!
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0,
     1 SIX=6.D0, NINE=9.0D0)
!
! ---------------------------------------------------------------------
!     STEP 1 - COMPUTING NECESSARY VARIABLES AND CONSTANTS
! ---------------------------------------------------------------------
! 
!     DEFINING INTERNAL STATE VARIABLE FOR FIRST ITERATION
      IF (KSTEP.EQ.1 .AND. KINC.EQ.1) THEN
            STATEV(1) = ONE
            STATEV(2) = ONE
            STATEV(3) = ONE
            STATEV(4) = ZERO
            STATEV(5) = ZERO
            STATEV(6) = ZERO
      END IF
!     
!     GET INTERNAL STATE-VARIABLE FROM PREVIOUS ITERATION
      BeOLD(1,1) = STATEV(1)
      BeOLD(1,2) = STATEV(4)
      BeOLD(1,3) = STATEV(5)
      BeOLD(2,1) = STATEV(4)
      BeOLD(2,2) = STATEV(2)
      BeOLD(2,3) = STATEV(6)
      BeOLD(3,1) = STATEV(5)
      BeOLD(3,2) = STATEV(6)
      BeOLD(3,3) = STATEV(3)
! 
!     ASSIGN PROPERTIES TO CONSTANTS
      MU          = PROPS(1)
      ALPHA       = PROPS(2)
      KELAS       = PROPS(3)
      MUVIS       = PROPS(4)
      ALPHAVIS    = PROPS(5)
      KVIS        = PROPS(6)
      ETADEV      = PROPS(7)
      ETAVOL      = PROPS(8)
!
!     DEFINE SECOND ORDER IDENTITY TENSOR
      IDT2 = ZERO
      DO I = 1,3
        IDT2(I,I) = ONE
      END DO
!
!     CALCULATE JACOBIAN OF DEFORMATION GRADIENT
      DET = DFGRD1(1,1)*DFGRD1(2,2)*DFGRD1(3,3)
     1     -DFGRD1(1,2)*DFGRD1(2,1)*DFGRD1(3,3)
      IF(NSHR.EQ.3) THEN
        DET=DET+DFGRD1(1,2)*DFGRD1(2,3)*DFGRD1(3,1)
     1         +DFGRD1(1,3)*DFGRD1(3,2)*DFGRD1(2,1)
     2         -DFGRD1(1,3)*DFGRD1(3,1)*DFGRD1(2,2)
     3         -DFGRD1(2,3)*DFGRD1(3,2)*DFGRD1(1,1)
      END IF
!
! ---------------------------------------------------------------------
!     STEP 2 - CALCULATION OF TRIAL LEFT CAUCHY-GREEN TENSOR
! ---------------------------------------------------------------------
!
      CALL M33INV(BeOLD, BeOLDINV, OK_FLAG)
      CiOLD = MATMUL(TRANSPOSE(DFGRD0), MATMUL(BeOLDINV, DFGRD0))
      CALL M33INV(CiOLD, CiOLDINV, OK_FLAG)
      BeTR  = MATMUL(DFGRD1, MATMUL(CiOLDINV, TRANSPOSE(DFGRD1)))
!
! ---------------------------------------------------------------------
!     STEP 3 - LOCAL NEWTON RHAPSON FOR NEQ EVOLUTION EQUATION
! ---------------------------------------------------------------------
!     IN:
!           BeTR
!           DTIME
!           MUVIS, ALPHAVIS, KVIS, ETADEV, ETAVOL
!
!     OUT:
!           Be, TAUNEQ
!           PVBeTR, PVTAU, PDTAU
!           CALG
! ---------------------------------------------------------------------
!
!     GET EIGEN VALUES AND EIGEN VECTORS OF BeTR
      DO I = 1,3
        DO J = 1,3
          BeTR_(I,J) = BeTR(I,J)
        END DO
      END DO
      CALL DSYEVJ3(BeTR_, PDBeTR, PVBeTR)
      PDBeTR = TRANSPOSE(PDBeTR)
!
!     CALCULATE EPSeTR
      DO I = 1,3
        EPSeTR(I)=ONE/TWO*LOG(PVBeTR(I))
      END DO
!
!     PERFORM LOCAL ITERATION
!     1. Initialize Iteration Variables
      DO I = 1,3
        PVBe(I) = PVBeTR(I)
        EPSe(I) = EPSeTR(I)
      END DO
!     2. Newton Rhapsons Iteration with MAXITER=200
      DO ITER = 1, 200 
!       - Calculating Jacobian
        Je = (PVBe(1)*PVBe(2)*PVBe(3))**(ONE/TWO)
!       - Calculating Deviatoric Principal Values of Be
        PVBeBAR(1) = Je**(-TWO/THREE)*PVBe(1)
        PVBeBAR(2) = Je**(-TWO/THREE)*PVBe(2)
        PVBeBAR(3) = Je**(-TWO/THREE)*PVBe(3)
!       - Calculating Principal Values of Deviatoric Kirchoff Stress
        DEVTAU(1) = MUVIS * ((TWO/THREE)*(PVBeBAR(1)**(ALPHAVIS/TWO))
     1                      -(ONE/THREE)*(PVBeBAR(2)**(ALPHAVIS/TWO))
     2                      -(ONE/THREE)*(PVBeBAR(3)**(ALPHAVIS/TWO)))
        DEVTAU(2) = MUVIS * ((TWO/THREE)*(PVBeBAR(2)**(ALPHAVIS/TWO))
     1                      -(ONE/THREE)*(PVBeBAR(3)**(ALPHAVIS/TWO))
     2                      -(ONE/THREE)*(PVBeBAR(1)**(ALPHAVIS/TWO)))
        DEVTAU(3) = MUVIS * ((TWO/THREE)*(PVBeBAR(3)**(ALPHAVIS/TWO))
     1                      -(ONE/THREE)*(PVBeBAR(1)**(ALPHAVIS/TWO))
     2                      -(ONE/THREE)*(PVBeBAR(2)**(ALPHAVIS/TWO)))
!       - Calculating Residual Vector
        DO N = 1, 3
          RESVEC(N) = EPSe(N) - EPSeTR(N) 
     1                + DTIME * (ONE/(TWO*ETADEV)*DEVTAU(N)
     2                          +KVIS/(THREE*FOUR*ETAVOL)*(Je*Je-ONE))
        END DO
!       - Calculating Norm of Residual Vector
        NORMRES = (RESVEC(1)**TWO 
     1            +RESVEC(2)**TWO
     2            +RESVEC(3)**TWO)**(ONE/TWO)
!       - Terminate if Norm of Residual Vector is less than Tolerance
        IF ((ITER.GT.1).AND.ABS(NORMRES).LT.RESTOL) THEN
            EXIT
        END IF
!       - Calculating dDEVTAU/dEPSe
        DDEVTAUDEPSe(1,1) = MUVIS * ALPHAVIS 
     1                   *((FOUR/NINE)*(PVBeBAR(1)**(ALPHAVIS/TWO))
     2                     +(ONE/NINE)*(PVBeBAR(2)**(ALPHAVIS/TWO)
     3                                 +PVBeBAR(3)**(ALPHAVIS/TWO)))
        DDEVTAUDEPSe(2,2) = MUVIS * ALPHAVIS 
     1                   *((FOUR/NINE)*(PVBeBAR(2)**(ALPHAVIS/TWO))
     2                     +(ONE/NINE)*(PVBeBAR(3)**(ALPHAVIS/TWO)
     3                                 +PVBeBAR(1)**(ALPHAVIS/TWO))) 
        DDEVTAUDEPSe(3,3) = MUVIS * ALPHAVIS 
     1                   *((FOUR/NINE)*(PVBeBAR(3)**(ALPHAVIS/TWO))
     2                     +(ONE/NINE)*(PVBeBAR(1)**(ALPHAVIS/TWO)
     3                                 +PVBeBAR(2)**(ALPHAVIS/TWO)))
        DDEVTAUDEPSe(1,2) = MUVIS * ALPHAVIS 
     1                   *((-TWO/NINE)*(PVBeBAR(1)**(ALPHAVIS/TWO)
     2                                 +PVBeBAR(2)**(ALPHAVIS/TWO))
     3                     +(ONE/NINE)*(PVBeBAR(3)**(ALPHAVIS/TWO)))
        DDEVTAUDEPSe(1,3) = MUVIS * ALPHAVIS 
     1                   *((-TWO/NINE)*(PVBeBAR(1)**(ALPHAVIS/TWO)
     2                                 +PVBeBAR(3)**(ALPHAVIS/TWO))
     3                     +(ONE/NINE)*(PVBeBAR(2)**(ALPHAVIS/TWO)))
        DDEVTAUDEPSe(2,3) = MUVIS * ALPHAVIS 
     1                   *((-TWO/NINE)*(PVBeBAR(2)**(ALPHAVIS/TWO)
     2                                 +PVBeBAR(3)**(ALPHAVIS/TWO))
     3                     +(ONE/NINE)*(PVBeBAR(1)**(ALPHAVIS/TWO)))
        DDEVTAUDEPSe(2,1) = DDEVTAUDEPSe(1,2)
        DDEVTAUDEPSe(3,1) = DDEVTAUDEPSe(1,3)
        DDEVTAUDEPSe(3,2) = DDEVTAUDEPSe(2,3)
!       - Calculating K Matrix for Newton Rhapson
        DO I = 1,3
          DO J = 1,3
            KMAT(I,J) = IDT2(I,J)
     1               + DTIME * ((ONE/(TWO*ETADEV))*DDEVTAUDEPSe(I,J) 
     2                         +(ONE/(THREE*ETAVOL)*KVIS*Je*Je))
          END DO
        END DO
!       - Calculating KINV
        CALL M33INV(KMAT, KINV, OK_FLAG)
!       - Calculating NEQ Strain Increment
        DELEPSe(1) = -(KINV(1,1)*RESVEC(1)
     1                +KINV(1,2)*RESVEC(2)
     2                +KINV(1,3)*RESVEC(3))
                      
        DELEPSe(2) = -(KINV(2,1)*RESVEC(1)
     1                +KINV(2,2)*RESVEC(2)
     2                +KINV(2,3)*RESVEC(3))
        
        DELEPSe(3) = -(KINV(3,1)*RESVEC(1)
     1                +KINV(3,2)*RESVEC(2)
     2                +KINV(3,3)*RESVEC(3))
!       - Calculating Updated NEQ Strain
        DO I = 1,3
        EPSe(I) = EPSe(I) + DELEPSe(I)
        END DO
!       - Updating Eigen Value of Be
        DO I = 1,3
          PVBe(I) = EXP(TWO*EPSe(I))
        END DO
      END DO
!     3. Non-Covergence Test for Newton Rhapson
      IF (ITER.GT.200) THEN
            WRITE(*,*) "LOCAL ITERATION DID NOT CONVERGE" 
      END IF
!
!     CALCULATE REQUIRED QUANTITIES
!     1. Be
      DO I=1,3
        DO J=1,3
            Be(J,I)=PVBe(1)*PDBeTR(1,J)*PDBeTR(1,I)
     1             +PVBe(2)*PDBeTR(2,J)*PDBeTR(2,I)
     2             +PVBe(3)*PDBeTR(3,J)*PDBeTR(3,I)
        END DO
      END DO
! 
!     2. PVTAU
      DO I=1,3
            PVTAU(I) = DEVTAU(I) + KVIS/TWO*(Je*Je-ONE)
      END DO
! 
!     3. PDTAU
      DO I=1,3
        DO J=1,3
          PDTAU(J,I) = PDBeTR(J,I)
        END DO
      END DO
! 
!     4. TAUNEQ
      DO I=1,3
        DO J=1,3
          TAUNEQ(J,I)=PVTAU(1)*PDTAU(1,J)*PDTAU(1,I)
     1               +PVTAU(2)*PDTAU(2,J)*PDTAU(2,I)
     2               +PVTAU(3)*PDTAU(3,J)*PDTAU(3,I)
        END DO
      END DO
! 
!     5. CALG
!     - Calculating Calculating dPVTAU/dEPSe
      DO I=1,3
        DO J=1,3
          DPVTAUDEPSe(J,I) = DDEVTAUDEPSe(J,I) + KVIS*Je*Je
        END DO
      END DO
!     - Calculating CALG
      CALG = MATMUL(DPVTAUDEPSe, KINV)
!     
! ---------------------------------------------------------------------
!     STEP 4 - UPDATE INTERNAL STATE VARIABLE - Be
! ---------------------------------------------------------------------
!     SET INTERNAL STATE-VARIABLE FROM FINAL LOCAL ITERATION
      STATEV(1) = Be(1,1)
      STATEV(2) = Be(2,2)
      STATEV(3) = Be(3,3)
      STATEV(4) = Be(1,2)
      STATEV(5) = Be(1,3)
      STATEV(6) = Be(2,3)
! ---------------------------------------------------------------------
!     STEP 5 - CALCULATION OF NEQ CAUCHY STRESS & ELASTICITY TENSOR
! ---------------------------------------------------------------------
!     CALCULATE NEQ CAUCHY STRESS TENSOR
      SIGMANEQ = ZERO
      DO J=1,3
        DO K=1,3
          SIGMANEQ(J,K) = TAUNEQ(J,K)/DET
        END DO
      END DO
!
!     CALCULATE MATERIAL TENSOR (L4NEQ) IN INTERMEDIATE CONFIGURATION 
!     EXPRESSED IN EIGEN BASIS 
      L4NEQ = ZERO
!     1. CALCULATE COEFFICIENTS FOR AABB
      DO J=1,3
        DO K=1,3
          L4NEQ(J,J,K,K) = (CALG(J,K) - PVTAU(J)*TWO*IDT2(J,K))
     1                    /(PVBeTR(J)*PVBeTR(K))
        END DO
      END DO 
!     2. CALCULATE COEFFICIENTS FOR ABAB / ABBA
      DO J=1,3
        DO K=1,3
          IF (J.EQ.K) THEN
            CYCLE
          END IF
!         - In case of Equal PVBeTR, Use Alternate Form
          IF (ABS(PVBeTR(J)-PVBeTR(K)).LT.EPS) THEN
            L4NEQ(J,K,J,K) = (L4NEQ(J,J,J,J)-L4NEQ(J,J,K,K))/TWO
          ELSE
            L4NEQ(J,K,J,K) = ((PVTAU(K)/PVBeTR(K))-(PVTAU(J)/PVBeTR(J)))
     1                      /(PVBeTR(K)-PVBeTR(J))
          END IF
          L4NEQ(J,K,K,J) = L4NEQ(J,K,J,K)
        END DO
      END DO
!     CALCULATE MATERIAL TENSOR (C4NEQ) IN CURRENT CONFIGURATION
!     EXPRESSED IN STANDARD BASIS
      C4NEQ = ZERO
      DO J=1,3
        DO K=1,3
          DO L=1,3
            DO M=1,3
              C4NEQ(J,K,L,M) = 
     1          L4NEQ(1,1,1,1)*PVBeTR(1)*PVBeTR(1)
     1                        *PDBeTR(1,J)*PDBeTR(1,K)
     1                        *PDBeTR(1,L)*PDBeTR(1,M)
     1         +L4NEQ(1,1,2,2)*PVBeTR(1)*PVBeTR(2)
     1                        *PDBeTR(1,J)*PDBeTR(1,K)
     1                        *PDBeTR(2,L)*PDBeTR(2,M)
     1         +L4NEQ(1,1,3,3)*PVBeTR(1)*PVBeTR(3)
     1                        *PDBeTR(1,J)*PDBeTR(1,K)
     1                        *PDBeTR(3,L)*PDBeTR(3,M)
     1         +L4NEQ(2,2,1,1)*PVBeTR(2)*PVBeTR(1)
     1                        *PDBeTR(2,J)*PDBeTR(2,K)
     1                        *PDBeTR(1,L)*PDBeTR(1,M)
     1         +L4NEQ(2,2,2,2)*PVBeTR(2)*PVBeTR(2)
     1                        *PDBeTR(2,J)*PDBeTR(2,K)
     1                        *PDBeTR(2,L)*PDBeTR(2,M)
     1         +L4NEQ(2,2,3,3)*PVBeTR(2)*PVBeTR(3)
     1                        *PDBeTR(2,J)*PDBeTR(2,K)
     1                        *PDBeTR(3,L)*PDBeTR(3,M)
     1         +L4NEQ(3,3,1,1)*PVBeTR(3)*PVBeTR(1)
     1                        *PDBeTR(3,J)*PDBeTR(3,K)
     1                        *PDBeTR(1,L)*PDBeTR(1,M)
     1         +L4NEQ(3,3,2,2)*PVBeTR(3)*PVBeTR(2)
     1                        *PDBeTR(3,J)*PDBeTR(3,K)
     1                        *PDBeTR(2,L)*PDBeTR(2,M)
     1         +L4NEQ(3,3,3,3)*PVBeTR(3)*PVBeTR(3)
     1                        *PDBeTR(3,J)*PDBeTR(3,K)
     1                        *PDBeTR(3,L)*PDBeTR(3,M)
     1         +L4NEQ(1,2,1,2)*PVBeTR(1)*PVBeTR(2)
     1                        *PDBeTR(1,J)*PDBeTR(2,K)
     1                        *PDBeTR(1,L)*PDBeTR(2,M)
     1         +L4NEQ(1,2,2,1)*PVBeTR(1)*PVBeTR(2)
     1                        *PDBeTR(1,J)*PDBeTR(2,K)
     1                        *PDBeTR(2,L)*PDBeTR(1,M)
     1         +L4NEQ(1,3,1,3)*PVBeTR(1)*PVBeTR(3)
     1                        *PDBeTR(1,J)*PDBeTR(3,K)
     1                        *PDBeTR(1,L)*PDBeTR(3,M)
     1         +L4NEQ(1,3,3,1)*PVBeTR(1)*PVBeTR(3)
     1                        *PDBeTR(1,J)*PDBeTR(3,K)
     1                        *PDBeTR(3,L)*PDBeTR(1,M)
     1         +L4NEQ(2,1,2,1)*PVBeTR(2)*PVBeTR(1)
     1                        *PDBeTR(2,J)*PDBeTR(1,K)
     1                        *PDBeTR(2,L)*PDBeTR(1,M)
     1         +L4NEQ(2,1,1,2)*PVBeTR(2)*PVBeTR(1)
     1                        *PDBeTR(2,J)*PDBeTR(1,K)
     1                        *PDBeTR(1,L)*PDBeTR(2,M)
     1         +L4NEQ(2,3,2,3)*PVBeTR(2)*PVBeTR(3)
     1                        *PDBeTR(2,J)*PDBeTR(3,K)
     1                        *PDBeTR(2,L)*PDBeTR(3,M)
     1         +L4NEQ(2,3,3,2)*PVBeTR(2)*PVBeTR(3)
     1                        *PDBeTR(2,J)*PDBeTR(3,K)
     1                        *PDBeTR(3,L)*PDBeTR(2,M)
     1         +L4NEQ(3,1,3,1)*PVBeTR(3)*PVBeTR(1)
     1                        *PDBeTR(3,J)*PDBeTR(1,K)
     1                        *PDBeTR(3,L)*PDBeTR(1,M)
     1         +L4NEQ(3,1,1,3)*PVBeTR(3)*PVBeTR(1)
     1                        *PDBeTR(3,J)*PDBeTR(1,K)
     1                        *PDBeTR(1,L)*PDBeTR(3,M)
     1         +L4NEQ(3,2,3,2)*PVBeTR(3)*PVBeTR(2)
     1                        *PDBeTR(3,J)*PDBeTR(2,K)
     1                        *PDBeTR(3,L)*PDBeTR(2,M)
     1         +L4NEQ(3,2,2,3)*PVBeTR(3)*PVBeTR(2)
     1                        *PDBeTR(3,J)*PDBeTR(2,K)
     1                        *PDBeTR(2,L)*PDBeTR(3,M)
            END DO
          END DO
        END DO
      END DO
!
!     CALCLUATE JAUMANN RATE FORM OF MATERIAL TENSOR
      C4NEQJ = ZERO
      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO L=1,3
              C4NEQJ(I,J,K,L) = C4NEQ(I,J,K,L)/DET
     1                        +(IDT2(I,K)*TAUNEQ(J,L)
     1                         +IDT2(I,L)*TAUNEQ(J,K)
     1                         +IDT2(J,K)*TAUNEQ(I,L)
     1                         +IDT2(J,L)*TAUNEQ(I,K))/(TWO*DET)        
            END DO
          END DO
        END DO
      END DO
!
!     GENERATE ABAQUS NEQ TANGENT STIFFNESS MATRIX (VOIGT NOTATION)
      CNEQ = ZERO
      CNEQ(1,1) = C4NEQJ(1,1,1,1)
      CNEQ(1,2) = C4NEQJ(1,1,2,2)
      CNEQ(1,3) = C4NEQJ(1,1,3,3)
      CNEQ(1,4) = C4NEQJ(1,1,1,2)
      CNEQ(1,5) = C4NEQJ(1,1,1,3)
      CNEQ(1,6) = C4NEQJ(1,1,2,3)
      CNEQ(2,1) = CNEQ(1,2)
      CNEQ(2,2) = C4NEQJ(2,2,2,2)
      CNEQ(2,3) = C4NEQJ(2,2,3,3)
      CNEQ(2,4) = C4NEQJ(2,2,1,2)
      CNEQ(2,5) = C4NEQJ(2,2,1,3)
      CNEQ(2,6) = C4NEQJ(2,2,2,3)
      CNEQ(3,1) = CNEQ(1,3)
      CNEQ(3,2) = CNEQ(2,3)
      CNEQ(3,3) = C4NEQJ(3,3,3,3)
      CNEQ(3,4) = C4NEQJ(3,3,1,2)
      CNEQ(3,5) = C4NEQJ(3,3,1,3)
      CNEQ(3,6) = C4NEQJ(3,3,2,3)
      CNEQ(4,1) = CNEQ(1,4)
      CNEQ(4,2) = CNEQ(2,4)
      CNEQ(4,3) = CNEQ(3,4)
      CNEQ(4,4) = C4NEQJ(1,2,1,2)
      CNEQ(4,5) = C4NEQJ(1,2,1,3)
      CNEQ(4,6) = C4NEQJ(1,2,2,3)
      CNEQ(5,1) = CNEQ(1,5)
      CNEQ(5,2) = CNEQ(2,5)
      CNEQ(5,3) = CNEQ(3,5)
      CNEQ(5,4) = CNEQ(4,5)
      CNEQ(5,5) = C4NEQJ(1,3,1,3)
      CNEQ(5,6) = C4NEQJ(1,3,2,3)
      CNEQ(6,1) = CNEQ(1,6)
      CNEQ(6,2) = CNEQ(2,6)
      CNEQ(6,3) = CNEQ(3,6)
      CNEQ(6,4) = CNEQ(4,6)
      CNEQ(6,5) = CNEQ(5,6)
      CNEQ(6,6) = C4NEQJ(2,3,2,3)
! ---------------------------------------------------------------------
!     STEP 6 - CALCULATION OF EQ CAUCHY STRESS & ELASTICITY TENSOR
! ---------------------------------------------------------------------
!     CALCULATE LEFT CAUCHY-GREEN DEFORMATION TENSOR
      BTOT = MATMUL(DFGRD1, TRANSPOSE(DFGRD1))
!
!     GET EIGEN VALUES AND EIGEN VECTORS OF BTOT
      DO I = 1,3
        DO J = 1,3
          BTOT_(I,J) = BTOT(I,J)
        END DO
      END DO     
      CALL DSYEVJ3(BTOT_, PDBTOT, PVBTOT)
      PDBTOT = TRANSPOSE(PDBTOT)
!
!     CALCULATE EIGEN VALUES OF DEVIATORIC LEFT CAUCHY-GREEN
!     DEFORMATION TENSOR
      DO J=1,3
        PVBBAR(J) = DET**(-TWO/THREE)*PVBTOT(J)
      END DO
!     
!     CALCULATING PRINCIPAL VALUES OF KIRCHOFF STRESS
      PVTAUEQ(1) = MU * ((TWO/THREE)*(PVBBAR(1)**(ALPHA/TWO))
     1                 -  (ONE/THREE)*(PVBBAR(2)**(ALPHA/TWO))
     2                 -  (ONE/THREE)*(PVBBAR(3)**(ALPHA/TWO)))
     3                 +   KELAS/TWO*(DET*DET - ONE)
      PVTAUEQ(2) = MU * ((TWO/THREE)*(PVBBAR(2)**(ALPHA/TWO))
     1                 -  (ONE/THREE)*(PVBBAR(3)**(ALPHA/TWO))
     2                 -  (ONE/THREE)*(PVBBAR(1)**(ALPHA/TWO)))
     3                 +   KELAS/TWO*(DET*DET - ONE)
      PVTAUEQ(3) = MU * ((TWO/THREE)*(PVBBAR(3)**(ALPHA/TWO))
     1                 -  (ONE/THREE)*(PVBBAR(1)**(ALPHA/TWO))
     2                 -  (ONE/THREE)*(PVBBAR(2)**(ALPHA/TWO)))
     3                 +   KELAS/TWO*(DET*DET - ONE)
!
!     CALCULATING EQ CAUCHY STRESS TENSOR
      DO J=1,3
        DO K=1,3
          SIGMAEQ(K,J)=(PVTAUEQ(1)*PDBTOT(1,K)*PDBTOT(1,J)
     1                 +PVTAUEQ(2)*PDBTOT(2,K)*PDBTOT(2,J)
     2                 +PVTAUEQ(3)*PDBTOT(3,K)*PDBTOT(3,J))/DET
        END DO
      END DO
!
!     CALCULATING CAB
      CAB(1,1) = MU * ALPHA 
     1          *((FOUR/NINE)*(PVBBAR(1)**(ALPHA/TWO))
     2            +(ONE/NINE)*(PVBBAR(2)**(ALPHA/TWO)
     3                        +PVBBAR(3)**(ALPHA/TWO)))
     4          + KELAS*DET*DET 
      CAB(2,2) = MU * ALPHA 
     1          *((FOUR/NINE)*(PVBBAR(2)**(ALPHA/TWO))
     2            +(ONE/NINE)*(PVBBAR(3)**(ALPHA/TWO)
     3                        +PVBBAR(1)**(ALPHA/TWO)))
     4          + KELAS*DET*DET  
      CAB(3,3) = MU * ALPHA 
     1          *((FOUR/NINE)*(PVBBAR(3)**(ALPHA/TWO))
     2            +(ONE/NINE)*(PVBBAR(1)**(ALPHA/TWO)
     3                        +PVBBAR(2)**(ALPHA/TWO)))
     4          + KELAS*DET*DET 
      CAB(1,2) = MU * ALPHA 
     1          *((-TWO/NINE)*(PVBBAR(1)**(ALPHA/TWO)
     2                        +PVBBAR(2)**(ALPHA/TWO))
     3            +(ONE/NINE)*(PVBBAR(3)**(ALPHA/TWO)))
     4          + KELAS*DET*DET 
      CAB(1,3) = MU * ALPHA 
     1          *((-TWO/NINE)*(PVBBAR(1)**(ALPHA/TWO)
     2                        +PVBBAR(3)**(ALPHA/TWO))
     3            +(ONE/NINE)*(PVBBAR(2)**(ALPHA/TWO)))
     4          + KELAS*DET*DET 
      CAB(2,3) = MU * ALPHA 
     1          *((-TWO/NINE)*(PVBBAR(2)**(ALPHA/TWO)
     2                        +PVBBAR(3)**(ALPHA/TWO))
     3            +(ONE/NINE)*(PVBBAR(1)**(ALPHA/TWO)))
     4          + KELAS*DET*DET 
      CAB(2,1) = CAB(1,2)
      CAB(3,1) = CAB(1,3)
      CAB(3,2) = CAB(2,3)
!
!     CALCULATING GAB 
      DO J=1,3
        DO K=1,3
          IF(J.EQ.K) THEN
            GAB(J,K) = MU * ALPHA * PVBBAR(J)**(ALPHA/TWO)
          ELSE
            IF (ABS(PVBTOT(J)-PVBTOT(K)).LT.EPS) THEN
              GAB(J,K) = MU * ALPHA 
     1                  *(TWO/THREE*PVBBAR(J)**(ALPHA/TWO)
     2                   +ONE/THREE*PVBBAR(K)**(ALPHA/TWO))
            ELSE
              GAB(J,K) = (PVTAUEQ(J)*PVBBAR(K) - PVTAUEQ(K)*PVBBAR(J))
     1                 / (PVBBAR(J) - PVBBAR(K))     
            END IF
          END IF
        END DO
      END DO
!
!     CALCULATE MATERIAL TENSOR COEFFICIENTS (C4EQMAT)
!     EXPRESSED IN EIGEN BASIS 
      C4EQMAT = ZERO
!     1. CALCULATE COEFFICIENTS FOR AABB
      DO J=1,3
        DO K=1,3
          C4EQMAT(J,J,K,K) = (CAB(J,K) - TWO*PVTAUEQ(J)*IDT2(J,K))
     1                       /DET
        END DO
      END DO 
!     2. CALCULATE COEFFICIENTS FOR ABAB / ABBA
      DO J=1,3
        DO K=1,3
          IF (J.EQ.K) THEN
            CYCLE
          END IF
          C4EQMAT(J,K,J,K) = GAB(J,K) / (TWO*DET)
          C4EQMAT(J,K,K,J) = C4EQMAT(J,K,J,K)
        END DO
      END DO
!     
!     CALCULATE MATERIAL TENSOR (C4EQ) IN CURRENT CONFIGURATION
!     EXPRESSED IN STANDARD BASIS
      C4EQ = ZERO
      DO J=1,3
        DO K=1,3
          DO L=1,3
            DO M=1,3
              C4EQ(J,K,L,M) = 
     1          C4EQMAT(1,1,1,1)*PDBTOT(1,J)*PDBTOT(1,K)
     1                          *PDBTOT(1,L)*PDBTOT(1,M)
     1         +C4EQMAT(1,1,2,2)*PDBTOT(1,J)*PDBTOT(1,K)
     1                          *PDBTOT(2,L)*PDBTOT(2,M)
     1         +C4EQMAT(1,1,3,3)*PDBTOT(1,J)*PDBTOT(1,K)
     1                          *PDBTOT(3,L)*PDBTOT(3,M)
     1         +C4EQMAT(2,2,1,1)*PDBTOT(2,J)*PDBTOT(2,K)
     1                          *PDBTOT(1,L)*PDBTOT(1,M)
     1         +C4EQMAT(2,2,2,2)*PDBTOT(2,J)*PDBTOT(2,K)
     1                          *PDBTOT(2,L)*PDBTOT(2,M)
     1         +C4EQMAT(2,2,3,3)*PDBTOT(2,J)*PDBTOT(2,K)
     1                          *PDBTOT(3,L)*PDBTOT(3,M)
     1         +C4EQMAT(3,3,1,1)*PDBTOT(3,J)*PDBTOT(3,K)
     1                          *PDBTOT(1,L)*PDBTOT(1,M)
     1         +C4EQMAT(3,3,2,2)*PDBTOT(3,J)*PDBTOT(3,K)
     1                          *PDBTOT(2,L)*PDBTOT(2,M)
     1         +C4EQMAT(3,3,3,3)*PDBTOT(3,J)*PDBTOT(3,K)
     1                          *PDBTOT(3,L)*PDBTOT(3,M)
     1         +C4EQMAT(1,2,1,2)*PDBTOT(1,J)*PDBTOT(2,K)
     1                          *PDBTOT(1,L)*PDBTOT(2,M)
     1         +C4EQMAT(1,2,2,1)*PDBTOT(1,J)*PDBTOT(2,K)
     1                          *PDBTOT(2,L)*PDBTOT(1,M)
     1         +C4EQMAT(1,3,1,3)*PDBTOT(1,J)*PDBTOT(3,K)
     1                          *PDBTOT(1,L)*PDBTOT(3,M)
     1         +C4EQMAT(1,3,3,1)*PDBTOT(1,J)*PDBTOT(3,K)
     1                          *PDBTOT(3,L)*PDBTOT(1,M)
     1         +C4EQMAT(2,1,2,1)*PDBTOT(2,J)*PDBTOT(1,K)
     1                          *PDBTOT(2,L)*PDBTOT(1,M)
     1         +C4EQMAT(2,1,1,2)*PDBTOT(2,J)*PDBTOT(1,K)
     1                          *PDBTOT(1,L)*PDBTOT(2,M)
     1         +C4EQMAT(2,3,2,3)*PDBTOT(2,J)*PDBTOT(3,K)
     1                          *PDBTOT(2,L)*PDBTOT(3,M)
     1         +C4EQMAT(2,3,3,2)*PDBTOT(2,J)*PDBTOT(3,K)
     1                          *PDBTOT(3,L)*PDBTOT(2,M)
     1         +C4EQMAT(3,1,3,1)*PDBTOT(3,J)*PDBTOT(1,K)
     1                          *PDBTOT(3,L)*PDBTOT(1,M)
     1         +C4EQMAT(3,1,1,3)*PDBTOT(3,J)*PDBTOT(1,K)
     1                          *PDBTOT(1,L)*PDBTOT(3,M)
     1         +C4EQMAT(3,2,3,2)*PDBTOT(3,J)*PDBTOT(2,K)
     1                          *PDBTOT(3,L)*PDBTOT(2,M)
     1         +C4EQMAT(3,2,2,3)*PDBTOT(3,J)*PDBTOT(2,K)
     1                          *PDBTOT(2,L)*PDBTOT(3,M)
            END DO
          END DO
        END DO
      END DO
!
!     CALCLUATE JAUMANN RATE FORM OF MATERIAL TENSOR
      C4EQJ = ZERO
      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO L=1,3
              C4EQJ(I,J,K,L) = C4EQ(I,J,K,L)/DET
     1                        +(IDT2(I,K)*SIGMAEQ(J,L)
     1                         +IDT2(I,L)*SIGMAEQ(J,K)
     1                         +IDT2(J,K)*SIGMAEQ(I,L)
     1                         +IDT2(J,L)*SIGMAEQ(I,K))/(TWO)        
            END DO
          END DO
        END DO
      END DO
!
!     GENERATE ABAQUS EQ TANGENT STIFFNESS MATRIX (VOIGT NOTATION)
      CEQ = ZERO
      CEQ(1,1) = C4EQJ(1,1,1,1)
      CEQ(1,2) = C4EQJ(1,1,2,2)
      CEQ(1,3) = C4EQJ(1,1,3,3)
      CEQ(1,4) = C4EQJ(1,1,1,2)
      CEQ(1,5) = C4EQJ(1,1,1,3)
      CEQ(1,6) = C4EQJ(1,1,2,3)
      CEQ(2,1) = CEQ(1,2)
      CEQ(2,2) = C4EQJ(2,2,2,2)
      CEQ(2,3) = C4EQJ(2,2,3,3)
      CEQ(2,4) = C4EQJ(2,2,1,2)
      CEQ(2,5) = C4EQJ(2,2,1,3)
      CEQ(2,6) = C4EQJ(2,2,2,3)
      CEQ(3,1) = CEQ(1,3)
      CEQ(3,2) = CEQ(2,3)
      CEQ(3,3) = C4EQJ(3,3,3,3)
      CEQ(3,4) = C4EQJ(3,3,1,2)
      CEQ(3,5) = C4EQJ(3,3,1,3)
      CEQ(3,6) = C4EQJ(3,3,2,3)
      CEQ(4,1) = CEQ(1,4)
      CEQ(4,2) = CEQ(2,4)
      CEQ(4,3) = CEQ(3,4)
      CEQ(4,4) = C4EQJ(1,2,1,2)
      CEQ(4,5) = C4EQJ(1,2,1,3)
      CEQ(4,6) = C4EQJ(1,2,2,3)
      CEQ(5,1) = CEQ(1,5)
      CEQ(5,2) = CEQ(2,5)
      CEQ(5,3) = CEQ(3,5)
      CEQ(5,4) = CEQ(4,5)
      CEQ(5,5) = C4EQJ(1,3,1,3)
      CEQ(5,6) = C4EQJ(1,3,2,3)
      CEQ(6,1) = CEQ(1,6)
      CEQ(6,2) = CEQ(2,6)
      CEQ(6,3) = CEQ(3,6)
      CEQ(6,4) = CEQ(4,6)
      CEQ(6,5) = CEQ(5,6)
      CEQ(6,6) = C4EQJ(2,3,2,3)
!
! ---------------------------------------------------------------------
!     STEP 7 - CALCULATION OF TOTAL CAUCHY STRESS & ELASTICITY TENSOR
! ---------------------------------------------------------------------
!     CALCULATE TOTAL CAUCHY STRESS TENSOR
      DO J=1,3
        DO K=1,3
          STRESSTOT(J,K) = SIGMANEQ(J,K) + SIGMAEQ(J,K)
        END DO
      END DO

!     VOIGT NOTATION FOR CAUCHY STRESS
      STRESS(1) = STRESSTOT(1,1)
      STRESS(2) = STRESSTOT(2,2)
      STRESS(3) = STRESSTOT(3,3)
      STRESS(4) = STRESSTOT(1,2)
      STRESS(5) = STRESSTOT(1,3)
      STRESS(6) = STRESSTOT(2,3)
!
!     CALCULATE TOTAL ABAQUS TANGENT STIFNESS MATRIX
      DO J=1,6
        DO K=1,6
          DDSDDE(J,K) = CNEQ(J,K) + CEQ(J,K)
        END DO
      END DO

      RETURN
      END SUBROUTINE UMAT

! ---------------------------------------------------------------------
!     UTILITY SUBROUTINES
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
!     M33INV  -  Compute the Inverse of a 3x3 Matrix.
!     SOURCE:  David G. Simpson - NASA Goddard Space Flight Center
!
!     A       - Input  - 3x3 Matrix to be Inverted
!     AINV    - Output - 3x3 Inverse of Matrix A
!     OK_FLAG - .TRUE. If A is non-singular and A is Inverted
! ---------------------------------------------------------------------
      SUBROUTINE M33INV (A, AINV, OK_FLAG)
            IMPLICIT NONE
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
            LOGICAL, INTENT(OUT) :: OK_FLAG

            DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-16
            DOUBLE PRECISION :: DET
            DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR
!
!
            DET = A(1,1) * A(2,2) * A(3,3)  
     1          - A(1,1) * A(2,3) * A(3,2)  
     2          - A(1,2) * A(2,1) * A(3,3)  
     3          + A(1,2) * A(2,3) * A(3,1)  
     4          + A(1,3) * A(2,1) * A(3,2)  
     5          - A(1,3) * A(2,2) * A(3,1)
!
            IF (ABS(DET) .LE. EPS) THEN
            AINV = 0.0D0
            OK_FLAG = .FALSE.
            RETURN
            END IF
!
            COFACTOR(1,1) = +(A(2,2) * A(3,3) - A(2,3) * A(3,2))
            COFACTOR(1,2) = -(A(2,1) * A(3,3) - A(2,3) * A(3,1))
            COFACTOR(1,3) = +(A(2,1) * A(3,2) - A(2,2) * A(3,1))
            COFACTOR(2,1) = -(A(1,2) * A(3,3) - A(1,3) * A(3,2))
            COFACTOR(2,2) = +(A(1,1) * A(3,3) - A(1,3) * A(3,1))
            COFACTOR(2,3) = -(A(1,1) * A(3,2) - A(1,2) * A(3,1))
            COFACTOR(3,1) = +(A(1,2) * A(2,3) - A(1,3) * A(2,2))
            COFACTOR(3,2) = -(A(1,1) * A(2,3) - A(1,3) * A(2,1))
            COFACTOR(3,3) = +(A(1,1) * A(2,2) - A(1,2) * A(2,1))
!
            AINV = TRANSPOSE(COFACTOR) / DET
!
            OK_FLAG = .TRUE.
            RETURN
!
      END SUBROUTINE M33INV

* ----------------------------------------------------------------------------
      SUBROUTINE DSYEVJ3(A, Q, W)
* ----------------------------------------------------------------------------
* Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
* matrix A using the Jacobi algorithm.
* The upper triangular part of A is destroyed during the calculation,
* the diagonal elements are read but not destroyed, and the lower
* triangular elements are not referenced at all.
* ----------------------------------------------------------------------------
* Parameters:
*   A: The symmetric input matrix
*   Q: Storage buffer for eigenvectors
*   W: Storage buffer for eigenvalues
* ----------------------------------------------------------------------------
*     .. Arguments ..
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION W(3)

*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 3 )
    
*     .. Local Variables ..
      DOUBLE PRECISION SD, SO
      DOUBLE PRECISION S, C, T
      DOUBLE PRECISION G, H, Z, THETA
      DOUBLE PRECISION THRESH
      INTEGER          I, X, Y, R

*     Initialize Q to the identitity matrix
*     --- This loop can be omitted if only the eigenvalues are desired ---
      DO 10 X = 1, N
        Q(X,X) = 1.0D0
        DO 11, Y = 1, X-1
          Q(X, Y) = 0.0D0
          Q(Y, X) = 0.0D0
   11   CONTINUE
   10 CONTINUE

*     Initialize W to diag(A)
      DO 20 X = 1, N
        W(X) = A(X, X)
   20 CONTINUE

*     Calculate SQR(tr(A))  
      SD = 0.0D0
      DO 30 X = 1, N
        SD = SD + ABS(W(X))
   30 CONTINUE
      SD = SD**2
 
*     Main iteration loop
      DO 40 I = 1, 50
*       Test for convergence
        SO = 0.0D0
        DO 50 X = 1, N
          DO 51 Y = X+1, N
            SO = SO + ABS(A(X, Y))
   51     CONTINUE
   50   CONTINUE
        IF (SO .EQ. 0.0D0) THEN
          RETURN
        END IF

        IF (I .LT. 4) THEN
          THRESH = 0.2D0 * SO / N**2
        ELSE
          THRESH = 0.0D0
        END IF

*       Do sweep
        DO 60 X = 1, N
          DO 61 Y = X+1, N
            G = 100.0D0 * ( ABS(A(X, Y)) )
            IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X))
     $                    .AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
              A(X, Y) = 0.0D0
            ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
*             Calculate Jacobi transformation
              H = W(Y) - W(X)
              IF ( ABS(H) + G .EQ. ABS(H) ) THEN
                T = A(X, Y) / H
              ELSE
                THETA = 0.5D0 * H / A(X, Y)
                IF (THETA .LT. 0.0D0) THEN
                  T = -1.0D0 / (SQRT(1.0D0 + THETA**2) - THETA)
                ELSE
                  T = 1.0D0 / (SQRT(1.0D0 + THETA**2) + THETA)
                END IF
              END IF

              C = 1.0D0 / SQRT( 1.0D0 + T**2 )
              S = T * C
              Z = T * A(X, Y)
              
*             Apply Jacobi transformation
              A(X, Y) = 0.0D0
              W(X)    = W(X) - Z
              W(Y)    = W(Y) + Z
              DO 70 R = 1, X-1
                T       = A(R, X)
                A(R, X) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   70         CONTINUE
              DO 80, R = X+1, Y-1
                T       = A(X, R)
                A(X, R) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   80         CONTINUE
              DO 90, R = Y+1, N
                T       = A(X, R)
                A(X, R) = C * T - S * A(Y, R)
                A(Y, R) = S * T + C * A(Y, R)
   90         CONTINUE

*             Update eigenvectors
*             --- This loop can be omitted if only the eigenvalues are desired ---
              DO 100, R = 1, N
                T       = Q(R, X)
                Q(R, X) = C * T - S * Q(R, Y)
                Q(R, Y) = S * T + C * Q(R, Y)
  100         CONTINUE
            END IF
   61     CONTINUE
   60   CONTINUE
   40 CONTINUE

      PRINT *, "DSYEVJ3: No convergence."
            
      END SUBROUTINE
* End of subroutine DSYEVJ3

