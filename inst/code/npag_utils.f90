module npag_utils
!
! npag_utils.f90 ! put this file _first_ in compilation command; the following must be created at first step of compilation
! -> npag_utils.mod ! zipped fortran interfaces
! -> npag_utils.o ! .o functions and subroutines
implicit none
! by default, this module should appear empty ...
private
! ... to access, put "USE npag_utils,only: list_of_what_you_want" at
! the head of the program unit (e.g. immediately after "subroutine name"
! and prior to "implicit ___" statements).  Where, "public" variables
! are what you can put in "list_of_what_you_want"
public SimConc, useranal, pause, shift, makevec, check_input_array_size, &
  cp_theta_to_rpar, cp_lrcs_to_rpar, expand_grid, practice_c_call, &
  verifyval, orderdelta,thesame, predlast3, do140, obs_is_BLQ, obs_is_missing, &
  print_matrix, i_cycle, i_do, i_jsub, i_ig, max_ODE_comps, max_ODE_params, &
  max_pop_rand_varbs, max_pop_varbs, max_pop_params, max_SS_doses, &
  max_covs, maxnumeq, max_meas_per_eqn, max_m_per_obs, max_obs, max_obs_dim, &
  max_doses, max_input_dim, max_RS_J, k_ig, k_jsub, k_sum_norm_err_sq, &
  i_dvode_reserved, k_dvode_reserved, k_sum_z_sq, k_prod_pr, k_resolve, &
  k_p_start, k_p_end, k_r_start, k_r_end, k_gamma, k_flat, k_sfac, k_ofac, &
  k_c0_base, k_c1_base, k_c2_base, k_c3_base, k_c4_base, k_c5_base, i_errmod, &
  i_is_log10, i_Nnormalobs, i_is_poisson_obs, i_Npoissonobs, i_misval, i_skip_ig, &
  comp_pred_amt, k_rtol, k_atol, i_mf, i_ndim, i_nparams, i_NBLQ
! ODE parameters
integer, parameter :: max_ODE_comps =      20	! check RWORK(1:LRW); LRW = 1002
integer, parameter :: k_dvode_reserved =   23   ! RPAR(1:23) are reserved by DVODE
integer, parameter :: i_dvode_reserved =   23   ! IPAR(1:23) are reserved by DVODE
integer, parameter :: max_ODE_params =     47   ! 32  ! NVAR+NRANFIX+NOFIX <= max_ODE_params , of types:
integer, parameter :: max_pop_rand_varbs = 30   ! 30  ! NVAR <= max_pop_rand_varbs; variables w/#support <= max_obs
integer, parameter :: max_pop_varbs =      20   ! NRANFIX <= max_pop_varbs; variables with #support = 1
integer, parameter :: max_pop_params =     20   ! NOFIX <= max_pop_params; #user-defined scalar pop parameters
integer, parameter :: max_SS_doses =      100   ! dim of SS calc = max_SS_doses X max_ODE_comps
! INPUT / OUTPUT parameters
integer, parameter :: max_input_dim =       7   ! #Input types, NDRUG <= max_input_dim
integer, parameter :: max_doses =        5000   ! ND <= max_doses per record
integer, parameter :: max_obs  =          800   ! = max #records; MAXSUB = 800
integer, parameter :: max_covs =           20   ! max #covariate measures per record; dim
integer, parameter :: max_obs_dim  =      150   ! = MAXOBDIM
integer, parameter :: maxnumeq =            7   ! #output equations ~ max_output_dim
integer, parameter :: max_meas_per_eqn =   99   ! #measurements per patient
integer, parameter :: obs_is_BLQ =        -88   ! observation is marked BLQ in datafile
integer, parameter :: obs_is_missing =    -99   ! observation is marked missing in datafile
! Parameters that are determined from above:
integer, parameter :: max_m_per_obs = max_meas_per_eqn * maxnumeq ! was hardcoded to 594
integer, parameter :: max_RS_J = 2*max_input_dim + max_covs           !  NI <= max_RS_J
! X(NBCOMP(I))=X(NBCOMP(I))+BS(KNS,I)*FA(I)
! RPAR and IPAR are concatenations of ODE parameters, following are
! shortcuts into RPAR (k) and ILIST (j) and IPAR (i)
integer, parameter :: k_p_start = k_dvode_reserved +  1
integer, parameter :: k_p_end = k_dvode_reserved + max_ODE_params
integer, parameter :: k_r_start = k_p_end + 1
integer, parameter :: k_r_end = k_p_end + max_RS_J
integer, parameter :: k_jsub = k_r_end + 1
integer, parameter :: k_ig = k_jsub + 1
integer, parameter :: k_gamma = k_ig + 2 ! this should be 96
integer, parameter :: k_flat = k_gamma + 1
integer, parameter :: k_c_start = k_flat + 1
integer, parameter :: k_c0_base = k_flat
integer, parameter :: k_c1_base = k_c0_base + maxnumeq
integer, parameter :: k_c2_base = k_c1_base + maxnumeq
integer, parameter :: k_c3_base = k_c2_base + maxnumeq
integer, parameter :: k_c4_base = k_c3_base + maxnumeq
integer, parameter :: k_c5_base = k_c4_base + maxnumeq
integer, parameter :: k_c_end = k_flat + maxnumeq * 6       ! maxnumeq * (C0,C1,C2,C3,C4,C5)
integer, parameter :: k_sfac = k_c_end + 1         ! SIGFAC
integer, parameter :: k_ofac = k_sfac + 1          ! OFAC
integer, parameter :: k_sum_z_sq = k_ofac + 1      ! = \sum ((obs - pred)/sigma)^2 -- for obs ~ Normal
integer, parameter :: k_prod_pr = k_sum_z_sq + 1   ! = \Prod pr(obs) -- for obs ~ Poisson
integer, parameter :: k_resolve = k_prod_pr + 1
integer, parameter :: k_rtol = k_resolve + 1
integer, parameter :: k_atol = k_rtol + 1
integer, parameter :: k_sum_norm_err_sq = k_atol + 1  ! = \sum ((obs - pred)/pred)^2 -- for all obs
! IPAR(110 + J) in {229, ~229} = flag for Poisson
! (note: set to 229 if C0 = C2 = C3 = 229 in model.txt)
! IPAR(120 + J) in {10, ~10} obs = log10(measurement)
! IPAR(23 + I) = INTLIST(I; 1:10); where,
! Values copied into INTLIST are ->
! (1) = AGE; (2) = ISEX; (3) = HEIGHT; (4) = IETHFLAG;
! (5) = /CNST2/NDRUG; (6) = /CNST2/NADD;
! (7) = /CNST/NI = 2*NDRUG+NADD; (8) = /CNST/ND;
! (9) = /CNST2/NUMEQT; (10) = /CNST2/M = NOBSER;
! Note: Some above are NOT integers!!!
integer, parameter :: i_jsub = i_dvode_reserved + 11 ! yes, 11, not 1
integer, parameter :: i_ig = i_jsub + 1
integer, parameter :: i_do = 36                ! IPAR(i_do) = DO# in 1000, 6000, 7000
integer, parameter :: i_cycle = 37
integer, parameter :: i_errmod = 38
integer, parameter :: i_misval = 39
integer, parameter :: i_Nnormalobs = 40
integer, parameter :: i_Npoissonobs = 41
integer, parameter :: i_mf = 42                ! ODE method flag
integer, parameter :: i_NBLQ = 43
integer, parameter :: i_ndim = 130             ! number of compartments in regression model
integer, parameter :: i_nparams = 131          ! total #parameters
integer, parameter :: i_skip_ig = 100          ! if ipar(i_skip_ig) = 0 set P(JSUB|IG)=0
integer, parameter :: i_is_poisson_obs = 110   ! 110 + eqno = 229 if Poisson
integer, parameter :: i_is_log10 = 120         ! 120 + eqno = 10, then obs recorded as log10(measurement)
double precision, dimension(max_m_per_obs,max_ODE_comps) :: comp_pred_amt ! ijob==4
! MODULE random is catted into MODULE npag_utils -----------------------
! IMPLICIT NONE
REAL, PRIVATE      :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0,   &
                       vsmall = TINY(1.0), vlarge = HUGE(1.0)
PRIVATE            :: integral
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)
!
! CONTAINS
! MODULE random is catted into MODULE npag_utils -----------------------
! Pmetrics erases all lines beginning with 'c' or "C", so
! the following line must be " contains", w/leading white-space
 contains
!
! subroutine SimConc
! subroutine useranal
! subroutine print_matrix
! subroutine pause
! subroutine do140
! subroutine cp_lrcs_to_rpar
! function check_input_array_size
! subroutine makevec
! subroutine practice_c_call
!
! #################################################################### !
! #################################################################### !
!
! SimConc is a wrapper around USERANAL and ANAL3.
! Based on (deprecated) NPAG:idm2*.f:FUNC2()
!
!  NOTE THAT THE DIMENSIONS RELATED TO THE NO. OF OUTPUT EQS. IN
!  YO, YT AND Y ARE CHANGED TO MAXNUMEQ (FROM 6). NUMEQT COULD NOT
!  BE USED BECAUSE THESE ARRAYS WERE NOT PASSED TO THIS ROUTINE AS
!  DUMMY ARGUMENTS.
!  THE 2ND DIMENSION OF YPRED IS CHANGED TO NUMEQT, SINCE IT IS PASSED
!  IN THE ARGUMENT LIST, AND CAN THEREFORE BE VARIABLY DIMENSIONED BY
!  NUMEQT.
!  NOTE THAT "7" IN THE ABOVE ARRAYS INDICATE THE NO. OF DRUGS ALLOWED.
! Error: $ operator is invalid for atomic vectors
!  above error prevents output report but PMreport(wd = ".") works.
	SUBROUTINE SimConc(ijob,JSUB,IG                                 &
         , M, YPRED, NUMEQT                                             &
         , NBCOMP                                                       &
         , NPP,ESTML,NDIM,MF,RTOL,ATOL                                  &
         , TIM,SIG,RS,BS                                                &
         , INTLIST,IPAR,RPAR,ERRFIL)
!  SimConc FINDS YPRED(I,J) = OUTPUT CONC. AT TIME I, I=1,M, and
!  OUTPUT Eqn. J, J=1,NOS, GIVEN PARAMETER VALUES IN P.
!       use npag_utils, only: useranal, pause                            &
!     &   , verifyval, shift, thesame, predlast3                         &
!     &   , maxnumeq, max_m_per_obs, max_SS_doses                        &
!     &   , max_ODE_params, max_doses, max_ODE_comps, max_RS_J           &
!     &   , max_input_dim, k_dvode_reserved, k_p_end, k_jsub, k_ig       &
!     &   , i_ig, i_jsub, i_dvode_reserved, comp_pred_amt
        IMPLICIT NONE
! INTLIST(1) = int(AGE)
! INTLIST(2) = ISEX
! INTLIST(3) = int(HEIGHT)
! INTLIST(4) = IETHFLAG
! INTLIST(5) = /CNST2/NDRUG
! INTLIST(6) = /CNST2/NADD = No. Additional Covariates
! INTLIST(7) = /CNST/NI = 2*NDRUG+NADD
! INTLIST(8) = /CNST/ND ; ND =  NDO ; #Dose Events
! INTLIST(9) = /CNST2/NOS (Called /CNST2/NUMEQTT in SUBROUTINE FILRED)
! INTLIST(10) = /CNST2/M = NOBSER
! ------ Argument List
      integer, intent(in) ::  ijob, JSUB, IG, M
      double precision, dimension(:,:), intent(inout) :: YPRED ! max_m_per_obs,NUMEQT
      integer, intent(in) ::  NUMEQT, NPP
      double precision, dimension(:), intent(in) :: ESTML ! max_ODE_params
      integer, intent(in) ::  NDIM,MF
      integer, dimension(:), intent(in) :: NBCOMP ! max_input_dim
      double precision, intent(in) :: RTOL
      double precision, dimension(:), intent(in)  :: ATOL ! max_ODE_comps
      double precision, dimension(:), intent(in) ::TIM     ! max_m_per_obs
      double precision, dimension(:), intent(inout) :: SIG  ! max_doses
      double precision, dimension(:,:), intent(inout) :: RS ! max_doses,max_RS_J
      double precision, dimension(:,:), intent(inout) :: BS ! max_doses,max_input_dim -- NOT USED
      integer, dimension(:),intent(in) :: INTLIST ! 128
      integer, dimension(:),intent(inout) :: IPAR ! 257
      double precision, dimension(:),intent(inout) :: RPAR ! 257
      CHARACTER ERRFIL*20
! ------ Local variables; remaining w/in SimConc
      integer JUMPTO45
      integer NI,NP
      integer NOS
      integer NTL, KNT, KNSNEW, KNS, J
      integer ISTEADY, ISKIPBOL, III
      integer IKNS, ID, I, KNSOLD
! Julian's test data
      double precision, dimension(max_ODE_comps) :: X_in
      double precision T_in
      integer TestData ! Set > 3, BELOW, to skip test
! ------ Local variables; passed to subroutines outside the !$omp do
      integer KNTM1,NN,NSET,ICONV,ND,NDRUG,NADD,N,ISAME,NDO,JJJ,KKK
      double precision T,TOUT
      double precision DOSEINT ! -sig(kns)
      double precision, dimension(max_RS_J) :: R
      double precision, dimension(max_ODE_comps) :: B
      double precision, dimension(max_ODE_comps) :: X
      double precision, dimension(max_SS_doses,max_ODE_comps) :: XSTORE
      double precision, dimension(max_ODE_comps) :: XPRED
      double precision, dimension(max_input_dim) :: TLAG, FA
      double precision, dimension(MAXNUMEQ) :: YT
      double precision, dimension(100) :: XVERIFY
!      double precision, dimension(max_m_per_obs) :: TIMO               ! TIM(KNT) Observation times
      double precision, dimension(max_doses) :: SIGO                   ! SIG(KNS) Stimulus/Dose times
      double precision, dimension(max_doses,max_RS_J) :: RSO           ! Rate info
      double precision, dimension(max_doses,max_input_dim) :: BSO      ! Bolus info
! (ijob == 3) AUC calculation in NPAG
!  M <- NUMT(JSUB)
!  YPRED <- YYPRED(71281,NUMEQT) dimensioned in main.
!  TIM <- TPRED(1:NUMT) i.e. 1:M
!
! replace above w/allocatables
!      double precision, allocatable :: R(:), B(:), X(:), XPRED(:),TLAG(:), FA(:), YT(:), XVERIFY(:)
!      double precision, allocatable ::  XSTORE(:,:)
!      double precision, allocatable :: TIMO(:), SIGO(:), RSO(:,:), BSO(:,:)
! $omp ThreadPrivate(RSO,BSO,SIGO) ! ,TIMO)
! $omp ThreadPrivate(KNTM1,NN,NSET,ICONV,ND,NDRUG,NADD,N,ISAME,NDO)
! $omp ThreadPrivate(TLAG,X,FA,R,B,YT,XSTORE,XPRED,T,TOUT,DOSEINT)
! $omp ThreadPrivate(XVERIFY)
!      allocate(R(max_RS_J), B(max_ODE_comps), X(max_ODE_comps), XSTORE(max_SS_doses,max_ODE_comps) &
!        , XPRED(max_ODE_comps), TLAG(max_input_dim), FA(max_input_dim), YT(MAXNUMEQ), XVERIFY(100))
!      allocate(TIMO(max_m_per_obs), SIGO(max_doses), RSO(max_doses,max_RS_J), BSO(max_doses,max_input_dim))
!       ijob \in 1,2,3: NPAG
!       if (IPAR(i_do)==1000) then ! NPAG
!         ijob = 1: FUNC
!       end if
!       if (IPAR(i_do)==6000) then ! NPAG post processing; Pop stats
!         ijob = 2: FUNC2 -- mean, median, mode
!         ijob = 3: FUNC3 -- AUC
!       end if
!       if (IPAR(i_do)==7000) then ! NPAG post processing; Bayes stats
!         ijob = 2: FUNC2
!         ijob = 3: FUNC3
!       end if
!       ijob = 4: SIMULATOR; FUNC2
!       ijob \in 5,6,7: IT2B
!       write (*,*) "In SimConc w/ijob=",ijob
! Initialize saved local variables to 0.D0 (note: modules save all local variables)
       YPRED = 0.D0  ! This is the return variable (an argument; not a local object)
       SIGO = 0.D0
       RSO = 0.d0
       BSO = 0.d0
!
       KNT = 0
       IKNS = 0
        KNSOLD = 0
        KNS = 0
! ----- Use INTLIST to initialize counters, if necessary
      JUMPTO45 = 0
      N = NDIM
      ND = intlist(8) ! ND = NDO, remove ND from code.
      NI = intlist(7)
      NP = NPP
      NDRUG = intlist(5)
      NADD = intlist(6)
      NOS = intlist(9)
! *****ODE CONSTANTS AND INITIALIZATION*****
      do III=1,max_ODE_comps
        X(III)=0.D0
        B(III)=0.D0 ! B(I) contains bolus info relevant to ANAL3, only.
      end do
      KNT=1
      KNS=1
!  NOTE THAT KNT IS THE RUNNING INDEX OF THE NEXT OBSERVATION TIME,
!  AND       KNS IS THE RUNNING INDEX OF THE NEXT DOSAGE TIME.
      T=0.0D0
      TOUT=0.D0
!  INITIALIZE ISKIPBOL = 0. SEE CODE BELOW. IT IS ONLY NEEDED FOR A
!  STEADY STATE DOSE SET WHICH HAS BOLUS DOSES.
      ISKIPBOL = 0
      DO I = 1,NDRUG
       R(2*I-1) = 0.D0   ! At T=0, IV rates = 0
      END DO
        DO I=1,ND
          DO J=1,NDRUG
            BSO(I,J)=0.D0 ! RS(I,2*J)
          END DO
        END DO
      if ( (ijob == 1) .or. (ijob == 5) )then ! NPAG func(ijob=1,...)
        DO I=1,ND ! DO I=1,INTLIST(8)
          DO J=1,NDRUG ! DO J = 1,INTLIST(5)
            BSO(I,J)=RS(I,2*J)
          END DO
        END DO
      end if
!  AS OF idm2x7.f, instead of R(1) = 0, the code has been changed to
!  set R(2*I-1) = 0, for I = 1,NDRUG. I.E., All IV rates for all NDRUG
!  drugs are initialized to be 0 ... in case the 1st obs. time is 0,
!  which means that OUTPUT is called before the R(I) are set below.
!  CALL SUBROUTINE GETFA IN npemdriv.f (THE FIRST TEMPLATE FILE TO
!  INCLUDE GETFA IS TSTMULTG.FOR) TO OBTAIN THE VALUE OF FA FOR EACH
!  OF THE NDRUG DRUGS.
!  AS OF idm2x12.f, BEFORE CALLING GETFA, MUST SET
!  THE R(.) IN CASE ANY OF THE FA(.) ARE FUNCTIONS OF THE
!  COVARIATES WHICH ARE ESTABLISHED FROM THE R(.) VALUES IN
!  GETFA.
      DO I=1,NI
       R(I)=RS(KNS,I) ! RS copied to RSO below
      END DO
! X, B are not initialized; and R, B, and INTLIST are not used; FA(I) <- ESTML(n)
      CALL GETFA(FA,X,ESTML,R,B,INTLIST)
!  NOTE THAT NBCOMP(I),I=1,NDRUG WAS SET IN SUBROUTINE SYMBOL AND
!  PASSED TO THIS ROUTINE VIA COMMON/BOLUSCOMP.
!  As of idm2x11.f, the code to save ND0, SIGO, RSO, is moved to before
!  the IF(N .EQ. 0) GO TO 75  statement. The reason is that before this
!  routine returns, ND, SIG, and RS are reset back to these values,
!  even if N = 0, and so they must be established at this time.
!  AS OF idm2x9.f, SAVE ND, SIG, AND RS WHETHER OR NOT NTL = 1, SINCE
!  IF THERE ARE STEADY STATE DOSE SETS, THE FIRST SIG(.) VALUE IN EACH
!  SET WILL BE CHANGED TO BE 0 BELOW.
      NDO = ND
      DO I=1,ND
        SIGO(I) = SIG(I)
        DO J=1,NI
          RSO(I,J) = RS(I,J)
        END DO
      END DO
!  IF N = 0, THE OUTPUT EQUATION(S) FOR Y ARE CODED EXPLICITLY INTO
!  SUBROUTINE OUTPUT, AND NO D.E. SOLUTIONS (VIA USERANAL/DIFFEQ) ARE
!  TO BE USED. IN THIS CASE, SKIP THE CODE REGARDING INITIAL CONDITIONS
!  OF THE COMPARTMENTS, SINCE THEY ARE IRRELEVANT (I.E., THE COMPARTMENT
!  AMOUNTS DON'T NEED TO BE INITIALIZED SINCE THEY WON'T BE UPDATED BY
!  INTEGRATING D.E.'S). IN FACT, COULD PROBABLY SKIP TIMELAGS TOO,
!  SINCE THEY CHANGE THE TIME THAT BOLUS DOSES ARE GIVEN, AND THIS
!  THEORETICALLY ONLY AFFECTS COMPARTMENT AMOUNTS (WHICH ARE NOT USED
!  IF N = 0), BUT JUST SKIP INITIAL CONDITIONS FOR NOW.
! IF(N .EQ. 0) GO TO 75
      IF(N /= 0) CALL GETIX(N,X,ESTML,R,B,INTLIST)
!
!  CALL SUBROUTINE GETIX IN npemdriv.f (THE FIRST TEMPLATE FILE TO
!  INCLUDE GETIX IS TSTMULTG.FOR) TO OBTAIN THE VALUE OF X (THE INITIAL
!  COMPARTMENT AMOUNT) FOR EACH OF THE N COMPARTMENTS.
!
!   75  Continue
!  CALL SUBROUTINE GETTLAG IN npemdriv.f (THE FIRST TEMPLATE FILE TO
!  INCLUDE GETTLAG IS TSTMULTG.FOR) TO OBTAIN THE VALUE OF THE TIMELAG
!  FOR EACH OF THE NDRUG DRUGS.
      CALL GETTLAG(TLAG,X,ESTML,R,B,INTLIST)
!  IF ANY TLAG(.) VALUES RETURN AS .NE. 0, THEN, CALL SUBROUTINE SHIFT
!  TO ADJUST THE DOSAGE REGIMEN APPROPRIATELY.
      NTL = 0
      DO ID = 1,NDRUG
        IF(TLAG(ID) /= 0) NTL = 1
      END DO
      IF(NTL == 1) THEN
!  STORE INCOMING VALUES IN ND, SIG, AND RS (WHICH CONTAINS BS VALUES)
!  SINCE THEY WILL BE CHANGED IN THE CALL TO SUBROUTINE SHIFT, WHICH
!  "SHIFTS" THE DOSAGE REGIMEN MATRIX TO ACCOUNT FOR THE TIMELAG
!  PARAMETER(S), TLAG(I). AT THE END OF THIS ROUTINE, THE VALUES IN ND,
!  SIG, AND RS WILL BE RESET TO THEIR INCOMING VALUES - TO BE READY FOR
!  THE NEXT CALL TO THIS ROUTINE WITH POSSIBLY DIFFERENT VALUES FOR
!  TLAG(I).
        CALL SHIFT(TLAG,ND,SIG,NDRUG,NADD,RS,INTLIST)
!  ESTABLISH THE VALUES IN EACH DRUG'S PO COLUMN TO THE CORRESPONDING
!  COLUMN IN ARRAY BS.
        DO I=1,ND
          DO J=1,NDRUG
            BS(I,J)=RS(I,2*J)  ! BSO replaced w/BS
          END DO
        END DO
      END IF ! (NTL == 1) CONDITION.
! *** KNT == KNS == 1; BEGIN INCREMENTING TIME
      IF(TIM(KNT) < SIG(KNS)) then ! IF(TIM(KNT).GE.SIG(KNS)) GO TO 12
        IF(TIM(KNT)== 0.0D0) then ! IF(TIM(KNT).NE.0.0D0) GO TO 45
!  THE ONLY WAY THE FOLLOWING CALL TO OUTPUT CAN OCCUR IS IF TIM(KNT)
!  = 0 --> OBTAIN YT = OUTPUT VALUE(S) AT TIME 0.0.
          CALL OUTPUT(0.D0,YT,X,RPAR,IPAR)
          DO 2000 I=1,NOS ! INTLIST(9)
            YPRED(KNT,I)=YT(I)
2000      end do
          KNT=KNT+1
        END IF
        JUMPTO45 = 1  ! GO TO 45
      ELSE
12      IF(TIM(KNT) == SIG(KNS)) then  ! IF(TIM(KNT).GT.SIG(KNS)) GO TO 13
          IF(TIM(KNT) == 0.0D0) then ! IF(TIM(KNT).NE.0.0D0) GO TO 45
!  THE ONLY WAY THE FOLLOWING CALL TO OUTPUT CAN OCCUR IS IF TIM(KNT)
!  = 0 --> OBTAIN YT = OUTPUT VALUE(S) AT TIME 0.0.
            CALL OUTPUT(0.D0,YT,X,RPAR,IPAR)
	    DO I=1,NOS
              YPRED(KNT,I)=YT(I)
2005        end do
            if (ijob == 4) then ! Simulator MONT114#7103,7115
              IF(N > 0) THEN
                DO III = 1,N
                  comp_pred_amt(KNT,III) = X(III)
                END DO
              ENDIF
              IF(N == -1) THEN
                DO III = 1,3
                  comp_pred_amt(KNT,III) = X(III)
                END DO
              ENDIF
           write (*,*) "In SimConc KNT and KNS=", KNT, KNS
            end if
            KNT=KNT+1
! Note: JUMPTO45 only if test on 13 is true
          else
            JUMPTO45 = 1
          end if
        END IF
13      IF(SIG(KNS) > 0.0D0) JUMPTO45 = 1
      END IF
!  CHECK TO SEE IF SIG(KNS) < 0. IF SO, IT MEANS THAT 100 STEADY STATE
!  DOSES SHOULD NOW BE APPLIED WITH AN INTERDOSE INTERVAL EQUAL TO
!  -SIG(KNS).
      ISTEADY = 0
      if (JUMPTO45 /= 1) then ! therefore, it is == 0
        IF(SIG(KNS) < 0.D0) THEN
          ISTEADY = 1
          NSET = 1
!  NOTE THAT ISTEADY = 1 TELLS THE PROGRAM BELOW TO PROCEED AS IF THE
!  DOSE TIME IS 0, AND START INTEGRATING THROUGH THE SET OF 100
!  DOSE SETS, ALL OF WHICH OCCUR BEFORE THE NEXT OBSERVATION TIME ...
!  BUT PAUSE AFTER THE END OF THE 5TH DOSE SET (NSET IS THE RUNNING NO.
!  OF THE CURRENT DOSE SETS THAT HAVE BEEN RUN) AND CALL SUBROUTINE
!  PREDLAST3 TO PREDICT THE STEADY STATE COMPARTMENT AMOUNTS AFTER THE
!  100 DOSE SETS (NOTE THAT THE COMPARTMENT AMOUNTS WILL HAVE TO BE
!  STORED AT THE END OF EACH OF THE STEADY STATE DOSE SETS AS THE LOGIC
!  OF PREDLAST3 REQUIRES).
!  IF "CONVERGENCE" IS ACHIEVED AT THAT POINT, ASSIGN THE COMPARTMENT
!  AMOUNTS TO BE THE PREDICTED AMOUNTS, AND ASSIGN KNS TO BE WHAT IT IS
!  WHEN THESE STEADY STATE DOSE SETS HAVE FINISHED. NOTE THAT THE END OF
!  THE 100TH DOSE SET WILL BE AT TIME 100*(-SIG(KNS)), SO KNS WILL BE
!  THE INDEX OF THE FIRST DOSE EVENT WHICH OCCURS AFTER THIS TIME.
!  IF "CONVERGENCE" IS NOT ACHIEVED, CONTINUE APPLYING THE LOGIC OF
!  PREDLAST3 UNTIL IT IS ACHIEVED, OR UNTIL THE 100 DOSE SETS ARE ALL
!  INTEGRATED THROUGH, WHICHEVER COMES FIRST.
          DOSEINT = -SIG(KNS)
!  RESET SIG(KNS) TO BE 0 SINCE THIS DOSE EVENT REPRESENTS THE START
!  OF 100 DOSE SETS THAT BEGIN AT TIME 0.
          SIG(KNS) = 0
        ENDIF
!  THE ABOVE ENDIF IS FOR THE  IF(SIG(KNS) .LT. 0.D0)  CONDITION.
        DO I=1,NI
          R(I)=RS(KNS,I)
        END DO
        if (NDRUG /= 0) then  ! IF(NDRUG .EQ. 0) GO TO 81
!  AS OF idm2x11.f: MUST CALL GETFA BEFORE EVERY TIME THAT
!  FA(.) ARE USED IN CASE THE EQUATION(S) FOR THE FA(.) ARE BASED
!  ON THE COVARIATES, WHICH CAN CHANGE DOSE TO DOSE.
          CALL GETFA(FA,X,ESTML,R,B,INTLIST)
          IF(N /= 0) then  !  IF(N .EQ. 0) GO TO 120
            DO I=1,NDRUG
              X(NBCOMP(I))=X(NBCOMP(I))+BS(KNS,I)*FA(I) ! BSO replaced w/BS
            END DO
!  NOTE THAT FA(I) IS THE FRACTION OF DRUG AVAILABLE FROM A BOLUS INPUT
!  FOR DRUG I INTO ITS ABSORPTIVE COMPARTMENT.
            ! GO TO 81
          else
120         DO I=1,NDRUG
              B(I)=BS(KNS,I)*FA(I) ! BSO replaced w/BS
            END DO
          end if
        end if
81      KNS=KNS+1
      END IF ! (JUMPTO45 /= 1) condition.
! ***** INTEGRATION OF EQUATIONS *****
      JUMPTO45 = 0
! 45    continue
!  45    IF(KNS.GT.ND) GO TO 15
      DO WHILE (KNT <= M)
      ! "IF(KNT .LE. M) GO TO 45" at end of the block is converted to DO WHILE
!      write (*,*) "at DO WHILE KNT <= M", KNT, M
!  DETERMINE IF, OBSER(ID=0), OR DOSE(ID=1), OR BOTH(ID=2).
!  Inserted label 46
46      IF(KNS > ND) then ! "GO TO 15" replaced with a copy of block 15, here
          ID=0
          TOUT=TIM(KNT)
          KNT=KNT+1
        else  ! IF(N .EQ. 0) GO TO 31
! CODE CHANGE BELOW FOR idm2x8.f.
          IF(TIM(KNT) == 0.D0 .AND. KNT > 1) THEN
!  AS OF idm2x7.f, A TIME RESET NO LONGER REQUIRES ALL INITIAL
!  COMPARTMENT AMOUNTS TO BE RESET TO 0. THIS IS BECAUSE A TIME RESET
!  NO LONGER HAS TO MEAN THAT AN "INFINITE" AMOUNT OF TIME HAS OCCURRED
!  WITH NO DOSING; IT CAN ALSO NOW MEAN THAT AN "INFINITE" AMOUNT OF
!  TIME HAS OCCURRED WITH UNKNOWN DOSING (IN THIS CASE, SUBROUTINE
!  GETIX WILL BE CALLED BELOW TO ESTABLISH INITIAL CONDITIONS FOR THIS
!  TIME PERIOD).
!  ADVANCE KNS TO THE NEXT VALUE THAT HAS SIG(KNS) .LE. 0. I.E., ONCE
!  TIMN(KNT) = 0, IT MEANS THAT WE ARE DONE WITH THE OUTPUT OBS.
!  TIMES IN THE PREVIOUS SECTION --> THERE IS NO POINT IN CONTINUING
!  TO INTEGRATE TILL THE END OF THE DOSES IN THE PREVIOUS SECTION
!  (IF THERE ARE ANY).
            DO IKNS = KNS,ND
              IF(SIG(IKNS) <= 0.D0) then ! GO TO 110
                kns = ikns
                exit ! This is command on line 110
              end if
            END DO ! Since sig(ikns) <= 0, the next block is skipped on exit.
            if (SIG(IKNS) > 0.D0) then
!  TO GET HERE MEANS THAT NO VALUE IN SIG(.) FROM KNS TO ND HAS A
!  VALUE .LE. 0, AND THIS IS AN ERROR. IT MEANS THAT THE PATIENT DATA
!  FILE HAS AN OBSERVATION TIME RESET ROW WITHOUT AN ACCOMPANYING
!  DOSE RESET ROW. TELL THE USER AND STOP.
!  REPLACE WRITING OF SIG() WITH XVERIFY (SEE LOGIC IN SUBROUTINE
!  VERIFYVAL.
              XVERIFY(1) = SIG(KNS)
              CALL VERIFYVAL(1,XVERIFY)
              WRITE(*,111) ND,KNS,XVERIFY(1)
              WRITE(25,111) ND,KNS,XVERIFY(1)
 111  FORMAT(//' THE CURRENT SUBJECT HAS AN OBSERVATION TIME RESET'/    &
     &' ROW WITHOUT AN ACCOMPANYING DOSE RESET ROW. THE PROGRAM NOW'/   &
     &' STOPS. '//                                                      &
     &' REVIEW YOUR PATIENT FILES AND CORRECT THE ERROR.'//             &
     &' NOTE THAT THE ',I4,' DOSE TIMES (POSSIBLY ALTERED BY TIMELAGS'/ &
     &' ARE THE FOLLOWING (AND THERE IS NO TIME .LE. 0 AFTER TIME'/     &
     &' NO. ',I4,' WHICH HAS CORRESPONDING TIME ',F15.4,'):')
              OPEN(42,FILE=ERRFIL)
              WRITE(42,111) ND,KNS,XVERIFY(1)
              DO I = 1,NDO
                WRITE(*,*) SIG(I)
                WRITE(25,*) SIG(I)
                WRITE(42,*) SIG(I)
              END DO
              CLOSE(42)
              CALL PAUSE
              STOP
            end if
  110       continue ! KNS = IKNS, here, I think!
!  THERE ARE TWO POSSIBILITES AT THIS POINT, EITHER SIG(KNS) = 0
!  OR SIG(KNS) < 0.
!  IF SIG(KNS) = 0, THIS REPRESENTS A TIME RESET (T WILL BE SET = 0
!  BELOW) WITH A SINGLE DOSE LINE TO START. IN THIS CASE, CALL GETIX
!  AGAIN (JUST AS WAS DONE NEAR THE TOP OF THIS ROUTINE) TO OBTAIN
!  INITIAL COMPARTMENT AMOUNTS. NOTE THAT BY DEFAULT, IN GETIX, ALL
!  COMPARTMENT AMOUNTS ARE SET = 0 (WHICH WOULD BE THE CASE IF IN THE
!  LONG TIME PERIOD BETWEEN THE LAST SET OF DOSES AND THIS NEW
!  BEGINNING, NO DOSES HAVE BEEN GIVEN). BUT THE USER MAY ALSO HAVE
!  CODED INTO GETIX EQUATIONS THAT SET ONE OR MORE OF THE X(I) TO
!  FUNCTIONS OF COVARIATE AND PARAMETER VALUES (WHICH WOULD BE THE
!  SITUATION IF AN UNKNOWN DOSING REGIMEN HAS TAKEN PLACE BUT IT
!  DOESN'T MATTER WHAT IT WAS BECAUSE THE PATIENT COMES TO A LAB AND
!  SIMPLY HAS HIS COMPARTMENT VALUES ESTABLISHED BEFORE CONTINUING
!  WITH THE OTHER VALUES IN HIS PATIENT FILE).
!  IF SIG(KNS) < 0, THIS REPRESENTS A TIME RESET WITH A STEADY STATE
!  SET OF 100 DOSES ABOUT TO BEGIN. IN THIS CASE, WE ASSUME THAT THE
!  PATIENT IS ABOUT TO GET 100 SETS OF DOSES SO THAT HIS COMPARTMENT
!  AMOUNTS WILL ACHIEVE STEADY STATE VALUES. THESE STEADY STATE VALUES
!  WILL BE ESTIMATED IN THE BLOCK OF CODE BELOW THAT STARTS WITH
!  IF(ISTEADY .EQ. 1). IN THIS CASE, WE WILL STILL CALL GETIX TO
!  MAKE SURE THAT ANY RESIDUAL COMPARTMENT AMOUNTS FROM A PREVIOUS
!  SET OF DOSES IS ZEROED OUT (OR SET = VALUES AS DETERMINED BY
!  SUBROUTINE GETIX).
!  AS OF idm2x13.f, BEFORE CALLING GETIX, MUST SET
!  THE R(.) IN CASE ANY OF THE INITIAL CONDITIONS FOR THE X(.)
!  ARE FUNCTIONS OF THE COVARIATES WHICH ARE ESTABLISHED FROM THE
!  R(.) VALUES IN GETFA.
            DO I=1,NI
              R(I)=RS(KNS,I)
            END DO
            CALL GETIX(N,X,ESTML,R,B,INTLIST)
!  MUST ALSO RESET T = 0 SINCE THE INTEGRATION WILL AGAIN START FROM
!  TIME 0.
            T = 0.D0
!  IF SIG(KNS) .LT. 0, THIS IS NOT ONLY A TIME RESET, IT IS THE
!  BEGINNING OF A STEADY STATE DOSE SET. IN THIS CASE, APPLY 100
!  STEADY STATE DOSES WITH AN INTERDOSE INTERVAL EQUAL TO -SIG(KNS).
            ISTEADY = 0
            IF(SIG(KNS) < 0.D0) THEN
              ISTEADY = 1
              NSET = 1
!  NOTE THAT ISTEADY = 1 TELLS THE PROGRAM BELOW TO PROCEED AS IF THE
!  DOSE TIME IS 0, AND START INTEGRATING THROUGH THE SET OF 100
!  DOSE SETS, ALL OF WHICH OCCUR BEFORE THE NEXT OBSERVATION TIME ...
!  BUT PAUSE AFTER THE END OF THE 5TH DOSE SET (NSET IS THE RUNNING NO.
!  OF THE CURRENT DOSE SETS THAT HAVE BEEN RUN) AND CALL SUBROUTINE
!  PREDLAST3 TO PREDICT THE STEADY STATE COMPARTMENT AMOUNTS AFTER THE
!  100 DOSE SETS (NOTE THAT THE COMPARTMENT AMOUNTS WILL HAVE TO BE
!  STORED AT THE END OF EACH OF THE STEADY STATE DOSE SETS AS THE LOGIC
!  OF PREDLAST3 REQUIRES).
!  IF "CONVERGENCE" IS ACHIEVED AT THAT POINT, ASSIGN THE COMPARTMENT
!  AMOUNTS TO BE THE PREDICTED AMOUNTS, AND ASSIGN KNS TO BE WHAT IT IS
!  WHEN THESE STEADY STATE DOSE SETS HAVE FINISHED. NOTE THAT THE END OF
!  THE 100TH DOSE SET WILL BE AT TIME 100*(-SIG(KNS)), SO KNS WILL BE
!  THE INDEX OF THE FIRST DOSE EVENT WHICH OCCURS AFTER THIS TIME.
!  IF "CONVERGENCE" IS NOT ACHIEVED, CONTINUE APPLYING THE LOGIC OF
!  PREDLAST3 UNTIL IT IS ACHIEVED, OR UNTIL THE max_SS_doses DOSE SETS ARE ALL
!  INTEGRATED THROUGH, WHICHEVER COMES FIRST.
              DOSEINT = -SIG(KNS)
!  RESET SIG(KNS) TO BE 0 SINCE THIS DOSE EVENT REPRESENTS THE START
!  OF max_SS_doses DOSE SETS THAT BEGIN AT TIME 0.
              SIG(KNS) = 0
            ENDIF
!  THE ABOVE ENDIF IS FOR THE  IF(SIG(KNS) .LT. 0.D0)  CONDITION.
          ENDIF
!  THE ABOVE ENDIF IS FOR THE IF(TIM(KNT) .EQ. 0.D0 .AND. KNT .GT. 1) CONDITION.
! IF(TIM(KNT) .NE .SIG(KNS)) GO TO 20
          IF(TIM(KNT) == SIG(KNS)) then
            ID=2
            TOUT=TIM(KNT)
            KNT=KNT+1
            KNS=KNS+1
            ! IF(N .EQ. 0) GO TO 31
            ! GO TO 30
          else
 20         IF(TIM(KNT) < SIG(KNS) .OR. SIG(KNS) <= 0) then
! 20        IF(TIM(KNT) .GT. SIG(KNS) .AND. SIG(KNS) .GT. 0) GO TO 25
15            ID=0
              TOUT=TIM(KNT)
              KNT=KNT+1
              ! IF(N .EQ. 0) GO TO 31
              ! GO TO 30
            else
25            ID=1
              TOUT=SIG(KNS)
              KNS=KNS+1
            end if
          end if
        end if !  46    IF(KNS.GT.ND) then ! GO TO 15
30    CONTINUE
!      write (*,*) "at if N \= 0", N
      IF (N /= 0) then  !  replaces IF(N .EQ. 0) GO TO 31
!--------------------------------------------------------------- 30 / 32
! --- wmy2018.06.15 These lines are lifted from USERANAL; they have to
! --- be here to make ANAL3() work.
! --- When you get a chance, go back to useranal and erase these lines
! --- there as those lines are now redundant.  Also, remove INTLIST
! --- from the USERANAL() arglist
        do III=1,max_ODE_params
          RPAR(k_dvode_reserved + III) = ESTML(III)
        end do
        do III = 1,max_RS_J
          RPAR(k_p_end + III) = R(III)
        end do
        RPAR(k_jsub) = dble(JSUB)
        RPAR(k_ig) = dble(IG)
        do III = 1,10
          IPAR(i_dvode_reserved + III) = INTLIST(III)
        end do
        IPAR(i_jsub) = JSUB
        IPAR(i_ig) = IG
!---------------------------------------------------------------
32      continue
        ! Test for Julian ---------------------
           TestData = 73 ! 0 to work
           if (TestData == 0) then
             if ((JSUB==3).and.(IG==2).and.(KNT==3).and.(ipar(i_cycle)==7)) then ! Choose these for Vori
               T_in = T
               do iii=1,max_ODE_comps
                 X_in(iii) = X(iii)
               end do
               TestData = 1
             end if
           end if
        ! Test for Julian ---------------------
        IF (N == -1) then
           CALL ANAL3(X,T,TOUT,RPAR,IPAR)
        else ! IF (N /= -1) then
           CALL USERANAL(JSUB,IG,X,T,TOUT                                 &
             , NDIM,MF,RTOL,ATOL,ESTML,R,INTLIST,IPAR,RPAR)
        endif
        ! write (*,*) "SimConc XP update Ret. w/",TOUT,X(1),X(2),X(3)
        ! Test for Julian ---------------------
           if (TestData == 1) then
             OPEN (21, FILE="ODEtest01", STATUS="REPLACE"                 &
             , FORM="UNFORMATTED", POSITION="REWIND", ACTION="WRITE")
             write (21) X_in,T_in,JSUB,IG,X,T,TOUT                        &
             , NDIM,MF,RTOL,ATOL,ESTML,R,INTLIST,IPAR,RPAR
             close(21)
             TestData = 2
             write (*,*) "NOT Stopping program after writing Julian's data"
             ! STOP ! there is no other way to preserve the new binary file
           end if
        ! End Test for Julian -----------------
        ! PK reqs X(.)>=0
             DO III = 1,max_ODE_comps
               if (X(III)<0.D0) X(III)=0.D0
             END DO
!-- Cycle past point, only if you are still optimizing
        if (ijob == 1) then
          if (IPAR(i_skip_ig) == 0) return
        end if
!  IF ISTEADY = 1, THIS IS INSIDE A STEADY STATE DOSE SET. CHECK TO SEE
!  IF TOUT IS A MULTIPLE OF DOSEINT. IF SO, RECORD THE COMPARTMENT
!  AMOUNTS. THEN, AFTER COMPARTMENT AMOUNTS HAVE BEEN STORED FOR AT
!  LEAST THE 1ST 5 MULTIPLES OF DOSEINT, STOP AND CALL SUBROUTINE
!  PREDLAST3 WHICH PREDICTS THE FINAL (STEADY STATE) COMPARTMENT AMOUNTS
!  AFTER THE LAST (100TH) DOSE SET.
!  IF PREDLAST3 HAS PREDICTED VALUES WHICH "CONVERGE", ASSIGN THE
!  PREDICTED VALUES TO X, INCREASE KNS TO BE THE INDEX OF THE FIRST DOSE
!  EVENT WHICH OCCURS AFTER THE STEADY STATE DOSE SET ENDS AND CONTINUE.
!  IF PREDLAST3 VALUES DON'T CONVERGE, CONTINUE THE PROCESS WITH
!  COMPARTMENT AMOUNTS FOR MULTIPLES 2 - 6 OF DOSEINT, TEST FOR
!  "CONVERGENCE", ETC. THIS PROCESS CONTINUES UNTIL "CONVERGENCE" IS
!  ACHIEVED FOR A SET OF 5 COMPARTMENT AMOUNTS (OR SETS OF AMOUNTS IF
!  NDRUG IS > 1), OR UNTIL ALL 100 DOSE SETS IN THE STEADY STATE REGIMEN
!  HAVE FINISHED.
        IF(ISTEADY == 1) THEN
!  THE NEXT DOSE SET END TIME IS DOSEINT*NSET. IF TOUT = DOSEINT*NSET,
!  STORE THE COMPARTMENT AMOUNTS. IF NSET .GE. 5, CALL PREDLAST3 AND
!  PROCEED AS INDICATED ABOVE.
          CALL THESAME(TOUT,DOSEINT*NSET,ISAME)
          IF(ISAME == 1) THEN
            NN = N
            IF(N == -1) NN = 3
            DO J = 1,NN
              XSTORE(NSET,J) = X(J)
            END DO
            IF(NSET >= 5) THEN
              CALL PREDLAST3(NN,NSET,XSTORE,XPRED,ICONV)
              IF(ICONV == 1) THEN
!  SINCE THE PREDICTED VALUES ARE CONSIDERED ACCURATE (I.E.,
!  "CONVERGENCE WAS ACHIEVED IN PREDLAST), RESET ISTEADY TO 0,
!  WHICH MEANS THAT THE STEADY STATE DOSES ARE FINISHED; ASSIGN THE
!  COMPARTMENT AMOUNTS TO BE THE PREDICTED VALUES; AND SET KNS TO THE
!  FIRST DOSE EVENT AFTER THE END OF THE STEADY STATE DOSE SET. ALSO,
!  SET T = THE ENDING TIME OF THE STEADY STATE DOSE SET = 100*DOSEINT,
!  SINCE THAT IS WHAT IT WOULD HAVE BEEN HAD ALL 100 DOSE SETS BEEN
!  RUN.
                ISTEADY = 0
                DO J = 1,NN
                  X(J) = XPRED(J)
                END DO
                T = 100.D0*DOSEINT
!  ADVANCE KNS TO BE THE FIRST DOSE PAST THE 100 DOSE SETS IN THIS
!  STEADY STATE SET. NOTE THAT THIS SET ENDS BEFORE 100*DOSEINT, SO
!  FIND THE FIRST SIG(.) THAT IS .GE. 100*DOSEINT, OR THAT IS = 0
!  (WHICH SIGNIFIES A TIME RESET) OR THAT IS < 0 (WHICH SIGNIFIES
!  ANOTHER STEADY STATE SET).
! ------------------------------- Replaced logic
!                DO I = KNS,ND
!                  IF(SIG(I) >= 100.D0*DOSEINT .OR. SIG(I) <= 0.D0) THEN
!                    KNSNEW = I
!                    GO TO 100
!                  ENDIF
!                END DO
!
!  TO GET HERE MEANS THAT THERE ARE NO DOSE TIMES PAST THE END OF THIS
!  STEADY STATE DOSE SET. IN THIS CASE, SET KNS TO ND+1.
!
!                KNS = ND+1
!                GO TO 200
!  100           KNS = KNSNEW
!  200           CONTINUE
! ----------------- Above is replaced w/logic below
                KNSOLD=KNS
                KNS = INTLIST(8) + 1
                DO I = KNSOLD,INTLIST(8)
                  IF(SIG(I) >= 100.D0*DOSEINT .OR. SIG(I) <= 0.D0) THEN
                    KNS = I
                    EXIT
                  ENDIF
                END DO
200             Continue
!  SET ISKIPBOL = 1 WHENEVER CONVERGENCE OCCURS IN
!  THE STEADY STATE DOSES SINCE IN THIS CASE, WE DON'T WANT TO
!  REAPPLY THE LAST BOLUS FROM THE STEADY STATE SET BELOW LABEL 83.
                ISKIPBOL = 1
              END IF  ! (ICONV .EQ. 1) CONDITION.
!  IF ICONV = 0, ISTEADY IS STILL = 1,
!  WHICH MEANS THAT THE ATTEMPT TO PREDICT THE FINAL (STEADY STATE)
!  COMPARTMENT AMOUNTS CONTINUES.
            END IF ! (NSET .GE. 5)  CONDITION.
!  SINCE ISAME = 1, THE END OF THE SET NO. NSET HAS OCCURRED -->
!  INCREASE NSET BY 1.
            NSET = NSET + 1
          END IF ! (ISAME .EQ. 1)  CONDITION.
        END IF ! (ISTEADY .EQ. 1)  CONDITION.
      end if ! (N /= 0) then {} else GO TO 31
31      CONTINUE
!  RECORD OBSERVATION AND SUPPLY NEW DOSE
      IF(ID /= 1) then ! else GO TO 35
	KNTM1=KNT-1
!  NOTE THAT THE TIME AT WHICH THE OUTPUT IS DESIRED IS TIM(KNTM1); THIS
!  IS CLEAR SINCE THE RETURNING VALUE(S) IN YT ARE PUT INTO ROW NO.
!  KNTM1 OF Y.
        CALL OUTPUT(TIM(KNTM1),YT,X,RPAR,IPAR)
        DO I=1,NOS
          YPRED(KNTM1,I)=YT(I)
2010    end do
! 20200713 Matches UseInlineSimConc=1
!        write (*,*) JSUB,IG,KNT,KNS,"YPRED(KNT,KNS)=",(YT(I),I=1,NOS)
        if (ijob == 4) then
          IF(N > 0) THEN
            DO III = 1,N
              comp_pred_amt(KNTM1,III) = X(III)
            END DO
          ENDIF
          IF(N == -1) THEN
            DO III = 1,3
              comp_pred_amt(KNTM1,III) = X(III)
            END DO
          ENDIF
        end if ! MONT114#7632,7642
! 55      IF(ID.EQ.0) GO TO 40
      end if ! (ID /= 1) condition.
  35  CONTINUE
 55   IF(ID /= 0) then ! else GO TO 40
        IF(NI /= 0) then ! else GO TO 83
          DO I=1,NI
            R(I)=RS(KNS-1,I)
          END DO
!  AS OF idm2x12.f: MUST CALL GETFA BEFORE EVERY TIME THAT
!  FA(.) ARE USED IN CASE THE EQUATION(S) FOR THE FA(.) ARE BASED
!  ON THE COVARIATES, WHICH CAN CHANGE DOSE TO DOSE.
          CALL GETFA(FA,X,ESTML,R,B,INTLIST)
        end if ! (NI /= 0) condition.
83      IF(NDRUG /= 0 .AND. N /= 0) then ! else GO TO 82
!  ADDING N .EQ. 0 TO ABOVE IF STATEMENT SHOWS CLEARLY THAT IF
!  N = 0 (IN WHICH CASE ANALYTIC SOLUTIONS ARE CODED DIRECTLY INTO
!  SUBROUTINE OUTPUT, WHICH MAKES THE COMPARTMENT AMOUNTS IRRELEVANT)
!  SETTING VALUES FOR THE COMPARTMENTS, X, IS UNNECESSARY.
!  IF ISKIPBOL = 1, DO NOT APPLY BOLUSES FROM DOSE KNS-1, SINCE THESE
!  BOLUSES WERE PART OF THE STEADY STATE DOSE SET WHICH ALREADY HAD
!  BOLUSES (EFFECTIVELY) APPLIED ABOVE WHERE "CONVERGENCE" OF THE
!  STEADY STATE DOSE SET WAS OBTAINED.
          IF(ISKIPBOL == 0) THEN
            DO I=1,NDRUG
              X(NBCOMP(I))=X(NBCOMP(I))+BS(KNS-1,I)*FA(I)
            END DO
          ENDIF
!  RESET ISKIPBOL = 0 HERE. IF IT IS NOW = 1, IT MEANS THAT
!  THE ABOVE APPLICATION OF BOLUSES WAS SKIPPED SINCE THERE HAS JUST
!  BEEN A STEADY STATE SET OF DOSES WHICH CONVERGED (AND WE DON'T
!  WANT THE LAST BOLUS DOSE REAPPLIED). BUT, GOING FORWARD, ISKIPBOL
!  SHOULD BE SET AGAIN TO 0 SO THE ABOVE APPLICATION OF BOLUSES WILL
!  OCCUR WHENEVER THERE IS A NEW BOLUS TO BE APPLIED.
          ISKIPBOL = 0
        end if ! (NDRUG /= 0 .OR. N /= 0) then ! else GO TO 82
82      CONTINUE
!  CHECK STOPPING TIME.
      end if ! (ID /= 0) then ! else GO TO 40
40    continue
      END DO ! IF(KNT .LE. M) GO TO 45 is converted to DO WHILE
! Endgame
      if (ijob==1) then
!        deallocate(R, B, X, XSTORE, XPRED, TLAG, FA, YT, XVERIFY)
!        deallocate(TIMO, SIGO, RSO, BSO)
        ! idm1::func() passed copies of TIM,SIG,RS, and BS to SimConc(),
        ! so there is no reason to preserve their values on return; SimConc()
        ! is used similar to a function, returning YPRED.
        return
      end if
!  AS OF idm2x9.f, RESTORE THE VALUES FOR ND, SIG, AND RS, IN CASE
!  THIS MODEL HAS TIME LAGS OR STEADY STATE DOSES - TO BE READY FOR THE
!  NEXT CALL TO THIS ROUTINE. Note: Observation time error should be included in
!  TIM prior to calling SimConc.
!
      ND = NDO
      DO I=1,ND
        SIG(I) = SIGO(I)
        DO J=1,NI
          RS(I,J) = RSO(I,J)
        END DO
      END DO
!  ESTABLISH THE VALUES IN EACH DRUG'S PO COLUMN TO THE CORRESPONDING
!  COLUMN IN ARRAY BS.
      DO I=1,ND
        DO J=1,NDRUG
          BS(I,J)=RS(I,2*J)
	END DO
      END DO
! DEALLOCATE previously saved (local) variables
!      deallocate(R, B, X, XSTORE, XPRED, TLAG, FA, YT, XVERIFY)
!      deallocate(TIMO, SIGO, RSO, BSO)
      RETURN
      END SUBROUTINE SimConc
! #################################################################### !
        SUBROUTINE USERANAL(JSUB,IG,X,TIN,TOUT,                         &
           NDIM,MF,RTOL,ATOL,P,R,INTLIST,IPAR,RPAR)
!  PURPOSE:
!       GIVEN X(:) at TIN, CALCULATE X(:) at TOUT, WHERE
!       X IS THE STATE VECTOR FOR THE MODEL UNDER CONSIDERATION (DEFINED
!       BY THE D.E.'S IN SUBROUTINE DIFFEQ). D.E.'S ARE SOLVED
!       USING THE LINPACK ROUTINE, DVODE.FOR (AND ASSOCIATED ROUTINES).
!       Edited for f90 by wmy.
!       THIS ROUTINE CALLS SUBROUTINE DVODE (VODE.FOR) ONCE FOR EACH
!       time POINT AT WHICH ANSWERS ARE DESIRED. NOTE THAT DVODE WILL CALL
!       SUBROUTINE DIFFEQ (SUPPLIED BY THE USER -- IT GIVES THE
!       DIFF. EQS. OF THE MODEL, XP(I)) AND, IF THE USER DESIRES,
!       SUBROUTINE JACOB (SUPPLIED BY THE USER -- IT GIVES THE
!       JACOBIAN OF PARTIAL DERIVATIVES, dXP(I)/dX(J)).
!
!       SUBROUTINES DIFFEQ AND JACOB ARE compiled into the NPAG optimization
!       program immediately prior to running -- see notes in Pmetrics documentation
!       describing the #DIFFeq block in the model file.
!  ARGUMENTS ON INPUT:
!           X - AN ARRAY OF DIMENSION max_ODE_comp. IN THE STANDARD 3-COMPARTMENT
!               MODEL,  X(1), X(2), X(3) SHOULD
!               BE SET TO THE AMOUNT OF DRUG IN THE ABSORBTION,
!               CENTRAL, AND PERIPHERAL COMPARTMENTS, RESPECTIVELY,
!               AT TIME T=TIN.
!         TIN - CURRENT VALUE OF TIME.
!        TOUT - TIME AT WHICH SOLUTION IS DESIRED.
!
!           P - Deprecated
!           R - Deprecated
!     INTLIST - Integer list from SUBROUTINE FILRED, mostly describing dose regimen
!        IPAR - Integer array passed throughout entire program
!        RPAR - Double array passed through entire program
!
!        NDIM = NO. OF STATES IN MODEL (.LE. 3 FOR NOW).
!        MF = METHOD FLAG.
!        RTOL = SCALAR RELATIVE TOLERANCE PARAMETER.
!        ATOL(I), I=1,NDIM = ABSOLUTE TOLERANCE PARAMETERS.
!  ARGUMENTS ON OUTPUT:
!           X - THE COMPARTMENT AMOUNTS AT T=TOUT.
!         TIN <- TOUT
!        use npag_utils, only: max_RS_J,max_ODE_params,max_ODE_comps
!     1   , k_dvode_reserved, k_p_end, k_ig, k_jsub, i_ig, i_ig
!     2   , i_jsub, i_skip_ig, k_dvode_reserved, i_dvode_reserved
        IMPLICIT none
! INPUT parameters ::
        integer, intent(in) :: JSUB, IG
        double precision, dimension(max_ODE_comps), intent(INOUT) :: X
        double precision,intent(inout) :: TIN, TOUT
        integer, intent(in) :: NDIM, MF
        double precision,intent(in) :: RTOL
        double precision,dimension(max_ODE_comps),intent(IN) :: ATOL
        real*8, dimension(max_ODE_params),intent(in) :: P
        real*8, dimension(max_RS_J),intent(inout) :: R
        integer, dimension(128),intent(in) :: INTLIST
        integer, dimension(257),intent(inout) :: IPAR
        double precision, dimension(257),intent(inout) :: RPAR
!  AS OF idm1x16.f, THE DIMENSION OF RWORK IS CHANGED FROM 300 TO
!  1002 TO ACCOMADATE AN NDIM (NEQ IN SUBROUTINE DVODE) OF UP TO max_ODE_comp. SO
!  CHANGE LRW BELOW TO 1002. SIMILARLY THE DIMENSION OF IWORK AND
!  LIW BELOW ARE CHANGED FROM 40 TO 50.
        EXTERNAL DIFFEQ,JACOB
        real*8, dimension(1002) ::  RWORK
        integer, dimension(50) :: IWORK
        integer :: ITOL, ITASK, ISTATE, IOPT, LRW, LIW, III
! wmy2017Oct23
! The following ThreadPrivates got the program "working" again i.e. the
! 3.5 week long search for the bug causing serial program to work
! perfectly, BUT parallel program throwing an Illegal Instruction : 4 error
! Is likely due to an inappropriate memory map imposed by OpenMP not
!
!$omp ThreadPrivate (RWORK, IWORK)
        save RWORK, IWORK
! This interface allowed the f77 version of useranal to safely call
! subroutine dvode--not sure if I still need it, but online advice
! suggests interfaces are safer than the external statement.
      INTERFACE
        SUBROUTINE DVODE (F, NEQ, Y, T, TOUT, ITOL, RtolIn,   &
     &    ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW,  &
     &    JAC, MF, RPAR, IPAR)
          EXTERNAL F
          integer, intent(IN) :: NEQ
          double precision, dimension(:), intent(INOUT) :: Y
          double precision, intent(INOUT) :: T, TOUT
          integer, intent(IN) :: ITOL
          double precision, intent(IN) :: RtolIn
          double precision, dimension(:), intent(IN) :: ATOL
          integer, intent(IN) :: ITASK
          integer, intent(INOUT) :: ISTATE
          integer, intent(IN) :: IOPT
          double precision, dimension(:), intent(INOUT) :: RWORK
          integer, intent(IN) :: LRW
          integer, dimension(:), intent(INOUT) :: IWORK
          integer, intent(IN) :: LIW
          EXTERNAL JAC
          integer, intent(IN) :: MF
          double precision, dimension(:), intent(INOUT) :: RPAR
          integer, dimension(:), intent(INOUT) :: IPAR
        end SUBROUTINE
      END INTERFACE
! -------- START OF DVODE DECLARATIONS SPECIFIC ------------
! wmy2017Oct06 -- See DVODE docs "Part iii." -- not sure if the
!   declaration is supposed to be here, in SUBROUTINE FUNC, or
!   maybe even as far back as SUBROUTINE NPAG, outside the
!   parallelized region.
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),           &
                      ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,    &
                      RC, RL1, TAU(13), TQ(5), TN, UROUND,              &
                      ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
                      L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,      &
                      LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
                      N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,    &
                      NSLP, NYH
      COMMON /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL, ETA, ETAMAX, &
       H, HMIN, HMXI, HNEW, HSCAL, PRL1, RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, L, LMAX, &
       LYH, LEWT, LACOR, LSAVF, LWM, LIWM, LOCJS, MAXORD, METH, MITER,  &
       MSBJ, MXHNIL, MXSTEP, N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT,   &
       NSLJ, NSLP, NYH
!
! Type declarations for labeled COMMON block DVOD02 --------------------
!
      DOUBLEPRECISION HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
! wmy2017Oct -- May also req. calls to DVSRCO -- see DVODE doc
!  "Interrupting and Restarting" -- which will store the above COMMON
!  variables in the following arrays
!      COMMON /DVOD01/ RVOD(48), IVOD(33)
!      COMMON /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
! wmy2017Oct06 -- I _assume_ this is required
!$omp Threadprivate(/DVOD01/,/DVOD02/)
! wmy2017Nov07 -- I _assume_ this is also required
      save /DVOD01/,/DVOD02/
! -------- END OF DVODE DECLARATIONS SPECIFIC -------------
!  THE LOGIC OF THIS CODE IS TAKEN FROM PROGRAM DESOLV3.FOR (4/28/96).
!  FOLLOWING VALUES ARE SUPPLIED TO SUBROUTINE DVODE:
!  DIFFEQ  = NAME OF SUBROUTINE COMPLETED BY USER WHICH GIVES THE D.E.'S
!            OF THE MODEL. IT MUST BE DECLARED EXTERNAL.
! TIN      = The initial value of the independent variable.
! TOUT   = First point where output is desired (.ne. TIN).
! ITOL   = 2 SINCE ATOL IS AN ARRAY.
! RTOL   = Relative tolerance parameter (scalar).
! ATOL   = Absolute tolerance parameter.
!          The estimated local error in X(i) will be controlled so as
!          to be roughly less (in magnitude) than
!             EWT(i) = RTOL*abs(X(i)) + ATOL(i)  SINCE ITOL = 2.
!          Thus the local error test passes if, in each component,
!          either the absolute error is less than ATOL (or ATOL(i)),
!          or the relative error is less than RTOL.
!          Use RTOL = 0.0 for pure absolute error control, and
!          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!          control.  Caution.. Actual (global) errors may exceed these
!          local tolerances, so choose them conservatively.
! ITASK  = 1 for normal computation of output values of X at t = TOUT.
! ISTATE = Integer flag (input and output).  Set ISTATE = 1.
! input used.
! RWORK  = Real work array of length at least..
!             20 + 16*NDIM                      for MF = 10,
!             22 +  9*NDIM + 2*NDIM**2           for MF = 21 or 22,
!             22 + 11*NDIM + (3*ML + 2*MU)*NDIM  for MF = 24 or 25.
!       ... I'LL USE AN ARRAY OF 300 (PLENTY FOR NDIM .LE. 8).
! LRW    = Declared length of RWORK (in user's DIMENSION statement).
! IWORK  = Integer work array of length at least..
!             30        for MF = 10,
!             30 + NDIM  for MF = 21, 22, 24, or 25.
!          If MF = 24 or 25, input in IWORK(1),IWORK(2) the lower
!          and upper half-bandwidths ML,MU.
!       ... I'LL USE AN ARRAY OF 40 (PLENTY FOR NDIM .LE. 8).
! LIW    = Declared length of IWORK (in user's DIMENSION).
! JACOB    = Name of subroutine COMPLETED BY USER for Jacobian matrix
!            (MF = 21 or 24). If used, this name must be declared
!            external.  If not used, pass a dummy name.
! MF     = Method flag.  Standard values are..
!          10 for nonstiff (Adams) method, no Jacobian used.
!          21 for stiff (BDF) method, user-supplied full Jacobian.
!          22 for stiff method, internally generated full Jacobian.
!          24 for stiff method, user-supplied banded Jacobian.
!          25 for stiff method, internally generated banded Jacobian.
! RPAR,IPAR = user-defined real and integer SCALARS OR arrays passed to
!             DIFFEQ AND JACOB.
! Note that the main program must declare arrays X, RWORK, IWORK,
! and possibly ATOL, RPAR, and IPAR.
!  THE FOLLOWING VALUES RETURN FROM CALLS TO SUBROUTINE DVODE.
!      X = Array of computed values of X vector (AT TIME TOUT).
!      T = Corresponding value of independent variable (normally TOUT).
! ISTATE = 2  if DVODE was successful, negative otherwise.
!          -1 means excess work done on this call. (Perhaps wrong MF.)
!          -2 means excess accuracy requested. (Tolerances too small.)
!          -3 means illegal input detected. (See printed message.)
!          -4 means repeated error test failures. (Check all input.)
!          -5 means repeated convergence failures. (Perhaps bad
!             Jacobian supplied or wrong choice of MF or tolerances.)
!          -6 means error weight became zero during problem. (Solution
!             component i vanished, and ATOL or ATOL(i) = 0.)
        ITOL=2
        ITASK=1
        ISTATE=1
        IOPT=0
        LRW=1002
        LIW=50
        do III=1,max_ODE_params
           RPAR(k_dvode_reserved + III) = P(III)
        end do
        do III = 1,max_RS_J
           RPAR(k_p_end + III) = R(III)
        end do
           RPAR(k_jsub) = dble(JSUB)
           RPAR(k_ig) = dble(IG)
        do III = 1,10
           IPAR(i_dvode_reserved + III) = INTLIST(III)
        end do
           IPAR(i_jsub) = JSUB
           IPAR(i_ig) = IG
        if ((TIN == TOUT).and.(TIN /= 0)) write (*,*) "WARNING: TIN=TOUT; JSUB=",JSUB
        CALL DVODE(DIFFEQ,NDIM,X,TIN,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,  &
                 IOPT,RWORK,LRW,IWORK,LIW,JACOB,MF,RPAR,IPAR)
        IF (ISTATE .LT. 0) THEN
!         WRITE(*,16) ISTATE
! 16      FORMAT(///' On return from DVODE, ISTATE =',I3)
         IPAR(i_skip_ig) = 0
        ENDIF
        TIN=TOUT
        RETURN
        END subroutine useranal
! #################################################################### !
      subroutine print_matrix(b,n,m)
        integer::n,m
        double precision::b(n,m) !n = # rows, m = # columns
        integer i,j
        do i=1,n; print '(11g20.12)',b(i,1:m); enddo
        ! do i=1,n ! works, too:
        !   write (*,*)(b(i,j),j=1,m)
        ! end do
      endsubroutine !2 decimal reals printed in fields of 6 spaces
! #################################################################### !
! Copied from npagranfix*.f
        SUBROUTINE PAUSE
!  THIS ROUTINE IS USED TO REPLACE A PAUSE STATEMENT, WHICH CAUSES
!  WARNINGS WHEN THIS PROGRAM IS COMPILED AND LINKED USING gfortran
!  (AND FORCES THE USER TO TYPE "go" INSTEAD OF SIMPLY HITTING THE
!  ENTER KEY).
        integer IKEY
        WRITE(*,1)
    1   FORMAT(' HIT ANY KEY TO CONTINUE: ')
        READ(*,*,ERR=10) IKEY
        IF(IKEY .EQ. 1) RETURN
   10   RETURN
        END subroutine
!------------------------------------------------------------------------------------------------------------------
! subroutine DO140
! Inputs observed YO(:,:) and predicted Y(:,:) values, where
! YO and Y are (1:NumOfOutEqns, 1:NumOfMeasurements)
! Outputs F(K) = zscore of the kth prediction
!    1:K ~ \series_i\series_j (outeq_i,measurement_j).
!    If measurement is "-99" or "-88" (BLQ or missing,
!    respectively) then F() <- 0.D0.
! RPAR(k_prod_pr), RPAR(k_sum_z_sq), RPAR(k_sfac), RPAR(k_ofac)
! RPAR(k_sum_sq_norm_err)
! Note: For all programs, DO140 is embedded in a JSUB/IG loop.
!   Usually, DO140 is called to support calculating Pr(JSUB|IG)
!    = \prod Pr(each measurement of JSUB | IG)
!    = Pr(all measurements ~ Poisson)
!     * Pr(all measurements ~ Normal)
!   As of 20190729 In main,
!    PYJGX(JSUB,IG) = Pr(\all Measurements ~ Poisson) * Pr(\all Measurements ~ Normal)
!     [ 10**RPAR(k_prod_pr) ]
!     * [ DEXP(-.5D0*RPAR(k_sum_z_sq))
!       / RPAR(k_sfac)/RPAR(k_ofac) ], where,
!    RPAR(k_prod_pr) = \sum log10(PoissonProb)
!    RPAR(k_sigfac) = \prod sd(YMEAN), \all measurements ~ Normal
!    RPAR(k_ofac) == Normalization constant for the R^N space, \Omega
! Note: ijob \in 1,2,3,4
!   1: NPAG or IT2B
!   2,3: NPAG post-proicessing statistics
!   4: Simulation
      subroutine DO140(ijob,JSUB,IG,INTLIST,IPAR,RPAR,F,ObsError,YO,Y,ERRFIL)
!      use npag_utils, only: k_sum_z_sq, k_prod_pr, k_sfac, k_ofac &
!        , k_resolve, k_gamma, k_flat &
!        , k_c0_base, k_c1_base, k_c2_base, k_c3_base, k_c4_base, k_c5_base &
!        , i_errmod, i_misval, i_Nnormalobs, i_Npoissonobs &
!        , i_is_log10, i_is_Poisson_obs
!        , max_m_per_obs, maxnumeq
      implicit none
!      integer k_sum_z_sq, k_prod_pr, k_sfac, k_ofac &
!        , k_resolve, k_gamma, k_flat &
!        , k_c0_base, k_c1_base, k_c2_base, k_c3_base, k_c4_base, k_c5_base &
!        , i_errmod, i_misval, i_Nnormalobs, i_Npoissonobs &
!        , i_is_log10, i_is_Poisson_obs
! Arguments
      integer, intent(in) :: ijob, JSUB, IG
      integer, dimension(1:), intent(in) :: INTLIST            ! dimension(128)
      integer, dimension(1:), intent(inout) :: IPAR            ! dimension(257)
      double precision, dimension(1:), intent(inout) :: RPAR   ! dimension(257)
      double precision, dimension(1:), intent(inout) :: F      ! dimension(3564)
      double precision, dimension(1:,1:), intent(inout) :: ObsError ! dimension(max_m_per_obs,maxnumeq)
      double precision, dimension(1:,1:), intent(inout) :: YO  ! dimension(max_m_per_obs,maxnumeq)
      double precision, dimension(1:,1:), intent(inout) :: Y   ! dimension(max_m_per_obs,maxnumeq)
      character*20 ERRFIL
! Local variables
      integer, parameter :: DEBUG =  0                         ! write debugging messages
      integer i,j,k                                            ! see notes above
      real pmean
      integer MISVAL, NNORMALOBS, NPOISSONOBS
      double precision YOBSERVED, YMEAN, ZSCORE, meanobsswap
      double precision, parameter :: OBSSWAPPOISSON = 0.5D0, OBSSWAPNORMAL = 0.0D0
      double precision poissonprob
      double precision, parameter :: UseNormal = 128.0D0, P_thresh = 32.0D0
!      double precision, parameter :: UseNormal = 1.0D12, P_thresh = 32.0D0
      double precision SIGFAC
      double precision ConvertSigmaToError
      logical :: PoissonInit = .TRUE.
      integer PoissonVariate
      integer TestData
      double precision, dimension(max_m_per_obs,maxnumeq) :: ObsError_in
      double precision, dimension(max_m_per_obs,maxnumeq) :: YO_in
! Block I. Initialization ----------------------------------------------
      k = 0
      MISVAL = 0
      NNORMALOBS = 0
      NPOISSONOBS = 0
      YOBSERVED = 0.D0
      YMEAN = 0.D0
      ZSCORE = 0.D0
      meanobsswap = 0.D0
      poissonprob = 0.D0
      SIGFAC = 1.D0
      TestData=0 ! if 1 then collect data and stop program
      if (TestData == 1) then
        if ((JSUB==3).and.(IG==3).and.(ipar(i_cycle)==10)) then ! Choose these for Poisson example
!          write (*,*) TestData, max_m_per_obs, maxnumeq, size(YO), size(ObsError)
          do i=1,max_m_per_obs
            do j=1,maxnumeq
              ObsError_in(i,j)=ObsError(i,j)
              YO_in(i,j)=YO(i,j)
            end do
          end do
        end if
      end if
! Reset return variables
      RPAR(k_sum_z_sq) = 0.D0
      RPAR(k_prod_pr) = 0.D0
      RPAR(k_sfac) = SIGFAC
      IPAR(i_misval) = MISVAL
      IPAR(i_NBLQ) = 0
      IPAR(i_Nnormalobs) = NNORMALOBS
      IPAR(i_Npoissonobs) = NPOISSONOBS
      if (DEBUG /= 0) then
        write (*,*) "In DO140(Debug/=0)", DEBUG, JSUB, IG
          ! 101 to check for C values
          ! 102, 103, ... 109
          ! 199
      endif
!      write (*,*) IG, "In DO140, r.v.s:",(RPAR(k_dvode_reserved + i), i = 1,5)
! Block II. Iterate through (Observation, Equation) ---------------------
!       DO 140 I=1,NOBSER
!         DO 140 J=1,NUMEQT
        DO J=1,INTLIST(9)     ! NUMEQT or number of output equations
            if (DEBUG == 101) then
              if (NINT(RPAR(J+k_c4_base)).eq.10) then  ! wmy: worked, 1/16/2019
                write (*,*) "Obs is log10",I,J
              endif
              if (NINT(RPAR(J+k_c5_base)).eq.229) then ! wmy: worked, 1/16/2019
                write (*,*) "Obs is Poisson",I,J
              endif
              write (*,*) "Checking C4 and C5", &
                RPAR(J+k_c4_base), RPAR(J+k_c5_base)   ! wmy: worked, 1/16/2019
            endif ! (DEBUG == 101)
           PoissonInit = .TRUE. ! ijob==4; each output equation has constant mu
          DO I=1,INTLIST(10)    ! Number of measurements, NOBSER
            k = (J-1)*INTLIST(10)+I ! k_0 = 1, k = k + 1
! YOBSERVED ~ BLQ -----------------------------------------------------
! 20200626 re: ijob==4: YO is 0.d0, Y is 2.0 -- the observations of 10.44 in the pt*.csv are lost
! fixed: set YO <- /OBSER/Y prior to call DO140
!            write (*,*) JSUB,IG,ijob,"on entry; eqn, m, obs, pred:",J,I,YO(I,J),Y(I,J)
! Output seems right, YO = Y = log10(.) ... CORRECT!
!
            if(int(YO(I,J)) == int(obs_is_BLQ)) then ! -88
               IPAR(i_NBLQ) = IPAR(i_NBLQ) + 1
! 2018Sep01 :: MN wants -99 to indicate BLQ, and -88 to indicate missing.
!    RPAR(130) is available, I think, to use for the Pr(obs<LQ). But we
!    also need a space for the LOQ -- maybe use RPAR(131) ?
! As of, 20200624 we have not coded BLQ in the datafile. All we have is
!    -99 and we use that for "missing"
               if (DEBUG == 199) then
                  write (*,*) JSUB,IG,J,I,"(JSUB,IG,EqnNo,m) is BLQ"
               endif
! Poisson
          ! Missing Poisson varbs < small# are automatically set to count = 0.
          ! If small# is too large, then pr of Yobs == 0 is too small and the
          ! program will fail b/c pr(JSUB|G), which is a product of all measures
          ! for the observation, will -> 0. So, we want to make sure that we
          ! only automatically set Yobs to 0 if the expected value is close to 0.
          ! For Yexp = 9, pr(yobs = 0) = 0.0001234.
          ! For Yexp = 2, pr(yobs = 0) = 0.1353353.
          ! For Yexp = 1, pr(yobs = 0) = 0.3678794.
          !
          ! Above strategy can not work: The YO that pass the sieve and are processed
          ! will have a lower likelihood than the YO that do not pass the sieve!
          ! Why? Because the ones that do not pass are not processed, but the
          ! ones that do, have prob < 1.D0. SOLUTION: You have to preprocess data
          ! and then process all obs (regardless of prediction) that for some \theta
          ! passed the sieve.  Alternatively, you can use the sieve below to remove
          ! some \theta from consideration. Do this either by setting prob above 1.1
          ! for passed \theta, or to zero for unpassed data.
          !
          !  if ( (YO(I,J) == -99) &
          !    .and. (IPAR(i_is_poisson_obs + J) == 229) &
          !    .and. (Y(I,J) <= 1) ) then
          !       YO(I,J) = 0.D0
          !  end if
              if (IPAR(i_is_poisson_obs + J) == 229) then
                if (ijob /= 4) then ! not simulating
                  if (ipar(i_is_log10 + J) /= 10) then
                    YO(I,J) = 0.D0 ! Nothing was seen; Normal variate = 0
                  ! else
                  !   do nothing, b/c YO = log10(0) -> large negative number
                  !   and BLQ = -88, a large negative number. But during analysis
                  !   check and reset YO to 0 there; b/c 10^0 = 1, can't do
                  !   it here.
                  end if
                ! else ! (ijob==4) and YO is the simulation template value
                !   do nothing; missing data is possible; but not BLQ
                end if
              else
              ! Obs is Normal -- carry the -88 forward.
              !
              !  A = <prob this obs is BLQ>
              !  A = 1/2 * ERF((LQ - 0)/sigma)
              !  RPAR(130) = RPAR(130) * A
                 YO(I,J) = 0.D0 ! A cluge b/c we don't know the LQ
              end if
            end if ! BLQ -88
! if YOBSERVED ~ MISSING -----------------------------------------------
            if (int(YO(I,J)) == obs_is_missing) then
              F(k) = 0.D0
              MISVAL = MISVAL+1
              YOBSERVED = obs_is_missing
              YMEAN = Y(I,J)
              ObsError(I,J) = 0.D0
              ZSCORE = 0.D0
              if (ijob == 4) then ! simulating; retain "-99"
                Y(I,J) = YO(I,J)
              end if
              if (DEBUG == 102) write (*,*) JSUB,IG,J,I,"is missing"
! else YOBSERVED is an acceptible measurement --------------------------
            else
! Assign YOBSERVED and YMEAN (predicted) for this measurement
!  If data.csv contains log10(obs) entries, then C4(J) should
!    be set to 10 in the model.txt file; also, SUBROUTINE OUTPUT
!    should convert predicted measurement to log10(pred).
!  Thus, at this point, calculating the standard deviations
!    requires conversion of Y and YO (both of which are log10)
!    to untransformed values:
              IF (ipar(i_is_log10 + J) == 10) THEN
                YMEAN = 10**Y(I,J)
                IF (ijob /=4 ) THEN
                  if (YO(I,J) < 0.d0) then ! obs is -88 (-99 was treated above)
                    YO(I,J) = 0.d0
                     ! This case covers (IPAR(i_is_poisson_obs + J) == 229)
                     ! And for now, the Normal case, too. But Normal obs should
                     ! really be treated w/an erf()
                  else
                    YOBSERVED = 10**(YO(I,J))
                  endif
                ELSE
                  YOBSERVED = YMEAN ! The simulated value
                ENDIF
              ELSE
                YMEAN = Y(I,J)
                IF (ijob /=4 ) THEN
                  YOBSERVED = YO(I,J)
                ELSE
                  YOBSERVED = YMEAN ! The simulated value
                ENDIF
              ENDIF
! ijob==4
! Y <- YPRED_SIMBIG, which is the OUTEQ() _after_ applying process errors
! YO <- YO_SIMBIG, and note: YO = Y (\= YPRED) prior to calling DO140.
! YO = Observed value from the SIMBIG templates, ZMQtemp.csv and XQZPJ001.ZMQ (but not abcde12.csv)
! YMEAN = linear measurement, so if OUTEQ() applied a LOG10 then YMEAN = 10^Y
! YOBSERVED = YMEAN ; and is correct to 17 digit precision
! for the Poisson test:
!   YPRED = Normal random variate, which will then be used as \lambda to
!     generate a Poisson variate, and
!   YO = 10^(mean of the simulated Normal distribution)
!              write (*,*) JSUB,"is log10 converted",Y(I,J),YO(I,J),YMEAN,YOBSERVED
! starts YOBSERVED ~ POISSON -------------------------------------------
! If the data is believed to best be described by a Poisson distribution,
!  the user should set C5 <- 229 (the extended ASCII code for lower
!  case lambda. Note that C4 and C5 are NOT recognized in the data.csv
!  file by the prep program, even though C4 and C5 are required to be in
!  the model.txt file bv pmetrics::NPrun.
! If YOBSERVED .ge. UseNormal, then skip to YOBSERVED ~ NORMAL. The
!  observation will be treated as normal with mean YMEAN and standard
!  deviation determined according to ierrmod
              IF ( (IPAR(i_is_poisson_obs + J) == 229) &
                .and. (YOBSERVED .lt. UseNormal) )  THEN
!               Poisson: p( YO(i,J) | Y(I,J) ) =  exp( -1 * Y(I,J) )
!                              * ( Y(I,J)^YO(I,J) / fact(YO(I,J)) ))
!                where, fact(int N) = gamma(real N+1), and
!                       mean = var = Y.
! if ( Yexp = 17 and Yobs = 0) then pr = 4.139938e-08 -> 0.D0
!               If RESOLVE is "high" then use max{YMEAN,YOBSERVED}
!               as surrogate for the mean.  For Poisson, var = mean,
!               therefore, choosing max{.} increases likelihood(IG|JSUB),
!               and increases chance for this support point to survive.
                if (RPAR(k_resolve) .gt. OBSSWAPPOISSON ) then
                  if (YOBSERVED .gt. YMEAN) then
                    meanobsswap = YOBSERVED
                    YOBSERVED = YMEAN
                    YMEAN = meanobsswap
                    write (*,*) JSUB,IG,I,J,"Swapped Poisson Obs (mu,obs)=",YMEAN,YOBSERVED
                  endif
                endif
                if (ijob == 4) then
                  ConvertSigmaToError = ObsError(I,J)
                  ObsError(I,J) = sqrt( YMEAN )
                  if (ObsError(I,J) < 1.d-8) ObsError(I,J) = 1.d-8
                  Y(I,J) = YMEAN + ConvertSigmaToError * ObsError(I,J)
                  IF(Y(I,J) < 0.0) Y(I,J) = 0.0
                  YOBSERVED = Y(I,J)
                ! YO is the template.csv file value (may be log10, or -99)
                ! YMEAN value that came out of OUTEQ(), or 10^value if value is log10
                ! YOBSERVED ~ N(YMEAN, unspec. obs. ERR)
                ! note: YO and YMEAN are ~ N(mean = PK model, var due to process noise)
                ! Y = Poisson(\lambda=YOBSERVED), below
                ! Now, convert YPRED_SIMBIG to a Poisson variate:
                  pmean = sngl(Y(I,J))
                  PoissonVariate = random_Poisson(pmean, PoissonInit)
                ! The initial Poisson validations have constant mean for each output
                ! equation. So set PoissonInit to FALSE here, and update to TRUE for
                ! each equation. but for PK models generally, the mean will be different
                ! every time the routine is called; so PoissonInit should always be true.
                !  PoissonInit = .FALSE. ! Only for constant valued mu
                  if (pmean > 1.0 .and. PoissonVariate > 10.0*pmean) then
                    do while (PoissonVariate > 10.0 * pmean)
                      PoissonVariate = random_Poisson(pmean, PoissonInit)
                      write (*,*) PoissonVariate, pmean, "stuck in loop for", JSUB
                    end do
                  end if
                  Y(I,J) = dble(PoissonVariate)
! debug 20200706
!                  if ( (pmean > 7).and.(pmean < 10) ) then
!                    write (*,*) "pmean=",pmean            &
!                    , "Poisson variate=",PoissonVariate   &
!                    , "Y(I,J)=",Y(I,J),I,J                &
!                    , "PoissonInit=",PoissonInit
!                  end if
!           if (JSUB == 64) then
!                write (*,*) JSUB,I,J,"Yobs, Ymu, Y, YO", &
!                     YOBSERVED,YMEAN,Y(I,J),YO(I,J),ConvertSigmaToError,sqrt(YMEAN)
!           endif
                else ! ijob /= 4: we are not simulating; we are optimizing.
                  ObsError(I,J) = rpar(k_gamma) * sqrt( YMEAN )
                endif
                if (DEBUG == 109) then
                  write (*,*) JSUB,IG,J,I,"(JSUB,IG,EqnNo,m,YO,mu,s)"   &
                  , YOBSERVED, YMEAN, ObsError(I,J), rpar(k_gamma)
                endif
          ! write (*,*) "PTST",Y(I,J),YO(I,J),YMEAN,YOBSERVED,ObsError(I,J),PoissonVariate
                IF (YMEAN > P_thresh) THEN ! Use Normal approx. to Poisson
                  NNORMALOBS = NNORMALOBS+1
                  SIGFAC = SIGFAC * ObsError(I,J)
                  RPAR(k_sum_z_sq) = RPAR(k_sum_z_sq) &
                   + ((YOBSERVED - YMEAN)/ObsError(I,J))**2
! 20200716 These seem correct
!              if ((IG == 786 ).and.(IPAR(i_cycle)==1)) then
!                write (*,*) IG,"sum z^2=",RPAR(k_sum_z_sq), YMEAN, YOBSERVED, ObsError(I,J)
!              endif
                ELSE !           Use Poisson
! wmy2019.12.16 if obs is true Poisson, then prob > 0.01 (always)
! wmy2019.12.16 if obs is true Poisson, then \sum prob <= 12
                  poissonprob = 0.00001 + 0.99999 * exp( -1.D0 * YMEAN ) &
                    * ( YMEAN**YOBSERVED ) / gamma( YOBSERVED + 1)
! wmy2019.01.29 -- if Pr -> 0, then don't count probability in product.
!   i.e. ignore this observation. Also, delete this support point from
!   potential list of active supports
                  if (poissonprob .gt. 0.0000001 ) then
                    NPOISSONOBS = NPOISSONOBS + 1
! wmy2019.02.05
!   prod of probabilities of Poisson sampled data almost always -> 0.D0.
!   so use RPAR(k_prod_pr) = \sum log10(poissonprob); thus, in main
!   the prod(probabilities) is 10^RPAR(k_prod_pr).
                    RPAR(k_prod_pr) = RPAR(k_prod_pr) + log10(PoissonProb)
                    if (DEBUG == 104) then
                      write (*,*) "JSUB,IG,M,EQN",JSUB,IG,I,J &
                      , "CAUGHT PoissonProb and SumLog10 =" &
                      , poissonprob,RPAR(k_prod_pr)
                    endif
                  else
                    IPAR(i_skip_ig) = 0
                    if (DEBUG == 105) then
                      write (*,*) "JSUB,IG,M,EQN",JSUB,IG,I,J &
                      , "is P w/mu,obs,err,Pr=,\Sum log10(Pr)" &
                      , YMEAN,YOBSERVED,ObsError(I,J) &
                      , poissonprob,RPAR(k_prod_pr)
                    endif
                  endif
                ENDIF
! ends YOBSERVED ~ Poisson ---------------------------------------------
! starts YOBSERVED ~ NORMAL --------------------------------------------
              ELSE
! ijob==4
! Y <- YPRED_SIMBIG, which is the OUTEQ() _after_ applying process errors
! YO <- YO_SIMBIG, and note: YO = Y (\= YPRED) prior to calling DO140.
! YO = Observed value from the SIMBIG templates, ZMQtemp.csv and XQZPJ001.ZMQ (but not abcde12.csv)
! YMEAN = linear measurement, so if OUTEQ() applied a LOG10 then YMEAN = 10^Y
! YOBSERVED = YMEAN ; and is correct to 17 digit precision
! for the Poisson test:
!   YPRED = Normal random variate, which will the be used as \lambda to
!     generate a Poisson variate, and
!   YO = 10^(mean of the simulated Normal distribution)
!              write (*,*) JSUB,"is log10 converted",Y(I,J),YO(I,J),YMEAN,YOBSERVED
                if ( RPAR(k_resolve) > OBSSWAPNORMAL ) then ! use YOBS as surrogate for mean
                  meanobsswap = YOBSERVED
                  YOBSERVED = YMEAN
                  YMEAN = meanobsswap
                endif
                if (ijob == 4) then
                  ConvertSigmaToError = ObsError(I,J)
                endif
!       if(ierrmod.eq.1) sig(i,j) = ...
                ObsError(I,J) = RPAR(J+k_c0_base) &
                + RPAR(J+k_c1_base)*YMEAN &
                + RPAR(J+k_c2_base)*YMEAN*YMEAN &
                + RPAR(J+k_c3_base)*YMEAN*YMEAN*YMEAN
!  NOTE THAT, THEORETICALLY, SIG SHOULD BE A CUBIC FNT. OF THE 'TRUE'
!  OBSERVED VALUES, NOT THE 'NOISY' OBSERVED VALUES or the predictions
!  (BUT THE 'TRUE' VALUES ARE UNKNOWN).
!       if(ierrmod.eq.2) sig(i,j) = sig(i,j)*gamma
!       if(ierrmod.eq.3) sig(i,j) = dsqrt(sig(i,j)**2 + gamma**2)
!       if(ierrmod.eq.4) sig(i,j) = gamma*flat
              if(IPAR(i_errmod).eq.2) ObsError(i,j) = &
                  ObsError(i,j)*RPAR(k_gamma)
              if(IPAR(i_errmod).eq.3) ObsError(i,j) = &
                  dsqrt(ObsError(i,j)**2 + RPAR(k_gamma)**2)
              if(IPAR(i_errmod).eq.4) ObsError(i,j) = &
                  RPAR(k_gamma)*RPAR(k_flat)
              IF((ObsError(I,J) == 0).and.(ijob/=4)) THEN
                WRITE(*,2345) JSUB
                WRITE(25,2345) JSUB
2345            FORMAT(//' A S.D. IS 0 FOR JSUB = ',I5,'. RERUN THE &
                           PROGRAM WITH C0 NOT = 0 FOR THIS SUBJECT, &
                           OR WITH THIS SUBJECT ELIMINATED.')
                CLOSE(27)
                CLOSE(25)
!  ABNORMAL TERMINATION; WRITE THE ERROR MESSAGE TO ERRFIL.
                OPEN(42,FILE=ERRFIL)
                WRITE(42,2345) JSUB
                CLOSE(42)
                CALL PAUSE
                STOP
              ENDIF
              IF(ObsError(I,J) .LT. 0) THEN
                WRITE(*,2346) JSUB,IG
                WRITE(25,2346) JSUB,IG
2346            FORMAT(//' A S.D. < 0 FOR JSUB,IG = ',I5,I6,'. &
                   RERUN THE PROGRAM WITH A BETTER CHOICE FOR THE &
                   ASSAY ERROR POLYNOMIAL COEFFICIENTS.')
                CLOSE(27)
                CLOSE(25)
! ABNORMAL TERMINATION; WRITE THE ERROR MESSAGE TO ERRFIL.
                OPEN(42,FILE=ERRFIL)
                WRITE(42,2346) JSUB,IG
                CLOSE(42)
                CALL PAUSE
                STOP
              ENDIF
!--------------  SIGFAC and NNORMALOBS
              if (ijob == 4) then
                YO(I,J) = Y(I,J)
                Y(I,J) = YMEAN + ObsError(I,J)*ConvertSigmaToError
                IF(Y(I,J) < 0.0) Y(I,J) = 0.0
                ! prevent catastrophic errors in case where Cn = 0
                if (ObsError(I,J) .lt. 0.D000001) ObsError(I,J) = 0.D000001
                YMEAN = Y(I,J)
                ! YMEAN and Y ~ PK model + Process Noise + Assay ERR
                ! YOBSERVED and YO contain only process noise
                ! For simulation, Y is returned as the predicted amount
                ! while YO is the template amount.  Beware that YO may be
                ! log10; but YOBSERVED is an absolute amount.
              endif
!
! 20200717
! MUST REDO CALC as product of each measurement ... separating out leads to
! huge numbers that are too difficult to deal
!
              SIGFAC=SIGFAC*ObsError(I,J)
              NNORMALOBS=NNORMALOBS+1
              RPAR(k_sum_z_sq) = RPAR(k_sum_z_sq)  &
                 + ((YOBSERVED-YMEAN)/ObsError(I,J))**2
            endif ! if(Poisson){...} else {Normal}
! ends YOBSERVED ~ NORMAL ----------------------------------------------
            ZSCORE=(YOBSERVED-YMEAN)/ObsError(I,J)
            F(k)=ZSCORE ! F((J-1)*INTLIST(10)+I)=ZSCORE
            RPAR(k_sum_norm_err_sq) = RPAR(k_sum_norm_err_sq) &
                 + ((YOBSERVED-YMEAN)/YMEAN)**2
            if (DEBUG == 107) then
              write (*,*) "JSUB,IG,...",JSUB,IG,I,J,YOBSERVED,YMEAN,zscore
            endif
          ENDIF ! if(YO(I,J) .EQ. obs_is_missing){...} else {YO ~ {Normal or Poisson}}
! Last checked on 202300717 Correct. Note Cn are read from csv file, and
! if not there, then read from the model file. So, be sure to set csv file
! constants to '.' if you want to use the model.txt #ERR block numbers.
          if ((IG == 1 ).and.(IPAR(i_cycle)==1)) then
           if (DEBUG==108) then
              write (*,*) "idm01",I,J,ZSCORE,RPAR(k_prod_pr) &
                , YMEAN, YOBSERVED, ObsError(I,J), YO(I,J), Y(I,J) &
                , NNORMALOBS, SIGFAC, RPAR(J+k_c0_base) &
                , RPAR(J+k_c1_base), RPAR(J+k_c2_base) &
                , RPAR(J+k_c3_base)
           endif
          endif
! 140 CONTINUE  ! is replaced by two END DO
          END DO ! number of measurements
        END DO  ! number of output equations
! Block III. Cleanup ---------------------------------------------------
! Final checks on NNORMALOBS, NPOISSONOBS, MISVAL, SIGFAC, and OFAC
! Copy values into appropriate IPAR and RPAR bins for communication
! back up to main.
! RPAR(k_prod_pr) = \prod Poisson Probs
! RPAR(k_sum_z_sq) = \sum Squared Normalized deviations for Y_obs ~ Normal
        RPAR(k_sfac) = SIGFAC
        RPAR(k_ofac) = 2.506628274631**NNORMALOBS ! = OFAC
        IPAR(i_misval) = MISVAL
        IPAR(i_Nnormalobs) = NNORMALOBS
        IPAR(i_Npoissonobs) = NPOISSONOBS
        ! Test for Julian --------------------
           if ((JSUB==3).and.(IG==3).and.(ipar(i_cycle)==10)) then ! Choose these for Poisson example
           if (TestData == 1) then
             OPEN (21, FILE="VARtest01", STATUS="REPLACE"                  &
             , FORM="UNFORMATTED", POSITION="REWIND", ACTION="WRITE")
             write (21) k_prod_pr, k_sum_z_sq, k_sfac, k_ofac              &
             , i_misval, i_Nnormalobs, i_Npoissonobs                       &
             , ijob, JSUB, IG, 128, INTLIST, 257, IPAR, 257, RPAR, 3564, F &
             , max_m_per_obs, maxnumeq, ObsError, YO, Y, ObsError_in, YO_in
             close(21)
        !     STOP ! there is no other way to preserve the new binary file
           end if
           end if
        ! End Test for Julian -----------------
! 20200714 -- This print is in idm1*
!       DO 140 I=1,NOBSER
!         DO 140 J=1,NUMEQT
!        write (*,*) IPAR(i_cycle),JSUB,IG,((ObsError(I,J),I=1,INTLIST(10)),J=1,INTLIST(9))
      RETURN
      end subroutine DO140
!------------------------------------------------------------------------------------------------------------------
      subroutine cp_theta_to_rpar(JSUB,IG,theta,rpar)
! Copies the IG-th support point to rpar(k_p_start:k_p_end)
! k_p_start = k_dvode_reserved + 1
      implicit none
      integer JSUB, IG, III
      double precision, intent(in), dimension(1:) ::  theta
      double precision, intent(inout), dimension(1:) ::  rpar
      do III=1,max_ODE_params
        RPAR(k_dvode_reserved + III) = theta(III)
      end do
      end subroutine cp_theta_to_rpar
!------------------------------------------------------------------------------------------------------------------
      subroutine cp_RS_J_to_rpar(JSUB,IG,RS_J,rpar)
! Copies the JSUB's input to rpar(k_r_start:k_r_end)
! k_r_start = k_p_end + 1
!   k_r_end = k_p_end + max_RS_J
      implicit none
      integer JSUB, IG, III
      double precision, intent(in), dimension(1:) ::  RS_J
      double precision, intent(inout), dimension(1:) ::  rpar
      do III=1,max_RS_J
        RPAR(k_p_end + III) = RS_J(III)
      end do
      end subroutine cp_RS_J_to_rpar
!------------------------------------------------------------------------------------------------------------------
      subroutine cp_lrcs_to_rpar(gamma,flat,numeqt,c0,c1,c2,c3,c4,c5,rpar)
      implicit none
! Copies regression error arguments to rpar, placing them immediately
! after IG (see k_xxx above).  The constants are stored in k_c_gamma,
! k_c_flat, and k_c_start = k_flat + 1 to k_c_end = k_flat + maxnumeq * 4,
! as cat(gamma,flat,C0,C1,C2,C3).
      double precision gamma, flat
      integer numeqt
      double precision c0(:), c1(:), c2(:), c3(:), c4(:), c5(:), rpar(:)
      integer III
!     integer Ic0base,Ic1base,Ic2base,Ic3base
! Verify that NUMEQT == 7
        if (numeqt.gt.maxnumeq) then
           write (*,*) "ERROR :: NUMEQT.gt.maxnumeq",   &
             NUMEQT
           return
        end if
        RPAR(k_gamma) = gamma
        RPAR(k_flat) = flat
!        write (*,*) k_p_start, k_p_end, k_r_start, k_r_end
!        write (*,*) k_jsub, k_ig, k_gamma, k_flat
!        write (*,*) k_c_start, k_c_end, k_sfac, k_ofac
!        Ic0base = k_flat
!        Ic1base = Ic0base + maxnumeq
!        Ic2base = Ic1base + maxnumeq
!        Ic3base = Ic2base + maxnumeq
!        write (*,*) "base:", Ic0base, Ic1base, Ic2base, Ic3base
        do III=1,maxnumeq
           if (III.le.numeqt) then
             RPAR(k_c0_base + III) = C0(III)
             RPAR(k_c1_base + III) = C1(III)
             RPAR(k_c2_base + III) = C2(III)
             RPAR(k_c3_base + III) = C3(III)
             RPAR(k_c4_base + III) = C4(III)
             RPAR(k_c5_base + III) = C5(III)
           else
             RPAR(k_c0_base + III) = 0.d0
             RPAR(k_c1_base + III) = 0.d0
             RPAR(k_c2_base + III) = 0.d0
             RPAR(k_c3_base + III) = 0.d0
             RPAR(k_c4_base + III) = 0.d0
             RPAR(k_c5_base + III) = 0.d0
           end if
!           write (*,*) "cp:", RPAR(k_c0_base + III), RPAR(k_c1_base + III), &
!             RPAR(k_c2_base + III), RPAR(k_c3_base + III)
        end do
      end subroutine cp_lrcs_to_rpar
! #################################################################### !
      logical function check_input_array_size(                &
       MAXSUB,MAXGRD,MAXDIM,MAXACT,NUMEQT,MAXOBS,WORK,WORKK,  &
       SPXGYJ,DXI,PYJGX,PYJGXX,DENSTOR,EXX,CORDEN,CORHOLD,    &
       YPREDPOP,YPREDPOPT,YPREDBAY,CORDLAST)
      implicit none
      integer maxsub,maxgrd,maxdim,maxact,numeqt,maxobs
      double precision, intent(in), dimension(1:) :: work,workk,spxgyj,dxi
      double precision, intent(in), dimension(1:,1:) :: pyjgx
      double precision, intent(in), dimension(1:) :: pyjgxx
      double precision, intent(in), dimension(1:,1:) :: denstor
      double precision, intent(in), dimension(1:,1:,1:) :: exx
      double precision, intent(in), dimension(1:,1:) :: corden,corhold
      double precision, intent(in), dimension(1:,1:,1:,1:) :: ypredpop,ypredpopt,ypredbay
      double precision, intent(in), dimension(1:,1:) :: cordlast
! Below are the declared dimensions in subroutine NPAG. We only
! neec to verify the six input integers are .le. the maximums
! declared in the module, above. But to be careful, this
! function will check the actual memory allocated by each array.
!
!        DIMENSION WORK(MAXGRD),WORKK(MAXGRD),
!     1  SPXGYJ(MAXGRD),DXI(MAXGRD),PYJGX(MAXSUB,MAXACT),
!     2  PYJGXX(MAXACT),DENSTOR(MAXGRD,4),
! double precision, dimension(MAXSUB,3,max_pop_rand_varbs) :: EXX
!     3  CORDEN(MAXGRD,MAXDIM+1),CORHOLD(MAXGRD,MAXDIM+1),
!     4  YPREDPOP(MAXSUB,NUMEQT,MAXOBS,3),
!     5  YPREDPOPT(MAXSUB,NUMEQT,7201,3),
!     6  YPREDBAY(MAXSUB,NUMEQT,MAXOBS,3),
!     7  CORDLAST(MAXGRD,MAXDIM+1)
        write (*,*) "MAXGRD =", MAXGRD
        write (*,*) "MAXSUB =", MAXSUB
        write (*,*) "MAXACT =", MAXACT
        write (*,*) "MAXDIM =", MAXDIM
        write (*,*) "NUMEQT =", NUMEQT
        write (*,*) "MAXOBS =", MAXOBS
        write (*,*) "size(work) =", shape(WORK)
        write (*,*) "size(workk) =", shape(WORKK)
        write (*,*) "size(spxgyj) =", shape(SPXGYJ)
        write (*,*) "size(dxi) =", shape(dxi)
        write (*,*) "size(pyjgx) =", shape(PYJGX)
        write (*,*) "size(pyjgxx) =", shape(PYJGXX)
        write (*,*) "size(denstor) =", shape(DENSTOR)
        write (*,*) "size(exx) =", shape(EXX)
        write (*,*) "size(corden) =", size(CORDEN)
        write (*,*) "size(corhold) =", size(CORHOLD)
        write (*,*) "shape(ypredpop) =", shape(YPREDPOP)
        write (*,*) "shape(ypredpopt) =", shape(YPREDPOPT)
        write (*,*) "shape(ypredbay) =", shape(YPREDBAY)
        write (*,*) "shape(cordlast) =", shape(CORDLAST)
        write (*,*) "size(shape(cordlast)) =", size(shape(CORDLAST))
        write (*,*) "size(cordlast,1) =", size(CORDLAST,1)
        write (*,*) "size(cordlast,2) =", size(CORDLAST,2)
        check_input_array_size = .true. ! Assume all array dimensions are .le. limit
      end function check_input_array_size
! ##################################################################### !
      subroutine makevec(NVAR,NOFIX,NRANFIX,IRAN,X,VALFIX,RANFIXEST,PX)
!  THIS ROUTINE, CALLED BY MAIN, INPUTS NVAR, NOFIX, NRANFIX, IRAN, X,
!  VALFIX, AND RANFIXEST, AND RETURNS PX(I) = A COMBINATION OF THE
!  VALUES IN X, VALFIX, AND RANFIXEST, IN THE PROPER ORDER (AS
!  DETERMINED BY IRAN).
!        IMPLICIT NONE
! wmy2018.9.22
! F77 parameter declarations replaced w/F90 assumed shape arrays
!       DIMENSION IRAN(max_ODE_params),X(max_pop_rand_varbs),
!         VALFIX(20),PX(max_ODE_params),RANFIXEST(20)
        integer nvar,nofix,nranfix
        integer, intent(in), dimension(1:) :: iran
        double precision, intent(in), dimension(1:) ::  x,valfix,ranfixest
        double precision, intent(out), dimension(1:) :: px
! ---
        integer NNNVAR, NNNFIX, NNNRANFIX, I
        I = size(iran)
        if (I .gt. max_ODE_params) then
          write (*,*) "in makevec() :: size(iran) > max_ODE_params", I, max_ODE_params
        end if
        I = size(X)
        if (I .gt. max_pop_rand_varbs) then
          write (*,*) "in makevec() :: size(X) > max_pop_rand_varbs", I, max_pop_rand_varbs
        end if
        I = size(valfix)
        if (I .gt. max_pop_params) then
          write (*,*) "in makevec() :: size(valfix) > max_pop_params", I, max_pop_params
        end if
        I = size(ranfixest)
        if (I .gt. max_pop_varbs) then
          write (*,*) "in makevec() :: size(ranfixest) > max_pop_params", I, max_pop_params
        end if
        if (NVAR+NOFIX+NRANFIX .gt. max_ODE_params) then
          write(*,*) "makevec() :: NVAR+NOFIX+NRANFIX .gt. max_ODE_params",    &
          NVAR+NOFIX+NRANFIX, max_ODE_params
        end if
        NNNVAR = 0
        NNNFIX = 0
        NNNRANFIX = 0
        DO I = 1,NVAR+NOFIX+NRANFIX
        IF(IRAN(I) .EQ. 1) THEN
          NNNVAR = NNNVAR+1
          PX(I) = X(NNNVAR)
        ENDIF
        IF(IRAN(I) .EQ. 0) THEN
          NNNFIX = NNNFIX+1
          PX(I) = VALFIX(NNNFIX)
        ENDIF
        IF(IRAN(I) .EQ. 2) THEN
          NNNRANFIX = NNNRANFIX+1
          PX(I) = RANFIXEST(NNNRANFIX)
        ENDIF
        END DO
        RETURN
      end subroutine
! ##################################################################### !
! expand.f90
!
! Expansion algorithm in NPAG
!
! Contains subroutine expand_grid(), function checkd()
!
!------------ end subroutine checkd() ----------------------------------
      integer function checkd(corden,new,nactveold,nvar,ab)
      implicit none ! real*8 (a-h,o-z)
      double precision, dimension(1:,1:), intent(in) :: corden
      integer, intent(in) :: new, nactveold, nvar
      double precision, dimension(1:,1:), intent(in) :: ab
      integer i, iclose, ibas
      double precision sum
        iclose=0
        do ibas=1,nactveold
          sum=0.
          do i=1,nvar
            sum=sum+abs(corden(new,i)-corden(ibas,i))/(ab(i,2)-ab(i,1))
          enddo
          if(sum.le.1.d-4) then
            iclose=1
          endif
        enddo
        checkd = iclose ! if good point then 0, else 1
      end function checkd
!------------ end subroutine checkd() ---------------------------
!------------ subroutine expand_grid() ---------------------------
      subroutine expand_grid(alg_type,isupres,nvar,nactve,ngridn, &
        resolve, corden, ab)
        implicit none
! Args
        integer, intent(in) :: alg_type, isupres, nvar
        integer, intent(inout) :: nactve, ngridn
        double precision, intent(in) :: resolve
        double precision, intent(inout), dimension(1:,1:) :: corden
        double precision, intent(in), dimension(1:,1:) :: ab
! Local variable declarations
        integer ntry, nactveold, new, ipoint, ivar, i, iclose
        double precision pcur, del
        IF(ISUPRES .EQ. 0) write(*,*) &
          'Number of active points =', nactve
! Add points near the current solution
        IF(ISUPRES .EQ. 0) &
          write(*,*) 'expanding current grid with new points'
        IF(ISUPRES .EQ. 0) write(*,5200) 100.*resolve
 5200     format(' current grid resolution = ',f8.3, '%')
        new=2*nvar+1
        nactveold=nactve
        do ipoint=1,nactveold
! first, divide current probability into 2*nvar+1 pieces
          pcur=corden(ipoint,nvar+1)/(2*nvar+1)
! update original point
          corden(ipoint,nvar+1)=pcur
          do ivar=1,nvar
	    del=(ab(ivar,2)-ab(ivar,1))*resolve
! create first new trial point at -eps in coordinate ivar
            do i=1,nvar
              corden(nactve+1,i)=corden(ipoint,i)
	    enddo
	    corden(nactve+1,ivar)=corden(nactve+1,ivar)-del
            corden(nactve+1,nvar+1)=pcur
            ntry=nactve+1
! icheck that new point is at least minimally distant from old points
            iclose = checkd(corden,ntry,nactve,nvar,ab)
! only keep trial lower point if it lies above lower bound and satisfies
! minimal distance requirement
            if(corden(nactve+1,ivar).ge.ab(ivar,1)) then
              if(iclose.eq.0) nactve=nactve+1
	    endif
! now create second trial point at +eps in coordinate ivar
            do i=1,nvar
              corden(nactve+1,i)=corden(ipoint,i)
	    enddo
	    corden(nactve+1,ivar)=corden(nactve+1,ivar)+del
            corden(nactve+1,nvar+1)=pcur
! only keep upper point if it lies below upper bound and
! satisfies distance requirement
            ntry=nactve+1
            iclose = checkd(corden,ntry,nactve,nvar,ab)
	    if(corden(nactve+1,ivar).le.ab(ivar,2)) then
	      if(iclose.eq.0) nactve=nactve+1
            endif
          enddo
!    above enddo for loop over ivar=1,nvar
        enddo
!    above enddo for loop over ipoint=1,nactveold
        IF(ISUPRES .EQ. 0) write (*,*) &
          'Number of active grid points after expansion =', nactve
        IF(ISUPRES .EQ. 0) write(*,*)
        ngridn=nactve
! end expansion
      return
    end subroutine expand_grid
!------------ end expand_grid() -------------------------
!  shift10.f                                               11/21/14
!  shift10 has the following changes from shift9:
!  It has the Threadprivate and Save statements to make it compatible
!  with the new npageng28.f program. These statements allow the
!  program to be run in parallel.
!-----------------------------------------------------------------------
!  shift9.f                                                9/28/12
!  shift9 has the following subtle change from shift8:
!  In step 4, the logic to assign the bolus time, BOL(I,IND,1) is
!  simplified in the case where a steady state dose set begins as a
!  time reset event. In this case, the bolus time will be TAU(I) only
!  if both TAU(I) and the bolus value (RR) are not 0. See the reason
!  in the code.
!-----------------------------------------------------------------------
!  shift8.f                                                9/20/12
!  shift8 has changes from shift7 in Step 4 to correct the code in the
!  case where bolus inputs are used in steady state dose sets. In
!  shift.f, a timelag for a bolus which was part of a steady state
!  dose set would not be applied properly. Now it will.
!-----------------------------------------------------------------------
!  shift7.f                                                11/6/11
!  shift7 differs from shift6 as follows:
!  1. The dimensions related to the no. of dose events are changed from
!  500 to 5000. This is needed as shift7 is compiled with idm1x7.f,
!  idm2x7.f, and idm3x7.f (part of the npageng16.f "engine"), which
!  accommodates steady state dose sets.
!  2. 3 lines testing for IF(SIG(IDOSE) .EQ. 0 .AND. IDOSE .GT. 1)
!  are replaced by 	  IF(SIG(IDOSE) .LE. 0 .AND. IDOSE .GT. 1)
!  since now a dose reset occurs when a dose time is 0 (a regular
!  time reset) or < 0 (a time reset occurring with a steady state
!  dose set indicator).
!-----------------------------------------------------------------------
!  SHIFT6.F                                                4/26/11
!  SHIFT5 HAS THE FOLLOWING CHANGES TO SHIFT5:
!  WT AND CCR ARE NO LONGER ASSUMED TO BE SPECIAL COVARIATES IN EACH
!  PATIENT'S WORKING COPY PATIENT DATA FILE. SO ALL DO LOOPS THAT
!  START WITH  DO I = 1, 2+NADD ARE CHANGED TO START WITH DO I = 1,NADD,
!  BUT ONLY IF NADD .GT. 0.
!-----------------------------------------------------------------------
!  SHIFT5.F							9/11/09
!  SHIFT5 HAS THE FOLLOWING CHANGES TO SHIFT4.F.
!  THE ARGUMENT LIST CONTAINS TAU(.) RATHER THAN NTLAG(.). THIS
!  MEANS THAT TAU(I) IS INPUT DIRECTLY AS THE TIMELAG FOR DRUG I.
!  I.E., IT NO LONGER HAS TO BE CALCULATED AS A FUNCTION OF THE
!  PARAMETER ARRAY, P. BECAUSE OF THIS, P IS REMOVED FROM THE ARGUMENT
!  LIST AND THE DIMENSION STATEMENT. ALSO, NTLAG IS REMOVED FROM
!  THT DIMENSION STATEMENT.
!  THE FIRST SET OF ID MODULES TO CALL SHIFT5.F ARE idm1x3.f,
!  idm2x3.f, AND idm3x3.f
!-----------------------------------------------------------------------
!  SHIFT4.FOR							9/1/09
!  SHIFT4 HAS THE FOLLOWING CHANGES FROM SHIFT3:
!  1. NTLAG(I) CAN NOW BE NEGATIVE. IF THIS OCCURS, IT MEANS THAT THE
!  TIMELAG PARAMETER FOR DRUG I WILL BE EXP(P(-NTLAG(I)).
!  2. A BUG IS CORRECTED RELATED TO TIME "RESETS". PREVIOUSLY, IF THE
!  USER HAD A TIME "RESET" IN HIS DOSAGE REGIMEN, THIS ROUTINE WOULD
!  NOT WORK. THE REASON IS THAT IN THE CODE BELOW, EACH NEXT TIME
!  FOR AN IV, COVARIATE, OR BOLUS IS COMPARED TO THE PREVIOUSLY
!  ESTABLISHED TIME IN THE DOSAGE ARRAY (TIMNXT) AND IS A CANDIDATE
!  TO BE THE NEXT TIMNXT IF IT IS .GE. TIMNXT. SO IF A TIME RESET
!  VALUE OF 0 OCCURS, IT WILL NEVER BE A CANDIATE SINCE IT IS NOT
!  .GE. THE LAST TIMNXT. TO FIX THIS, AND MAKE SURE THAT A TIME
!  RESET VALUE OF 0 IS INCLUDED IN THE ADJUSTED DOSAGE BLOCK, THE
!  CODE WILL ADD TO EACH IV, BOLUS, AND COVARIATE ARRAY AN EXTRA
!  LINE WHEN A TIME RESET OCCURS. THIS LINE WILL HAVE A TIME OF
!  1.D19 (I.E., A LARGE VALUE WHICH REPRSENTS INFINITY); AND IT
!  WILL BE FOLLOWED BY A LINE WITH THE ADJUSTED RESET TIME (0 FOR
! AND COVARIATES, AND 0 + TAU(I) FOR BOLI.
!-----------------------------------------------------------------------
!  SHIFT3.FOR							5-23-02
!  SHIFT3 HAS MAJOR CHANGES FROM SHIFT2 TO ALLOW FOR MULTIPLE TIMELAGS,
!  ONE POTENTIALLY FOR EACH BOLUS INPUT OF UP TO NTAU DRUGS.
        SUBROUTINE SHIFT(TAU,ND,SIG,NDRUG,NADD,RS,INTLIST)
! wmy20190318 Moved into npag_utils last week.
!        use npag_utils, only: max_doses,max_RS_J,max_covs,max_input_dim
        IMPLICIT NONE
! Arg list
        double precision, intent(IN), dimension(1:)  :: tau
        integer, intent(INOUT) :: ND
        double precision, intent(INOUT), dimension(1:) :: SIG      ! dimension(max_doses)
        integer, intent(IN) :: NDRUG, NADD
        double precision, intent(INOUT), dimension(1:,1:) :: RS    ! dimension(max_doses,max_RS_J)
        integer, intent(IN), dimension(1:) :: INTLIST           ! dimension(128)
! Local Variables
        integer, dimension(max_input_dim) :: INDIV, INDBOL
        integer, dimension(max_covs) :: INDCOV
        double precision, dimension(max_RS_J) :: TIMCAN
        integer I,IDOSE,IND,ISAME1,ISAME2,ISAME3,NI
        double precision RR,TIMNXT,VALAST
! NEW PARALLEL CODE BELOW AS OF npageng28.f.
        double precision, dimension(max_input_dim,max_doses,2), save :: & ! f90 line continuation
     &    XIV, BOL
        double precision, dimension(max_covs,max_doses,2), save :: COV
!$omp   Threadprivate(XIV,BOL,COV)
!  INPUT ARE:
!  TAU(I) =  THE VALUE OF THE TIMELAG FOR DRUG I.
!  ND = ORIGINAL NO. OF DOSE EVENTS.
!  SIG(I) = TIME FOR ITH DOSE EVENT IN THE ORIGINAL DOSAGE REGIMEN,
!           I=1,ND.
!  NDRUG = NO. OF DRUGS (EACH HAS AN IV, FOLLOWED BY A BOLUS COLUMN).
!  NADD = NO. OF ADDITIONAL COVARIATES (EACH IS IN ITS OWN COLUMN
!         FOLLOWING THE IV/BOLUS COLUMNS.
!  RS(I,J) = "RATE" J FOR THE ITH DOSE EVENT IN THE ORIGINAL DOSAGE
!            REGIMEN; J=1,NI, I=1,ND, WHERE NI = 2*NDRUG + NADD
!            BECAUSE THE "RATES" CONTAIN, IN ORDER, 2 ENTRIES FOR
!            EACH DRUG (1 FOR THE IV AND 1 FOR THE BOLUS) AND 1 EACH
!            FOR THE NADD ADDITIONAL COVARIATES.
!  OUTPUT ARE:
!  ND, SIG, RS, AS ABOVE, EXCEPT FOR THE ALTERED DOSAGE REGIMEN.
!-----------------------------------------------------------------------
!  SHIFT2.FOR							11-16-99
!  SHIFT2 HAS THE FOLLOWING CHANGE FROM SHIFT. AT THE END OF THE
!  FORMATION OF ARRAY XMAT, ALL ROWS WHICH HAVE 0 BOLUS INPUT AND THE
!  SAME OTHER DATA VALUES (EXCEPT TIME) AS THE PREVIOUS ROW ARE NOT
!  USED IN THE NEW ARRAY XMAT2 WHICH HAS ONLY NON-REDUNDANT ROWS.
!  THIS, THEORETICALLY, SHOULDN'T HAVE ANY EFFECT ON CALCULATIONS, BUT
!  NUMERICALLY IT DOES SINCE WHEN THE DVODE ROUTINE SOLVES D.E.'S, IT
!  INTEGRATES OVER DIFFERENT INTERVALS IF EXTRA DOSAGE LINES ARE
!  INCLUDED.
!  EX: TIME   IV   BOLUS	TIME   IV   BOLUS
!       0    100     0		 0    100     0
!       5    100   1000		 2    100   1000
!  NOTE THAT BOTH ABOVE CASES SHOULD GIVE THE SAME RESULTS IF THERE IS
!  A TIME-LAG = 3 IN THE 2ND CASE. BUT, AS THE CODE IS WRITTEN IN
!  SHIFT.FOR, THE 2ND CASE WOULD TRANSLATE TO THE FOLLOWING:
!	 TIME   IV   BOLUS
!         0    100     0
!         2    100     0
!         5    100   1000
!  ... AND THIS WOULD MEAN THAT THE 1ST INTEGRATION BY DVODE WOULD END
!      AT T = 2, RATHER THAN 5 (OR, E.G., 3 IF 3 WAS THE
!      FIRST OBSERVATION TIME). THIS CREATES NUMERICAL DIFFERENCES DUE
!      TO SMALL ROUNDOFF ERRORS WHICH CAN GROW SIGNIFICANTLY.
!-----------------------------------------------------------------------
!  SHIFT.FOR							7-27-99
!  SHIFT.FOR IS A MODULE WHICH INCLUDES SUBROUTINE SHIFT. SHIFT WILL BE
!  CALLED BY ROUTINES OF THE "BIG" NPEM AND IT2B PROGRAMS WHICH HAVE
!  SUBROUTINES FUNC, FUNC1, FUNC2, OR FUNC3 IN THEM.
!  SHIFT INPUTS THE DOSAGE REGIMEN VIA THE INPUT ARGUMENTS (SEE BELOW),
!  AND RETURNS AN ALTERED DOSAGE REGIMEN, WHICH HAS EACH BOLUS INPUT
!  TIME INCREASED BY THE INPUT VALUE OF TAU (THE TIME LAG). NOTE THAT
!  EACH ROW WITH A NON-0 BOLUS INPUT VALUE WILL RESULT IN A NEW ROW IN
!  THE DOSAGE REGIMEN.
!-----------------------------------------------------------------------
!  PROCEDURE FOR THE DOSAGE REGIMEN MODIFICATION:
!  1. ESTABLISH TAU(I) AS THE TIMELAG FOR DRUG I'S BOLUS COLUMN.
!     NO. AS OF SHIFT5.F, THIS VALUE IS INPUT AS AN ARGUMENT.
!  2. ESTABLISH THE IV VALUES AND TIMES INTO XIV(I,J,K). IN PARTICULAR,
!     XIV(I,J,2) IS THE JTH IV VALUE FOR DRUG I, AND XIV(I,J,1) IS THE
!     TIME THIS IV VALUE FIRST OCCURRED. SET THE LAST TIME TO 1.D29 AS
!     AN INDICATOR THAT THERE ARE NO MORE ENTRIES IN THE ARRAY.
!  3. ESTABLISH THE COVARIATE VALUES AND TIMES INTO COV(I,J,K). IN
!     PARTICULAR, COV(I,J,2) IS THE JTH VALUE FOR COVARIATE I, AND
!     COV(I,J,1) IS THE TIME THIS COV VALUE FIRST OCCURRED. SET THE
!     LAST TIME TO 1.D29 AS AN INDICATOR THAT THERE ARE NO MORE ENTRIES
!     IN THE ARRAY.
!  4. ESTABLISH THE BOLUS VALUES AND TIMES INTO BOL(I,J,K).
!     IN PARTICULAR, BOL(I,J,2) IS THE JTH BOLUS VALUE FOR DRUG I, AND
!     BOL(I,J,1) IS THE TIME THIS BOLUS OCCURRED. THE TIMES FOR EACH
!     BOLUS VALUE ARE THOSE ADJUSTED TIMES FROM THE ASSOCIATED TIMELAGS
!     TAU(I),I=1,NDRUG, FROM STEP 1. SET THE LAST TIME TO 1.D29 AS AN
!     INDICATOR THAT THERE ARE NO MORE ENTRIES IN THE ARRAY.
!  5. REASSIGN THE VALUES IN IV, BOL, AND COV TO THE APPROPRIATE ENTRIES
!     OF RS, KEEPING TRACK OF THE RUNNING INDEX, ND, OF DOSE EVENTS. IF
!     ND EXCEEDS 5000, STOP THE PROGRAM WITH A MESSAGE TO THE USER. ALSO
!     REASSIGN THE CORRESPONDING TIME VALUES TO ARRAY SIG.
!  STEP 1.
! TO DO. AS OF SHIFT5.F, TAU(I), I=1,NDRUG, IS INPUT AS
!  AN ARGUMENT TO THIS ROUTINE.
! wmy2018.07.20 -- Changed tau(7) to tau(:)
        I = size(tau)
        IDOSE = size(cov,3)
        IND = size(rs,2)
        if (I /= max_input_dim) then
          write (*,*) "ERR: In SR SHIFT(), size(tlag) \= max_input_dim"
          write (*,*) "size: tau,IDOSE_3, rs_2 =", I, IDOSE, IND
        endif
!  STEP 2:
!  ESTABLISH THE IV VALUES AND TIMES INTO XIV(I,J,K). IN PARTICULAR,
!  XIV(I,J,2) IS THE JTH IV VALUE FOR DRUG I, AND XIV(I,J,1) IS THE
!  TIME THIS IV VALUE FIRST OCCURRED.
	DO I = 1,NDRUG
!  ESTABLISH XIV(I,J,K) FOR DRUG I'S IV. PRESET THE LAST VALUE TO
!  -99 SO THAT THE FIRST VALUE WILL BE DIFFERENT AND THEREFORE ENGAGE
!  THE LOGIC (WHICH ONLY WRITES A ROW INTO THE ARRAY IF THE VALUE IS
!  DIFFERENT THAN THE PREVIOUS VALUE).
!*** MODIFICATION IN SHIFT4.F: IF A TIME RESET OCCURS (I.E., A
!    SIG(IDOSE) = 0, WHERE IDOSE > 1), IT WILL BE HANDLED BY ASSIGNING
!    AN EXTRA TIME VALUE OF 1.D19 (I.E., A LARGE VALUE REPRESENTING
!    TIME = INFINITY) TO THE IV TIME ARRAY. THEN THE REST OF THE
!    THE IV TIME ARRAY WILL BE ESTABLISHED WITH THE REST OF THE VALUES
!    IN SIG, STARTING, OF COURSE, WITH THE TIME RESET VALUE OF 0.
!    THE SAME LOGIC WILL APPLY TO THE COVARIATES AND THE BOLI.
!  NOTE THAT IND WILL BE THE RUNNING INDEX OF THE LATEST ENTRY INTO
!  THE ARRAY. PLACE 1.D29 INTO THE LAST TIME ENTRY OF EACH SUB-ARRAY
!  AS AN INDICATOR THAT THERE ARE NO MORE ENTRIES.
	 XIV(I,1,1) = 1.D29
	 IND = 0
	 VALAST = -99.D0
!  FOR DRUG I, THE IV VALUE IS IN COLUMN 2*I-1 OF ARRAY RS.
	DO IDOSE = 1,ND
	  RR = RS(IDOSE,2*I-1)
!*** MODIFICATION IN SHIFT7.F: A TIME RESET IS NOW DESIGNATED BY A
!  SIG(IDOSE) .LE. 0, RATHER THAN JUST .EQ. 0 (SINCE A STEADY STATE
!  DOSE INDICATOR HAS A NEGATIVE DOSE TIME).
	  IF(SIG(IDOSE) .LE. 0 .AND. IDOSE .GT. 1) THEN
!  THIS REPRESENTS A TIME "RESET". IN THIS CASE, AS INDICATED ABOVE,
!  PUT IN AN EXTRA ROW FOR THE IV REPRESENTING A VERY LARGE TIME
!  AND THE SAME IV VALUE AS THE PREVIOUS VALUE. THEN PUT IN THE
!  LINE REPRESENTING THE RESET TIME OF 0.
	    IND = IND + 1
	    XIV(I,IND,1) = 1.D19
	    XIV(I,IND,2) = XIV(I,IND-1,2)
	    IND = IND + 1
!*** MODIFICATION IN SHIFT7.F. SET THE NEXT XIV(I,IND,1) TO BE
!  SIG(IDOSE), NOT 0, SINCE SIG(IDOSE) MAY BE < 0 (SINCE A STEADY STATE
!  DOSE INDICATOR HAS A NEGATIVE DOSE TIME).
	    XIV(I,IND,1) = SIG(IDOSE)
	    XIV(I,IND,2) = RR
	    XIV(I,IND+1,1) = 1.D29
	    VALAST = RR
	    GO TO 200
	  ENDIF
!  TO GET HERE, THIS DOSE LINE DOES NOT REPRESENT A TIME RESET.
	  IF(RR .NE. VALAST) THEN
         IND = IND + 1
	   XIV(I,IND,1) = SIG(IDOSE)
	   XIV(I,IND,2) = RR
	   XIV(I,IND+1,1) = 1.D29
	   VALAST = RR
	  ENDIF
  200     CONTINUE
	 END DO
!  THE ABOVE END DO IS FOR THE  DO IDOSE = 1,ND  LOOP.
	END DO
!  THE ABOVE END DO IS FOR THE 	DO I = 1,NDRUG  LOOP.
!  STEP 3:
!  ESTABLISH THE COVARIATE VALUES AND TIMES INTO COV(I,J,K). IN
!  PARTICULAR, COV(I,J,2) IS THE JTH VALUE FOR COVARIATE I, AND
!  COV(I,J,1) IS THE TIME THIS COV VALUE FIRST OCCURRED. SET THE
!  LAST TIME TO 1.D29 AS AN INDICATOR THAT THERE ARE NO MORE ENTRIES
!  IN THE ARRAY.
        IF(NADD .GT. 0) THEN
	DO I = 1, NADD
!  ESTABLISH COV(I,J,K) FOR COVARIATE NO. I.
!  PRESET THE LAST VALUE TO -99 SO THAT THE FIRST VALUE WILL BE
!  DIFFERENT AND THEREFORE ENGAGE THE LOGIC (WHICH ONLY WRITES A ROW
!  INTO THE ARRAY IF THE VALUE IS DIFFERENT THAN THE PREVIOUS VALUE).
!  NOTE THAT IND WILL BE THE RUNNING INDEX OF THE LATEST ENTRY INTO THE
!  ARRAY. PLACE 1.D29 INTO THE LAST TIME ENTRY OF EACH SUB-ARRAY AS AN
!  INDICATOR THAT THERE ARE NO MORE ENTRIES.
	 COV(I,1,1) = 1.D29
	 IND = 0
	 VALAST = -99.D0
!  FOR COVARIATE I, THE VALUE IS IN COLUMN 2*NDRUG+I OF ARRAY RS.
	 DO IDOSE = 1,ND
	  RR = RS(IDOSE,2*NDRUG+I)
!*** MODIFICATION IN SHIFT7.F: A TIME RESET IS NOW DESIGNATED BY A
!  SIG(IDOSE) .LE. 0, RATHER THAN JUST .EQ. 0 (SINCE A STEADY STATE
!  DOSE INDICATOR HAS A NEGATIVE DOSE TIME).
	  IF(SIG(IDOSE) .LE. 0 .AND. IDOSE .GT. 1) THEN
!  THIS REPRESENTS A TIME "RESET". IN THIS CASE, AS INDICATED ABOVE,
!  PUT IN AN EXTRA ROW FOR THE COVARIATE REPRESENTING A VERY LARGE TIME
!  AND THE SAME COV VALUE AS THE PREVIOUS VALUE. THEN PUT IN THE
!  LINE REPRESENTING THE RESET TIME OF 0.
	    IND = IND + 1
	    COV(I,IND,1) = 1.D19
	    COV(I,IND,2) = COV(I,IND-1,2)
	    IND = IND + 1
!*** MODIFICATION IN SHIFT7.F. SET THE NEXT COV(I,IND,1) TO BE
!  SIG(IDOSE), NOT 0, SINCE SIG(IDOSE) MAY BE < 0 (SINCE A STEADY STATE
!  DOSE INDICATOR HAS A NEGATIVE DOSE TIME).
	    COV(I,IND,1) = SIG(IDOSE)
	    COV(I,IND,2) = RR
	    COV(I,IND+1,1) = 1.D29
	    VALAST = RR
	    GO TO 300
	  ENDIF
!  TO GET HERE, THIS DOSE LINE DOES NOT REPRESENT A TIME RESET.
	  IF(RR .NE. VALAST) THEN
           IND = IND + 1
	   COV(I,IND,1) = SIG(IDOSE)
	   COV(I,IND,2) = RR
	   COV(I,IND+1,1) = 1.D29
	   VALAST = RR
	  ENDIF
  300     CONTINUE
	 END DO
!  THE ABOVE END DO IS FOR THE   DO IDOSE = 1,ND  LOOP.
	END DO
!  THE ABOVE END DO IS FOR THE  DO I = 1, NADD  LOOP.
        ENDIF
!  THE ABOVE ENDIF IS FOR THE   IF(NADD .GT. 0)  CONDITION.
!  STEP 4:
!  ESTABLISH THE BOLUS VALUES AND TIMES INTO BOL(I,J,K). IN PARTICULAR,
!  BOL(I,J,2) IS THE JTH BOLUS VALUE FOR DRUG I, AND BOL(I,J,1) IS THE
!  ADJUSTED (USING THE ASSOCIATED TIMELAGS TAU(I),I=1,NDRUG) TIME THIS
!  BOLUS OCCURRED.
	DO I = 1,NDRUG
!  ESTABLISH BOL(I,J,K) FOR DRUG I'S BOLUS. EACH ARRAY IS FILLED ONLY
!  WITH NON-0 BOLUS VALUES. NOTE THAT IND WILL BE THE RUNNING INDEX OF
!  THE LATEST ENTRY INTO THE ARRAY. PLACE 1.D29 INTO THE LAST TIME ENTRY
!  OF EACH SUB-ARRAY AS AN INDICATOR THAT THERE ARE NO MORE ENTRIES.
	 BOL(I,1,1) = 1.D29
	 IND = 0
!  FOR DRUG I, THE BOLUS VALUE IS IN COLUMN 2*I OF ARRAY RS.
	 DO IDOSE = 1,ND
	  RR = RS(IDOSE,2*I)
!*** MODIFICATION IN SHIFT7.F: A TIME RESET IS NOW DESIGNATED BY A
!  SIG(IDOSE) .LE. 0, RATHER THAN JUST .EQ. 0 (SINCE A STEADY STATE
!  DOSE INDICATOR HAS A NEGATIVE DOSE TIME).
	  IF(SIG(IDOSE) .LE. 0 .AND. IDOSE .GT. 1) THEN
!  THIS REPRESENTS A TIME "RESET". IN THIS CASE, AS INDICATED ABOVE,
!  PUT IN AN EXTRA ROW FOR THE BOLUS REPRESENTING A VERY LARGE TIME
!  AND AN ACCOMPANYING BOLUS VALUE OF 0. THEN PUT IN THE
!  LINE REPRESENTING THE RESET TIME OF 0 + THE TIMELAG ... IF
!  RR .NE. 0.
	    IND = IND + 1
	    BOL(I,IND,1) = 1.D19
	    BOL(I,IND,2) = 0.D0
	    IND = IND + 1
!*** THE FOLLOWING CODE IS CHANGED IN SHIFT8.F. NOW BOLUS VALUES
!  WORK PROPERLY EVEN WITH TIMELAGS. AND AN ADDITIONAL SUBTLE CHANGE
!  WAS ADDED IN shift9.f (SEE THE COMMENTS AT THE TOP OF shift9.f),
!  AND THE EXTRA COMMENTS BELOW.
!  LOGIC IS NOW AS FOLLOWS:
!  IF SIG(IDOSE) = 0, THIS IS A TIME RESET WHICH IS NOT THE START OF
!     A STEADY STATE DOSE SET. IN THIS CASE, A BOLUS WITH A TIMELAG OF
!     TAU(I) WILL OCCUR AT SIG(IDOSE) + TAU(I) = TAU(I).
!  IF SIG(IDOSE) < 0, THIS IS A TIME RESET WHICH IS THE START OF A
!     STEADY STATE DOSE SET. IN THIS CASE:
!     THE BOLUS TIME WILL BE TAU(I) ONLY IF BOTH TAU(I) AND RR
!     ARE NOT 0. OTHERWISE, IT WILL BE SIG(IDOSE).
!     REASON: IF RR = 0, THERE IS NO BOLUS TO BE GIVEN, SO IT WOULD
!     BE SILLY TO INCLUDE AN EXTRA LINE IN THE DOSAGE REGIMEN WITH
!     A 0 BOLUS (AND IT WOULD VERY SLIGHTLY CHANGE THE RESULTS SINCE
!     THE NUMERICAL INTEGRATION THEN HAS TO INTEGRATE THROUGH AN EXTRA
!     TIME). IN AN EXAMPLE (REMARK 4.b IN NPAG109.EXP, THIS CHANGED THE
!     VALUES IN THE LOG-LIKELIHOODS OUT IN THE 13TH DIGIT, BUT SOME
!     VALUES IN THE DENSITY FILE WERE CHANGED IN THE 4TH DIGIT).
!     ALSO, IF TAU(I) = 0, THE BOLUS HAS NO TIMELAG AND THEREFORE
!     OCCURS AT SIG(IDOSE).
!  THE FOLLOWING EXAMPLE SHOWS WHY A NON-0 BOLUS IN A STEADY STATE DOSE
!  SET, WITH TAU(I) .NE. 0, MUST BE GIVEN AT TAU(I) AND NOT
!  SIG(IDOSE) + TAU(I).
!  EX: IF SIG(IDOSE) = -12, IT MEANS THAT A STEADY STATE DOSE SET IS
!      STARTING WITH AN INTERDOSE INTERVAL OF 12 HOURS. SO, IF A
!      BOLUS WITH A TLAG OF 1.5 HOURS IS GIVEN, ITS TIME MUST BE
!      1.5, NOT -12 + 1.5 = -10.5. REASON: AFTER THE SIG(IDOSE) OF
!      -12 IS CONVERTED IN SUBROUTINE FUNC2 TO 0, THE 1.5 WILL CORRECTLY
!      INDICATE THAT THE BOLUS IS GIVEN 1.5 HOURS AFTER THE START OF THE
!      STEADY STATE DOSE SET. ALSO, A TIME OF -10.5 WOULD COMPLETELY
!      SCREW UP THE FUNC2 LOGIC WHICH WOULD INTERPRET IT AS THE START
!      OF ANOTHER STEADY STATE DOSE SEST.
!      ON THE OTHER HAND, IF A DRUG HAS A TAU(I) = 0, IT CANNOT SHOW
!      UP AS OCCURRING AT TAU(I) = 0 SINCE THIS WILL COMPLETELY SCREW
!      UP FUNC2'S LOGIC, WHICH WILL INTERPRET THE TIME OF 0 AS A
!      TIME RESET EVENT. IN THIS CASE, THE BOLUS OCCURS AT THE START OF
!      THE STEADY STATE DOSE SET, I.E., AT SIG(IDOSE) = -12, WHICH WILL
!      BE CONVERTED TO 0 BY FUNC2).
      CALL THESAME(SIG(IDOSE),0.D0,ISAME1)
      CALL THESAME(TAU(I),0.D0,ISAME2)
      CALL THESAME(RR,0.D0,ISAME3)
      IF(ISAME1 .EQ. 1) BOL(I,IND,1) = TAU(I)
!  NOTE THAT, TECHNICALLY, WE SHOULD SET BOL(I,IND,1) = SIG(IDOSE) = 0
!  IF RR = 0, SINCE THERE IS NO REASON TO HAVE AN EXTRA LINE IN THE
!  DOSAGE REGIMEN FOR A 0 BOLUS ... BUT CHANGING THIS WOULD CHANGE
!  VERY SLIGHTLY THE RESULTS IN A 0 BOLUS CASE SINCE THERE WOULD BE ONE
!  LESS DOSAGE LINE FOR THE NUMERICAL INTEGRATOR TO INTEGRATE THROUGH,
!  SO THE CODE WILL BE LEFT AS IS, FOR CONSISTENCY SAKE.
      IF(ISAME1 .EQ. 0) THEN
       BOL(I,IND,1) = SIG(IDOSE)
       IF(ISAME2 .EQ. 0 .AND. ISAME3 .EQ. 0) BOL(I,IND,1) = TAU(I)
      ENDIF
	    BOL(I,IND,2) = RR
	    BOL(I,IND+1,1) = 1.D29
	    VALAST = RR
	    GO TO 400
	  ENDIF
!  TO GET HERE, THIS DOSE LINE DOES NOT REPRESENT A TIME RESET.
	  IF(RR .NE. 0.D0) THEN
           IND = IND + 1
!  *** CHANGE FOR SHIFT8.F.
!  NOW BOLUS VALUES CAN OCCUR IN STEADY STATE DOSES. AND IF THEY DO,
!  THE FIRST ONE MUST OCCUR AT TIME TAU(I), NOT SIG(IDOSE) + TAU(I)
!  AS THE FOLLOWING EXAMPLE ILLUSTRATES:
!  EX: SIG(1) = -12 INDICATING THAT THE STEADY STATE DOSE SET HAS
!      AN INTERDOSE INTERVAL OF 12 HOURS. TAU(1) = 1.5 -->
!      DRUG 1 HAS A TIMELAG OF 1.5 HOURS. SO, IF THE FIRST BOLUS TIME IS
!      SET =  SIG(1) + TAU(1) = -12 + 1.5 = -10.5, THIS WILL SCREW
!      UP THE FUNC2 LOGIC SINCE IN THAT CODE, THE FIRST TIME OF
!      -12 WILL BE RESET TO BE 0, AND THIS WILL BE FOLLOWED BY -10.5,
!      WHICH WILL LOOK LIKE THE START OF ANOTHER STEADY STATE DOSE
!      SET. INSTEAD, SET FIRST BOLUS TIME = TAU(1) = 1.5, WHICH IS
!      CORRECT SINCE IT OCCURS 1.5 HOURS AFTER THE STEADY STATE DOSE
!      STARTS.
         IF(SIG(IDOSE) .GE. 0.D0) BOL(I,IND,1) = SIG(IDOSE) + TAU(I)
         IF(SIG(IDOSE) .LT. 0.D0) BOL(I,IND,1) = TAU(I)
	   BOL(I,IND,2) = RR
	   BOL(I,IND+1,1) = 1.D29
	  ENDIF
  400     CONTINUE
	 END DO
!  THE ABOVE END DO IS FOR THE  DO IDOSE = 1,ND  LOOP.
	END DO
!  THE ABOVE END DO IS FOR THE  DO I = 1,NDRUG  LOOP.
!  STEP 5:
!  REASSIGN THE VALUES IN IV, BOL, AND COV TO THE APPROPRIATE ENTRIES
!  OF RS, KEEPING TRACK OF THE RUNNING INDEX, ND, OF DOSE EVENTS. IF
!  ND EXCEEDS 5000, STOP THE PROGRAM WITH A MESSAGE TO THE USER. ALSO,
!  REASSIGN THE CORRESPONDING TIME VALUES TO ARRAY SIG.
	NI = 2*NDRUG + NADD
	ND = 0
!  GO THROUGH THE ARRAYS IV, BOL, AND COV TO DETERMINE THE NEXT
!  LOWEST DOSE TIME. PUT THIS VALUE INTO RS, ALONG WITH THE
!  CORRESPONDING VALUES FOR THE IV'S, THE BOLI, AND THE COVARIATES.
!  IN THE LOOP BELOW, IT IS NECESSARY TO KNOW TO WHAT POINT IN THE
!  IV, BOL, AND COV ARRAYS THE TIMES AND VALUES HAVE ALREADY BEEN
!  STORED INTO RS. THESE INDICES ARE INDIV(I), I=1,NDRUG; INDBOL(I),
!  I=1,NDRUG; AND INDCOV(I), I=1,NADD, RESPECTIVELY. E.G.,
!  INDIV(2) = 4 MEANS THAT ALL VALUES IN THE IV, BOL, AND COV ARRAYS,
!  THROUGH THE 4TH TIME FOR IV DRUG 2 (I.E., THROUGH TIME = XIV(2,4,1))
!  HAVE BEEN OR ARE ABOUT TO BE STORED INTO THE RS ARRAY.
!  SO PRESET ALL THESE INDEX INDICATORS = 1, AND INITIALIZE THE
!  CURRENT DOSE TIME TO A NEGATIVE NO. SO THAT THE FIRST TIME
!  THROUGH THE FOLLOWING LOOP WILL ENGAGE THE LOGIC.
	DO I = 1,NDRUG
	 INDIV(I) = 1
	 INDBOL(I) = 1
	END DO
        IF(NADD .GT. 0) THEN
         DO I = 1,NADD
          INDCOV(I) = 1
         END DO
        ENDIF
	TIMNXT = -9999999.D0
  100   CONTINUE
!  FIND THE NEXT LOWEST TIME AMONG THE IV, BOL, AND COV ARRAYS.
!  ESTABLISH INTO TIMCAN(J) THE CANDIDATES FOR THE NEXT DOSE TIME
!  (AND CORRESPONDING VALUES FOR THE IV'S, BOLI, AND COVARIATES) TO
!  BE PUT INTO RS.
        DO I = 1,NDRUG
	 IF(XIV(I,INDIV(I),1) .GT. TIMNXT) TIMCAN(I)=XIV(I,INDIV(I),1)
	 IF(XIV(I,INDIV(I),1) .EQ. TIMNXT) TIMCAN(I)=XIV(I,INDIV(I)+1,1)
	END DO
        DO I = 1,NDRUG
	 IF(BOL(I,INDBOL(I),1) .GT. TIMNXT) TIMCAN(NDRUG+I) =           &
     &    BOL(I,INDBOL(I),1)
	 IF(BOL(I,INDBOL(I),1) .EQ. TIMNXT) TIMCAN(NDRUG+I) =           &
     &    BOL(I,INDBOL(I)+1,1)
	END DO
        IF(NADD .GT. 0) THEN
         DO I = 1,NADD
          IF(COV(I,INDCOV(I),1) .GT. TIMNXT) TIMCAN(2*NDRUG+I) =        &
     &     COV(I,INDCOV(I),1)
          IF(COV(I,INDCOV(I),1) .EQ. TIMNXT) TIMCAN(2*NDRUG+I) =        &
     &     COV(I,INDCOV(I)+1,1)
         END DO
        ENDIF
!  FIND THE NEXT TIMNXT, THE MINIMUM VALUE AMONG THE NI ENTRIES IN
!  TIMCAN. TIMNXT WILL BE THE NEXT TIME TO BE PUT INTO ARRAY RS (ALONG
!  WITH ALL THE CORRESPONDING IV'S, BOLI, AND COVARIATE VALUES). IF
!  TIMNXT = 1.D29, IT IS BECAUSE THERE ARE NO FURTHER VALUES TO BE PUT
!  INTO RS (I.E, THE PROCESS IS FINISHED).
	TIMNXT = TIMCAN(1)
	DO I = 2,NI
	 IF(TIMCAN(I) .LT. TIMNXT) TIMNXT = TIMCAN(I)
	END DO
	IF(TIMNXT .EQ. 1.D29) RETURN
!  SINCE TIMNXT < 1.D29, THERE ARE MORE VALUES TO BE PUT INTO RS.
!  GO THROUGH ALL THE SUBARRAYS AND PUT IN VALUES AS FOLLOWS. IF THE
!  CURRENT TIME FOR AN IV, BOLUS, OR COVARIATE IS THE SAME AS TIMNXT,
!  PUT THE CORRESPONDING IV, BOLUS, OR COVARIATE VALUE INTO RS, AND
!  INCREASE THE INDEX FOR THAT SUB-ARRAY TO THE NEXT VALUE. IF THE
!  CURRENT TIME FOR AN IV OR A COVARIATE IS .GT. TIMNXT, PUT THE IV OR
!  COVARIATE VALUE FROM THE PREVIOUS ROW INTO RS, AND LEAVE THE INDEX
!  UNCHANGED. IF THE CURRENT TIME FOR A BOLUS IS .GT. TIMNXT, PUT 0.0
!  INTO RS (I.E., BOLUS VALUES ARE INSTANTANEOUS, WHEREAS IV AND
!  COVARIATE VALUES CONTINUE UNTIL CHANGED), AND LEAVE THE INDEX
!  UNCHANGED.
!  TEST FOR TIMNXT = 1.D19, WHICH INDICATES A TIME RESET.
	IF(TIMNXT .EQ. 1.D19) THEN
!  TIMNXT = 1.D19 MEANS THAT THE NEXT TIME IN EACH ARRAY IS THE
!  TIME AT OR AFTER THE RESET. SO INCRASE ALL THE ARRAY INDICES BY
!  1, RESET TIMNXT TO A NEGATIVE NO. AND RETURN TO LABEL 100.
       DO I = 1,NDRUG
	  INDIV(I) = INDIV(I) + 1
	  INDBOL(I) = INDBOL(I) + 1
	 END DO
        IF(NADD .GT. 0) THEN
         DO I = 1,NADD
          INDCOV(I) = INDCOV(I) + 1
         END DO
        ENDIF
	 TIMNXT = -9999999.D0
	 GO TO 100
	ENDIF
	ND = ND+1
	IF(ND .GT. max_doses) THEN
!  IF ND > max_doses, STOP WITH A MESSAGE TO THE USER THAT THE
!  PROGRAM ONLY ALLOWS A TOTAL OF 5000 DOSE EVENTS.
   10	 WRITE(*,1) ND
    1    FORMAT(/' THE NUMBER OF DOSE EVENTS, AFTER TAKING INTO'/       &
     &' ACCOUNT DIFFERING TIMES DUE TO TIMELAGS IS ',I6,', MORE THAN'/  &
     &' THE ALLOWABLE MAXIMUM OF 5000. THE PROGRAM IS STOPPING. PLEASE'/&
     &' RERUN WITH PATIENTS HAVING FEWER DOSE EVENTS, OR WITH FEWER'/   &
     &' TIMELAG VALUES SELECTED AS FIXED OR RANDOM PARAMETERS.'//)
	 STOP
	ENDIF
!  ND .LE. max_doses, SO CONTINUE. FOR THIS DOSE EVENT, PUT IN THE CURRENT
!  TIME, AND THE CORRESPONDING IV, BOLUS, AND COVARIATE VALUES.
	SIG(ND) = TIMNXT
        DO I = 1,NDRUG
	 IF(TIMNXT .LT. XIV(I,INDIV(I),1)) THEN
	  RS(ND,2*I-1) = RS(ND-1,2*I-1)
	 ENDIF
	 IF(TIMNXT .EQ. XIV(I,INDIV(I),1)) THEN
	  RS(ND,2*I-1) = XIV(I,INDIV(I),2)
	  INDIV(I) = INDIV(I) + 1
	 ENDIF
	 IF(TIMNXT .LT. BOL(I,INDBOL(I),1)) THEN
	  RS(ND,2*I) = 0.D0
	 ENDIF
	 IF(TIMNXT .EQ. BOL(I,INDBOL(I),1)) THEN
	  RS(ND,2*I) = BOL(I,INDBOL(I),2)
	  INDBOL(I) = INDBOL(I) + 1
	 ENDIF
	END DO
        IF(NADD .GT. 0) THEN
         DO I = 1,NADD
          IF(TIMNXT .LT. COV(I,INDCOV(I),1))                            &
     &     RS(ND,2*NDRUG+I) = RS(ND-1,2*NDRUG+I)
          IF(TIMNXT .EQ. COV(I,INDCOV(I),1)) THEN
           RS(ND,2*NDRUG+I) = COV(I,INDCOV(I),2)
           INDCOV(I) = INDCOV(I) + 1
          ENDIF
         END DO
        ENDIF
	GO TO 100
        return
        END subroutine shift
! ##################################################################### !
!
! wmy20190315 Translated to F90 and Moved into npag_utils.f90.
!
!  All numbers written out in F or G format are now tested to see if
!  they are inside [-1.D-99, 1.D-99]. If so, they are changed to be
!  0. This makes text I/O consistent across machines and for small or
!  big numbers.
      SUBROUTINE VERIFYVAL(N,X)
      IMPLICIT NONE
      integer, intent(in) :: N
      double precision, dimension(1:), intent(inout) :: X
!  THIS ROUTINE INPUTS X(I),I=1,N.
!  ON OUTPUT, EACH X(.) WHICH IS INSIDE [-1.D-99, 1.D-99] IS REPLACED
!  BY 0. THIS PREVENTS THIS VALUE FROM BEING WRITTEN OUT IMPROPERLY,
!  E.G., AS .934-106, RATHER THAN .934E-106.
!  ANY X(.) VALUE NOT INSIDE THE ABOVE RANGE WILL BE UNCHANGED ON
!  OUTPUT.
      integer I
      DO I = 1,N
       IF(X(I) .GE. -1.D-99 .AND. X(I) .LE. 1.D-99) X(I) = 0.D0
      END DO
      RETURN
      END subroutine verifyval
! ##################################################################### !
      SUBROUTINE ORDERDELTA(NDRUG,DELTAIV,NDELTA,ORDELT)
!
! wmy20190318 Validated on NPrun(30,...), SIMrun(1,...), and NPrun(1,...)
!
!      IMPLICIT REAL*8(A-H,O-Z)
!
! wmy20190316
      implicit none
!      DIMENSION DELTAIV(7),ORDELT(7),X(7)
!
! wmy20190317 In arrays above, 7 = max_input_dim.  In future, this will
!   be larger; but no error checking on length is done.
      integer, intent(in) ::  NDRUG
      double precision, intent(in), dimension(1:) :: DELTAIV
      integer, intent(inout) :: NDELTA
      double precision, intent(inout), dimension(1:) :: ORDELT
!
! wmy20190317 -- I find it difficult to believe there is not a (better) canned
!  fortran routine that orders arrays; but for now, use this subroutine.
!
!  SUBROUTINE ORDERDELTA IS CALLED BY NEWWORK1 TO OBTAIN NDELTA, THE NO.
!  OF UNIQUE NON-0 VALUES IN THE DELTAIV(.) ARRAY. THEN THE ORDERED SET
!  OF THESE NDELTA VALUES IS PUT INTO ORDELT(.). NOTE THAT
!  NDELTA WILL BE 0 IF ALL THE PARTICIPATING DRUGS ARE BOLUSES SINCE
!  THEY WOULDN'T NEED AN ENDING TIME THEN.
      integer IDRUG, IDRUGNEW, ICOMP, ISAME
      double precision X(NDRUG), VALUE
!  FIRST STORE ALL THE VALUES IN DELTAIV INTO X SO THAT DELTAIV WILL
!  NOT BE CHANGED.
      DO IDRUG = 1,NDRUG
       X(IDRUG) = DELTAIV(IDRUG)
      END DO
!  THE LOGIC OF THIS ROUTINE IS BASED ON \PERSONAL\FINANCE\ORDER.FOR.
!  TO DO THIS, EACH VALUE IN X(.) WILL BE COMPARED TO THE
!  PREVIOUS ONE. IF IT IS < THE PREVIOUS ONE, THE VALUE WILL EXCHANGE
!  PLACES WITH THE PREVIOUS ONE, AND THE TESTING WILL CONTINUE. THE
!  TESTING WILL STOP FOR A VALUE WHEN IT IS COMPARED TO A PREVIOUS
!  VALUE WHICH IS .LE. ITS VALUE.
      DO IDRUG = 2, NDRUG
!  COMPARE VALUE FOR IDRUG WITH EACH PREVIOUS VALUE, AND HAVE IT
!  EXCHANGE PLACES WITH THAT VALUE, UNTIL IT REACHES ONE WHICH HAS A
!  SMALLER VALUE. FIRST SET IDRUGNEW = IDRUG; AFTER THE FOLLOWING
!  CODE, IDRUGNEW WILL BE THE INDEX NO. FOR VALUE AT THE OLD IDRUG
!  POSITION.
       IDRUGNEW = IDRUG
       ICOMP = IDRUG
  110  ICOMP = ICOMP - 1
!  NOW COMPARE VALUE IN LOCATION ICOMP WITH THE VALUE IN LOCATION
!  IDRUGNEW. IF THE LATTER IS .LT. THE FORMER, INTERCHANGE THE RECORDS.
       IF(X(IDRUGNEW) .LT. X(ICOMP)) THEN
        VALUE = X(IDRUGNEW)
        X(IDRUGNEW) = X(ICOMP)
        X(ICOMP) = VALUE
        IDRUGNEW = ICOMP
!  IF IDRUGNEW = 1, IT HAS BEEN CHECKED AGAINST ALL RECORDS (AND IS
!  THE SMALLEST VALUE); IF IS IS > 1, CONTINUE THE PROCESS.
        IF(IDRUGNEW .EQ. 1) GO TO 150
        IF(IDRUGNEW .GT. 1) GO TO 110
       ENDIF
!  THE ABOVE ENDIF IS FOR THE
!   IF(X(IDRUGNEW) .LT. X(ICOMP))  CONDITION.
  150 END DO
!  THE ABOVE END DO IS FOR THE  DO IDRUG = 2, NDRUG LOOP.
!  NOW THE NDRUG VALUES ARE ORDERED, FROM SMALL TO LARGE IN X.
!  REWRITE THEM INTO ORDELT, BUT PUT ONLY THE NON-0 AND
!  UNIQUE VALUES INTO ORDELT, AND KEEP TRACK OF NOW MANY OF THESE
!  UNIQUE NON O VALUES THERE ARE - IT WILL BE NDELTA AT THE END OF
!  THE FOLLOWING LOOP.
      NDELTA = 0
      DO IDRUG = 1,NDRUG
!  FOR THIS VALUE TO BE COUNTED, IT CANNOT = THE PREVIOUS VALUE, AND
!  IT CANNOT = 0.
       IF(IDRUG .EQ. 1 .AND. X(IDRUG) .GT. 0) THEN
        NDELTA = NDELTA + 1
        ORDELT(NDELTA) = X(IDRUG)
       ENDIF
       IF(IDRUG .GE. 2) THEN
        CALL THESAME(X(IDRUG),X(IDRUG-1),ISAME)
        IF(ISAME .EQ. 0) THEN
         NDELTA = NDELTA + 1
         ORDELT(NDELTA) = X(IDRUG)
        ENDIF
       ENDIF
      END DO
!  THE ABOVE END DO IS FOR THE  DO IDRUG = 1,NDRUG  LOOP.
      RETURN
      END subroutine orderdelta
! ##################################################################### !
      SUBROUTINE THESAME(X1,X2,ISAME)
!
! wmy20190318 Validated on Pmetrics examples NPrun(30,...), NPrun(1,...),
!  and SIMrun(1,...).
!
!        IMPLICIT REAL*8(A-H,O-Z)
      implicit none
      double precision, intent(in) ::  X1, X2
      integer, intent(inout) :: ISAME
! wmy20190317  This subroutine should be a function.  The point is to
!  ensure consistency, So the threshold should be a common parameter.
!
!  THIS ROUTINE CHECKS TO SEE IF X1 AND X2 ARE VIRTUALLY THE SAME
!  VALUES (I.E., IF THEY ARE WITHIN 1.D-10 OF EACH OTHER). IF SO,
!  ISAME RETURNS AS 1; IF NOT ISAME RETURNS AS 0.
      double precision XDEL
        ISAME = 0
        XDEL = DABS(X1-X2)
        IF(XDEL .LE. 1.D-10) ISAME = 1
        RETURN
        END subroutine thesame
! ##################################################################### !
      SUBROUTINE PREDLAST3(NN,NSET,XSTORE,XPRED,ICONV)
!      use npag_utils, only: thesame
!      IMPLICIT REAL*8(A-H,O-Z)
       implicit none
!      DIMENSION XSTORE(100,20),XPRED(20),COMP(5,20)
      integer, intent(in) ::  NN, NSET
      double precision, dimension(1:,1:), intent(in) :: XSTORE ! dimension(max_SS_doses,max_ODE_comps)
      double precision, dimension(1:), intent(inout) :: XPRED  ! dimension(max_ODE_comps)
      integer, intent(inout) :: ICONV
!
! TODO
! 1) NN is the number of compartments used in current model; verify that
!    NN .le. max_ODE_comps, or (even more cautious) verify size of
!    appropriate diemnsion of input arrays is .le. max_ODE_comps
! 2) Likewise, verify that NSET .le. max_SS_doses
! 3) Exit w/error if the above inequalities not true; or procede w/no more
!    error checking needed.
! wmy20190814 4) Finally, for either approach above, set dimensions above to assumed shape.
!
! Locals
      double precision, dimension(5,max_ODE_comps) :: COMP
      double precision TOL1, TOL2, A1, A2, A3, DEL1, DEL2
      double precision F, PRED1, PRED2, PRED3, DEN, PREDNEG
      integer I, II, J, IN
      integer ISAME1, ISAMEF1, ISAME2, ISAMEF2, ISAME3, ISAMEF3
      integer ISAMETOT, ISAMEFTOT, ISAMEDEN, ISAME
!  NOTE THAT AS OF MONT109.FOR, THE DIMENSIONS OF 6 IN XSTORE, XPRED,
!  AND COMP HAVE BEEN CHANGED TO 20, WHICH IS WHAT THEY SHOULD HAVE BEEN
!  ALL ALONG (SEE SUBROUTINE FUNC2). Also, the dimension 5 in COMP is
!  a hardcoded computational parameter: PREDLAST3 is called after every
!  SS dose _after_ the 5th dose.
!  THIS SUBROUTINE IS CALLED BY SUBROUTINE FUNC WITH NSET SETS OF NN
!  COMPARTMENT  VALUES IN XSTORE. USE THE LAST 5 SETS OF VALUES TO
!  PREDICT THE FINAL (STEADY STATE) COMPARTMENT AMOUNTS AFTER THE LAST
!  (100TH) DOSE SET. NOTE THAT AS OF MONT110.FOR, THERE WILL ALSO BE AN
!  ADDITIONAL 101ST DOSE SET APPLIED. BUT THE PREDICTION SHOULD STILL
!  BE AFTER THE 100TH SET.
!  IF THESE VALUES "CONVERGE", SET ICONV = 1, AND WRITE THE PREDICTED
!  VALUES INTO XPRED. IF THEY DON'T CONVERGE, SET ICONV = 0.
!  TOL1 AND TOL2 ARE, FOR NOW, HARDCODED TO BE .0005.
        TOL1 = .0005D0
        TOL2 = .0005D0
!  THE LAST 5 SETS OF VALUES ARE IN XSTORE(NSET-4:NSET,.). PUT THESE
!  VALUES INTO COMP(.,.).
      II = 0
      DO I = NSET-4,NSET
       II = II+1
       DO J = 1,NN
        COMP(II,J) = XSTORE(I,J)
       END DO
      END DO
!  FOR EACH COMPARTMENT AMOUNT, SEE IF THE FINAL STEADY STATE COMP.
!  AMOUNT CAN BE PREDICTED ACCURATELY.
      DO IN = 1,NN
       A1 = COMP(1,IN)
       A2 = COMP(2,IN)
       A3 = COMP(3,IN)
       DEL1 = A2 - A1
       DEL2 = A3 - A2
!  TEST FOR DEL1 = 0. IF SO, SEE ISAMETOT BELOW.
       CALL THESAME(DEL1,0.D0,ISAME1)
       IF(ISAME1 .EQ. 0) THEN
        F = DEL2/DEL1
!  THE UNDERLYING ASSUMPTION IS THAT THE RATIO F = DEL2/DEL1
!  IS CONTANT BETWEEN CONSECUTIVE OUTPUT DIFFERENCES (At SS) . IF SO, THEN
!  THE STEADY STATE VALUE WILL BE A1 + DEL1/(1 - F) (SEE SS.EXP
!  IN \ALAN3\STEADYSTATE). CALCULATE THIS VALUE AND CALL IT PRED1.
!  BUT, IF DEL2 = DEL1, THEN F = 1. IN THIS CASE, CAN'T DO THE FOLLOWING
!  CALCULATION FOR PRED1, AND WE WOULDN'T WANT TO DO IT SINCE
!  DEL2 = DEL1 --> A2 - A1 = A3 - A2 --> A1, A2, AND A3 ARE IN AN
!  ARITHMETIC PROGRESSION --> THERE OBVIOUSLY CAN BE NO CONVERGENCE
!  SINCE, AFTER 100 DOSES, THE VALUE WOULD JUST A1 + 99*DEL1 ...
!  UNLESS DEL1 = 0, IN WHICH CASE THE VALUE WOULD CONVERGE TO A1.
!  IN THIS CASE SET ISAMEF1 = 1, AND SKIP CALC. OF PRED1. AND THEN
!  SEE THE LOGIC RELATED TO ISAMEF1 BELOW.
        CALL THESAME(F,1.D0,ISAMEF1)
        IF(ISAMEF1 .EQ. 0) PRED1 = A1 + DEL1/(1.D0 - F)
       ENDIF
!  SIMILARLY, CALCULATE PRED2 (BASED ON (A2,A3,A4)) AND PRED3 (BASED
!  ON (A3,A4,A5).
       A1 = COMP(2,IN)
       A2 = COMP(3,IN)
       A3 = COMP(4,IN)
       DEL1 = A2 - A1
       DEL2 = A3 - A2
!  TEST FOR DEL1 = 0. IF SO, SEE ISAMETOT BELOW.
       CALL THESAME(DEL1,0.D0,ISAME2)
       IF(ISAME2 .EQ. 0) THEN
        F = DEL2/DEL1
        CALL THESAME(F,1.D0,ISAMEF2)
        IF(ISAMEF2 .EQ. 0) PRED2 = A1 + DEL1/(1.D0 - F)
       ENDIF
       A1 = COMP(3,IN)
       A2 = COMP(4,IN)
       A3 = COMP(5,IN)
       DEL1 = A2 - A1
       DEL2 = A3 - A2
!  TEST FOR DEL1 = 0. IF SO, SEE ISAMETOT BELOW.
       CALL THESAME(DEL1,0.D0,ISAME3)
       IF(ISAME3 .EQ. 0) THEN
        F = DEL2/DEL1
        CALL THESAME(F,1.D0,ISAMEF3)
        IF(ISAMEF3 .EQ. 0) PRED3 = A1 + DEL1/(1.D0 - F)
       ENDIF
!  ASSUMING A NEGATIVE EXPONENTIAL PATTERN FIT (SEE SS.EXP IN
!  \ALAN3\STEADYSTATE OR HOME NOTES, PG.2 ON 9/11/11 FOR DETAILS) ON
! (PRED1,PRED2,PRED3), CALCULATE PREDNEG.
!  BUT ONLY DO THIS CALCULATION, AND THE SUBSEQUENT
!  CONVERGENCE DETERMINATION IF ISAME1 = ISAME2 = ISAME3 = 0, AND
!  ISAMEF1 = ISAMEF2 = ISAMEF3 = 0. OTHERWISE, AT LEAST ONE OF THE
!  PREDICTED VALUES ABOVE WAS NOT CALCULATED.
       ISAMETOT = ISAME1 + ISAME2 + ISAME3
       ISAMEFTOT = ISAMEF1 + ISAMEF2 + ISAMEF3
       IF(ISAMETOT .EQ. 0 .AND. ISAMEFTOT .EQ. 0) THEN
! EDITED CODE BELOW FOR MONT103.FOR.
!  IF PRED1 + PRED3 - 2*PRED2 = 0, PREDNEG (SEE BELOW) CANNOT BE
!  CALCULATED. IN THIS CASE, PRED2 - PRED1 = PRED3 - PRED2 -->
!  THE SEQUENCE (PRED1, PRED2, PRED3) IS LINEAR, WHICH CANNOT BE
!  MODELED WITH AN EXPONENTIAL FIT (SEE COMMENTS ABOVE). SO, IF THIS
!  HAPPENS, CONVERGENCE WILL BE SATISFIED IF THESE 3 VALUES ARE
!  VIRTUALLY THE SAME - I.E., ONLY THE REQUIREMENT INVOLVING TOL1
!  WILL BE NEEDED FOR CONVERGENCE (RECALL THE ONLY REASON FOR THE
!  EXTRA NEGATIVE EXPONENTIAL FIT, AND THE CALCULATION OF PREDNEG IS FOR
!  THOSE CASES WHERE PRED1, PRED2, AND PRED3 ARE NOT ALL VIRTUALLY THE
!  SAME VALUE).
        DEN = PRED1+PRED3-2.D0*PRED2
        CALL THESAME(DEN,0.D0,ISAMEDEN)
        IF(ISAMEDEN .EQ. 0) &
         PREDNEG = (PRED1*PRED3 - PRED2*PRED2)/DEN
!  NOW CHECK FOR CONVERGENCE, WHICH HAS BEEN OBTAINED IF
!  |PRED3/PRED2 - 1| < TOL1 AND |PREDNEG/PRED3 - 1| < TOL2.
        ICONV = 1
        IF(DABS(PRED3/PRED2 - 1.D0) .GE. TOL1) ICONV = 0
        IF(ISAMEDEN .EQ. 0 .AND. DABS(PREDNEG/PRED3 - 1.D0) .GE. TOL2) &
         ICONV = 0
!  IF ICONV = 1 FOR THIS COMPARTMENT, IN, STORE THE PREDICTED AMOUNT,
!  AND CONTINUE TO THE NEXT COMPARTMENT. NOTE BELOW THAT
!  NON-CONVERGENCE IN ANY COMPARTMENT ENDS THE PROCESS SINCE TO
!  CONVERGE, ALL COMPARTMENT PREDICTIONS MUST CONVERGE.
        IF(ICONV .EQ. 1 .AND. ISAMEDEN .EQ. 1) XPRED(IN) = PRED3
        IF(ICONV .EQ. 1 .AND. ISAMEDEN .EQ. 0) XPRED(IN) = PREDNEG
! EDITED CODE ABOVE FOR MONT103.FOR.
       ENDIF
!  THE ABOVE ENDIF IS FOR THE  IF(ISAMETOT .EQ. 0 .AND. ISAMEFTOT .EQ.0)
!  CONDITION.
!  IF ISAMETOT .GT. 0, THERE ARE TWO POSSIBILITIES (AND NOTE THAT IT
!  DOSEN'T MATTER WHAT ISAMEFTOT IS IN THIS CASE):
!   ISAMETOT = 3, IN WHICH CASE COMP(1:4,IN) ARE ALL THE SAME.
!   ISAMETOT = 1 OR 2, IN WHICH CASE SOME OF THE COMP(1:4,IN) VALUES
!     ARE THE SAME, AND SOME ARE NOT.
!  IN THE FORMER CASE, VERIFY THAT COMP(5,IN) IS THE SAME VALUE AS
!  THE COMP(1:4,IN). IF SO, SET THE PREDICTED VALUE = THIS VALUE
!  (I.E., THE PREDICTED VALUE FOR A CONSTANT FUNCTION IS THE
!  CONSTANT VALUE), AND SET ICONV = 1. OTHERWISE, SET ICONV = 0
!  SINCE THERE IS NO WAY TO FIT 4 VALUES WHICH ARE THE SAME AND ONE
!  WHICH IS NOT USING A NEGATIVE EXPONENTIAL FUNCTION.
!  IN THE LATTER CASE, SINCE SOME OF THE COMP(1:4,IN) VALUES ARE THE
!  SAME, AND SOME ARE NOT, SET ICONV = 0 FOR THE SAME REASON AS
!  STATED IN THE PREVIOUS PARAGRAPH.
       IF(ISAMETOT .EQ. 3) THEN
        CALL THESAME(COMP(5,IN),COMP(1,IN),ISAME)
        IF(ISAME .EQ. 1) THEN
         ICONV = 1
         XPRED(IN) = COMP(1,IN)
        ENDIF
        IF(ISAME .EQ. 0) ICONV = 0
       ENDIF
       IF(ISAMETOT .EQ. 1 .OR. ISAMETOT .EQ. 2) ICONV = 0
!  IF ICONV = 0, CONVERGENCE WAS NOT ACHIEVED.
       IF(ICONV .EQ. 0) RETURN
      END DO
!  THE ABOVE END DO IS FOR THE  DO IN = 1,NN  LOOP.
!  TO GET TO THIS POINT, ALL COMPARTMENT AMOUNTS HAVE CONVERGED, AND
!  THEIR PREDICTED AMOUNTS HAVE BEEN STORED INTO XPRED(IN),IN=1,NN.
      RETURN
      END subroutine PREDLAST3
! ##################################################################### !
      subroutine practice_c_call(NP,PVALS)
!
! 20200730 This is an F90 envelope around a void C-function, that can
!   be called from the F77 routines in Pmetrics, in particular, IT2B.
!   This routine takes in the NP-length double vector PVALS, and copies
!   PVALS to a local allocatable array, ArgForPkSim. ArgsForPkSim is
!   passed to a C function, void r8_abs(int N, double *x), in c_utils.c
!   *x is operated on, and "returned".
        use, intrinsic :: iso_c_binding
        implicit none
        !  TIMESTAMP is provided by the C library, and so the following
        !  INTERFACE block must set up the interface.
        interface
          subroutine timestamp ( ) bind ( c )
            use iso_c_binding
          end subroutine timestamp
        end interface
        interface
          subroutine r8_abs ( n, x ) bind ( c )
            use iso_c_binding
            integer ( c_int ) , VALUE :: n
            real (  c_double ) :: x( * )
          end subroutine r8_abs
        end interface
        integer NP
        double precision, dimension ( 1: ) :: PVALS
        integer iii
        double precision, dimension(1) :: r8sample
        double precision, dimension ( : ), allocatable :: ArgsForPkSim
        allocate ( ArgsForPkSim ( NP ) )
        call timestamp ( )
        r8sample( 1 ) = -3.14159d12
        call r8_abs ( 1, r8sample )
        write (*,*) "absolute value of r8sample is", r8sample( 1 )
        do iii = 1,NP
           ArgsForPkSim ( iii ) = PVALS ( iii )
        end do
        call r8_abs ( NP, ArgsForPkSim )
        write (*,*)  "ret. fr. r8_abs", ( ArgsForPkSim ( iii ), iii=1,NP )
        write (*,*) "PVALS^0.5", (sqrt( PVALS ( iii )), iii = 1,NP)
        ! 20200730 -- the above writes match to the 16 digit precision printout
        deallocate ( ArgsForPkSim )
      end subroutine
! ##################################################################### !
! end module ! npag_utils.f90
! MODULE random
! A module for random number generation from the following distributions:
!
!     Distribution                    Function/subroutine name
!
!     Normal (Gaussian)               random_normal
!     Gamma                           random_gamma
!     Chi-squared                     random_chisq
!     Exponential                     random_exponential
!     Weibull                         random_Weibull
!     Beta                            random_beta
!     t                               random_t
!     Multivariate normal             random_mvnorm
!     Generalized inverse Gaussian    random_inv_gauss
!     Poisson                         random_Poisson
!     Binomial                        random_binomial1   *
!                                     random_binomial2   *
!     Negative binomial               random_neg_binomial
!     von Mises                       random_von_Mises
!     Cauchy                          random_Cauchy
!
!  Generate a random ordering of the integers 1 .. N
!                                     random_order
!     Initialize (seed) the uniform random number generator for ANY compiler
!                                     seed_random_number
!     Lognormal - see note below.
!  ** Two functions are provided for the binomial distribution.
!  If the parameter values remain constant, it is recommended that the
!  first function is used (random_binomial1).   If one or both of the
!  parameters change, use the second function (random_binomial2).
! The compilers own random number generator, SUBROUTINE RANDOM_NUMBER(r),
! is used to provide a source of uniformly distributed random numbers.
! N.B. At this stage, only one random number is generated at each call to
!      one of the functions above.
! The module uses the following functions which are included here:
! bin_prob to calculate a single binomial probability
! lngamma  to calculate the logarithm to base e of the gamma function
! Some of the code is adapted from Dagpunar's book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
!
! In most of Dagpunar's routines, there is a test to see whether the value
! of one or two floating-point parameters has changed since the last call.
! These tests have been replaced by using a logical variable FIRST.
! This should be set to .TRUE. on the first call using new values of the
! parameters, and .FALSE. if the parameter values are the same as for the
! previous call.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lognormal distribution
! If X has a lognormal distribution, then log(X) is normally distributed.
! Here the logarithm is the natural logarithm, that is to base e, sometimes
! denoted as ln.  To generate random variates from this distribution, generate
! a random deviate from the normal distribution with mean and variance equal
! to the mean and variance of the logarithms of X, then take its exponential.
! Relationship between the mean & variance of log(X) and the mean & variance
! of X, when X has a lognormal distribution.
! Let m = mean of log(X), and s^2 = variance of log(X)
! Then
! mean of X     = exp(m + 0.5s^2)
! variance of X = (mean(X))^2.[exp(s^2) - 1]
! In the reverse direction (rarely used)
! variance of log(X) = log[1 + var(X)/(mean(X))^2]
! mean of log(X)     = log(mean(X) - 0.5var(log(X))
! N.B. The above formulae relate to population parameters; they will only be
!      approximate if applied to sample values.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Version 1.13, 2 October 2000
! Changes from version 1.01
! 1. The random_order, random_Poisson & random_binomial routines have been
!    replaced with more efficient routines.
! 2. A routine, seed_random_number, has been added to seed the uniform random
!    number generator.   This requires input of the required number of seeds
!    for the particular compiler from a specified I/O unit such as a keyboard.
! 3. Made compatible with Lahey's ELF90.
! 4. Marsaglia & Tsang algorithm used for random_gamma when shape parameter > 1.
! 5. INTENT for array f corrected in random_mvnorm.
!     Author: Alan Miller
!     e-mail: amiller @ bigpond.net.au
! MODULE random is catted into MODULE npag_utils -----------------------
! IMPLICIT NONE
! REAL, PRIVATE      :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0,   &
!                       vsmall = TINY(1.0), vlarge = HUGE(1.0)
! PRIVATE            :: integral
! INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)
!
! CONTAINS
! MODULE random is catted into MODULE npag_utils -----------------------
FUNCTION random_normal() RESULT(fn_val)
! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.
!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.
REAL :: fn_val
!     Local variables
REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,    &
            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
!     Generate P = (u,v) uniform in rectangle enclosing acceptance region
DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156 * (v - half)
!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)
!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO
!     Return ratio of P's coordinates as the normal deviate
fn_val = v/u
RETURN
END FUNCTION random_normal
FUNCTION random_gamma(s, first) RESULT(fn_val)
! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
!     FUNCTION GENERATES A RANDOM GAMMA VARIATE.
!     CALLS EITHER random_gamma1 (S > 1.0)
!     OR random_exponential (S = 1.0)
!     OR random_gamma2 (S < 1.0).
!     S = SHAPE PARAMETER OF DISTRIBUTION (0 < REAL).
REAL, INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val
IF (s <= zero) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE MUST BE POSITIVE'
  STOP
END IF
IF (s > one) THEN
  fn_val = random_gamma1(s, first)
ELSE IF (s < one) THEN
  fn_val = random_gamma2(s, first)
ELSE
  fn_val = random_exponential()
END IF
RETURN
END FUNCTION random_gamma
FUNCTION random_gamma1(s, first) RESULT(fn_val)
! Uses the algorithm in
! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
! gamma variables', Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.
! Generates a random gamma deviate for shape parameter s >= 1.
REAL, INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val
! Local variables
REAL, SAVE  :: c, d
REAL        :: u, v, x
IF (first) THEN
  d = s - one/3.
  c = one/SQRT(9.0*d)
END IF
! Start of main loop
DO
! Generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.
  DO
    x = random_normal()
    v = (one + c*x)**3
    IF (v > zero) EXIT
  END DO
! Generate uniform variable U
  CALL RANDOM_NUMBER(u)
  IF (u < one - 0.0331*x**4) THEN
    fn_val = d*v
    EXIT
  ELSE IF (LOG(u) < half*x**2 + d*(one - v + LOG(v))) THEN
    fn_val = d*v
    EXIT
  END IF
END DO
RETURN
END FUNCTION random_gamma1
FUNCTION random_gamma2(s, first) RESULT(fn_val)
! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO
! GAMMA2**(S-1) * EXP(-GAMMA2),
! USING A SWITCHING METHOD.
!    S = SHAPE PARAMETER OF DISTRIBUTION
!          (REAL < 1.0)
REAL, INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val
!     Local variables
REAL       :: r, x, w
REAL, SAVE :: a, p, c, uf, vr, d
IF (s <= zero .OR. s >= one) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE OUTSIDE PERMITTED RANGE'
  STOP
END IF
IF (first) THEN                        ! Initialization, if necessary
  a = one - s
  p = a/(a + s*EXP(-a))
  IF (s < vsmall) THEN
    WRITE(*, *) 'SHAPE PARAMETER VALUE TOO SMALL'
    STOP
  END IF
  c = one/s
  uf = p*(vsmall/a)**s
  vr = one - vsmall
  d = a*LOG(a)
END IF
DO
  CALL RANDOM_NUMBER(r)
  IF (r >= vr) THEN
    CYCLE
  ELSE IF (r > p) THEN
    x = a - LOG((one - r)/(one - p))
    w = a*LOG(x)-d
  ELSE IF (r > uf) THEN
    x = a*(r/p)**c
    w = x
  ELSE
    fn_val = zero
    RETURN
  END IF
  CALL RANDOM_NUMBER(r)
  IF (one-r <= w .AND. r > zero) THEN
    IF (r*(w + one) >= one) CYCLE
    IF (-LOG(r) <= w) CYCLE
  END IF
  EXIT
END DO
fn_val = x
RETURN
END FUNCTION random_gamma2
FUNCTION random_chisq(ndf, first) RESULT(fn_val)
!     Generates a random variate from the chi-squared distribution with
!     ndf degrees of freedom
INTEGER, INTENT(IN) :: ndf
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val
fn_val = two * random_gamma(half*ndf, first)
RETURN
END FUNCTION random_chisq
FUNCTION random_exponential() RESULT(fn_val)
! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
! TO EXP(-random_exponential), USING INVERSION.
REAL  :: fn_val
!     Local variable
REAL  :: r
DO
  CALL RANDOM_NUMBER(r)
  IF (r > zero) EXIT
END DO
fn_val = -LOG(r)
RETURN
END FUNCTION random_exponential
FUNCTION random_Weibull(a) RESULT(fn_val)
!     Generates a random variate from the Weibull distribution with
!     probability density:
!                      a
!               a-1  -x
!     f(x) = a.x    e
REAL, INTENT(IN) :: a
REAL             :: fn_val
!     For speed, there is no checking that a is not zero or very small.
fn_val = random_exponential() ** (one/a)
RETURN
END FUNCTION random_Weibull
FUNCTION random_beta(aa, bb, first) RESULT(fn_val)
! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
! FUNCTION GENERATES A RANDOM VARIATE IN [0,1]
! FROM A BETA DISTRIBUTION WITH DENSITY
! PROPORTIONAL TO BETA**(AA-1) * (1-BETA)**(BB-1).
! USING CHENG'S LOG LOGISTIC METHOD.
!     AA = SHAPE PARAMETER FROM DISTRIBUTION (0 < REAL)
!     BB = SHAPE PARAMETER FROM DISTRIBUTION (0 < REAL)
REAL, INTENT(IN)    :: aa, bb
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val
!     Local variables
REAL, PARAMETER  :: aln4 = 1.3862944
REAL             :: a, b, g, r, s, x, y, z
REAL, SAVE       :: d, f, h, t, c
LOGICAL, SAVE    :: swap
IF (aa <= zero .OR. bb <= zero) THEN
  WRITE(*, *) 'IMPERMISSIBLE SHAPE PARAMETER VALUE(S)'
  STOP
END IF
IF (first) THEN                        ! Initialization, if necessary
  a = aa
  b = bb
  swap = b > a
  IF (swap) THEN
    g = b
    b = a
    a = g
  END IF
  d = a/b
  f = a+b
  IF (b > one) THEN
    h = SQRT((two*a*b - f)/(f - two))
    t = one
  ELSE
    h = b
    t = one/(one + (a/(vlarge*b))**b)
  END IF
  c = a+h
END IF
DO
  CALL RANDOM_NUMBER(r)
  CALL RANDOM_NUMBER(x)
  s = r*r*x
  IF (r < vsmall .OR. s <= zero) CYCLE
  IF (r < t) THEN
    x = LOG(r/(one - r))/h
    y = d*EXP(x)
    z = c*x + f*LOG((one + d)/(one + y)) - aln4
    IF (s - one > z) THEN
      IF (s - s*z > one) CYCLE
      IF (LOG(s) > z) CYCLE
    END IF
    fn_val = y/(one + y)
  ELSE
    IF (4.0*s > (one + one/d)**f) CYCLE
    fn_val = one
  END IF
  EXIT
END DO
IF (swap) fn_val = one - fn_val
RETURN
END FUNCTION random_beta
FUNCTION random_t(m) RESULT(fn_val)
! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
! FUNCTION GENERATES A RANDOM VARIATE FROM A
! T DISTRIBUTION USING KINDERMAN AND MONAHAN'S RATIO METHOD.
!     M = DEGREES OF FREEDOM OF DISTRIBUTION
!           (1 <= 1NTEGER)
INTEGER, INTENT(IN) :: m
REAL                :: fn_val
!     Local variables
REAL, SAVE      :: s, c, a, f, g
REAL            :: r, x, v
REAL, PARAMETER :: three = 3.0, four = 4.0, quart = 0.25,   &
                   five = 5.0, sixteen = 16.0
INTEGER         :: mm = 0
IF (m < 1) THEN
  WRITE(*, *) 'IMPERMISSIBLE DEGREES OF FREEDOM'
  STOP
END IF
IF (m /= mm) THEN                    ! Initialization, if necessary
  s = m
  c = -quart*(s + one)
  a = four/(one + one/s)**c
  f = sixteen/a
  IF (m > 1) THEN
    g = s - one
    g = ((s + one)/g)**c*SQRT((s+s)/g)
  ELSE
    g = one
  END IF
  mm = m
END IF
DO
  CALL RANDOM_NUMBER(r)
  IF (r <= zero) CYCLE
  CALL RANDOM_NUMBER(v)
  x = (two*v - one)*g/r
  v = x*x
  IF (v > five - a*r) THEN
    IF (m >= 1 .AND. r*(v + three) > f) CYCLE
    IF (r > (one + v/s)**c) CYCLE
  END IF
  EXIT
END DO
fn_val = x
RETURN
END FUNCTION random_t
SUBROUTINE random_mvnorm(n, h, d, f, first, x, ier)
! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
! N.B. An extra argument, ier, has been added to Dagpunar's routine
!     SUBROUTINE GENERATES AN N VARIATE RANDOM NORMAL
!     VECTOR USING A CHOLESKY DECOMPOSITION.
! ARGUMENTS:
!        N = NUMBER OF VARIATES IN VECTOR
!           (INPUT,INTEGER >= 1)
!     H(J) = J'TH ELEMENT OF VECTOR OF MEANS
!           (INPUT,REAL)
!     X(J) = J'TH ELEMENT OF DELIVERED VECTOR
!           (OUTPUT,REAL)
!
!    D(J*(J-1)/2+I) = (I,J)'TH ELEMENT OF VARIANCE MATRIX (J> = I)
!            (INPUT,REAL)
!    F((J-1)*(2*N-J)/2+I) = (I,J)'TH ELEMENT OF LOWER TRIANGULAR
!           DECOMPOSITION OF VARIANCE MATRIX (J <= I)
!            (OUTPUT,REAL)
!    FIRST = .TRUE. IF THIS IS THE FIRST CALL OF THE ROUTINE
!    OR IF THE DISTRIBUTION HAS CHANGED SINCE THE LAST CALL OF THE ROUTINE.
!    OTHERWISE SET TO .FALSE.
!            (INPUT,LOGICAL)
!    ier = 1 if the input covariance matrix is not +ve definite
!        = 0 otherwise
INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN)      :: h(:), d(:)   ! d(n*(n+1)/2)
REAL, INTENT(IN OUT)  :: f(:)         ! f(n*(n+1)/2)
REAL, INTENT(OUT)     :: x(:)
LOGICAL, INTENT(IN)   :: first
INTEGER, INTENT(OUT)  :: ier
!     Local variables
INTEGER       :: j, i, m
REAL          :: y, v
INTEGER, SAVE :: n2
IF (n < 1) THEN
  WRITE(*, *) 'SIZE OF VECTOR IS NON POSITIVE'
  STOP
END IF
ier = 0
IF (first) THEN                        ! Initialization, if necessary
  n2 = 2*n
  IF (d(1) < zero) THEN
    ier = 1
    RETURN
  END IF
  f(1) = SQRT(d(1))
  y = one/f(1)
  DO j = 2,n
    f(j) = d(1+j*(j-1)/2) * y
  END DO
  DO i = 2,n
    v = d(i*(i-1)/2+i)
    DO m = 1,i-1
      v = v - f((m-1)*(n2-m)/2+i)**2
    END DO
    IF (v < zero) THEN
      ier = 1
      RETURN
    END IF
    v = SQRT(v)
    y = one/v
    f((i-1)*(n2-i)/2+i) = v
    DO j = i+1,n
      v = d(j*(j-1)/2+i)
      DO m = 1,i-1
        v = v - f((m-1)*(n2-m)/2+i)*f((m-1)*(n2-m)/2 + j)
      END DO ! m = 1,i-1
      f((i-1)*(n2-i)/2 + j) = v*y
    END DO ! j = i+1,n
  END DO ! i = 2,n
END IF
x(1:n) = h(1:n)
DO j = 1,n
  y = random_normal()
  DO i = j,n
    x(i) = x(i) + f((j-1)*(n2-j)/2 + i) * y
  END DO ! i = j,n
END DO ! j = 1,n
RETURN
END SUBROUTINE random_mvnorm
FUNCTION random_inv_gauss(h, b, first) RESULT(fn_val)
! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY] FROM
! A REPARAMETERISED GENERALISED INVERSE GAUSSIAN (GIG) DISTRIBUTION
! WITH DENSITY PROPORTIONAL TO  GIG**(H-1) * EXP(-0.5*B*(GIG+1/GIG))
! USING A RATIO METHOD.
!     H = PARAMETER OF DISTRIBUTION (0 <= REAL)
!     B = PARAMETER OF DISTRIBUTION (0 < REAL)
REAL, INTENT(IN)    :: h, b
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val
!     Local variables
REAL            :: ym, xm, r, w, r1, r2, x
REAL, SAVE      :: a, c, d, e
REAL, PARAMETER :: quart = 0.25
IF (h < zero .OR. b <= zero) THEN
  WRITE(*, *) 'IMPERMISSIBLE DISTRIBUTION PARAMETER VALUES'
  STOP
END IF
IF (first) THEN                        ! Initialization, if necessary
  IF (h > quart*b*SQRT(vlarge)) THEN
    WRITE(*, *) 'THE RATIO H:B IS TOO SMALL'
    STOP
  END IF
  e = b*b
  d = h + one
  ym = (-d + SQRT(d*d + e))/b
  IF (ym < vsmall) THEN
    WRITE(*, *) 'THE VALUE OF B IS TOO SMALL'
    STOP
  END IF
  d = h - one
  xm = (d + SQRT(d*d + e))/b
  d = half*d
  e = -quart*b
  r = xm + one/xm
  w = xm*ym
  a = w**(-half*h) * SQRT(xm/ym) * EXP(-e*(r - ym - one/ym))
  IF (a < vsmall) THEN
    WRITE(*, *) 'THE VALUE OF H IS TOO LARGE'
    STOP
  END IF
  c = -d*LOG(xm) - e*r
END IF
DO
  CALL RANDOM_NUMBER(r1)
  IF (r1 <= zero) CYCLE
  CALL RANDOM_NUMBER(r2)
  x = a*r2/r1
  IF (x <= zero) CYCLE
  IF (LOG(r1) < d*LOG(x) + e*(x + one/x) + c) EXIT
END DO
fn_val = x
RETURN
END FUNCTION random_inv_gauss
FUNCTION random_Poisson(mu, first) RESULT(ival)
!**********************************************************************
!     Translated to Fortran 90 by Alan Miller from:
!                           RANLIB
!
!     Library of Fortran Routines for Random Number Generation
!
!                    Compiled and Written by:
!
!                         Barry W. Brown
!                          James Lovato
!
!             Department of Biomathematics, Box 237
!             The University of Texas, M.D. Anderson Cancer Center
!             1515 Holcombe Boulevard
!             Houston, TX      77030
!
! This work was supported by grant CA-16672 from the National Cancer Institute.
!                    GENerate POIsson random deviate
!                            Function
! Generates a single random deviate from a Poisson distribution with mean mu.
!                            Arguments
!     mu --> The mean of the Poisson distribution from which
!            a random deviate is to be generated.
!                              REAL mu
!                              Method
!     For details see:
!               Ahrens, J.H. and Dieter, U.
!               Computer Generation of Poisson Deviates
!               From Modified Normal Distributions.
!               ACM Trans. Math. Software, 8, 2
!               (June 1982),163-179
!     TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
!     COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL
!     SEPARATION OF CASES A AND B
!     .. Scalar Arguments ..
REAL, INTENT(IN)    :: mu
LOGICAL, INTENT(IN) :: first
INTEGER             :: ival
!     ..
!     .. Local Scalars ..
REAL          :: b1, b2, c, c0, c1, c2, c3, del, difmuk, e, fk, fx, fy, g,  &
                 omega, px, py, t, u, v, x, xx
REAL, SAVE    :: s, d, p, q, p0
INTEGER       :: j, k, kflag
LOGICAL, SAVE :: full_init
INTEGER, SAVE :: l, m
!     ..
!     .. Local Arrays ..
REAL, SAVE    :: pp(35)
!     ..
!     .. Data statements ..
REAL, PARAMETER :: a0 = -.5, a1 = .3333333, a2 = -.2500068, a3 = .2000118,  &
                   a4 = -.1661269, a5 = .1421878, a6 = -.1384794,   &
                   a7 = .1250060
REAL, PARAMETER :: fact(10) = (/ 1., 1., 2., 6., 24., 120., 720., 5040.,  &
                                 40320., 362880. /)
!     ..
!     .. Executable Statements ..
IF (mu > 10.0) THEN
!     C A S E  A. (RECALCULATION OF S, D, L IF MU HAS CHANGED)
  IF (first) THEN
    s = SQRT(mu)
    d = 6.0*mu*mu
!             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
!             PROBABILITIES FK WHENEVER K >= M(MU). L=IFIX(MU-1.1484)
!             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .
    l = mu - 1.1484
    full_init = .false.
  END IF
!     STEP N. NORMAL SAMPLE - random_normal() FOR STANDARD NORMAL DEVIATE
  g = mu + s*random_normal()
  IF (g > 0.0) THEN
    ival = g
!     STEP I. IMMEDIATE ACCEPTANCE IF ival IS LARGE ENOUGH
    IF (ival>=l) RETURN
!     STEP S. SQUEEZE ACCEPTANCE - SAMPLE U
    fk = ival
    difmuk = mu - fk
    CALL RANDOM_NUMBER(u)
    IF (d*u >= difmuk*difmuk*difmuk) RETURN
  END IF
!     STEP P. PREPARATIONS FOR STEPS Q AND H.
!             (RECALCULATIONS OF PARAMETERS IF NECESSARY)
!             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
!             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
!             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
!             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.
  IF (.NOT. full_init) THEN
    omega = .3989423/s
    b1 = .4166667E-1/mu
    b2 = .3*b1*b1
    c3 = .1428571*b1*b2
    c2 = b2 - 15.*c3
    c1 = b1 - 6.*b2 + 45.*c3
    c0 = 1. - b1 + 3.*b2 - 15.*c3
    c = .1069/mu
    full_init = .true.
  END IF
  IF (g < 0.0) GO TO 50
!             'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)
  kflag = 0
  GO TO 70
!     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)
  40 IF (fy-u*fy <= py*EXP(px-fx)) RETURN
!     STEP E. EXPONENTIAL SAMPLE - random_exponential() FOR STANDARD EXPONENTIAL
!             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
!             (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)
  50 e = random_exponential()
  CALL RANDOM_NUMBER(u)
  u = u + u - one
  t = 1.8 + SIGN(e, u)
  IF (t <= (-.6744)) GO TO 50
  ival = mu + s*t
  fk = ival
  difmuk = mu - fk
!             'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)
  kflag = 1
  GO TO 70
!     STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)
  60 IF (c*ABS(u) > py*EXP(px+e) - fy*EXP(fx+e)) GO TO 50
  RETURN
!     STEP F. 'SUBROUTINE' F. CALCULATION OF PX, PY, FX, FY.
!             CASE ival < 10 USES FACTORIALS FROM TABLE FACT
  70 IF (ival>=10) GO TO 80
  px = -mu
  py = mu**ival/fact(ival+1)
  GO TO 110
!             CASE ival >= 10 USES POLYNOMIAL APPROXIMATION
!             A0-A7 FOR ACCURACY WHEN ADVISABLE
!             .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)
  80 del = .8333333E-1/fk
  del = del - 4.8*del*del*del
  v = difmuk/fk
  IF (ABS(v)>0.25) THEN
    px = fk*LOG(one + v) - difmuk - del
  ELSE
    px = fk*v*v* (((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0) - del
  END IF
  py = .3989423/SQRT(fk)
  110 x = (half - difmuk)/s
  xx = x*x
  fx = -half*xx
  fy = omega* (((c3*xx + c2)*xx + c1)*xx + c0)
  IF (kflag <= 0) GO TO 40
  GO TO 60
!---------------------------------------------------------------------------
!     C A S E  B.    mu < 10
!     START NEW TABLE AND CALCULATE P0 IF NECESSARY
ELSE
  IF (first) THEN
    m = MAX(1, INT(mu))
    l = 0
    p = EXP(-mu)
    q = p
    p0 = p
  END IF
!     STEP U. UNIFORM SAMPLE FOR INVERSION METHOD
  DO
    CALL RANDOM_NUMBER(u)
    ival = 0
    IF (u <= p0) RETURN
!     STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
!             PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
!             (0.458=PP(9) FOR MU=10)
    IF (l == 0) GO TO 150
    j = 1
    IF (u > 0.458) j = MIN(l, m)
    DO k = j, l
      IF (u <= pp(k)) GO TO 180
    END DO
    IF (l == 35) CYCLE
!     STEP C. CREATION OF NEW POISSON PROBABILITIES P
!             AND THEIR CUMULATIVES Q=PP(K)
    150 l = l + 1
    DO k = l, 35
      p = p*mu / k
      q = q + p
      pp(k) = q
      IF (u <= q) GO TO 170
    END DO
    l = 35
  END DO
  170 l = k
  180 ival = k
  RETURN
END IF
RETURN
END FUNCTION random_Poisson
FUNCTION random_binomial1(n, p, first) RESULT(ival)
! FUNCTION GENERATES A RANDOM BINOMIAL VARIATE USING C.D.Kemp's method.
! This algorithm is suitable when many random variates are required
! with the SAME parameter values for n & p.
!    P = BERNOULLI SUCCESS PROBABILITY
!           (0 <= REAL <= 1)
!    N = NUMBER OF BERNOULLI TRIALS
!           (1 <= INTEGER)
!    FIRST = .TRUE. for the first call using the current parameter values
!          = .FALSE. if the values of (n,p) are unchanged from last call
! Reference: Kemp, C.D. (1986). `A modal method for generating binomial
!            variables', Commun. Statist. - Theor. Meth. 15(3), 805-813.
INTEGER, INTENT(IN) :: n
REAL, INTENT(IN)    :: p
LOGICAL, INTENT(IN) :: first
INTEGER             :: ival
!     Local variables
INTEGER         :: ru, rd
INTEGER, SAVE   :: r0
REAL            :: u, pd, pu
REAL, SAVE      :: odds_ratio, p_r
REAL, PARAMETER :: zero = 0.0, one = 1.0
IF (first) THEN
  r0 = (n+1)*p
  p_r = bin_prob(n, p, r0)
  odds_ratio = p / (one - p)
END IF
u = u - p_r
IF (u < zero) THEN
  ival = r0
  RETURN
END IF
pu = p_r
ru = r0
pd = p_r
rd = r0
DO
  rd = rd - 1
  IF (rd >= 0) THEN
    pd = pd * (rd+1) / (odds_ratio * (n-rd))
    u = u - pd
    IF (u < zero) THEN
      ival = rd
      RETURN
    END IF
  END IF
  ru = ru + 1
  IF (ru <= n) THEN
    pu = pu * (n-ru+1) * odds_ratio / ru
    u = u - pu
    IF (u < zero) THEN
      ival = ru
      RETURN
    END IF
  END IF
END DO
!     This point should not be reached, but just in case:
ival = r0
RETURN
END FUNCTION random_binomial1
FUNCTION bin_prob(n, p, r) RESULT(fn_val)
!     Calculate a binomial probability
INTEGER, INTENT(IN) :: n, r
REAL, INTENT(IN)    :: p
REAL                :: fn_val
!     Local variable
REAL                :: one = 1.0
fn_val = EXP( lngamma(DBLE(n+1)) - lngamma(DBLE(r+1)) - lngamma(DBLE(n-r+1)) &
              + r*LOG(p) + (n-r)*LOG(one - p) )
RETURN
END FUNCTION bin_prob
FUNCTION lngamma(x) RESULT(fn_val)
! Logarithm to base e of the gamma function.
!
! Accurate to about 1.e-14.
! Programmer: Alan Miller
! Latest revision of Fortran 77 version - 28 February 1988
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val
!       Local variables
REAL (dp) :: a1 = -4.166666666554424D-02, a2 = 2.430554511376954D-03,  &
             a3 = -7.685928044064347D-04, a4 = 5.660478426014386D-04,  &
             temp, arg, product, lnrt2pi = 9.189385332046727D-1,       &
             pi = 3.141592653589793D0
LOGICAL   :: reflect
!       lngamma is not defined if x = 0 or a negative integer.
IF (x > 0.d0) GO TO 10
IF (x /= INT(x)) GO TO 10
fn_val = 0.d0
RETURN
!       If x < 0, use the reflection formula:
!               gamma(x) * gamma(1-x) = pi * cosec(pi.x)
10 reflect = (x < 0.d0)
IF (reflect) THEN
  arg = 1.d0 - x
ELSE
  arg = x
END IF
!       Increase the argument, if necessary, to make it > 10.
product = 1.d0
20 IF (arg <= 10.d0) THEN
  product = product * arg
  arg = arg + 1.d0
  GO TO 20
END IF
!  Use a polynomial approximation to Stirling's formula.
!  N.B. The real Stirling's formula is used here, not the simpler, but less
!       accurate formula given by De Moivre in a letter to Stirling, which
!       is the one usually quoted.
arg = arg - 0.5D0
temp = 1.d0/arg**2
fn_val = lnrt2pi + arg * (LOG(arg) - 1.d0 + &
                  (((a4*temp + a3)*temp + a2)*temp + a1)*temp) - LOG(product)
IF (reflect) THEN
  temp = SIN(pi * x)
  fn_val = LOG(pi/temp) - fn_val
END IF
RETURN
END FUNCTION lngamma
FUNCTION random_binomial2(n, pp, first) RESULT(ival)
!**********************************************************************
!     Translated to Fortran 90 by Alan Miller from:
!                              RANLIB
!
!     Library of Fortran Routines for Random Number Generation
!
!                      Compiled and Written by:
!
!                           Barry W. Brown
!                            James Lovato
!
!               Department of Biomathematics, Box 237
!               The University of Texas, M.D. Anderson Cancer Center
!               1515 Holcombe Boulevard
!               Houston, TX      77030
!
! This work was supported by grant CA-16672 from the National Cancer Institute.
!                    GENerate BINomial random deviate
!                              Function
!     Generates a single random deviate from a binomial
!     distribution whose number of trials is N and whose
!     probability of an event in each trial is P.
!                              Arguments
!     N  --> The number of trials in the binomial distribution
!            from which a random deviate is to be generated.
!                              INTEGER N
!     P  --> The probability of an event in each trial of the
!            binomial distribution from which a random deviate
!            is to be generated.
!                              REAL P
!     FIRST --> Set FIRST = .TRUE. for the first call to perform initialization
!               the set FIRST = .FALSE. for further calls using the same pair
!               of parameter values (N, P).
!                              LOGICAL FIRST
!     random_binomial2 <-- A random deviate yielding the number of events
!                from N independent trials, each of which has
!                a probability of event P.
!                              INTEGER random_binomial
!                              Method
!     This is algorithm BTPE from:
!         Kachitvichyanukul, V. and Schmeiser, B. W.
!         Binomial Random Variate Generation.
!         Communications of the ACM, 31, 2 (February, 1988) 216.
!**********************************************************************
!*****DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY
!     ..
!     .. Scalar Arguments ..
REAL, INTENT(IN)    :: pp
INTEGER, INTENT(IN) :: n
LOGICAL, INTENT(IN) :: first
INTEGER             :: ival
!     ..
!     .. Local Scalars ..
REAL            :: alv, amaxp, f, f1, f2, u, v, w, w2, x, x1, x2, ynorm, z, z2
REAL, PARAMETER :: zero = 0.0, half = 0.5, one = 1.0
INTEGER         :: i, ix, ix1, k, mp
INTEGER, SAVE   :: m
REAL, SAVE      :: p, q, xnp, ffm, fm, xnpq, p1, xm, xl, xr, c, al, xll,  &
                   xlr, p2, p3, p4, qn, r, g
!     ..
!     .. Executable Statements ..
!*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE
IF (first) THEN
  p = MIN(pp, one-pp)
  q = one - p
  xnp = n * p
END IF
IF (xnp > 30.) THEN
  IF (first) THEN
    ffm = xnp + p
    m = ffm
    fm = m
    xnpq = xnp * q
    p1 = INT(2.195*SQRT(xnpq) - 4.6*q) + half
    xm = fm + half
    xl = xm - p1
    xr = xm + p1
    c = 0.134 + 20.5 / (15.3 + fm)
    al = (ffm-xl) / (ffm - xl*p)
    xll = al * (one + half*al)
    al = (xr - ffm) / (xr*q)
    xlr = al * (one + half*al)
    p2 = p1 * (one + c + c)
    p3 = p2 + c / xll
    p4 = p3 + c / xlr
  END IF
!*****GENERATE VARIATE, Binomial mean at least 30.
  20 CALL RANDOM_NUMBER(u)
  u = u * p4
  CALL RANDOM_NUMBER(v)
!     TRIANGULAR REGION
  IF (u <= p1) THEN
    ix = xm - p1 * v + u
    GO TO 110
  END IF
!     PARALLELOGRAM REGION
  IF (u <= p2) THEN
    x = xl + (u-p1) / c
    v = v * c + one - ABS(xm-x) / p1
    IF (v > one .OR. v <= zero) GO TO 20
    ix = x
  ELSE
!     LEFT TAIL
    IF (u <= p3) THEN
      ix = xl + LOG(v) / xll
      IF (ix < 0) GO TO 20
      v = v * (u-p2) * xll
    ELSE
!     RIGHT TAIL
      ix = xr - LOG(v) / xlr
      IF (ix > n) GO TO 20
      v = v * (u-p3) * xlr
    END IF
  END IF
!*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST
  k = ABS(ix-m)
  IF (k <= 20 .OR. k >= xnpq/2-1) THEN
!     EXPLICIT EVALUATION
    f = one
    r = p / q
    g = (n+1) * r
    IF (m < ix) THEN
      mp = m + 1
      DO i = mp, ix
        f = f * (g/i-r)
      END DO
    ELSE IF (m > ix) THEN
      ix1 = ix + 1
      DO i = ix1, m
        f = f / (g/i-r)
      END DO
    END IF
    IF (v > f) THEN
      GO TO 20
    ELSE
      GO TO 110
    END IF
  END IF
!     SQUEEZING USING UPPER AND LOWER BOUNDS ON LOG(F(X))
  amaxp = (k/xnpq) * ((k*(k/3. + .625) + .1666666666666)/xnpq + half)
  ynorm = -k * k / (2.*xnpq)
  alv = LOG(v)
  IF (alv<ynorm - amaxp) GO TO 110
  IF (alv>ynorm + amaxp) GO TO 20
!     STIRLING'S (actually de Moivre's) FORMULA TO MACHINE ACCURACY FOR
!     THE FINAL ACCEPTANCE/REJECTION TEST
  x1 = ix + 1
  f1 = fm + one
  z = n + 1 - fm
  w = n - ix + one
  z2 = z * z
  x2 = x1 * x1
  f2 = f1 * f1
  w2 = w * w
  IF (alv - (xm*LOG(f1/x1) + (n-m+half)*LOG(z/w) + (ix-m)*LOG(w*p/(x1*q)) +    &
      (13860.-(462.-(132.-(99.-140./f2)/f2)/f2)/f2)/f1/166320. +               &
      (13860.-(462.-(132.-(99.-140./z2)/z2)/z2)/z2)/z/166320. +                &
      (13860.-(462.-(132.-(99.-140./x2)/x2)/x2)/x2)/x1/166320. +               &
      (13860.-(462.-(132.-(99.-140./w2)/w2)/w2)/w2)/w/166320.) > zero) THEN
    GO TO 20
  ELSE
    GO TO 110
  END IF
ELSE
!     INVERSE CDF LOGIC FOR MEAN LESS THAN 30
  IF (first) THEN
    qn = q ** n
    r = p / q
    g = r * (n+1)
  END IF
  90 ix = 0
  f = qn
  CALL RANDOM_NUMBER(u)
  100 IF (u >= f) THEN
    IF (ix > 110) GO TO 90
    u = u - f
    ix = ix + 1
    f = f * (g/ix - r)
    GO TO 100
  END IF
END IF
110 IF (pp > half) ix = n - ix
ival = ix
RETURN
END FUNCTION random_binomial2
FUNCTION random_neg_binomial(sk, p) RESULT(ival)
! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
! FUNCTION GENERATES A RANDOM NEGATIVE BINOMIAL VARIATE USING UNSTORED
! INVERSION AND/OR THE REPRODUCTIVE PROPERTY.
!    SK = NUMBER OF FAILURES REQUIRED (Dagpunar's words!)
!       = the `power' parameter of the negative binomial
!           (0 < REAL)
!    P = BERNOULLI SUCCESS PROBABILITY
!           (0 < REAL < 1)
! THE PARAMETER H IS SET SO THAT UNSTORED INVERSION ONLY IS USED WHEN P <= H,
! OTHERWISE A COMBINATION OF UNSTORED INVERSION AND
! THE REPRODUCTIVE PROPERTY IS USED.
REAL, INTENT(IN)   :: sk, p
INTEGER            :: ival
!     Local variables
! THE PARAMETER ULN = -LOG(MACHINE'S SMALLEST REAL NUMBER).
REAL, PARAMETER    :: h = 0.7
REAL               :: q, x, st, uln, v, r, s, y, g
INTEGER            :: k, i, n
IF (sk <= zero .OR. p <= zero .OR. p >= one) THEN
  WRITE(*, *) 'IMPERMISSIBLE DISTRIBUTION PARAMETER VALUES'
  STOP
END IF
q = one - p
x = zero
st = sk
IF (p > h) THEN
  v = one/LOG(p)
  k = st
  DO i = 1,k
    DO
      CALL RANDOM_NUMBER(r)
      IF (r > zero) EXIT
    END DO
    n = v*LOG(r)
    x = x + n
  END DO
  st = st - k
END IF
s = zero
uln = -LOG(vsmall)
IF (st > -uln/LOG(q)) THEN
  WRITE(*, *) ' P IS TOO LARGE FOR THIS VALUE OF SK'
  STOP
END IF
y = q**st
g = st
DO
  IF (y > r) EXIT
  r = r - y
  s = s + one
  y = y*p*g/s
  g = g + one
END DO
ival = x + s + half
RETURN
END FUNCTION random_neg_binomial
FUNCTION random_von_Mises(k, first) RESULT(fn_val)
!     Algorithm VMD from:
!     Dagpunar, J.S. (1990) `Sampling from the von Mises distribution via a
!     comparison of random numbers', J. of Appl. Statist., 17, 165-168.
!     Fortran 90 code by Alan Miller
!     CSIRO Division of Mathematical & Information Sciences
!     Arguments:
!     k (real)        parameter of the von Mises distribution.
!     first (logical) set to .TRUE. the first time that the function
!                     is called, or the first time with a new value
!                     for k.   When first = .TRUE., the function sets
!                     up starting values and may be very much slower.
REAL, INTENT(IN)     :: k
LOGICAL, INTENT(IN)  :: first
REAL                 :: fn_val
!     Local variables
INTEGER          :: j, n
INTEGER, SAVE    :: nk
REAL, PARAMETER  :: pi = 3.14159265
REAL, SAVE       :: p(20), theta(0:20)
REAL             :: sump, r, th, lambda, rlast
REAL (dp)        :: dk
IF (first) THEN                        ! Initialization, if necessary
  IF (k < zero) THEN
    WRITE(*, *) '** Error: argument k for random_von_Mises = ', k
    RETURN
  END IF
  nk = k + k + one
  IF (nk > 20) THEN
    WRITE(*, *) '** Error: argument k for random_von_Mises = ', k
    RETURN
  END IF
  dk = k
  theta(0) = zero
  IF (k > half) THEN
!     Set up array p of probabilities.
    sump = zero
    DO j = 1, nk
      IF (j < nk) THEN
        theta(j) = ACOS(one - j/k)
      ELSE
        theta(nk) = pi
      END IF
!     Numerical integration of e^[k.cos(x)] from theta(j-1) to theta(j)
      CALL integral(theta(j-1), theta(j), p(j), dk)
      sump = sump + p(j)
    END DO
    p(1:nk) = p(1:nk) / sump
  ELSE
    p(1) = one
    theta(1) = pi
  END IF                         ! if k > 0.5
END IF                           ! if first
DO j = 1, nk
  r = r - p(j)
  IF (r < zero) EXIT
END DO
r = -r/p(j)
DO
  th = theta(j-1) + r*(theta(j) - theta(j-1))
  lambda = k - j + one - k*COS(th)
  n = 1
  rlast = lambda
  DO
    CALL RANDOM_NUMBER(r)
    IF (r > rlast) EXIT
    n = n + 1
    rlast = r
  END DO
  IF (n .NE. 2*(n/2)) EXIT         ! is n even?
  CALL RANDOM_NUMBER(r)
END DO
fn_val = SIGN(th, (r - rlast)/(one - rlast) - half)
RETURN
END FUNCTION random_von_Mises
SUBROUTINE integral(a, b, result, dk)
!     Gaussian integration of exp(k.cosx) from a to b.
REAL (dp), INTENT(IN) :: dk
REAL, INTENT(IN)      :: a, b
REAL, INTENT(OUT)     :: result
!     Local variables
REAL (dp)  :: xmid, range, x1, x2,                                    &
  x(3) = (/0.238619186083197_dp, 0.661209386466265_dp, 0.932469514203152_dp/), &
  w(3) = (/0.467913934572691_dp, 0.360761573048139_dp, 0.171324492379170_dp/)
INTEGER    :: i
xmid = (a + b)/2._dp
range = (b - a)/2._dp
result = 0._dp
DO i = 1, 3
  x1 = xmid + x(i)*range
  x2 = xmid - x(i)*range
  result = result + w(i)*(EXP(dk*COS(x1)) + EXP(dk*COS(x2)))
END DO
result = result * range
RETURN
END SUBROUTINE integral
FUNCTION random_Cauchy() RESULT(fn_val)
!     Generate a random deviate from the standard Cauchy distribution
REAL     :: fn_val
!     Local variables
REAL     :: v(2)
DO
  CALL RANDOM_NUMBER(v)
  v = two*(v - half)
  IF (ABS(v(2)) < vsmall) CYCLE               ! Test for zero
  IF (v(1)**2 + v(2)**2 < one) EXIT
END DO
fn_val = v(1) / v(2)
RETURN
END FUNCTION random_Cauchy
SUBROUTINE random_order(order, n)
!     Generate a random ordering of the integers 1 ... n.
INTEGER, INTENT(IN)  :: n
INTEGER, INTENT(OUT) :: order(n)
!     Local variables
INTEGER :: i, j, k
REAL    :: wk
DO i = 1, n
  order(i) = i
END DO
!     Starting at the end, swap the current last indicator with one
!     randomly chosen from those preceeding it.
DO i = n, 2, -1
  CALL RANDOM_NUMBER(wk)
  j = 1 + i * wk
  IF (j < i) THEN
    k = order(i)
    order(i) = order(j)
    order(j) = k
  END IF
END DO
RETURN
END SUBROUTINE random_order
SUBROUTINE seed_random_number(iounit)
INTEGER, INTENT(IN)  :: iounit
! Local variables
INTEGER              :: k
INTEGER, ALLOCATABLE :: seed(:)
ALLOCATE( seed(k) )
WRITE(*, '(a, i2, a)')' Enter ', k, ' integers for random no. seeds: '
READ(*, *) seed
WRITE(iounit, '(a, (7i10))') ' Random no. seeds: ', seed
DEALLOCATE( seed )
RETURN
END SUBROUTINE seed_random_number
! END MODULE random ! random.f90 is cat-ted (in)to npag_utils.f90 to support DO140
END MODULE npag_utils
