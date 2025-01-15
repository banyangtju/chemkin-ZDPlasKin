program discomb
!  the module used
    use ZDPlasKin
	
!  define scaler field
    IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
    double precision :: gas_pressure,gas_temperature,gas_density,MOL_C2H4,MOL_N2,MOL_O2,MOL_TOT,MASS_TOT,&
    gas_density_C2H4,gas_density_O2,gas_density_N2,X_tot,t,e_density,reduced_field,source,SPEC_HEAT_RATIO,&
    Y_CH4,Y_O2,Y_N2
    integer :: i, j, idata_max,w,numberdivide
    double precision :: time,time1, dtime1,time2,dtim2, EN, NE,density_tot,density_my(species_max),dtime,&
    power,densityTot,eleDiffV,powerSum,energySum,energy
    double precision, allocatable :: data_in(:,:)
!   integer, parameter          :: ipulses = 10
    parameter (LENIWK=500000, LENRWK=500000, LENCWK=10000, NK=5, NLMAX=55) ! LENIMC=100000,LENRMC=100000)
    parameter (LIN=5, LOUT=6,LINKCK=25,LINKMC=26)
    parameter (KMAX=200, ITOL=1, IOPT=1, RTOL=1.0E-12)
    parameter (ITASK=1, ATOL=1.0E-16,elementaryCharge=1.60217662d-19)
    CHARACTER*16   KSYM(KMAX), CWORK(LENCWK),specName(200)
!    real, dimension (:), allocatable :: specNmae
    DIMENSION IWORK(LENIWK), RWORK(LENRWK), X(KMAX), Z(KMAX)
    LOGICAL  IERR, ZD_init
    EXTERNAL FUN 
    COMMON /ICONS/ KK, NP, NWT, NH, NWDOT, MM
    COMMON /ICONS_specName/ specName
    COMMON /ICONS_ZD_init/ ZD_init
    
    DATA X/KMAX*0.0/, KSYM/KMAX*' '/

    character(*), parameter :: file_out = 'out.dat' ! ,file_out_power ='power.dat'
    OPEN (unit=LINKCK, file = 'chem.bin',form ='UNFORMATTED')

! ------------------------------------------------- start initialize ZDPlasKin ----------------------------------------------- !
    ZD_init = .FALSE.
! ------------------------------------------------- end initialize ZDPlasKIN ----------------------------------------------- !

! ------------------------------------------------- start initialize CHEMKIN ----------------------------------------------- !
    write(*,*)'-----------------------------start initialize CHEMKIN-----------------------------'
    call CKLEN (LINKCK, LOUT, LENI, LENR, LENC)
    CALL CKINIT (LENIWK, LENRWK, LENCWK, LINKCK, LOUT, IWORK, RWORK, CWORK)
!   CLOSE (LINKCK)
    CALL CKINDX (IWORK, RWORK, MM, KK, II, NFIT)
!   call MCINIT (LINKMC, LOUT, LENIMC, LENRMC, IWORK(50000), RWORK(50000))
!   CLOSE (LINKMC)
    NEQ   = KK + 1 
    LRW   = 22 + 9*NEQ + 2*NEQ**2  
    NVODE = LENR + 1 
    NP    = NVODE + LRW 
    NWT   = NP + 1
    NH    = NWT  + KK
    NWDOT = NH   + KK 
    NTOT  = NWDOT+ KK - 1
    LIW   = 30 + NEQ
    IVODE = LENI + 1
    ITOT  = IVODE + LIW - 1
    MF = 22
    ISTATE = 1

! IOPT = 1
    DO i=5,10
        RWORK(NVODE+i-1) = 0.d0
        IWORK(IVODE+i-1) = 0
    ENDDO
    RWORK(NVODE+6-1) = 1E-9 ! HMAX
    IWORK(IVODE+6-1) = 50000

    CALL CKSYMS (CWORK, LOUT, KSYM, IERR)
    write(*,*)'ksym = ',KSYM
    specName(1:KK) = KSYM(1:KK)

    IF (IERR) STOP
    
    IF (KK.GT.KMAX .OR. LENRWK.LT.NTOT .OR. LENIWK.LT.ITOT) THEN
        IF (KK .GT. KMAX)  WRITE (LOUT, *)' Error...KMAX too small...must be at least ', KK
        IF (LENRWK .LT. NTOT) WRITE (LOUT, *)' Error...LENRWK too small...must be at least', NTOT
        IF (LENIWK .LT. ITOT) WRITE (LOUT, *)' Error...LENIWK too small...must be at least', ITOT
        STOP
    ENDIF
    
    CALL CKWT (IWORK, RWORK, RWORK(NWT))  
    CALL CKRP (IWORK, RWORK, RU, RUC, PATM)
! ------------------------------------------------- end initialize CHEMKIN ----------------------------------------------- !


! ------------------------------------------------- open file for output and write headers ----------------------------- !
    ! open(11,file=file_out_power)   ! write to power.dat
    ! write(11,7130)

    open(21,file=file_out)   !  write to out.dat
    WRITE (21, 7100) (KSYM(K)(1:10), K=1,KK)
! ------------------------------------------------- input initial value ----------------------------------------------- !

    write(*,*)'-----------------------------start load initial temperature pressure-----------------------------'
    
    gas_pressure    = 1.0 ! atm
    RWORK(NP) = gas_pressure*PATM
    gas_temperature = 900.0
    
    ER = 0.5
    a = 2.0/ER

    CALL CKSNUM ('CH4', 1, LOUT, KSYM, KK, KNUM, NVAL, VAL, IERR)
    X(KNUM) = 1.0/(1.0+a+a*3.76)
   
    CALL CKSNUM ('O2', 1, LOUT, KSYM, KK, KNUM, NVAL, VAL, IERR)
    X(KNUM) = a/(1.0+a+a*3.76)

    CALL CKSNUM ('N2', 1, LOUT, KSYM, KK, KNUM, NVAL, VAL, IERR)
    X(KNUM) = a*3.76/(1.0+a+a*3.76) - 1E-8
    
    CALL CKSNUM ('E', 1, LOUT, KSYM, KK, KNUM, NVAL, VAL, IERR)
    X(KNUM) = 1E-8

    CALL CKXTY (X, IWORK, RWORK, Z(2))
    Z(1) = gas_temperature

! ------------------------------------------------- time integration ------------------------------------------------- !
! +++++++++++++++++++++++++++++++++++++++++++++++++ time integration +++++++++++++++++++++++++++++++++++++++++++++++++ !
! ************************************************* time integration ************************************************* !
! ------------------------------------------------- time integration ------------------------------------------------- !

    j=0
    w=1
    time =0.d0
    dtime = 1.0d-9
    endtime = 2E-2
    ipulsesTot = 200
    ipulse = 0
    frequency = 100.D4 ! 1 MHz
    disPulseWidth = 10.D-9 ! 10 ns
    
    dtimePulseWidth = 1E-9
    numberdivide = 10 
    dtimeIntervel = 1E-8
    numberdivideTot = 10+99 ! 10*1E-9 + 99*1E-8 = 1/1E6

    numbDisCount = 1

do while(time < endtime)

    CALL CKYTX (Z(2), IWORK, RWORK, X)
    
	if (ipulse < ipulsesTot .AND. numbDisCount < numberdivide) then
        dtime = dtimePulseWidth
        ISTATE = 1
        RWORK(NVODE+6-1) = dtimePulseWidth ! HMAX
    endif
    
    if (ipulse < ipulsesTot .AND. numbDisCount == numberdivide) then
        dtime = dtimeIntervel
        ISTATE = 1
        RWORK(NVODE+6-1) = 0 ! HMAX
    endif

    t = time
    TOUT = t + dtime
    
    350 CONTINUE
    CALL DVODE(FUN, NEQ, Z, t, TOUT, ITOL, RTOL, ATOL, ITASK,ISTATE, IOPT, &
    RWORK(NVODE), LRW, IWORK(IVODE),LIW, JAC, MF, RWORK, IWORK)
    
    IF (ISTATE .LE. -2) THEN
        IF (ISTATE .EQ. -1) THEN
            ISTATE = 2
            GO TO 350
        ELSE
            write (LOUT,*) ' ISTATE= ',ISTATE
            STOP
        ENDIF
    ENDIF
!    
    time=TOUT
    
    DO K=2,KK+1
        IF(Z(k) < 1.D-20)THEN
            Z(k) = 0.D0
        ENDIF
    ENDDO
    
    CALL CKYTX (Z(2), IWORK, RWORK, X)

    j = j+1
    
    numbDisCount = numbDisCount + 1
    if (numbDisCount == numberdivideTot) then
        numbDisCount = 1
        ipulse = ipulse + 1
    endif
    
    if( mod(j,100) == 0 )THEN
        write (21, 7105) time, Z(1), (X(K), K=1,KK)
    endif
    
    write(*,*)'Z(1) = ',Z(1),' for numbDisCount = ',numbDisCount,' time = ',time,' dtime = ',dtime
enddo

    write(*,*)'-----------------------------All END-----------------------------'


7100 FORMAT (8X, 'T(SEC)', 10X, 'TMP(K)', 6X, 200(1X,A19))
7105 FORMAT (200E20.10)
7130 FORMAT (8X, 'T(SEC)', 10X, 'TMP(K)', 10X, 'power', 10X, 'powerSUM', 10X, 'energy', 10X, 'energySUM')
7135 FORMAT (20E16.8)
END program discomb
!
!
SUBROUTINE get_EvN (TIME, EvN)
!
    IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N) 
    COMMON /ICONS/ KK, NP, NWT, NH, NWDOT, MM
    double precision :: EvN, disInterval
    
    real, parameter    :: frequency = 100.D4 ! 100 kHz
    real, parameter    :: disPulseWidth = 10.D-9 ! 10 ns
    real, parameter    :: EvNmax = 100.d0 ! 100 Td
    integer, parameter   :: ipulsesTot = 200
    
    ipulse = FLOOR(real(TIME)*frequency)
    if (ipulse < ipulsesTot .AND. (real(TIME)-ipulse/frequency) <= disPulseWidth) then
        EvN = EvNmax
    else
        EvN = 0.d0
    endif
	
    RETURN
END

SUBROUTINE FUN (N, TIME, Z, ZP, RPAR, IPAR)
!
use ZDPlasKin
    IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N) 
    COMMON /ICONS/ KK, NP, NWT, NH, NWDOT, MM
    CHARACTER*16   KSYM(KK), specName(200)
    LOGICAL  ZD_init 
    COMMON /ICONS_specName/ specName
    COMMON /ICONS_ZD_init/ ZD_init

    DIMENSION Z(*), ZP(*), RPAR(*), IPAR(*), X(KK)
    double precision :: density_my(species_max), ZDrates(species_max), ZDsource,&
                        ZDsource_T,ND_E,elasPower, EvN

    if (.NOT. ZD_init) then
        call ZDPlasKin_init()
        ZD_init = .true.
    endif

!    call ZDPlasKin_reset()

    call get_EvN (TIME, EvN)
	
	if( EvN > 1e-10 ) then
		call ZDPlasKin_set_conditions(REDUCED_FIELD=EvN)
		call ZDPlasKin_set_density('E',LDENS_CONST=.true.)
		
		gas_density = RPAR(NP) * 6.02D23 / (Z(1) * RU)

		call CKRHOY(RPAR(NP),Z(1),Z(2),IPAR,RPAR,RHO)
		DO k = 1,species_max
			! CALL CKSNUM (species_name(k), 1, LOUT, specName, KK, KNUM, NVAL, VAL, IERR)
			density_my(k)=Z(k+1)*RHO*6.02d23/RPAR(k+MM)
			call ZDPlasKin_set_density(trim(species_name(k)),density_my(k))
		ENDDO
		call ZDPlasKin_set_conditions(GAS_TEMPERATURE=Z(1))

		CALL ZDPlasKin_get_conditions(ELEC_DRIFT_VELOCITY = eleDiffV)
		CALL ZDPlasKin_get_density_total(ALL_SPECIES = densityTot)
		disPower = 2.2E10
		density_E = disPower/(1.60217662d-19 * eleDiffV * EvN * 1.D-17 * densityTot*1E6)
		call ZDPlasKin_set_density("e",density_E)
		
		! call ZDPlasKin_get_rates(SOURCE_TERMS_MATRIX = ZDrates(1:species_max,1:reactions_max))
		! call ZDPlasKin_get_rates(MEAN_SOURCE_TERMS = ZDrates(1:species_max))
		call ZDPlasKin_get_rates(SOURCE_TERMS = ZDrates(1:species_max))
		
		DO k = 1,species_max
			! CALL CKSNUM (species_name(k), 1, LOUT, specName, KK, KNUM, NVAL, VAL, IERR)
			ZDrates(k) = ZDrates(k) / RHO / 6.02d23 * RPAR(k+MM) ! 1/s
		ENDDO
		
		call ZDPlasKin_get_conditions(ELEC_POWER_ELASTIC_N = elasPower)
		call ZDPlasKin_get_density('E',ND_E)
		ZDsource_T = eV_to_K * elasPower * ND_E * ZDPlasKin_cfg(13) ! K/s

		! call ZDPlasKin_get_conditions(ELEC_POWER_N = dis_Power)
		! call ZDPlasKin_get_density_total(ALL_SPECIES = densityTot)
		! dis_Power = dis_Power* densityTot * 1.60217662d-19* density(species_electrons)* 1E6 ! J/(m3*s)

    endif
    if ( EvN < 1E-10 ) then
        ZDrates = 0.d0
    endif

    CALL CKRHOY (RPAR(NP), Z(1), Z(2), IPAR, RPAR, RHO)
    CALL CKCPBS (Z(1), Z(2), IPAR, RPAR, CPB)
    CALL CKWYP  (RPAR(NP), Z(1), Z(2), IPAR, RPAR, RPAR(NWDOT))
    CALL CKHMS  (Z(1), IPAR, RPAR, RPAR(NH))

!     Form governing equation
!
!     Variables in Z are:  Z(1)   = T
!                          Z(K+1) = Y(K)
!     Call CHEMKIN subroutines

    SUM = 0.0
    
    DO K = 1, KK
        ! CALL CKSNUM (species_name(k), 1, LOUT, specName, KK, KNUM, NVAL, VAL, IERR)
        ZDsource = ZDrates(k)
        
        H    = RPAR(NH    + k - 1) ! erg/g
        WDOT = RPAR(NWDOT + k - 1) ! mole/(cm3*s)
        WT   = RPAR(NWT   + k - 1) ! g/mol
        ZP(k+1) = WDOT * WT / RHO + ZDsource
        SUM = SUM + H * WDOT * WT
    ENDDO
    
!    ZP(1) = -SUM / (RHO*CPB) + ZDsource_T
    ZP(1) = -SUM / (RHO*CPB)
!    ZP(2) = 0.d0
    RETURN
END