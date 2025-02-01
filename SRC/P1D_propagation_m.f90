MODULE P1D_propagation_m
  USE QDUtil_m 
  USE P1D_lanczos_m 
  IMPLICIT NONE


  PRIVATE
  PUBLIC :: P1D_propagation_t, P1D_Propagation


  TYPE :: P1D_propagation_t   
    real(kind=Rkind)  :: t0
    real(kind=Rkind)  :: tf   
    real(kind=Rkind)  :: dt  
    character(len=10) :: Propa_name
    real(kind=Rkind)  :: epsilon
    integer           :: niopsit
    integer           :: nioprop

  END TYPE


  CONTAINS


  SUBROUTINE P1D_Read_propagation_parameter(propa, nio, args)                       ! nio le fichier dans lequel prendre les valeurs
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(P1D_propagation_t), intent(inout) :: propa
    integer,                 intent(in)    :: nio
    character(len=12),       intent(in)    :: args(:)

    real(kind=Rkind)                       :: t0, tf, dt, epsilon
    character(len=10)                      :: name_propa
    integer                                :: err_io, niopsit, nioprop
    
    NAMELIST /PROPAGATION_PARAMETER/ t0, tf, dt, name_propa, epsilon                ! indique la liste des éléments qu'il va trouver dans nml (dans l'odre) et déclare une namelist

    OPEN(NEWUNIT = niopsit, FILE = 'OUT/P1D_App_Propagation_1D_'//trim(args(1))//'_psit.log', &
       & FORM = 'formatted', ACTION = 'write', POSITION = 'rewind')
    OPEN(NEWUNIT = nioprop, FILE = 'OUT/P1D_App_Propagation_1D_'//trim(args(1))//'_propagation.log', &
       & FORM = 'formatted', ACTION = 'write', POSITION = 'rewind')

    !---------------------Initialization to default values---------------------
    t0 = ZERO
    tf = ONE
    dt = ONETENTH
    name_propa = 'spectral'
    epsilon = ZERO

    !----------------------------Reading of the nml----------------------------
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) 'PROPAGATION_PARAMETER :'
    READ(nio, nml = PROPAGATION_PARAMETER, iostat = err_io)                         ! assigne à la liste annoncée les valeurs lues dans la nml du fichier nio

    WRITE(out_unit, nml = PROPAGATION_PARAMETER)
    
    !---------------------Control of potential reading errors------------------
    IF(err_io < 0) THEN
      WRITE(out_unit, *) 'err in READ PROPAGATION_PARAMETER'
      WRITE(out_unit, *) 'err_io = ', err_io
      STOP 'check your propagation parameters'
    ELSE IF( err_io > 0) THEN
      WRITE(out_unit, *) 'err in READ PROPAGATION_PARAMETER'
      WRITE(out_unit, *) 'err_io = ', err_io
      STOP 'check propagation parameters'
    END IF

    !--Affectation of the read values to those of the P1D_basis_t type object-
    Propa%t0 = t0
    Propa%tf = tf
    Propa%dt = dt
    Propa%propa_name = trim(name_propa)
    Propa%niopsit = niopsit
    Propa%nioprop = nioprop

    IF(Propa%propa_name == 'SIL') THEN
      Propa%epsilon = epsilon
    END IF   

  END SUBROUTINE


  SUBROUTINE P1D_Marche_spec(psidt, psi, REigVal, REigVec , t)
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind = Rkind), intent(inout) :: psidt(:)
    complex(kind = Rkind), intent(in)    :: psi(:)
    real(kind = Rkind),    intent(in)    :: REigVal(:), REigVec(:,:)  
    real(kind = Rkind),    intent(in)    :: t

    complex(kind = Rkind), allocatable   :: c(:)                                    ! ci = <REigvec(:,i)|psi0(:)>

    psidt(:) = CZERO                                                                ! initialise vecteur nul
    ALLOCATE(c(size(psi)))

    c(:) = matmul(transpose(REigVec), psi)*exp(-EYE*REigVal(:)*t)
   
    psidt = matmul(REigVec, c)
    
    DEALLOCATE(c)

  END SUBROUTINE


  SUBROUTINE P1D_Marche_euler(psidt, psi, H, dt)
    USE P1D_wavepacket_m
    USE QDUtil_m
    USE P1D_operators_m
    IMPLICIT NONE

    complex(kind = Rkind), intent(inout) :: psidt(:)
    complex(kind = Rkind), intent(in)    :: psi(:)
    real(kind = Rkind),    intent(in)    :: H(:,:)
    real(kind = Rkind),    intent(in)    :: dt

    complex(kind = Rkind), allocatable   :: psi_prime(:)

    psidt(:) = CZERO
    ALLOCATE(psi_prime(size(psi)))

    CALL P1D_Derivee_t(psi_prime, psi, H)

    psidt(:) = psi(:) + psi_prime(:)*dt

    DEALLOCATE(psi_prime)

  END SUBROUTINE


  SUBROUTINE P1D_Marche_taylor_4(psidt, psi, H, dt)
    USE P1D_wavepacket_m
    USE QDUtil_m
    USE P1D_operators_m
    IMPLICIT NONE

    complex(kind = Rkind), intent(inout) :: psidt(:)
    complex(kind = Rkind), intent(in)    :: psi(:)
    real(kind = Rkind),    intent(in)    :: H(:,:)
    real(kind = Rkind),    intent(in)    :: dt

    complex(kind = Rkind), allocatable   :: derivee_a(:), derivee_b(:)
    integer                              :: fact, i
     
    ALLOCATE(derivee_a(size(psi)))
    derivee_a(:) = psi(:)
    psidt(:) = psi(:)
    fact = 1
    
    DO i = 1, 4, 1                                                                  ! ajoute termes d'ordre 1 à 4
      ALLOCATE(derivee_b(size(psi)))
      fact = i*fact
      
      CALL P1D_Derivee_t(derivee_b, derivee_a, H) 
      psidt(:) = psidt(:) + (derivee_b(:)*dt**i)/fact
      derivee_a(:) = derivee_b(:)
      
      DEALLOCATE(derivee_b)
    END DO

    DEALLOCATE(derivee_a)

  END SUBROUTINE


  SUBROUTINE P1D_Marche_RK4(psidt, psi, H, dt)
    USE P1D_wavepacket_m
    USE QDUtil_m
    USE P1D_operators_m
    IMPLICIT NONE

    complex(kind = Rkind), intent(inout) :: psidt(:)
    complex(kind = Rkind), intent(in)    :: psi(:)
    real(kind = Rkind),    intent(in)    :: H(:,:)
    real(kind = Rkind),    intent(in)    :: dt

    complex(kind = Rkind), allocatable   :: k1(:), k2(:), k3(:), k4(:), psi_a(:)
    integer                              :: nb
    
    nb = size(psi(:))
    ALLOCATE(K1(nb), K2(nb), K3(nb), K4(nb))
    
    ALLOCATE(psi_a(nb))

    CALL P1D_Derivee_t(k1, psi, H)
    psi_a(:) = psi(:) + k1(:)*dt*HALF
    CALL P1D_Derivee_t(k2, psi_a, H)
    psi_a(:) = psi(:) + k2(:)*dt*HALF
    CALL P1D_Derivee_t(k3, psi, H)
    psi_a(:) = psi(:) + k3(:)*dt
    CALL P1D_Derivee_t(k4, psi, H)
    psidt(:) = psi(:) + (dt/SIX)*(k1(:) + TWO*k2(:) + TWO*k3(:) + k4(:))
    
    DEALLOCATE(k1, k2, k3, k4, psi_a) 

  END SUBROUTINE
    
    
  SUBROUTINE P1D_Marche_SIL(psidt, psi, H, Propa)                                   ! SIL = Short Iterative Lanczos
    USE QDUtil_m
    USE P1D_lanczos_m
    USE P1D_operators_m
    IMPLICIT NONE

    complex(kind = Rkind),   intent(inout) :: psidt(:)
    complex(kind = Rkind),   intent(in)    :: psi(:)
    real(kind = Rkind),      intent(in)    :: H(:,:)
    TYPE(P1D_propagation_t), intent(in)    :: Propa 

    complex(kind = Rkind), allocatable     :: K(:,:), triband_H(:,:), psi1(:), psi2(:), psi_diff(:)
    complex(kind = Rkind), allocatable     :: Valp(:), Vecp(:,:), Vec_B(:,:)
    real(kind=Rkind),      allocatable     :: ValpH(:), VecpH(:,:)
    real(kind=Rkind)                       :: eps_temp, E0, E1, E2, E3
    integer                                :: m, nb

    psidt(:) = CZERO  
    nb = size(psi)

    !------------------------------First iteration-----------------------------
    m = 5
    ALLOCATE(Valp(m), Vecp(m,m), ValpH(nb), VecpH(nb,nb), triband_H(m,m), Vec_B(nb, m))
    ALLOCATE(psi1(nb), psi2(nb), psi_diff(nb))

    CALL P1D_Base_krylov(K, H, psi, m)
    CALL P1D_Construct_H_tribande(triband_H, K, H)
    CALL diagonalization(triband_H, Valp, Vecp) 
    CALL P1D_KrylovTOBasis(Vec_B, Vecp, K)
    CALL P1D_Construct_psi_approx(psi1, psi, Vec_B, Valp, Propa%dt)
   
    DEALLOCATE(Valp, Vecp, triband_H, Vec_B)

    DO m = 6, nb
      WRITE(out_unit,*) 
      WRITE(out_unit,*) '**********************************************************&
                       & *********************************************************', &
                       & m,'ieme iteration'
      ALLOCATE(Valp(m), Vecp(m,m), triband_H(m,m), Vec_B(nb, m))

      CALL P1D_Base_krylov_plus(K, H)
      CALL P1D_Construct_H_tribande(triband_H, K, H)
      CALL diagonalization(triband_H, Valp, Vecp)  
      CALL P1D_KrylovTOBasis(Vec_B, Vecp, K)
      
      WRITE(out_unit,*) 
      CALL P1D_Calc_energie(E0, H, Vec_B(:,1))
      CALL P1D_Calc_energie(E1, H, CMPLX(vecpH(:,1), kind = Rkind))
      WRITE(out_unit,*) '<Vec_B(:,1)|H|Vec_B(:,1)> = ', E0, 'ET Valp(1) = ', Valp(1), '<VecpH(:,1)|H|VecpH(:,1)>', E1

      CALL P1D_Construct_psi_approx(psi2, psi, Vec_B, Valp, Propa%dt)
      
      psi_diff(:) = psi2(:) - psi1(:)
      eps_temp = sqrt(dot_product(psi_diff, psi_diff))
      WRITE(out_unit,*) 
      WRITE(out_unit,*) 'eps_temp ', sqrt(real(dot_product(psi_diff, psi_diff), kind=Rkind)), &
                                   & sqrt(real(dot_product(psi1, psi1), kind=Rkind)), &
                                   & sqrt(real(dot_product(psi2, psi2), kind=Rkind))

      CALL P1D_Calc_energie(E2, H, psi1)
      CALL P1D_Calc_energie(E3, H, psi2)
      WRITE(out_unit,*) '<psi1|H|psi1> = ', E2, '<psi2|H|psi2> = ', E3


      IF(eps_temp < Propa%epsilon )  THEN 
        psidt(:) = psi2(:)
        WRITE(out_unit,*) 
        WRITE(out_unit,*) 'La précision est atteinte pour m = ', m
        EXIT

      ELSE 
        psi1(:) = psi2(:)
      END IF 

      DEALLOCATE(Valp, Vecp, triband_H, Vec_B)
    END DO 

  END SUBROUTINE


  SUBROUTINE P1D_Propagation(psif, psi0, H, Basis, Op, nio, args)
    USE P1D_wavepacket_m
    USE P1D_basis_m
    USE P1D_operators_m
    USE QDUtil_m
    USE P1D_lanczos_m
    IMPLICIT NONE

    complex(kind = Rkind), intent(inout) :: psif(:) 
    real(kind = Rkind),    intent(in)    :: H(:,:)
    TYPE(P1D_basis_t),     intent(in)    :: Basis
    complex(kind = Rkind), intent(in)    :: psi0(:)                                 ! def sur la base des sinus
    character(len=12),     intent(in)    :: args(:)
    integer,               intent(in)    :: nio

    TYPE(P1D_operators_t)                :: Op
    TYPE(P1D_propagation_t)              :: propa
    complex(kind=Rkind),   allocatable   :: mat_X(:,:), mat_Px(:,:) 
    complex(kind = Rkind), allocatable   :: psidt(:), psi(:)                        ! psidt la valeur à psi(t+dt) et psi sert à la stocker pour l'iteration suivante
    real(kind = Rkind)                   :: t, E, N, X0, P0, A_t
    integer                              :: nt, it, ib

    !------------------------------Initialization------------------------------
    CALL P1D_Read_propagation_parameter(propa, nio, args)
    
    nt = (propa%tf - propa%t0) / Propa%dt

    CALL P1D_Construct_Px(mat_Px, Basis) 
    CALL P1D_Construct_matX(mat_X, Basis)
    ALLOCATE(psi(Basis%nb))
    ALLOCATE(psidt(Basis%nb))

    psidt(:) = CZERO                                                                ! initialise vecteurs nuls
    psi(:)   = CZERO                                                                ! initialise vecteurs nuls

    psi(:) = psi0(:)

    IF(Propa%propa_name == 'spectral') THEN
      IF(allocated(Op%Reigval)) THEN
        WRITE(out_unit,*) 'REigval previously allocated'
      ELSE
        ALLOCATE(Op%REigVal(Basis%nb))
        ALLOCATE(Op%REigVec(Basis%nb,Basis%nb))
        CALL diagonalization(H,Op%REigVal,Op%REigVec)
      END IF
    END IF

    !--------------------Compute for t = t0 : initial state--------------------
    CALL P1D_Calc_energie(E, H, psi)   
    CALL P1D_Calc_norme_basis(N, psi)  
    CALL P1D_Position_moyenne(X0, mat_X, psi)
    CALL P1D_Impulsion_moyenne(P0, mat_Px, psi) 
    CALL P1D_Autocorrelation(A_t, psi, psi0)
    WRITE(out_unit,*) 'POSITION MOYENNE INITIALE : <Psi|X|Psi> = ', X0
    WRITE(out_unit,*) 'IMPULSION MOYENNE INITIALE: <Psi|Px|Psi> = ', P0
    WRITE(out_unit,*) 'ENERGIE MOYENNE INITIALE: <Psi|H|Psi>[Ha] = ', E
    WRITE(out_unit,*) 'ENERGIE MOYENNE INITIALE: <Psi|H|Psi>[cm-1] = ', E*219474.63_Rkind
    WRITE(out_unit,*) 'ENERGIE MOYENNE INITIALE: <Psi|H|Psi>[eV] = ', E*27.2114_Rkind
    WRITE(propa%nioprop,*) '--t-------------------------POSITION MOYENNE----------IMPULSION MOYENNE-----&
    &----ENERGIE MOYENNE----------NORME--------------------<Psi(0)|Psi(t)>' 
    WRITE(propa%nioprop,*) Propa%t0, X0, P0, E, N, A_t
    WRITE(propa%niopsit,*) '---------ib-----------------------------------psi_ib(t)'
    WRITE(propa%niopsit,*) 'psi_b(t  = ', propa%t0, ')'
    WRITE(out_unit,*) 'Psi0_b :'
    WRITE(out_unit,*) '.........ib,................Psi0(ib)'
    DO ib = 1, Basis%nb
      WRITE(out_unit,*) ib, psi0(ib)
      WRITE(propa%niopsit,*) ib, psi0(ib)
    END DO
    WRITE(propa%niopsit,*)

    !---------------------------Iterative propagation--------------------------
    DO it = 1, nt, 1  
      t = propa%t0 + propa%dt*it  
      
      CALL P1D_Marche_general(psidt, psi, H, Propa, Op%REigVal, Op%REigVec, psi0, t)

      psi(:) = psidt(:)

      CALL P1D_Calc_energie(E, H, psi)                                              ! vérifie que l'énergie se conserve
      CALL P1D_Calc_norme_basis(N, psi)                                             ! vérifie que la norme se conserve
      CALL P1D_Position_moyenne(X0, mat_X, psi)
      CALL P1D_Impulsion_moyenne(P0, mat_Px, psi) 
      CALL P1D_Autocorrelation(A_t, psi, psi0)   
      WRITE(propa%nioprop,*) t, X0, P0, E, N, A_t
      WRITE(propa%niopsit,*) '*************************************************'
      WRITE(propa%niopsit,*) 'psi_b(t  = ', t, ')'
      DO ib = 1, Basis%nb
        WRITE(propa%niopsit,*) ib, psi(ib)
      END DO
      WRITE(propa%niopsit,*)

    END DO

    psif(:) = psi(:)

    IF(Propa%propa_name == 'spectral') THEN
      DEALLOCATE(Op%REigVal, Op%REigVec)
    END IF
    DEALLOCATE(psidt, psi) 

  END SUBROUTINE


  SUBROUTINE P1D_Marche_general(psidt, psi, H, Propa, REigVal, REigVec, psi0, t) 
    USE QDUtil_m
    USE P1D_lanczos_m
    IMPLICIT NONE

    complex(kind = Rkind),   intent(inout) :: psidt(:)
    complex(kind = Rkind),   intent(in)    :: psi(:)
    real(kind = Rkind),      intent(in)    :: H(:,:), t, REigVal(:), REigVec(:,:)
    TYPE(P1D_propagation_t), intent(in)    :: Propa
    complex(kind = Rkind),   intent(in)    :: psi0(:)                               ! def sur la base des sinus

    integer                                :: nb

    SELECT CASE (propa%propa_name)

    CASE ('spectral')
      CALL P1D_Marche_spec(psidt, psi0, REigVal, REigVec, t)
      
    CASE ('euler')
      CALL P1D_Marche_euler(psidt, psi, H, Propa%dt)

    CASE ('Taylor_4')
      CALL P1D_Marche_taylor_4(psidt, psi, H, Propa%dt)

    CASE ('RK4')
      CALL P1D_Marche_RK4(psidt, psi, H, Propa%dt)  

    CASE ('SIL')
      CALL P1D_Marche_SIL(psidt, psi, H, Propa)    

    CASE DEFAULT
      STOP 'marche name unknown'

    END SELECT

  END SUBROUTINE


  SUBROUTINE P1D_Autocorrelation(A_t, psit, psi0)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind = Rkind),    intent(inout) :: A_t
    complex(kind = Rkind), intent(in)    :: psit(:), psi0(:)
  
    A_t = dot_product(psi0, psit)
  
  END SUBROUTINE
  

END MODULE