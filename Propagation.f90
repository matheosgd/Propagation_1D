MODULE Propagation_m
  USE QDUtil_m 
  USE Lanczos_m 
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Propa_t, propagation

  TYPE :: Propa_t   
    real(kind=Rkind)  :: t0
    real(kind=Rkind)  :: tf   
    real(kind=Rkind)  :: dt  
    character(len=10) :: Propa_name
    real(kind=Rkind)  :: epsilon
    integer           :: niopsit
    integer           :: nioprop
  END TYPE



CONTAINS



  SUBROUTINE read_propagation_parameter(propa, nio)          !nio le fichier dans lequel prendre les valeurs
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(Propa_t), intent(inout) :: propa
    integer, intent(in)          :: nio
    real(kind=Rkind)             :: t0, tf, dt, epsilon
    character(len=10)            :: name_propa
    integer                      :: err_io, niopsit, nioprop
    
    NAMELIST /PROPAGATION_PARAMETER/ t0, tf, dt, name_propa, epsilon

    OPEN(newunit = niopsit, file = 'psit')
    OPEN(newunit = nioprop, file = 'prop')

!Initialisation valeurs par défaut
    t0 = ZERO
    tf = ONE
    dt = ONETENTH
    name_propa = 'spectral'
    epsilon = ZERO

!Lecture
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) 'PROPAGATION_PARAMETER :'
    READ(nio, nml = PROPAGATION_PARAMETER, iostat = err_io)

    WRITE(out_unit, nml = PROPAGATION_PARAMETER)
    
!Contrôle des erreures de lectures
    IF(err_io < 0) THEN
      WRITE(out_unit, *) 'err in READ PROPAGATION_PARAMETER'
      WRITE(out_unit, *) 'err_io = ', err_io
      STOP 'check your propagation parameters'
    ELSE IF( err_io > 0) THEN
      WRITE(out_unit, *) 'err in READ PROPAGATION_PARAMETER'
      WRITE(out_unit, *) 'err_io = ', err_io
      STOP 'check propagation parameters'
    END IF

!Assignation valeurs nml
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



  SUBROUTINE marche_spec(psidt, psi, REigVal, REigVec , t)
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind = Rkind), intent(in)    :: psi(:)
    real(kind = Rkind), intent(in)       :: REigVal(:), REigVec(:,:)  
    real(kind = Rkind), intent(in)       :: t
    complex(kind = Rkind), intent(inout) :: psidt(:)
    complex(kind = Rkind), allocatable   :: c(:) !c(i) = <REigvec(:,i)|psi0(:)>

    psidt(:) = CZERO
    ALLOCATE(c(size(psi)))

    c(:) = matmul(transpose(REigVec), psi)*exp(-EYE*REigVal(:)*t)
   
    psidt = matmul(REigVec, c)
    
    DEALLOCATE(c)
  END SUBROUTINE



  SUBROUTINE marche_euler(psidt, psi, H, dt)
    USE Psi_m
    USE QDUtil_m
    USE Operateurs_m
    IMPLICIT NONE

    complex(kind = Rkind), intent(inout) :: psidt(:)
    complex(kind = Rkind), intent(in)    :: psi(:)
    real(kind = Rkind), intent(in)       :: H(:,:)
    real(kind = Rkind), intent(in)       :: dt
    complex(kind = Rkind), allocatable   :: psi_prime(:)

    psidt(:) = CZERO
    ALLOCATE(psi_prime(size(psi)))

    CALL derivee_t(psi_prime, psi, H)

    psidt(:) = psi(:) + psi_prime(:)*dt

    DEALLOCATE(psi_prime)
  END SUBROUTINE



  SUBROUTINE marche_Taylor_4(psidt, psi, H, dt)
    USE Psi_m
    USE QDUtil_m
    USE Operateurs_m
    IMPLICIT NONE

    complex(kind = Rkind), intent(inout) :: psidt(:)
    complex(kind = Rkind), intent(in)    :: psi(:)
    real(kind = Rkind), intent(in)       :: H(:,:)
    real(kind = Rkind), intent(in)       :: dt
    integer                              :: fact, i
    complex(kind = Rkind), allocatable   :: derivee_a(:), derivee_b(:)
     
    ALLOCATE(derivee_a(size(psi)))
    derivee_a(:) = psi(:)
    psidt(:) = psi(:)
    fact = 1
    
    DO i = 1, 4, 1
      ALLOCATE(derivee_b(size(psi)))
      fact = i*fact
      
      CALL derivee_t(derivee_b, derivee_a, H) 
      psidt(:) = psidt(:) + (derivee_b(:)*dt**i)/fact
      derivee_a(:) = derivee_b(:)
      
      DEALLOCATE(derivee_b)
    END DO

    DEALLOCATE(derivee_a)
  END SUBROUTINE



  SUBROUTINE marche_RK4(psidt, psi, H, dt)
    USE Psi_m
    USE QDUtil_m
    USE Operateurs_m
    IMPLICIT NONE

    complex(kind = Rkind), intent(in)    :: psi(:)
    real(kind = Rkind), intent(in)       :: H(:,:)
    real(kind = Rkind), intent(in)       :: dt
    complex(kind = Rkind), intent(inout) :: psidt(:)
    complex(kind = Rkind), allocatable   :: k1(:), k2(:), k3(:), k4(:), psi_a(:)
    integer                               :: nb
    
    nb = size(psi(:))
    ALLOCATE(K1(nb), K2(nb), K3(nb), K4(nb))
    
    ALLOCATE(psi_a(nb))

    CALL derivee_t(k1, psi, H)
    psi_a(:) = psi(:) + k1(:)*dt*HALF
    CALL derivee_t(k2, psi_a, H)
    psi_a(:) = psi(:) + k2(:)*dt*HALF
    CALL derivee_t(k3, psi, H)
    psi_a(:) = psi(:) + k3(:)*dt
    CALL derivee_t(k4, psi, H)
    psidt(:) = psi(:) + (dt/SIX)*(k1(:) + TWO*k2(:) + TWO*k3(:) + k4(:))
    
    DEALLOCATE(k1, k2, k3, k4, psi_a)  
  END SUBROUTINE    
    

    
  SUBROUTINE marche_SIL(psidt, psi, H, Propa) !SIL = Short Iterative Lanczos
    USE QDUtil_m
    USE Lanczos_m
    USE Operateurs_m
    IMPLICIT NONE

    complex(kind = Rkind), intent(inout) :: psidt(:)
    complex(kind = Rkind), intent(in)    :: psi(:)
    real(kind = Rkind), intent(in)       :: H(:,:)
    TYPE(Propa_t), intent(in)            :: Propa 

!Variables locales----------------------------------------------------------------
    complex(kind = Rkind), allocatable   :: K(:,:), triband_H(:,:), psi1(:), psi2(:), psi_diff(:)
    complex(kind = Rkind), allocatable   :: Valp(:), Vecp(:,:), Vec_B(:,:)
    real(kind=Rkind), allocatable        :: ValpH(:), VecpH(:,:)
    integer                              :: m, nb
    real(kind=Rkind)                     :: eps_temp, E0, E1, E2, E3

    psidt(:) = CZERO  
    nb = size(psi)

!Premiere itération---------------------------------------------------    
    m = 5
    ALLOCATE(Valp(m), Vecp(m,m), ValpH(nb), VecpH(nb,nb), triband_H(m,m), Vec_B(nb, m))
    ALLOCATE(psi1(nb), psi2(nb), psi_diff(nb))

    CALL Base_Krylov(K, H, psi, m)
    CALL Construct_H_tribande(triband_H, K, H)
    CALL diagonalization(triband_H, Valp, Vecp) 
    CALL Krylov_TO_Basis(Vec_B, Vecp, K)
    CALL construct_psi_approx(psi1, psi, Vec_B, Valp, Propa%dt)
   
    DEALLOCATE(Valp, Vecp, triband_H, Vec_B)

    DO m = 6, nb
      WRITE(out_unit,*) 
      WRITE(out_unit,*) '**********************************************************&
      &*********************************************************',m,'ieme iteration'
      ALLOCATE(Valp(m), Vecp(m,m), triband_H(m,m), Vec_B(nb, m))

      CALL Base_Krylov_plus(K, H)
      CALL Construct_H_tribande(triband_H, K, H)
      CALL diagonalization(triband_H, Valp, Vecp)  
      CALL Krylov_TO_Basis(Vec_B, Vecp, K)
      
      WRITE(out_unit,*) 
      CALL calc_energie(E0, H, Vec_B(:,1))
      CALL calc_energie(E1, H, CMPLX(vecpH(:,1), kind = Rkind))
      WRITE(out_unit,*) '<Vec_B(:,1)|H|Vec_B(:,1)> = ', E0, 'ET Valp(1) = ', Valp(1), '<VecpH(:,1)|H|VecpH(:,1)>', E1

      CALL construct_psi_approx(psi2, psi, Vec_B, Valp, Propa%dt)
      
      psi_diff(:) = psi2(:) - psi1(:)
      eps_temp = sqrt(dot_product(psi_diff, psi_diff))
      WRITE(out_unit,*) 
      WRITE(out_unit,*) 'eps_temp ', sqrt(real(dot_product(psi_diff, psi_diff), kind=Rkind)),&
      & sqrt(real(dot_product(psi1, psi1), kind=Rkind)), sqrt(real(dot_product(psi2, psi2), kind=Rkind))

      CALL calc_energie(E2, H, psi1)
      CALL calc_energie(E3, H, psi2)
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



  SUBROUTINE propagation(psif, psi0, H, Basis, Op, nio)
    USE Psi_m
    USE Basis_m
    USE Operateurs_m
    USE QDUtil_m
    USE Lanczos_m
    IMPLICIT NONE

    complex(kind = Rkind), intent(inout) :: psif(:) 
    real(kind = Rkind), intent(in)       :: H(:,:)
    TYPE(Basis_t), intent(in)            :: Basis
    TYPE(Operateur_t)                    :: Op
    complex(kind = Rkind), intent(in)    :: psi0(:)          !defini sur la base des sinus
    integer, intent(in)                  :: nio
!Variables locales---------------------------------------------------------
    TYPE(Propa_t)                        :: propa
    integer                              :: nt, it, ib
    complex(kind=Rkind), allocatable     :: mat_X(:,:), mat_Px(:,:) 
    complex(kind = Rkind), allocatable   :: psidt(:), psi(:)
    real(kind = Rkind)                   :: t, E, N, X0, P0, A_t

!Initialisation------------------------------------------------------------
    CALL read_propagation_parameter(propa, nio)
    
    nt = (propa%tf - propa%t0) / Propa%dt

    CALL Construct_Px(mat_Px, Basis) 
    CALL construct_matX(mat_X, Basis)
    ALLOCATE(psi(Basis%nb))
    ALLOCATE(psidt(Basis%nb))

    psidt(:) = CZERO
    psi(:)   = CZERO

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

!Calcul t = t0 : état initial----------------------------------------------
    CALL calc_energie(E, H, psi)   
    CALL calc_norme_basis(N, psi)  
    CALL position_moyenne(X0, mat_X, psi)
    CALL impulsion_moyenne(P0, mat_Px, psi) 
    CALL autocorrelation(A_t, psi, psi0)
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

!Propagation itérative-----------------------------------------------------
    DO it = 1, nt, 1  
      t = propa%t0 + propa%dt*it  
      
      CALL marche_general(psidt, psi, H, Propa, Op%REigVal, Op%REigVec, psi0, t)

      psi(:) = psidt(:)

      CALL calc_energie(E, H, psi)
      CALL calc_norme_basis(N, psi)
      CALL position_moyenne(X0, mat_X, psi)
      CALL impulsion_moyenne(P0, mat_Px, psi) 
      CALL autocorrelation(A_t, psi, psi0)   
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



  SUBROUTINE marche_general(psidt, psi, H, Propa, REigVal, REigVec, psi0, t) 
    USE QDUtil_m
    USE Lanczos_m
    IMPLICIT NONE

    complex(kind = Rkind), intent(inout) :: psidt(:)
    complex(kind = Rkind), intent(in)    :: psi(:)
    real(kind = Rkind), intent(in)       :: H(:,:), t, REigVal(:), REigVec(:,:)
    TYPE(Propa_t), intent(in)            :: Propa
    complex(kind = Rkind), intent(in)    :: psi0(:)          !defini sur la base des sinus
 !Variables locales
    integer                              :: nb

    SELECT CASE (propa%propa_name)

    CASE ('spectral')
      CALL marche_spec(psidt, psi0, REigVal, REigVec, t)
      
    CASE ('euler')
      CALL marche_euler(psidt, psi, H, Propa%dt)

    CASE ('Taylor_4')
      CALL marche_Taylor_4(psidt, psi, H, Propa%dt)

    CASE ('RK4')
      CALL marche_RK4(psidt, psi, H, Propa%dt)  

    CASE ('SIL')
      CALL marche_SIL(psidt, psi, H, Propa)    

    CASE DEFAULT
      STOP 'marche name unknown'

  END SELECT

  END SUBROUTINE 



  SUBROUTINE autocorrelation(A_t, psit, psi0)
    USE QDUtil_m
    IMPLICIT NONE
  
    real(kind = Rkind), intent(inout)              :: A_t
    complex(kind = Rkind), intent(in)              :: psit(:), psi0(:)
  
    A_t = dot_product(psi0, psit)
  
  END SUBROUTINE
  
  

END MODULE