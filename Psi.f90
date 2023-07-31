MODULE Psi_m
  USE QDUtil_m
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: psi_t, init_psi, derivee_t, calc_Norme_basis

  TYPE :: psi_t   
    real(kind=Rkind) :: sigma
    real(kind=Rkind) :: Q0
    real(kind=Rkind) :: P0  
    real(kind=Rkind) :: phasis
  END TYPE



CONTAINS


  SUBROUTINE read_psi(psi_parameters, nio)            !nio le fichier dans lequel prendre les valeurs
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(Psi_t),    intent(inout)   :: psi_parameters
    integer, intent(in)             :: nio
    real(kind=Rkind)                :: sigma, Q0, P0, phasis
    integer                         :: err_io

    NAMELIST /PSI0_PARAMETERS/ sigma, Q0, P0, phasis

!Initialisation valeurs par défaut
    sigma = sqrt(TWO)
    Q0 = ZERO
    P0 = ZERO
    phasis = ZERO  

!Lecture
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) 'CONSTRUCT PSI0_PARAMETERS :'
    READ(nio, nml = PSI0_PARAMETERS, iostat = err_io)
    WRITE(out_unit, nml = PSI0_PARAMETERS)

!Contrôle des erreures de lectures
    IF(err_io < 0) THEN
      WRITE(*, *) 'err in read PSI0_PARAMETERS'
      STOP 'check psi0 parameters data'
    ELSE IF( err_io > 0) THEN
      WRITE(*, *) 'err in read PSI0_PARAMETERS'
      STOP 'check psi0 parameters data'
    END IF

!Assignation valeurs nml
    psi_parameters%sigma = sigma
    psi_parameters%Q0 = Q0
    psi_parameters%P0 = P0
    psi_parameters%phasis = phasis


  END SUBROUTINE



  SUBROUTINE init_Psi(Psi_b, Basis, nio)
    USE Basis_m
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: Psi_b(:)
    TYPE(Basis_t), intent(in)          :: Basis
    TYPE(Psi_t)                        :: psi0_parameters
    integer, intent(in)                :: nio
    complex(kind=Rkind), allocatable   :: Psi_g(:)
    integer                            :: iq
    real(kind=Rkind)                   :: N_g, N_b, Norme

    CALL read_psi(psi0_parameters, nio)

    Norme = sqrt( sqrt(TWO) / (psi0_parameters%sigma*sqrt(PI)) ) !facteur de normalisation déterminé analytiquement (intégrales de Gausse)

    ALLOCATE(Psi_g(Basis%nq))
    Psi_g(:) = exp(-((Basis%x(:)-psi0_parameters%Q0)/psi0_parameters%sigma )**2 )
    Psi_g(:) = Psi_g(:)*exp(EYE*psi0_parameters%P0*(Basis%x(:)-psi0_parameters%Q0))
    Psi_g(:) = Psi_g(:)*exp(EYE*psi0_parameters%phasis)
    Psi_g(:) = Norme*Psi_g(:)

    N_g = sqrt(sum( abs(Psi_g(:) )**2 * Basis%w(:) ))
    WRITE(*,*) '<Psi_g|Psi_g> = ', N_g

!Normalisation de Psi_g ---------------------------------------------
    Psi_g(:) = Psi_g(:) / N_g
    N_g = sqrt(sum( abs(Psi_g(:) )**2 * Basis%w(:) ))
    WRITE(*,*) 'NORMALISATION : <Psi_g|Psi_g> = ', N_g

!Construction de Psi_b ----------------------------------------------
    CALL GridTOBasis(Psi_b, Psi_g, Basis)
    N_b = sqrt( sum( abs(Psi_b(:) )**2 ) )
    WRITE(*,*) 'PROJECTION SUR LA BASE : <Psi_b|Psi_b> = ', N_b

!Normalisation de Psi_b ---------------------------------------------
    Psi_b(:) =Psi_b(:) / N_b
    N_b = sqrt(sum( abs(Psi_b(:) )**2))
    WRITE(*,*) 'RENORMALISATION : <Psi_b|Psi_b> = ', N_b
     
  deallocate(Psi_g)
  END SUBROUTINE



  SUBROUTINE derivee_t(Psi1_b, Psi_b, H)
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout)              :: Psi1_b(:)
    complex(kind=Rkind), intent(in)                 :: Psi_b(:)
    real(kind = Rkind), intent(in)                  :: H(:,:)

    Psi1_b = -EYE * matmul(H, Psi_b)  ! par l'(ESDT)

  END SUBROUTINE 



  SUBROUTINE calc_Norme_basis(N, psi)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout)              :: N
    complex(kind=Rkind), intent(in)              :: psi(:)


    N = sqrt(dot_product(psi, psi))
    
  END SUBROUTINE 



END MODULE