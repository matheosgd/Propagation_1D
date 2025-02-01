MODULE P1D_wavepacket_m
  USE QDUtil_m
  IMPLICIT NONE


  PRIVATE
  PUBLIC :: P1D_wavepacket_t, P1D_Init_psi, P1D_Derivee_t, P1D_Calc_norme_basis


  TYPE :: P1D_wavepacket_t   
    real(kind=Rkind) :: sigma
    real(kind=Rkind) :: Q0
    real(kind=Rkind) :: P0  
    real(kind=Rkind) :: phasis

  END TYPE


  CONTAINS


  SUBROUTINE P1D_Read_psi(psi_parameters, nio)                                      ! nio le fichier dans lequel prendre les valeurs
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(P1D_wavepacket_t), intent(inout)   :: psi_parameters
    integer,                intent(in)      :: nio

    real(kind=Rkind)                        :: sigma, Q0, P0, phasis
    integer                                 :: err_io

    NAMELIST /PSI0_PARAMETERS/ sigma, Q0, P0, phasis                                ! indique la liste des éléments qu'il va trouver dans nml (dans l'odre) et déclare une namelist
   
    !---------------------Initialization to default values---------------------
    sigma = sqrt(TWO)
    Q0 = ZERO
    P0 = ZERO
    phasis = ZERO  

    !----------------------------Reading of the nml----------------------------
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) 'CONSTRUCT PSI0_PARAMETERS :'
    READ(nio, nml = PSI0_PARAMETERS, iostat = err_io)                               ! assigne à la liste annoncée les valeurs lues dans la nml du fichier nio
    WRITE(out_unit, nml = PSI0_PARAMETERS)

    !---------------------Control of potential reading errors------------------
    IF(err_io < 0) THEN
      WRITE(*, *) 'err in read PSI0_PARAMETERS'
      STOP 'check psi0 parameters data'
    ELSE IF( err_io > 0) THEN
      WRITE(*, *) 'err in read PSI0_PARAMETERS'
      STOP 'check psi0 parameters data'
    END IF

    !--Affectation of the read values to those of the P1D_basis_t type object-
    psi_parameters%sigma = sigma
    psi_parameters%Q0 = Q0
    psi_parameters%P0 = P0
    psi_parameters%phasis = phasis

  END SUBROUTINE


  SUBROUTINE P1D_Init_psi(Psi_b, Basis, nio)                                        ! psi the wavepacket
    USE P1D_basis_m
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: Psi_b(:)
    TYPE(P1D_basis_t),   intent(in)    :: Basis
    integer,             intent(in)    :: nio

    TYPE(P1D_wavepacket_t)             :: psi0_parameters
    complex(kind=Rkind), allocatable   :: Psi_g(:)
    real(kind=Rkind)                   :: N_g, N_b, Norme
    integer                            :: iq


    CALL P1D_Read_psi(psi0_parameters, nio)

    Norme = sqrt( sqrt(TWO) / (psi0_parameters%sigma*sqrt(PI)) )                    ! facteur de normalisation déterminé analytiquement

    ALLOCATE(Psi_g(Basis%nq))
    Psi_g(:) = exp(-((Basis%x(:)-psi0_parameters%Q0)/psi0_parameters%sigma )**2 )
    Psi_g(:) = Psi_g(:)*exp(EYE*psi0_parameters%P0*(Basis%x(:)-psi0_parameters%Q0))
    Psi_g(:) = Psi_g(:)*exp(EYE*psi0_parameters%phasis)
    Psi_g(:) = Norme*Psi_g(:)                                                       ! def sur la grille

    N_g = sqrt(sum( abs(Psi_g(:) )**2 * Basis%w(:) ))                               ! calcule norme Psi sur grille
    WRITE(*,*) '<Psi_g|Psi_g> = ', N_g

  !-----------------------------Psi_g normalisation----------------------------
    Psi_g(:) = Psi_g(:) / N_g                                                       ! normalise Psi sur la grille 
    N_g = sqrt(sum( abs(Psi_g(:) )**2 * Basis%w(:) ))                               ! recalcule norme
    WRITE(*,*) 'NORMALISATION : <Psi_g|Psi_g> = ', N_g

  !-----------------------------Psi_b construction-----------------------------
    CALL P1D_GridTOBasis(Psi_b, Psi_g, Basis)                                       ! projette sur la base
    N_b = sqrt( sum( abs(Psi_b(:) )**2 ) )
    WRITE(*,*) 'PROJECTION SUR LA BASE : <Psi_b|Psi_b> = ', N_b

  !-----------------------------Psi_b normalisation----------------------------
    Psi_b(:) =Psi_b(:) / N_b                                                        ! normalise Psi sur la base
    N_b = sqrt(sum( abs(Psi_b(:) )**2))
    WRITE(*,*) 'RENORMALISATION : <Psi_b|Psi_b> = ', N_b
     
    DEALLOCATE(Psi_g)

  END SUBROUTINE


  SUBROUTINE P1D_Derivee_t(Psi1_b, Psi_b, H)
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: Psi1_b(:)
    complex(kind=Rkind), intent(in)    :: Psi_b(:)
    real(kind = Rkind),  intent(in)    :: H(:,:)

    Psi1_b = -EYE * matmul(H, Psi_b)                                           ! par l'(ESDT)

  END SUBROUTINE


  SUBROUTINE P1D_Calc_norme_basis(N, psi)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: N
    complex(kind=Rkind), intent(in)    :: psi(:)


    N = sqrt(dot_product(psi, psi))
    !WRITE(*,*) 'NORME : <Psi|Psi> = ', N

  END SUBROUTINE


END MODULE