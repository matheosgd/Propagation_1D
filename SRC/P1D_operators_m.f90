MODULE P1D_operators_m
  USE QDUtil_m
  IMPLICIT NONE


  PUBLIC
  PRIVATE :: P1D_Potentiel, P1D_Read_operator_parameter


  TYPE :: P1D_operators_t
    !----------------------------General parameters----------------------------
    character(len=20)               :: pot_name
    real(kind = Rkind)              :: mass
    real(kind = Rkind), allocatable :: V(:)                                         ! grid representation of the potential                          

    !---------------------Parameters for harmonic potential--------------------
    real(kind = Rkind)              :: raideur
    real(kind = Rkind)              :: x0

    !----------------------Parameters for Morse potential----------------------
    real(kind = Rkind)              :: De
    real(kind = Rkind)              :: alpha
    real(kind = Rkind)              :: Re

    !----------------------Parameters for stairs potential----------------------
    real(kind = Rkind)              :: V_0
    real(kind = Rkind)              :: V_1
    real(kind = Rkind)              :: V_2
    real(kind = Rkind)              :: x1
    real(kind = Rkind)              :: x2

    !---------------Parameters for the H spectral representation---------------
    real(kind = Rkind), allocatable :: REigval(:)
    real(kind = Rkind), allocatable :: REigvec(:,:)
    integer                         :: niopot

  END TYPE
  
    
  CONTAINS
 
  
  FUNCTION P1D_Potentiel(x, Op) RESULT(V_x)
    USE QDUtil_m
    IMPLICIT NONE
  
    TYPE(P1D_operators_t), intent(in) :: Op
    real(kind = Rkind)                :: x
    real(kind = Rkind)                :: V_x
  
    SELECT CASE (Op%pot_name)
      CASE ('harmonique')
        v_x = HALF*Op%raideur*(x-Op%x0)**2

      CASE ('Morse')
        v_x = Op%De * (ONE - exp(-Op%alpha * (x - Op%Re)))**2

      CASE ('Double_puit')
        v_x =  (x**2 - ONE)**2

      CASE ('escalier')
        IF(x < Op%x1) THEN
          v_x = Op%V_0
        ELSE IF(x < Op%x2) THEN
          v_x = Op%V_1
        ELSE
          v_x = Op%V_2
        END IF

      CASE DEFAULT
        STOP 'pot_name unknown at evaluating on the grid'

    END SELECT
  
  END


  SUBROUTINE P1D_Read_operator_parameter(Op, nio, args)
    USE QDUtil_m
    IMPLICIT NONE
    
    type(P1D_operators_t), intent(inout) :: Op
    integer,               intent(in)    :: nio
    character(len=12),     intent(in)    :: args(:)

    character(len=20)                    :: pot_name
    integer                              :: niopot, err_io
    
    NAMELIST /POTENTIAL_NAME/ pot_name                                              ! indique la liste des éléments qu'il va trouver dans nml (dans l'odre) et déclare une namelist
    
    OPEN(NEWUNIT = niopot, FILE = 'OUT/P1D_App_Propagation_1D_'//trim(args(1))//'_potential.log', &
       & FORM = 'formatted', ACTION = 'write', POSITION = 'rewind')

    !----------------------Initialization to default value---------------------
    pot_name = '0'

    !----------------------------------Reading---------------------------------
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) 'OPERATEUR_PARAMETER :'
    READ(nio, nml = POTENTIAL_NAME, iostat = err_io)                                ! assigne à la liste annoncée les valeurs lues dans la nml du fichier nio

    WRITE(out_unit, nml = POTENTIAL_NAME)

    !---------------------Control of potential reading errors------------------
    IF(err_io < 0) THEN
      WRITE(out_unit, *) 'err in READ POTENTIAL_NAME'
      STOP 'check your potential name'
    ELSE IF( err_io > 0) THEN
      WRITE(out_unit, *) 'err in READ POTENTIAL_NAME'
      STOP 'check potential name'
    END IF
    
    !--Affectation of the read values to those of the operators_t type object--
    Op%pot_name = trim(pot_name)
    Op%niopot = niopot
    IF (Op%pot_name == 'harmonique') THEN
      CALL P1D_Read_operator_parameter_harmonique(Op, nio)
    ELSE IF (Op%pot_name == 'Morse') THEN
      CALL P1D_Read_operator_parameter_Morse(Op, nio)
    ELSE IF (Op%pot_name == 'Double_puit') THEN
      CALL P1D_Read_operator_parameter_double_puit(Op, nio)
    ELSE IF (Op%pot_name == 'escalier') THEN
      CALL P1D_Read_operator_parameter_escalier(Op, nio)
    ELSE
      STOP 'unknown pot name at assignating in Operateur_t from the namelist'
    END IF 
    
  END SUBROUTINE


  SUBROUTINE P1D_Read_operator_parameter_harmonique(Op, nio)
    USE QDUtil_m
    IMPLICIT NONE
    
    type(P1D_operators_t), intent(inout) :: Op
    integer, intent(in)              :: nio

    real(kind = Rkind)               :: mass, raideur, x0
    integer                          :: err_io
    
    NAMELIST /OPERATEUR_PARAMETER_HARMONIQUE/ mass, raideur, x0                     ! indique la liste des éléments qu'il va trouver dans nml (dans l'odre) et déclare une namelist

    !----------------------Initialization to default value---------------------
    mass = ONE
    raideur = ONE
    x0 = ZERO
    
    !----------------------------------Reading---------------------------------
    READ(nio, nml = OPERATEUR_PARAMETER_HARMONIQUE, iostat = err_io)                ! assigne à la liste annoncée les valeurs lues dans la nml du fichier nio

    WRITE(out_unit, nml = OPERATEUR_PARAMETER_HARMONIQUE)

    !---------------------Control of potential reading errors------------------
    IF(err_io < 0) THEN
      WRITE(out_unit, *) 'err in READ OPERATEUR_PARAMETER_H'
      STOP 'check your operateur parameters (h)'
    ELSE IF( err_io > 0) THEN
      WRITE(out_unit, *) 'err in READ OPERATEUR_PARAMETER_H'
      STOP 'check operateur parameters (h)'
    END IF
    
    !--Affectation of the read values to those of the operators_t type object--
    Op%mass = mass
    Op%raideur = raideur
    Op%x0 = x0
    
  END SUBROUTINE


  SUBROUTINE P1D_Read_operator_parameter_morse(Op, nio)
    USE QDUtil_m
    IMPLICIT NONE
    
    type(P1D_operators_t), intent(inout) :: Op
    integer,               intent(in)    :: nio

    real(kind = Rkind)                   :: mass, De, alpha, Re
    integer                              :: err_io
    
    NAMELIST /OPERATEUR_PARAMETER_MORSE/ mass, De, alpha, Re                        ! indique la liste des éléments qu'il va trouver dans nml (dans l'odre) et déclare une namelist
        
    !----------------------Initialization to default value---------------------
    De = ONE
    alpha = ONE
    Re = ZERO
    
    !----------------------------------Reading---------------------------------
    READ(nio, nml = OPERATEUR_PARAMETER_MORSE, iostat = err_io)  !Assigne à la liste annoncée les valeurs lues dans la nml du fichier nio

    WRITE(out_unit, nml = OPERATEUR_PARAMETER_MORSE)

    !---------------------Control of potential reading errors------------------
    IF(err_io < 0) THEN
      WRITE(out_unit, *) 'err in READ OPERATEUR_PARAMETER_M'
      STOP 'check your operateur parameters (m)'
    ELSE IF( err_io > 0) THEN
      WRITE(out_unit, *) 'err in READ OPERATEUR_PARAMETER_M'
      STOP 'check operateur parameters (m)'
    END IF
    
    !--Affectation of the read values to those of the operators_t type object--
    Op%mass = mass
    Op%De = De
    Op%alpha = alpha
    Op%Re = Re
    
  END SUBROUTINE
 
  
  SUBROUTINE P1D_Read_operator_parameter_double_puit(Op, nio)
    USE QDUtil_m
    IMPLICIT NONE
    
    type(P1D_operators_t), intent(inout) :: Op
    integer,               intent(in)    :: nio

    real(kind = Rkind)                   :: mass
    integer                              :: err_io
    
    NAMELIST /OPERATEUR_PARAMETER_DOUBLE_PUIT/ mass                                 ! indique la liste des éléments qu'il va trouver dans nml (dans l'odre) et déclare une namelist

    !----------------------Initialization to default value---------------------
    mass = ONE
    
    !----------------------------------Reading---------------------------------
    READ(nio, nml = OPERATEUR_PARAMETER_DOUBLE_PUIT, iostat = err_io)               ! assigne à la liste annoncée les valeurs lues dans la nml du fichier nio

    WRITE(out_unit, nml = OPERATEUR_PARAMETER_DOUBLE_PUIT)

    !---------------------Control of potential reading errors------------------
    IF(err_io < 0) THEN
      WRITE(out_unit, *) 'err in READ OPERATEUR_PARAMETER_dp', err_io
      STOP 'check your operateur parameters (dp)'
    ELSE IF( err_io > 0) THEN
      WRITE(out_unit, *) 'err in READ OPERATEUR_PARAMETER_dp', err_io
      STOP 'check operateur parameters (dp)'
    END IF
    
    !--Affectation of the read values to those of the operators_t type object--
    Op%mass = mass

  END SUBROUTINE


  SUBROUTINE P1D_Read_operator_parameter_escalier(Op, nio)
    USE QDUtil_m
    IMPLICIT NONE
    
    type(P1D_operators_t), intent(inout) :: Op
    integer,               intent(in)    :: nio

    real(kind = Rkind)                   :: mass, V_0, V_1, V_2, x1, x2
    integer                              :: err_io
    
    NAMELIST /OPERATEUR_PARAMETER_ESCALIER/ mass, V_0, V_1, V_2, x1, x2             ! indique la liste des éléments qu'il va trouver dans nml (dans l'odre) et déclare une namelist

    !----------------------Initialization to default value---------------------
    mass = ONE
    V_0 = ZERO
    V_1 = ONE
    V_2 = TWO
    x1 = ONE
    x2 = TWO

    !----------------------------------Reading---------------------------------
    READ(nio, nml = OPERATEUR_PARAMETER_ESCALIER, iostat = err_io)                  ! assigne à la liste annoncée les valeurs lues dans la nml du fichier nio

    WRITE(out_unit, nml = OPERATEUR_PARAMETER_ESCALIER)

    !---------------------Control of potential reading errors------------------
    IF(err_io < 0) THEN
      WRITE(out_unit, *) 'err in READ OPERATEUR_PARAMETER_e'
      STOP 'check your operateur parameters (e)'
    ELSE IF( err_io > 0) THEN
      WRITE(out_unit, *) 'err in READ OPERATEUR_PARAMETER_e'
      STOP 'check operateur parameters (e)'
    END IF
    
    !--Affectation of the read values to those of the operators_t type object--
    Op%mass = mass
    Op%V_0 = V_0
    Op%V_1 = V_1
    Op%V_2 = V_2
    Op%x1 = x1
    Op%x2 = x2
    
  END SUBROUTINE


  SUBROUTINE P1D_Construct_operator_grid(Op, Basis, nio, args)
    USE P1D_basis_m
    USE QDUtil_m
    IMPLICIT NONE
  
    TYPE(P1D_operators_t), intent(inout) :: Op
    TYPE(P1D_basis_t),     intent(in)    :: Basis
    character(len=12),     intent(in)    :: args(:)

    integer                              :: iq, nio

    ALLOCATE(Op%V(Basis%nq))

    !-------------------------------Construction-------------------------------
    CALL P1D_Read_operator_parameter(Op, nio, args)

    DO iq = 1, Basis%nq
      Op%V(iq) = P1D_Potentiel(Basis%x(iq), Op)                                     ! /!\ colonne peut-être à transposer pour pdts matriciels
    END DO
    
    !----------------------------------Writing---------------------------------
    WRITE(Op%niopot,*) 'GRILLE - POTENTIEL V(x) :'
    WRITE(Op%niopot,*) '---------iq----------------x(iq)---------------------w(iq)-----------------V(iq)[Ha]&
    &---------------------[cm-1]---------------------[eV]'    
    DO iq = 1, Basis%nq
      WRITE(Op%niopot,*) iq, Basis%x(iq), Basis%w(iq),Op%V(iq), Op%V(iq)*219474.63_Rkind, Op%V(iq)*27.2114_Rkind
    END DO
    
  END SUBROUTINE


  SUBROUTINE P1D_Construct_H(H, Basis, Op)
    USE P1D_basis_m
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), allocatable, intent(inout) :: H(:,:)
    TYPE(P1D_basis_t),             intent(in)    :: Basis
    TYPE(P1D_operators_t),         intent(in)    :: Op

    real(kind=Rkind), allocatable                :: G(:), d0bgw(:,:)                ! G_ib(:)=Ĥ|bi(:)>=Ĥ|bi>(:)
    integer                                      :: ib

    !CALL P1D_Construct_operator_grid(Basis, Op)                                    ! déjà appelé dans programme principal

    ALLOCATE(H(Basis%nb,Basis%nb))
    ALLOCATE(G(Basis%nq))
    ALLOCATE(d0bgw(Basis%nb,Basis%nq))                                              ! inutile : alloue automatiquement quand transpose ou pdt mat 

    !-------------------------------Construction-------------------------------
    d0bgw = transpose(Basis%d0gb)
    DO ib = 1, Basis%nb
      d0bgw(ib,:) = d0bgw(ib,:)*Basis%w(:)                                          ! inutile (?) calc G(:)*w(:) après (?) 
    END DO                                                                          ! dim(dbgw) = nb*nq
    
    DO ib = 1, Basis%nb
      G(:) = -(HALF/Op%mass)*Basis%d2gb(:,ib) + Op%V(:)*Basis%d0gb(:,ib)            ! dim(G_ib) = nq*1
      G(:) =  G(:)*Basis%w(:)
      H(:,ib) = matmul(transpose(Basis%d0gb), G)
    END DO

    DEALLOCATE(G, d0bgw)

  END SUBROUTINE


  SUBROUTINE P1D_Calc_energie(E, H, B)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: E
    real(kind=Rkind),    intent(in)    :: H(:,:)
    complex(kind=Rkind), intent(in)    :: B(:)                                      ! sur la base des sinus

    E = dot_product(B, matmul(H, B)) 
    !WRITE(out_unit,*) 'ENERGIE MOYENNE : <Psi|H|Psi> = ', E

  END SUBROUTINE


  SUBROUTINE P1D_Construct_matX(mat_X, Basis)                                       ! mat de l'opérateur position
    USE P1D_basis_m
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), allocatable, intent(inout) :: mat_X(:,:)
    TYPE(P1D_basis_t),                intent(in)    :: Basis
    real(kind=Rkind),    allocatable                :: B(:), d0bgw(:,:)             ! B_ib(:)=X|bi(:)>=X|bi>(:)

    integer                                         :: ib

    ALLOCATE(mat_X(Basis%nb,Basis%nb))
    ALLOCATE(B(Basis%nq))

    !-------------------------------Construction-------------------------------
    d0bgw = transpose(Basis%d0gb)
    DO ib = 1, Basis%nb
      d0bgw(ib,:) = d0bgw(ib,:)*Basis%w(:)
    END DO                                                                          ! dim(d0bgw) = nb*nq

    DO ib = 1, Basis%nb
      B(:) = Basis%x(:)*Basis%d0gb(:,ib)
      mat_X(:,ib) = matmul(d0bgw, B)
    END DO

    DEALLOCATE(B, d0bgw)

  END SUBROUTINE


  SUBROUTINE P1D_Position_moyenne(pos, mat_X, B)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: pos
    complex(kind=Rkind), intent(in)    :: mat_X(:,:)
    complex(kind=Rkind), intent(in)    :: B(:)

    pos = dot_product(B, matmul(mat_X, B))
    !WRITE(out_unit,*) 'POSITION MOYENNE : <Psi|X|Psi> = ', pos

  END SUBROUTINE


  SUBROUTINE P1D_Construct_Px(mat_Px, Basis)                                        ! mat de l'opérateur impulsion
    USE P1D_basis_m
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), allocatable, intent(inout) :: mat_Px(:,:)
    TYPE(P1D_basis_t),                intent(in)    :: Basis

    complex(kind=Rkind), allocatable                :: G(:)                         ! G_ib(:)=Px|bi(:)>=Px|bi>(:)
    integer                                         :: ib, i

    ALLOCATE(mat_Px(Basis%nb,Basis%nb))
    ALLOCATE(G(Basis%nq))

    !-------------------------------Construction-------------------------------
    DO ib = 1, Basis%nb              
      G(:) = -EYE * Basis%d1gb(:,ib)               
      G(:) = G(:) * Basis%w(:)                              
      mat_Px(:,ib) = matmul(transpose(Basis%d0gb), G)              
    END DO           
      
    DEALLOCATE(G)

  END SUBROUTINE


  SUBROUTINE P1D_Impulsion_moyenne(px, mat_Px, B)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind),    intent(inout) :: px
    complex(kind=Rkind), intent(in)    :: mat_Px(:,:)
    complex(kind=Rkind), intent(in)    :: B(:)

    px = dot_product(B, matmul(mat_Px, B))
    !WRITE(out_unit,*) 'IMPULSION MOYENNE : <Psi|Px|Psi> = ', px

  END SUBROUTINE


END MODULE