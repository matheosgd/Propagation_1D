MODULE P1D_basis_m
  USE QDUtil_m
  IMPLICIT NONE


  PUBLIC
  PRIVATE :: P1D_Read_basis


  TYPE :: P1D_basis_t
    character(len=:), allocatable :: basis_name   
    integer                       :: nb = 3                                         ! nombre de fonctions de base (initialise à 3)
    integer                       :: nq = 3                                         ! nombre de points de grille (initialise à 3)
    real(kind=Rkind)              :: a = ZERO                                       ! borne inferieure (initialise à 0)
    real(kind=Rkind)              :: b = ZERO                                       ! borne supérieure (initialise à 0)
    real(kind=Rkind), allocatable :: w(:)                                           ! poids  
    real(kind=Rkind), allocatable :: x(:)                                           ! liste des différents points de la grille
    real(kind=Rkind), allocatable :: d0gb(:, :)                                     ! fonction  initiale de la base = matrice def par ||b1(x1)  b2(x1)  b3(x1)| ->(fns de base en col))
                                                                                    ! les fonctions de base appliquées en les pts    v|b1(x2)  b2(x2)  b3(x2)|
                                                                                    ! de la grille                (points grille en l)|b1(x3)  b2(x3)  b3(x3)|
    real(kind=Rkind), allocatable :: d1gb(:, :)                                     ! derivee premiere   
    real(kind=Rkind), allocatable :: d2gb(:, :)                                     ! derivee seconde de la base 
    real(kind=Rkind)              :: scaleQ                                         ! scaleQ facteur d'echelle
    real(kind=Rkind)              :: x0 = ZERO                                      ! centre de la base (si on a pas une base de sinus => xo!= a)

  END TYPE


  CONTAINS


  SUBROUTINE P1D_Read_basis(Basis, nio)                                             ! nio le fichier dans lequel prendre les valeurs
    USE QDUtil_m
    IMPLICIT NONE
    
    TYPE(P1D_basis_t), intent(inout) :: Basis   
    integer,           intent(in)    :: nio

    character(len=10)                :: basis_name   
    real(kind=Rkind)                 :: a, b, x0, scaleQ                            ! déclare deux réels a,b. Pas a priori = a,b du TYPE dérivé
    integer                          :: nb, nq, err_io                              ! err_io cf contrôle erreures

    NAMELIST /BASIS_INFO/ basis_name, nb, nq, a, b, x0, scaleQ                      ! indique la liste des éléments qu'il va trouver dans nml (dans l'odre) et déclare une namelist

    !---------------------Initialization to default values---------------------
    Basis_name = '0'
    nb = 3
    nq = 3
    a = ZERO
    b = TEN
    x0 = ZERO
    scaleQ = ONE
 
    !----------------------------Reading of the nml----------------------------
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) 'CONSTRUCT BASIS :'
    READ(nio, nml = BASIS_INFO, iostat = err_io)                                    ! assigne à la liste annoncée les valeurs lues dans la nml du fichier nio

    WRITE(out_unit, nml = BASIS_INFO)
    
    !---------------------Control of potential reading errors------------------
    IF(err_io < 0) THEN
      WRITE(out_unit, *) 'err in READ BASIS_INFO'
      STOP 'check basis data'
    ELSE IF( err_io > 0) THEN
      WRITE(out_unit, *) 'err in READ BASIS_INFO'
      STOP 'check basis data'
    END IF

    !--Affectation of the read values to those of the P1D_basis_t type object-
    Basis%basis_name = trim(basis_name)
    Basis%nb = nb
    Basis%nq = nq
    Basis%a = a
    Basis%b = b
    Basis%scaleQ = scaleQ                                                           ! utilise valeurs par défaut (cf la nml) ...
    Basis%x0 = x0

    IF(Basis%basis_name == 'boxAB') then                                            ! ... sauf si la base est de TYPE boxAB auquel cas on calcule scaleQ et x0 EX a et b.
      Basis%scaleQ = PI / (Basis%b-Basis%a)
      Basis%x0 = Basis%a
    END IF

  END SUBROUTINE


  SUBROUTINE P1D_Construct_basis(Basis, nio)                                        ! il reste à construire w(:), x(:), d0,1,2gb
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(P1D_basis_t), intent(inout) :: Basis
    integer,           intent(in)    :: nio

    integer                          :: ib, iq
    real(kind=Rkind)                 :: dx

    CALL P1D_Read_basis(Basis, nio)                                                 ! pour avoir a,b,nb,nq assignés à Basis%a,b,nb,nq 

    !----------------------------Grid initialization---------------------------
    ALLOCATE(Basis%w(Basis%nq))
    ALLOCATE(Basis%x(Basis%nq))
    ALLOCATE(Basis%d0gb(Basis%nq,Basis%nb))
    ALLOCATE(Basis%d1gb(Basis%nq,Basis%nb))
    ALLOCATE(Basis%d2gb(Basis%nq,Basis%nb))

    dx = PI/Basis%nq
    Basis%w(:) = dx                                                                 ! here we decide constant step size
    
    Basis%x = [(dx*(iq- HALF), iq = 1, Basis%nq)]                                   ! generation of the grid /!\ ligne peut-être à transposer pour pdts matriciels

    DO ib = 1, Basis%nb, 1
      Basis%d0gb(:,ib) = sin(ib*Basis%x(:))*sqrt(TWO*PI**(-1))
      Basis%d1gb(:,ib) = ib*cos(ib*Basis%x(:))*sqrt(TWO*PI**(-1))
      Basis%d2gb(:,ib) = -ib**2*sin(ib*Basis%x(:))*sqrt(TWO*PI**(-1))
    END DO

    IF (Basis%nb == Basis%nq) then
      Basis%d0gb(:,Basis%nb) = Basis%d0gb(:,Basis%nb)/sqrt(TWO)
      Basis%d1gb(:,Basis%nb) = Basis%d1gb(:,Basis%nb)/sqrt(TWO)
      Basis%d2gb(:,Basis%nb) = Basis%d2gb(:,Basis%nb)/sqrt(TWO)
    END IF

    !----------------------------Change in variables---------------------------
    Basis%x(:) = Basis%x(:)/Basis%scaleQ + Basis%x0
    Basis%w(:) = Basis%w(:) / Basis%scaleQ
    Basis%d0gb(:,:) = Basis%d0gb(:,:)* sqrt(Basis%scaleQ)
    Basis%d1gb(:,:) = Basis%d1gb(:,:)* sqrt(Basis%scaleQ)*Basis%scaleQ 
    Basis%d2gb(:,:) = Basis%d2gb(:,:)* sqrt(Basis%scaleQ)*Basis%scaleQ*Basis%scaleQ

  END SUBROUTINE


  SUBROUTINE P1D_Check_ortho(Basis, niomat)
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(P1D_basis_t), intent(inout) :: Basis
    integer,           intent(in)    :: niomat

    real(kind=Rkind), allocatable    :: S0(:,:), S1(:,:), S2(:,:), d0bgw(:,:)
    integer                          :: ib

    ALLOCATE(S0(Basis%nb,Basis%nb))
    ALLOCATE(S1(Basis%nb,Basis%nb))
    ALLOCATE(S2(Basis%nb,Basis%nb))
    ALLOCATE(d0bgw(Basis%nb,Basis%nq))

    d0bgw = transpose(Basis%d0gb)

    DO ib = 1, Basis%nb
      d0bgw(ib,:) = d0bgw(ib,:)*Basis%w(:)
    END DO

    S0 = matmul(d0bgw, Basis%d0gb)
    S1 = matmul(d0bgw, Basis%d1gb)
    S2 = matmul(d0bgw, Basis%d2gb)
 
    WRITE(niomat,*) '***************************************************'

    WRITE(niomat,*) 'S0 = '
    CALL WRITE_Mat(S0,niomat,5,info='<d0gb|d0gb>') 
    WRITE(niomat,*)

    WRITE(niomat,*) '***************************************************'

    WRITE(niomat,*) 'S1 = '
    CALL WRITE_Mat(S1,niomat,5,info='<d0gb|d1gb>')
    WRITE(niomat,*)
    
    WRITE(niomat,*) '***************************************************'

    WRITE(niomat,*) 'S2 = '
    CALL WRITE_Mat(S2,niomat,5,info='<d0gb|d2gb>')
    WRITE(niomat,*)

    DEALLOCATE(S0, S1, S2, d0bgw)

  END SUBROUTINE


  SUBROUTINE P1D_GridTOBasis(B, G, Basis)
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: B(:)
    TYPE(P1D_basis_t),   intent(in)    :: Basis
    complex(kind=Rkind), intent(in)    :: G(:)

    real(kind=Rkind), allocatable      :: d0bgw(:,:)
    integer                            :: ib

    ALLOCATE(d0bgw(Basis%nb, Basis%nq))

    d0bgw = transpose(Basis%d0gb)
    DO ib = 1, Basis%nb
      d0bgw(ib,:) = d0bgw(ib,:)*Basis%w(:)
    END DO  

    B(:) = matmul(d0bgw, G)

    DEALLOCATE(d0bgw)

  END SUBROUTINE


  SUBROUTINE P1D_BasisTOGrid_cplx(G, B, Basis)
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: G(:)                                      ! la projection de B(:) sur la grille 
    TYPE(P1D_basis_t),   intent(in)    :: Basis                                     ! la base des sinus
    complex(kind=Rkind), intent(in)    :: B(:)                                      ! vecteur def sur la base des sinus

    G(:) = matmul(Basis%d0gb, B)
  
  END SUBROUTINE


  SUBROUTINE P1D_BasisTOGrid_real(G, B, Basis)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout) :: G(:)                                         ! la projection de B(:) sur la grille 
    TYPE(P1D_basis_t), intent(in)   :: Basis                                        ! la base des sinus
    real(kind=Rkind), intent(in)    :: B(:)                                         ! vecteur def sur la base des sinus

    G(:) = matmul(Basis%d0gb, B)
  
  END SUBROUTINE


END MODULE
