MODULE Basis_m
  USE QDUtil_m
  IMPLICIT NONE

  PUBLIC
  PRIVATE :: read_Basis

  TYPE :: Basis_t
    character(len=:), allocatable     :: basis_name   
    integer                           :: nb = 3     !nombre de fonctions de base (initialise à 3)
    integer                           :: nq = 3     !nombre de points de grille (initialise à 3)
    real(kind=Rkind)                  :: a = ZERO   !borne inferieure grille (initialise à 0)
    real(kind=Rkind)                  :: b = ZERO   !borne supérieure grille (initialise à 0)
    real(kind=Rkind), allocatable     :: w(:)       !poids  
    real(kind=Rkind), allocatable     :: x(:)       !liste des différents points de la grille
    real(kind=Rkind), allocatable     :: d0gb(:, :) !fonction  initiale de la base = matrice def par ||b1(x1)  b2(x1)  b3(x1)| ->(fonctions de base en col))
                                                    !                                                v|b1(x2)  b2(x2)  b3(x2)|
                                                    !                             (points grille en l)|b1(x3)  b2(x3)  b3(x3)|
    real(kind=Rkind), allocatable     :: d1gb(:, :) !derivee premiere   
    real(kind=Rkind), allocatable     :: d2gb(:, :) !derivee seconde des fonctions de base 
    real(kind=Rkind)                  :: scaleQ     !facteur d'echelle du changement de variables
    real(kind=Rkind)                  :: x0 = ZERO  !centre de la base
  END TYPE



CONTAINS



  SUBROUTINE read_Basis(Basis, nio)                           !nio le fichier de données
    USE QDUtil_m
    IMPLICIT NONE
    
    TYPE(Basis_t)  , intent(inout) :: Basis   
    integer        , intent(in)    :: nio
    character(len=10)              :: basis_name   
    real(kind=Rkind)               :: a, b, x0, scaleQ
    integer                        :: nb, nq, err_io           !err_io cf contrôle erreures

    NAMELIST /BASIS_INFO/ basis_name, nb, nq, a, b, x0, scaleQ

!Initialise des valeurs par défaut
    Basis_name = '0'
    nb = 3
    nq = 3
    a = ZERO
    b = TEN
    x0 = ZERO
    scaleQ = ONE
 
!Lecture
    WRITE(out_unit,*) '*******************************************************'
    WRITE(out_unit,*) 'CONSTRUCT BASIS :'
    READ(nio, nml = BASIS_INFO, iostat = err_io)

    WRITE(out_unit, nml = BASIS_INFO)
    
!Contrôle des erreures de lectures
    IF(err_io < 0) then
      WRITE(out_unit, *) 'err in READ BASIS_INFO'
      STOP 'check basis data'
    ELSE IF( err_io > 0) then
      WRITE(out_unit, *) 'err in READ BASIS_INFO'
      STOP 'check basis data'
    END IF

!Assigne les variables de la liste à celles de l'objet Basis du TYPE Basis_t
    Basis%basis_name = trim(basis_name)
    Basis%nb = nb
    Basis%nq = nq
    Basis%a = a
    Basis%b = b
    Basis%scaleQ = scaleQ                 !utilise valeurs par défaut (cf la nml) ...
    Basis%x0 = x0

    IF(Basis%basis_name == 'boxAB') then  !... sauf si la base est de TYPE boxAB auquel cas on calcule scaleQ et x0 EX a et b.
      Basis%scaleQ = PI / (Basis%b-Basis%a)
      Basis%x0 = Basis%a
    END IF

  END SUBROUTINE



  SUBROUTINE construct_basis(Basis, nio)
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(Basis_t)  , intent(inout) :: Basis
    integer, intent(in)            :: nio
    integer                        :: ib, iq
    real(kind=Rkind)               :: dx

    call read_Basis(Basis, nio)

! Allocation-----------------------------------------------------------
    ALLOCATE(Basis%w(Basis%nq))
    ALLOCATE(Basis%x(Basis%nq))
    ALLOCATE(Basis%d0gb(Basis%nq,Basis%nb))
    ALLOCATE(Basis%d1gb(Basis%nq,Basis%nb))
    ALLOCATE(Basis%d2gb(Basis%nq,Basis%nb))

! Construction---------------------------------------------------------
    dx = PI/Basis%nq
    Basis%w(:) = dx !ici les poids sont constants
    
    Basis%x = [(dx*(iq- HALF), iq = 1, Basis%nq)] !/!\ ligne peut-être à transposer pour pdts matriciels 
                                                

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
!Changement de variable---------------------------------------------

    Basis%x(:) = Basis%x(:)/Basis%scaleQ + Basis%x0
    Basis%w(:) = Basis%w(:) / Basis%scaleQ
    Basis%d0gb(:,:) = Basis%d0gb(:,:)* sqrt(Basis%scaleQ)
    Basis%d1gb(:,:) = Basis%d1gb(:,:)* sqrt(Basis%scaleQ)*Basis%scaleQ 
    Basis%d2gb(:,:) = Basis%d2gb(:,:)* sqrt(Basis%scaleQ)*Basis%scaleQ*Basis%scaleQ

  END SUBROUTINE   



  SUBROUTINE check_ortho(Basis, niomat)
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(Basis_t)  , intent(inout) :: Basis
    integer, intent(in)            :: niomat           
    real(kind=Rkind), allocatable  :: S0(:,:), S1(:,:), S2(:,:), d0bgw(:,:)
    integer                        :: ib

! Allocation-----------------------------------------------------------
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



  SUBROUTINE GridTOBasis(B, G, Basis)
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(Basis_t), intent(in)                :: Basis
    complex(kind=Rkind), intent(in)          :: G(:)
    complex(kind=Rkind), intent(inout)       :: B(:)
    real(kind=Rkind), allocatable            :: d0bgw(:,:)
    integer                                  :: ib

    ALLOCATE(d0bgw(Basis%nb, Basis%nq))

    d0bgw = transpose(Basis%d0gb)
    DO ib = 1, Basis%nb
      d0bgw(ib,:) = d0bgw(ib,:)*Basis%w(:)
    END DO  

    B(:) = matmul(d0bgw, G)

    DEALLOCATE(d0bgw)
  END SUBROUTINE



  SUBROUTINE BasisTOGrid_cplx(G, B, Basis)
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(Basis_t), intent(in)                :: Basis  !la base des sinus
    complex(kind=Rkind), intent(inout)       :: G(:)   !la projection de B(:) sur la grille 
    complex(kind=Rkind), intent(in)          :: B(:)   !vecteur complexe def sur la base des sinus

    G(:) = matmul(Basis%d0gb, B)
  
  END SUBROUTINE



  SUBROUTINE BasisTOGrid_real(G, B, Basis)
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(Basis_t), intent(in)             :: Basis  !la base des sinus
    real(kind=Rkind), intent(inout)       :: G(:)   !la projection de B(:) sur la grille 
    real(kind=Rkind), intent(in)          :: B(:)   !vecteur reel def sur la base des sinus

    G(:) = matmul(Basis%d0gb, B)
  
  END SUBROUTINE



END MODULE
