PROGRAM TD_Matheo
  USE QDUtil_m
  USE Basis_m
  USE Operateurs_m
  USE Psi_m
  USE Lanczos_m
  USE Propagation_m
  IMPLICIT NONE
  
  real(kind=Rkind)                 :: E = ZERO
  TYPE(Basis_t)                    :: Basis
  TYPE(Operateur_t)                :: Op
  real(kind=Rkind), allocatable    :: H(:,:)
  complex(kind=Rkind), allocatable :: mat_X(:,:), mat_Px(:,:)
  complex(kind=Rkind), allocatable :: psi0(:), psif(:)
  integer                          :: niomat
!TRACES PSIF
  integer                          :: niovec, i
  complex(kind=Rkind), allocatable :: psifgrid(:)
!DEUXIEME, TROISIEME ET QUATRIEME PROPAGATIONS
  complex(kind=Rkind), allocatable :: psif2(:), psif3(:), psif4(:)
  integer                          :: niovec

  OPEN(newunit = niomat, file = 'mat')

!CONSTRUCTION BASE--------------------------------------------  
  CALL construct_basis(Basis=Basis, nio=in_unit )

!CONSTRUCTION POTENTIEL---------------------------------------  
  CALL construct_operator_grid(Op, Basis, nio=in_unit)
 
!CONSTRUCTION HAMILTONIEN-------------------------------------    
  CALL Construct_H(H, Basis, Op)

!CONSTRUCTION MAT OP POSITION---------------------------------
  CALL construct_matX(mat_X, Basis)

!CONSTRUCTION MAT OP IMPULSION--------------------------------
  CALL Construct_Px(mat_Px, Basis) 

!CONSTRUCTION PAQUET D'ONDE-----------------------------------
  ALLOCATE(psi0(Basis%nb)) 
  CALL init_Psi(psi0, Basis, in_unit)

!PREMIERE PROPAGATION DU PAQUET D ONDE------------------------
  ALLOCATE(psif(Basis%nb))
  CALL propagation(psif, psi0, H, Basis, Op, in_unit) 
  WRITE(out_unit,*) '*******************************************************'

!DEUXIEME PROPAGATION DU PAQUET D ONDE------------------------
  ALLOCATE(psif2(Basis%nb))
  CALL propagation(psif2, psi0, H, Basis, Op, in_unit) 
  WRITE(out_unit,*) '*******************************************************'

!TROISIEME PROPAGATION DU PAQUET D ONDE-----------------------
  ALLOCATE(psif3(Basis%nb))
  CALL propagation(psif3, psi0, H, Basis, Op, in_unit) 
  WRITE(out_unit,*) '*******************************************************'

!QUATRIEME PROPAGATION DU PAQUET D ONDE-----------------------
  ALLOCATE(psif4(Basis%nb))
  CALL propagation(psif4, psi0, H, Basis, Op, in_unit) 
  WRITE(out_unit,*) '*******************************************************'

!TRACES PSIF--------------------------------------------------
  OPEN(newunit = niovec, file = 'psif')
  ALLOCATE(psifgrid(Basis%nq)) 
  CALL BasisTOGrid_cplx(psifgrid, psif, Basis)
  DO i = 1, Basis%nq
    WRITE(niovec, *) Basis%x(i), (ABS(psifgrid(i))**2)*219474.631443_Rkind
  END DO

!AFFICHAGE EIGVAL ET TROIS PREMIERS EIGVEC DANS RESULTS----------------------
  ALLOCATE(Op%REigval(Basis%nb))
  ALLOCATE(Op%REigvec(Basis%nb,Basis%nb))
  CALL diagonalization(H,Op%REigval,Op%REigvec)
  !EIGVAL-----------------------------------------------------
  WRITE(out_unit,*) 'VALEURS PROPRES'
  CALL WRITE_Vec(Op%Reigval, out_unit, 5, info = 'VP[Ha]')
  WRITE(out_unit,*)
  CALL WRITE_Vec(Op%Reigval*219474.631443_Rkind, out_unit, 5, info = 'VP[cm-1]')
  !EIGVEC-----------------------------------------------------
  WRITE(out_unit,*) '*******************************************************'
  WRITE(out_unit,*) 'TROIS PREMIERS VECTEURS PROPRES SUR LA BASE DES SINUS'
  DO i = 1, 3
    WRITE(out_unit, *) Op%Reigvec(:,i)
  END DO
  
END PROGRAM 