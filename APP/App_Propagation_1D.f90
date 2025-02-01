PROGRAM TD_Matheo
  USE QDUtil_m
  USE P1D_basis_m
  USE P1D_operators_m
  USE P1D_wavepacket_m
  USE P1D_lanczos_m
  USE P1D_propagation_m
  !USE P1D_perturbations_m
  IMPLICIT NONE


  !-----------------Variables for the recovery of the arguments----------------
  integer                          :: num_args, ix                                  ! the number of arguments passed and a loop increment dedicated to this process
  character(len=12), allocatable   :: args(:)                                       ! the list that will hold the aforementioned arguments

  !--------------------Variables for the system description--------------------
  TYPE(P1D_basis_t)                :: Basis
  TYPE(P1D_operators_t)            :: Op
  real(kind=Rkind), allocatable    :: H(:,:)
  complex(kind=Rkind), allocatable :: mat_X(:,:), mat_Px(:,:)
  complex(kind=Rkind), allocatable :: psi0(:)                                       ! the initial wavepacket

  !-------------Variables for the animated plot of the propagation-------------
  complex(kind=Rkind), allocatable :: psif1(:),  psif2(:),  psif3(:),  psif4(:),  psif5(:),  psif6(:),  &
                                    & psif7(:),  psif8(:),  psif9(:),  psif10(:), psif11(:), psif12(:), &
                                    & psif13(:), psif14(:), psif15(:), psifgrid(:)  ! the final wavepacket after successive propagation times
  complex(kind=Rkind), allocatable :: psi0grid(:),   psif1grid(:),  psif2grid(:),  &
                                    & psif3grid(:),  psif4grid(:),  psif5grid(:),  &
                                    & psif6grid(:),  psif7grid(:),  psif8grid(:),  &
                                    & psif9grid(:),  psif10grid(:), psif11grid(:), &
                                    & psif12grid(:), psif13grid(:), psif14grid(:), &
                                    & psif15grid(:)                                 ! the projection of the latter on the grid (to be plottable)

  !-------------Variables for the plot chosen potential eigenstates------------ 
  real(kind=Rkind), allocatable    :: psigrid0(:), psigrid1(:), psigrid2(:)         ! their grid projection

  !------------------------------Utility variables-----------------------------
  integer                          :: niomat, niovec, i                             ! two file Fortran units, and a loop increment


!-----------------------------PRELIMINARY PROCESSES-----------------------------
  !--------------------Recovery of the command-line arguments-------------------
  num_args = command_argument_count()
  ALLOCATE(args(num_args))

  IF (num_args < 1) THEN
    STOP "No argument has been provided, the writing of the result files will not be realised. &
        & cf. the script run.sh or the manual when it will be available."
  ELSE
    DO ix = 1, num_args
      CALL get_command_argument(ix,args(ix))
      WRITE(out_unit,*) "The", ix, "^th argument has been read as : ", args(ix)
    END DO
  END IF

  !--------------opening (and creation if needed) of a result file--------------
  OPEN(NEWUNIT = niomat, FILE = 'OUT/P1D_App_Propagation_1D_'//trim(args(1))//'_matrix.log', & 
     & FORM = 'formatted', ACTION = 'write', POSITION = 'rewind')


!-----------------------------SYSTEM INITIALIZATION----------------------------
  !-----------------------------Basis construction-----------------------------
  CALL P1D_Construct_basis(Basis=Basis, nio=in_unit )                               ! cf. run.sh : in_unit = data
  
  !---------------------------Potential construction---------------------------
  CALL P1D_Construct_operator_grid(Op, Basis, nio=in_unit, args=args)
 
  !--------------------------Hamiltonian construction--------------------------
  CALL P1D_Construct_H(H, Basis, Op)

  !-------------------Position operator's matrix construction------------------
  CALL P1D_Construct_matX(mat_X, Basis)

  !-------------------Momentum operator's matrix construction------------------
  CALL P1D_Construct_Px(mat_Px, Basis) 

  !---------------------------Wavepacket construction--------------------------
  ALLOCATE(psi0(Basis%nb)) 
  CALL P1D_Init_Psi(psi0, Basis, in_unit)


!------------------------------SYSTEM PROPAGATION------------------------------
  !-----------------------Wavepacket's first propagation-----------------------
  ALLOCATE(psif1(Basis%nb))
  CALL P1D_Propagation(psif1, psi0, H, Basis, Op, in_unit, args) 
  WRITE(out_unit,*) '*******************************************************'

  !-----------------------Wavepacket's second propagation----------------------
  ALLOCATE(psif2(Basis%nb))
  CALL P1D_Propagation(psif2, psi0, H, Basis, Op, in_unit, args) 
  WRITE(out_unit,*) '*******************************************************'

  !-----------------------Wavepacket's third propagation-----------------------
  ALLOCATE(psif3(Basis%nb))
  CALL P1D_Propagation(psif3, psi0, H, Basis, Op, in_unit, args) 
  WRITE(out_unit,*) '*******************************************************'

  !-----------------------Wavepacket's fourth propagation----------------------
  ALLOCATE(psif4(Basis%nb))
  CALL P1D_Propagation(psif4, psi0, H, Basis, Op, in_unit, args) 
  WRITE(out_unit,*) '*******************************************************'

  !-----------------------Wavepacket's fifth propagation-----------------------
  ALLOCATE(psif5(Basis%nb))
  CALL P1D_Propagation(psif5, psi0, H, Basis, Op, in_unit, args) 
  WRITE(out_unit,*) '*******************************************************'

  !-----------------------Wavepacket's sixth propagation-----------------------
  ALLOCATE(psif6(Basis%nb))
  CALL P1D_Propagation(psif6, psi0, H, Basis, Op, in_unit, args) 
  WRITE(out_unit,*) '****************************************************'

  !----------------------Wavepacket's seventh propagation----------------------
  ALLOCATE(psif7(Basis%nb))
  CALL P1D_Propagation(psif7, psi0, H, Basis, Op, in_unit, args) 
  WRITE(out_unit,*) '****************************************************'

  !-----------------------Wavepacket's eigth propagation-----------------------
  ALLOCATE(psif8(Basis%nb))
  CALL P1D_Propagation(psif8, psi0, H, Basis, Op, in_unit, args) 
  WRITE(out_unit,*) '****************************************************'

  !-----------------------Wavepacket's ninth propagation-----------------------
  ALLOCATE(psif9(Basis%nb))
  CALL P1D_Propagation(psif9, psi0, H, Basis, Op, in_unit, args) 
  WRITE(out_unit,*) '****************************************************'

  !-----------------------Wavepacket's tenth propagation-----------------------
  ALLOCATE(psif10(Basis%nb))
  CALL P1D_Propagation(psif10, psi0, H, Basis, Op, in_unit, args) 
  WRITE(out_unit,*) '****************************************************'

  !----------------------Wavepacket's eleventh propagation---------------------
  ALLOCATE(psif11(Basis%nb))
  CALL P1D_Propagation(psif11, psi0, H, Basis, Op, in_unit, args) 
  WRITE(out_unit,*) '****************************************************'

  !----------------------Wavepacket's twelvth propagation----------------------
  ALLOCATE(psif12(Basis%nb))
  CALL P1D_Propagation(psif12, psi0, H, Basis, Op, in_unit, args) 
  WRITE(out_unit,*) '****************************************************'

  !---------------------Wavepacket's thirteenth propagation--------------------
  ALLOCATE(psif13(Basis%nb))
  CALL P1D_Propagation(psif13, psi0, H, Basis, Op, in_unit, args) 
  WRITE(out_unit,*) '****************************************************'

  !---------------------Wavepacket's fourteenth propagation--------------------
  ALLOCATE(psif14(Basis%nb))
  CALL P1D_Propagation(psif14, psi0, H, Basis, Op, in_unit, args) 
  WRITE(out_unit,*) '****************************************************'

  !---------------------Wavepacket's fifteenth propagation---------------------
  ALLOCATE(psif15(Basis%nb))
  CALL P1D_Propagation(psif15, psi0, H, Basis, Op, in_unit, args) 
  WRITE(out_unit,*) '****************************************************'

  !-------------------------------Animation plot-------------------------------
  OPEN(NEWUNIT = niovec, FILE = 'OUT/P1D_App_Propagation_1D_'//trim(args(1))//'_psif.log', &
     & FORM = 'formatted', ACTION = 'write', POSITION = 'rewind')

  ALLOCATE(psi0grid(Basis%nq)) 
  ALLOCATE(psif1grid(Basis%nq)) 
  ALLOCATE(psif2grid(Basis%nq)) 
  ALLOCATE(psif3grid(Basis%nq)) 
  ALLOCATE(psif4grid(Basis%nq)) 
  ALLOCATE(psif5grid(Basis%nq)) 
  ALLOCATE(psif6grid(Basis%nq)) 
  ALLOCATE(psif7grid(Basis%nq)) 
  ALLOCATE(psif8grid(Basis%nq)) 
  ALLOCATE(psif9grid(Basis%nq)) 
  ALLOCATE(psif10grid(Basis%nq))
  ALLOCATE(psif11grid(Basis%nq)) 
  ALLOCATE(psif12grid(Basis%nq)) 
  ALLOCATE(psif13grid(Basis%nq)) 
  ALLOCATE(psif14grid(Basis%nq)) 
  ALLOCATE(psif15grid(Basis%nq)) 

  CALL P1D_BasisTOGrid_cplx(psi0grid,   psi0,   Basis)
  CALL P1D_BasisTOGrid_cplx(psif1grid,  psif1,  Basis)
  CALL P1D_BasisTOGrid_cplx(psif2grid,  psif2,  Basis)
  CALL P1D_BasisTOGrid_cplx(psif3grid,  psif3,  Basis)
  CALL P1D_BasisTOGrid_cplx(psif4grid,  psif4,  Basis)
  CALL P1D_BasisTOGrid_cplx(psif5grid,  psif5,  Basis)
  CALL P1D_BasisTOGrid_cplx(psif6grid,  psif6,  Basis)
  CALL P1D_BasisTOGrid_cplx(psif7grid,  psif7,  Basis)
  CALL P1D_BasisTOGrid_cplx(psif8grid,  psif8,  Basis)
  CALL P1D_BasisTOGrid_cplx(psif9grid,  psif9,  Basis)
  CALL P1D_BasisTOGrid_cplx(psif10grid, psif10, Basis)
  CALL P1D_BasisTOGrid_cplx(psif11grid, psif11, Basis)
  CALL P1D_BasisTOGrid_cplx(psif12grid, psif12, Basis)
  CALL P1D_BasisTOGrid_cplx(psif13grid, psif13, Basis)
  CALL P1D_BasisTOGrid_cplx(psif14grid, psif14, Basis)
  CALL P1D_BasisTOGrid_cplx(psif15grid, psif15, Basis)

  DO i = 1, Basis%nq
    WRITE(niovec, *) Basis%x(i), (ABS(psi0grid(i))**2)   * 219474.631443_Rkind, &
                               & (ABS(psif1grid(i))**2)  * 219474.631443_Rkind, &
                               & (ABS(psif2grid(i))**2)  * 219474.631443_Rkind, &
                               & (ABS(psif3grid(i))**2)  * 219474.631443_Rkind, &
                               & (ABS(psif4grid(i))**2)  * 219474.631443_Rkind, &
                               & (ABS(psif5grid(i))**2)  * 219474.631443_Rkind, &
                               & (ABS(psif6grid(i))**2)  * 219474.631443_Rkind, &
                               & (ABS(psif7grid(i))**2)  * 219474.631443_Rkind, &
                               & (ABS(psif8grid(i))**2)  * 219474.631443_Rkind, &
                               & (ABS(psif9grid(i))**2)  * 219474.631443_Rkind, &
                               & (ABS(psif10grid(i))**2) * 219474.631443_Rkind, &
                               & (ABS(psif11grid(i))**2) * 219474.631443_Rkind, &
                               & (ABS(psif12grid(i))**2) * 219474.631443_Rkind, &
                               & (ABS(psif13grid(i))**2) * 219474.631443_Rkind, &
                               & (ABS(psif14grid(i))**2) * 219474.631443_Rkind, &
                               & (ABS(psif15grid(i))**2) * 219474.631443_Rkind
  END DO


!----------------------PLOT OF THE POTENTIAL'S EIGENSTATES---------------------
  ALLOCATE(Op%REigval(Basis%nb))
  ALLOCATE(Op%REigvec(Basis%nb,Basis%nb))
  CALL diagonalization(H,Op%REigval,Op%REigvec)

  !---------------------------Writing the eigenvalues--------------------------
  WRITE(out_unit,*) 'EIGENVALUES'
  CALL WRITE_Vec(Op%Reigval, out_unit, 5, info = 'VP[Ha]')
  WRITE(out_unit,*)
  CALL WRITE_Vec(Op%Reigval*219474.631443_Rkind, out_unit, 5, info = 'VP[cm-1]')
  WRITE(out_unit,*) '*******************************************************'

  !--------------------------Writing the eigenvectors--------------------------
  WRITE(out_unit,*) 'EIGENVECTORS'
  CALL WRITE_Mat(Op%Reigvec, out_unit, 5, info = 'eigenvectors')
 
  !--------------------------Plot of the three firsts--------------------------
  ALLOCATE(psigrid0(Basis%nq)) 
  ALLOCATE(psigrid1(Basis%nq)) 
  ALLOCATE(psigrid2(Basis%nq))

  CALL P1D_BasisTOGrid_real(psigrid0, Op%REigvec(:,1), Basis)
  CALL P1D_BasisTOGrid_real(psigrid1, Op%REigvec(:,2), Basis)
  CALL P1D_BasisTOGrid_real(psigrid2, Op%REigvec(:,3), Basis)

  DO i = 1, Basis%nq
    WRITE(niomat, *) Basis%x(i), psigrid0(i)*219474.631443_Rkind, &
                   & psigrid1(i)*219474.631443_Rkind, psigrid2(i)*219474.631443_Rkind
  END DO

  
END PROGRAM 