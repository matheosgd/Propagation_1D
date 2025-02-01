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


!-----------------------------SYSTEM INITIALIZATION----------------------------
  WRITE(out_unit, *) "Nothing so far !"


END PROGRAM 