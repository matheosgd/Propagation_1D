MODULE Perturbation_m
  USE QDUtil_m
  IMPLICIT NONE


  CONTAINS


  SUBROUTINE P1D_Correction_energie_ordre_1(eps_n_1, W, phi_n_0)
    real(kind=Rkind),    intent(inout) :: eps_n_1                                   ! correction de l'énergie de l'état n NON DEGENERE à l'ordre 1
    complex(kind=Rkind), intent(in)    :: W(:,:)                                    ! mat de la perturbation sur la base des sinus
    real(kind=Rkind),    intent(in)    :: phi_n_0(:)                                ! nieme vecteur propre de l'H d'ordre 0 sur la base des sinus

    eps_n_1 = dot_product(phi_n_0, matmul(W, phi_n_0))

  END SUBROUTINE

    
END MODULE
