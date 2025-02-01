MODULE P1D_lanczos_m
  IMPLICIT NONE
  

  PUBLIC 
  PRIVATE :: P1D_Gram_schmidt
  

  CONTAINS


  SUBROUTINE P1D_Base_krylov(K, H, psi, m)
    USE QDUtil_m
    IMPLICIT NONE
  
    complex(kind=Rkind), allocatable, intent(inout) :: K(:,:)
    real(kind=Rkind),                 intent(in)    :: H(:,:)
    complex(kind=Rkind),              intent(in)    :: psi(:)                       ! sur la base des sinus => dim = nb
    integer,                          intent(in)    :: m

    complex(kind=Rkind), allocatable                :: Q(:,:), triband_H(:,:)
    real(kind=Rkind),    allocatable                :: S(:,:), identite(:,:)
    integer                                         :: im
  
    !-----------------------Creation of the Krylov basis-----------------------
    ALLOCATE(Q(size(psi), m))

    Q(:,:) = CZERO
    Q(:,1) = Psi(:)
  
    DO im = 2, m
      Q(:,im) = matmul(H, Q(:, im-1))
    END DO
  
    !----------------------------Orthonormalization----------------------------
    CALL P1D_Gram_schmidt(K, Q)
   
    !--------------------------------Check ortho-------------------------------
    ALLOCATE(identite(m,m))
    identite(:,:) = CZERO

    DO im = 1, m
      identite(im,im) = 1
    END DO

    ALLOCATE(S(m,m))
    S = matmul(conjg(transpose(K)), K)

    WRITE(out_unit, *) 'Smax = ', maxval(abs(S(:,:))-identite(:,:))
    
  END SUBROUTINE

 
  SUBROUTINE P1D_Base_krylov_plus(K, H)
    USE QDUtil_m
    IMPLICIT NONE
  
    complex(kind=Rkind), allocatable, intent(inout) :: K(:,:)
    real(kind=Rkind),                 intent(in)    :: H(:,:)

    complex(kind=Rkind), allocatable                :: Ktemp(:,:), V(:)
    real(kind=Rkind),    allocatable                :: S(:,:), identite(:,:)
    integer                                         :: m, nb, i

    !------------------------------Initialization------------------------------
    nb = size(K, dim = 1)
    m = size(K, dim = 2)
    Ktemp = K
    ALLOCATE(V(nb))
    DEALLOCATE(K)
    ALLOCATE(K(nb, m+1))
    ALLOCATE(identite(m+1,m+1))
    identite(:,:) = CZERO
    DO i = 1, m+1
      identite(i,i) = 1
    END DO
    
    !-------------------------------Construction-------------------------------
    V(:) = matmul(H,Ktemp(:,m)) 

    K(:,m+1) = V(:)
    DO i = 1, m
      K(:,m+1) = K(:,m+1) - dot_product(Ktemp(:,i), V(:))*Ktemp(:,i)
    END DO
    K(:,m+1) = K(:,m+1) / sqrt(dot_product(K(:,m+1),K(:,m+1)))

    K(:,1:m) = Ktemp
    DEALLOCATE(Ktemp, V)

    !--------------------------------Check ortho-------------------------------
    ALLOCATE(S(m+1,m+1))

    S = matmul(conjg(transpose(K)), K)
    WRITE(out_unit, *) 'Smax = ', maxval(abs(S(:,:))-identite(:,:))
    
  END SUBROUTINE


  SUBROUTINE P1D_Gram_schmidt(K, Q)
    USE QDUtil_m
    IMPLICIT NONE
    
    complex(kind=Rkind), allocatable, intent(inout) :: K(:,:)
    complex(kind=Rkind),              intent(in)    :: Q(:,:)

    complex(kind=Rkind), allocatable                :: v(:)
    integer                                         :: im, jm, n, m

    n = size(Q, dim=1)
    m = size(Q, dim=2)    
    ALLOCATE(K(n,m))
    ALLOCATE(v(n))
    K(:,:) = CZERO
    v(:) = CZERO
   
    K(:,1) = Q(:,1) / sqrt(dot_product(Q(:,1),Q(:,1)))
    
    DO im = 2, m
      v(:) = Q(:, im)
      DO jm = 1, im-1
        v(:) = v(:) - dot_product(Q(:,im), K(:,jm)) * K(:,jm)
      END DO
      K(:, im) = v(:) / sqrt(dot_product(v(:),v(:)))
    END DO
            
  END SUBROUTINE
  
  
  SUBROUTINE P1D_Construct_H_tribande(triband_H, K, H)
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: triband_H(:,:)
    complex(kind=Rkind), intent(in)    :: K(:,:)
    real(kind=Rkind),    intent(in)    :: H(:,:)

    integer                            :: ib, m

    !--------------------------------Allocation--------------------------------
    m = size(K, dim = 2)    
   
    !-------------------------------Construction-------------------------------
    triband_H = matmul(conjg(transpose(K)), matmul(H, K))
    !WRITE(out_unit, *) '------------------------------------'
    !CALL WRITE_Mat(triband_H,out_unit,5,info='triband_H') !ça prenait trop de place dans results
      
  END SUBROUTINE


  SUBROUTINE P1D_KrylovTOBasis(Vec_B, Vec_K, K)
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: Vec_B(:,:)                                ! mat des vecteurs colonnes def sur a base des sinus
    complex(kind=Rkind), intent(in)    :: Vec_K(:,:)                                ! mat des vecteurs colonnes def sur a base de Krylov
    complex(kind=Rkind), intent(in)    :: K(:,:)

    integer                            :: i, j, m

    m = size(K, dim = 2)
    Vec_B(:,:) = CZERO
    Vec_B = matmul(K, Vec_K)

  END SUBROUTINE


  SUBROUTINE P1D_Construct_psi_approx(Psidt, psi, Vec_B, Valp, dt)
    USE QDUtil_m
    IMPLICIT NONE

    complex(kind=Rkind), intent(inout) :: Psidt(:)
    complex(kind=Rkind), intent(in)    :: Vec_B(:,:), Valp(:), psi(:)
    real(kind=Rkind),    intent(in)    :: dt

    complex(kind=Rkind), allocatable   :: c1(:), c2(:)
    integer                            :: m, i
   
    m = size(Vec_B, dim = 2)
    ALLOCATE(c1(m), c2(m))

    c1 = matmul(conjg(transpose(Vec_B)), Psi)
    c2(:) = c1(:)*exp(-EYE*dt*Valp(:))
  
    WRITE(out_unit, *) 'coeff m-1 = ', c2(m-1), 'coeff m = ', c2(m)

    Psidt = matmul(Vec_B, c2)

  END SUBROUTINE
  

END MODULE
