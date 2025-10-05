!> Gauss-Jordan elimination
!> Ax=b
subroutine Gauss_Jordan(A,b,tol,x,error)
    use vectors_m
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> square matrix
    real(kind=8), intent(in) :: b(:) !> vector
    real(kind=8), intent(in) :: tol !> tolerance for solution
    real(kind=8), intent(out) :: x(:) !> solution of linear system
    integer(kind=4), intent(out), optional :: error !> Error flag: 0->No error, 1->Singular equations
    
    real(kind=8), parameter :: epsilon=1d-16
    real(kind=8) :: aux,factor
    real(kind=8), allocatable :: A_aux(:,:)
    integer(kind=4) :: n,icol,jcol,irow,jrow,ipeak
    n=size(A,1)
    if (n/=size(b)) then
        error stop "Dimension error in Gauss-Jordan"
    end if
    A_aux=A
    x=b
    !> Process n times to get all equations
    main_loop: do icol = 1,n
        !> Find peak pivot for column icol
        irow=icol
        ipeak = 1
        max_pivot: do jrow = 2,n
            if (abs(A_aux(jrow,icol))>abs(A_aux(ipeak,icol))) then
                ipeak = jrow
            end if
        end do max_pivot
        !> Check for singular equations
        if (abs(A_aux(ipeak,icol))<epsilon) then
            if (present(error)) then
                error=1
            end if
            return
        end if
        !> Otherwise, if ipeak /= irow, swap equations irow & ipeak
        if (ipeak/=irow .and. icol<n) then
            do jcol = 1,n
                aux=A_aux(ipeak,jcol)
                A_aux(ipeak,jcol)=A_aux(irow,jcol)
                A_aux(irow,jcol)=aux
            end do
            aux=x(ipeak)
            x(ipeak)=x(irow)
            x(irow)=aux
        end if
        !> Multiply equation irow by -A_aux(jrow,icol)/A_aux(irow,icol) and add it to eqn jrow
        do jrow = 1,n
            if (jrow/=irow) then
                factor=A_aux(jrow,icol)/A_aux(irow,icol)
                do jcol = 1,n
                    A_aux(jrow,jcol)=A_aux(jrow,jcol)-A_aux(irow,jcol)*factor
                end do
                x(jrow)=x(jrow)-x(irow)*factor
            end if
        end do
    end do main_loop
    !> End of main loop over all equations. All off-diagonal terms are now zero.
    !> To get the final answer, we must divide each equation by the coefficient of its diagonal term
    do irow = 1,n
        x(irow)=x(irow)/A_aux(irow,irow)
        A_aux(irow,irow)=1d0
    end do
    if (present(error)) then
        error=0
    end if
    if (inf_norm_vec_real(matmul(A,x)-b)>=tol) then
        print *, inf_norm_vec_real(matmul(A,x)-b)
        error stop "Wrong solution in Gauss-Jordan"
    end if
end subroutine

