!----------------------------------------------------
!
!some useful functions
!Zhang Zeying
!
!---------------------------------------------------




subroutine invr(n,matrix)
    use compara,only:acc
    implicit none
    integer info,lwork,lda
    integer,intent(in):: n
    integer ::ipiv(n)
    real(acc)  :: work(n)
    real(acc),intent(in) ::matrix(n,n)
    lwork=n*n
    lda=n
    call dgetrf(n,n,matrix,lda,ipiv,info)
    if (info.ne.0) then
        stop 'zgetrf error'
    endif

    call dgetri(n,matrix,lda,ipiv,work,lwork,info)
    if (info.ne.0) then
        stop 'zgetri error'
    endif
end subroutine invr


subroutine inv(n,matrix)
    use compara,only:acc
    implicit none
    integer info,lwork,lda
    integer,intent(in):: n
    integer ::ipiv(n)
    complex(acc)  :: work(n)
    complex(acc),intent(in) ::matrix(n,n)
    lwork=n*n
    lda=n
    call zgetrf(n,n,matrix,lda,ipiv,info)
    if (info.ne.0) then
        stop 'zgetrf error'
    endif

    call zgetri(n,matrix,lda,ipiv,work,lwork,info)
    if (info.ne.0) then
        stop 'zgetri error'
    endif
end subroutine inv


subroutine inv2(n,matrix)
    ! from wann_tools
    use compara,only:acc
    implicit none
    integer,parameter :: dp=8
    integer :: i,info
    integer,intent(in):: n
    integer,allocatable   :: ipiv(:)
    complex(dp),parameter :: zone=(1.0_dp,0.0_dp)
    complex(dp),intent(inout):: matrix(n,n)
    complex(dp),allocatable :: eye(:,:)
    allocate(ipiv(n))
    allocate(eye(n,n))
    ipiv=0

    eye=(0.0_dp, 0.0_dp)
    do i=1,n
      eye(i,i)= zone
    enddo
    call zgesv(n,n,matrix,n,ipiv,eye,n,info)
    matrix=eye
    if (info.ne.0) then  
            stop 'Inverse Matrix Failed'
    endif
end subroutine inv2



subroutine norm(n,matrix,mnorm)
    use compara,only:acc
    implicit none
    integer :: i,j
    integer,intent(in):: n
    complex(acc),intent(in) ::matrix(n,n)
    real(acc),intent(out):: mnorm
    mnorm=0_acc
    do i=1,n
      do j=1,n
        mnorm=mnorm+zabs(matrix(i,j))
      enddo
    enddo
end subroutine norm


subroutine eigensystem(jobz,n,matrix,eig)
    !for symmetric matrix
    use compara,only:acc
    implicit none
    character*1, intent(in) :: jobz ! N only eigenvalues, V eigenvalues and eigenvectors
    character*1 :: uplo
    integer,   intent(in) :: n
    complex(acc),intent(inout) :: matrix(n,n)
    integer   :: lda
    real(acc), intent(inout) :: eig(n)
    complex(acc), allocatable ::  work(:)
    integer   :: lwork
    real(acc), allocatable ::  rwork(:)
    integer :: info
    allocate(rwork(16*n))
    allocate( work(16*n))
    uplo='U'
    lda=n
    eig=0.0_acc
    work= 0_acc
    lwork=16*n

    call zheev(jobz, uplo,n,matrix,lda,eig,work,lwork,rwork,info)
    if (info.ne.0) then
            stop 'Digonalized Matrix Failed'
    endif
end subroutine eigensystem



subroutine checkhermit(n,matrix)
    use compara,only:acc
    use constants,only:eps5
    implicit none
    integer,   intent(in) :: n
    complex(acc),intent(in) :: matrix(n,n)
    integer                 :: i,j
    do i=1,n
        do j=i,n
            if (real(matrix(i,j)-real(matrix(j,i)))>eps5 .and.&
                aimag(matrix(i,j)+aimag(matrix(j,i)))>eps5) then
                print '(4F12.6)',matrix(i,j),matrix(j,i)
                print *,'not Hermite'
                return
            endif
        enddo
    enddo
    print *,'Hermite'
end subroutine checkhermit



subroutine eigensystem_p(jobz,n,matrix,eig)
    !remove the random phase,for symmetric matrix
    use compara,only:acc,pi
    use constants,only:zzero,zi,zero
    implicit none
    character*1, intent(in) :: jobz ! N only eigenvalues, V eigenvalues and eigenvectors
    character*1 :: uplo
    integer,   intent(in) :: n
    complex(acc),intent(inout) :: matrix(n,n)
    integer   :: lda
    real(acc), intent(inout) :: eig(n)
    complex(acc), allocatable ::  work(:)
    integer   :: lwork
    real(acc), allocatable ::  rwork(:)
    integer ::   lrwork
    integer, allocatable ::  iwork(:)
    integer ::  liwork
    integer :: info
    integer :: test,i,j
    complex(acc) :: phase
    real(acc) :: rphase
    test=16
    allocate(rwork(test*n))
    allocate( work(test*n))
    allocate( iwork(test*n))
    uplo='U'
    lda=n
    eig=0.0_acc
    work= 0_acc
    iwork= 0_acc
    lwork=test*n
    lrwork=test*n*n
    liwork=test*n
    phase=zzero
    call zheevd(jobz, uplo,n,matrix,lda,eig,work,lwork,rwork,lrwork,iwork,liwork,info)
    if (info.ne.0) then
        stop 'Digonalized Matrix Failed'
    endif
    do i=1,n
        matrix(:,i)=matrix(:,i)/matrix(n,i)
    enddo
    do i=1,n
    phase=zzero
        do j=1,n
            phase=phase+matrix(j,i)*conjg(matrix(j,i))
        enddo
        rphase=sqrt(real(phase))
        matrix(:,i)=matrix(:,i)/rphase
        !print *,i,j,rphase
    enddo
  !  do i=1,n
  !  phase=zzero
  !      do j=1,n
  !          phase=phase+matrix(j,i)*conjg(matrix(j,i))
  !      enddo
  !      rphase=real(phase)
  !      print *,rphase
  !  enddo

end subroutine eigensystem_p

subroutine eigensystem_non(n,matrix,eval,evecr)
    !for non symmetric matrix 
    use compara,only:acc
    implicit none
    integer,   intent(in)       :: n
    complex(acc),intent(in)     :: matrix(n,n)
    complex(acc)                :: a(n,n)
    integer                     :: lda
    complex(acc), intent(inout)    :: eval(n)
    complex(acc)                :: evecl(n,n),evecr(n,n)
    integer                     :: ldvl,ldvr
    complex(acc), allocatable   :: work(:)
    integer                     :: lwork
    real(acc), allocatable      :: rwork(:)
    integer                     :: info
    lwork=16*n
    allocate(rwork(lwork))
    allocate( work(lwork))
    lda=n
    work= 0_acc
    ldvl=n
    ldvr=n
    a=matrix
    call zgeev('N', 'V', n, a, lda, eval, evecl, ldvl, evecr, ldvr, work, lwork, rwork, info)
    if (info.ne.0) then
            stop 'Digonalized Matrix Failed'
    endif

end subroutine eigensystem_non

subroutine sortmirreigsys(n,eval,evec)
    !form zsteqr.f
    !sort the mirr eigvals by imagpart
    use compara,only:acc
    implicit none
    integer,   intent(in)       :: n
    complex(acc), intent(inout) :: eval(n)
    complex(acc)                :: p
    complex(acc),intent(inout)  :: evec(n,n)
    integer                     :: i,j,k
    integer                     :: ii
    if (n.eq.1)then
        return
    endif
    do ii=2,n
        i=ii-1
        k=i
        p=eval(i)
        do j=ii,n
            if (aimag(eval(j)).lt.(aimag(p))) then
                k=j
                p=eval(j)
            endif
        enddo
        if (k.ne.i)then
            eval(k)=eval(i)
            eval(i)=p
            call zswap(n,evec(1,i),1,evec(1,k),1)
        endif
    enddo
end subroutine sortmirreigsys

subroutine getnocc(n,hk,nocc)
    !get num of eigval lower than 0 
    use compara,only:acc
    use constants,only:zero3,zero
    implicit none
    integer,   intent(in)       :: n
    complex(acc), intent(inout)    :: hk(n,n)
    integer,    intent(out)    :: nocc
    real(acc)                   :: eig(n)
    integer                     :: i
    call eigensystem('N',n,hk,eig)
    nocc=0
    do i=1,n
        if (eig(i).le.zero) then
            nocc=nocc+1
        endif
    enddo
!    print *,nocc
end subroutine getnocc



subroutine eigensystem_select(jobz,n,matrix,eig)
    ! not avilible yet
    use compara,only:acc
    character*1, intent(in) :: jobz
    character*1  :: ran
    character*1 :: uplo
    integer,   intent(in) :: n
    complex(acc),intent(inout) :: matrix(n,n)
    integer   :: lda
    real(acc) ::  vl,vu
    integer :: il,iu
    real(acc) :: abstol
    integer :: m
    real(acc), intent(inout) :: eig(n)
    complex(acc), allocatable ::  z(:,:)
    integer :: ldz
    complex(acc), allocatable ::  work(:)
    integer   :: lwork
    real(acc), allocatable ::  rwork(:)
    integer,allocatable     :: ifail(:),iwork(:)
    integer :: info
    allocate(rwork(16*n))
    allocate( work(16*n))
    rwork=0_acc
    work=0_acc
    call zheevx('V','A','U',n,matrix,n,vl,vu,il,iu,abstol,&
            m,eigs,z,ldz,work,lwork,rwork,iwork, ifail,info)
    if (info.ne.0) then
            stop 'Digonalized Matrix Failed'
    endif
end subroutine eigensystem_select


subroutine trace(n,matrix,tr)
    use compara, only: acc
    implicit none
    integer :: i
    integer, intent(in) :: n
    complex(acc), intent(in) :: matrix(n,n)
    complex(acc), intent(out) :: tr
    tr=0_acc
    do i=1,n
        tr= tr + matrix(i,i)
    enddo
end subroutine trace



subroutine detz(n,matrix,det)
    use compara, only: acc
    use constants,only:zi
    implicit none
    complex(acc), intent(inout) :: matrix(n, n)
    complex(acc), intent(out) :: det
    integer :: i,info, n
    integer :: ipiv(n)
!    allocate(ipiv(n))
    call zgetrf(n, n, matrix, n, ipiv, info)
    if(info.ne. 0) then
        stop 'LU factorization error'
    endif
    det=1.0_acc + 0*zi
    do i=1,n
        if (ipiv(i) .ne. i) then
            det = -det*matrix(i,i)
        else
            det = det*matrix(i,i)
        endif
    enddo
    !print '(8F8.2)',matrix
   ! print '(8F8.2)',det
end subroutine detz



subroutine nullspace(n,mat,dimnul,nulspace)
    use constants,only:eps5
    implicit none
    integer,parameter :: dp=8
    integer,intent(in) :: n
    real(dp),intent(in) :: mat(n,n)
    real(dp)  :: a(n,n)
    real(dp)  :: s(n)
    real(dp), allocatable :: work(:)
    real(dp) :: u(n,n)
    real(dp)  :: vt(n,n)
    integer :: ldvt
    integer :: info, lwork,lw
    !real(dp), intent(out), allocatable :: nulspace(:,:)
    real(dp), intent(out) :: nulspace(n,n)
    integer :: i
    integer, intent(out) :: dimnul
    a=mat

!    print '(3F8.2)',a
    lw=1000
    allocate(work(lw))
    ldvt=n
    call dgesvd('A', 'A', n, n, a, n, s, u, n, vt, ldvt, work, -1, info)
    lwork = min(lw,int(work(1)))
!    print *,info
    !deallocate(work)
    !allocate(work(lwork))
    call dgesvd('A', 'A', n, n, a, n, s, u, n, vt, ldvt, work, lwork, info)
!    print *,info

    dimnul=0
    do i=1,n
        if (dabs(s(i))<eps5) then
            dimnul=dimnul+1
        endif
    enddo
!    print *,dimnul
!    allocate(nulspace(3,3))
    dimnul=0
    nulspace=0.0_dp
    do i=1,n
!        print *,s(i,i)
        !if (s(i,i).eq.0.0_dp) then
        if (dabs(s(i))<eps5) then
            dimnul=dimnul+1
            nulspace(dimnul,:)=vt(i,:)
        endif
    enddo
    deallocate(work)


 !   print '(3F22.12)',s
    !print *,dimnul
!    print '(3F8.2)',vt
    
    
end subroutine nullspace





subroutine kproduct(ra,rb,rab,a,b,ab)
    !from https://rosettacode.org/wiki/Kronecker_product#Fortran
    use compara, only: acc
    implicit none
    integer, intent(in) :: ra,rb,rab
    complex(acc),intent(in)        :: a(ra,ra),b(rb,rb)
    complex(acc)                :: ab(rab,rab)
    integer r,c,i,j
    r = 0
    do i = 1,ra
        c = 0
        do j = 1,ra
            ab(r + 1:r + rb,c + 1:c + rb) = a(i,j)*b
            c = c + rb
        end do
        r = r + rb
    end do
end subroutine kproduct

subroutine pfaffian(n,matrix)
    use compara, only: acc
    implicit none
    integer,intent(in)          :: n
    complex(acc),intent(in)     :: matrix(n,n)




end subroutine pfaffian



























