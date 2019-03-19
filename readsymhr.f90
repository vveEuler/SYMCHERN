!------------------------------------------------------------
!
!file io 
!
!------------------------------------------------------------


subroutine readsym
    use compara
    use para,only:ws
    use constants,only:zi,zero,zzero
    implicit none
    integer :: i,j
    real(acc), allocatable :: symre(:, :)
    real(acc), allocatable :: symim(:, :)

    open(symfile, file='sym.out')
    read(symfile, *)
    read(symfile, *)nbands

    allocate(wcc(3,nbands))
    !allocate(wcc(nbands,3))
    allocate(sym(nbands,nbands))
    allocate(symre(nbands,nbands))
    allocate(symim(nbands,nbands))

    sym=zzero
    symre=zero
    symim=zero
    wcc=zero
    centering=zero


    read(symfile, *)
    do i=1,3
      read(symfile, *)centering(i,:)
    enddo

    read(symfile, *)
    do i=1,nbands
      read(symfile, *)wcc(:,i)
    enddo

    read(symfile, *)
    do i=1,3
      read(symfile, *)lacsym(i,:)
    !  print *,lacsym(:,i)
    enddo

    read(symfile, *)
    read(symfile, *)trans(:)

    read(symfile, *)
    do i=1,nbands
      !read(symfile, *) (symre(i,j),j=1,nbands)
      read(symfile, *) symre(i,:)
      !print '(14F8.2)',symre(:,i)
    enddo

    read(symfile, *)
    do i=1,nbands
      read(symfile, *) symim(i,:)
    enddo

    sym=symre+zi*symim




end subroutine readsym



subroutine readhr
    use compara
    use para
    use constants,only:zero,zzero
    implicit none
    integer :: tmp
    integer :: i,j,r
    integer :: m,n,r1,r2,r3
    real(acc) :: re, im
    !real(acc) :: wsmn(nbands,nbands,3,nrpts)
    real(acc) :: pij(3,nbands,nbands)


    open(hrfile, file='wannier90_hr.dat')
    read(hrfile, *)
    read(hrfile, *)
    read(hrfile, *)nrpts
    read(hrfile, *)(tmp, i=1, nrpts)

    allocate(ws(3,nrpts))
    allocate(hmnr(nbands,nbands,nrpts))

    hmnr=zzero
    pij=zero
    do r=1, nrpts
      do i=1, nbands
        do j=1, nbands
          read(hrfile,*)r1,r2,r3,m,n,re,im
         ! hmnr(n,m,r)=dcmplx(re , im)
          hmnr(m,n,r)=dcmplx(re , im)
        enddo
      enddo
      ws(:,r)=(/r1,r2,r3/)
     ! print *, (/r1,r2,r3/)
    enddo

    allocate(wsmn(nbands,nbands,3,nrpts))
    wsmn=zero


    do i=1,nbands
      do j=1,nbands
        pij(:,i,j)=-wcc(:,i)+wcc(:,j)
      enddo
    enddo
!    pij=zero
!    print '(3F8.2)',pij

    do r=1,nrpts
      do i=1,nbands
        do j=1,nbands
          !wsmn(i,j,:,r)=ws(:,r)+pij(:,i,j)
          wsmn(i,j,:,r)=ws(:,r)+pij(:,i,j)
          !wsmn(j,i,:,r)=ws(:,r)-pij(:,j,i)
         ! print *,r,i,j
        enddo
      enddo
    enddo
!    print *, hmnr(:,:,1)
end subroutine readhr


subroutine readspin
    use compara
    use constants,only:zero,zzero
!    use para
    implicit none
    integer :: tmp
    integer :: i,j,r
    integer :: m,n,r1,r2,r3
    real(acc) :: re, im
    !real(acc) :: wsmn(nbands,nbands,3,nrpts)
!    real(acc) :: pij(3,nbands,nbands)


    open(hrfile, file='wannier90_spin.dat')
    read(hrfile, *)
    read(hrfile, *)
    read(hrfile, *)nrpts
    read(hrfile, *)(tmp, i=1, nrpts)

    !allocate(ws(3,nrpts))
    allocate(hmnr(nbands,nbands,nrpts))

    hmnr=zero
!    pij=zero
    do r=1, nrpts
      do i=1, nbands
        do j=1, nbands
          read(hrfile,*)r1,r2,r3,m,n,re,im
          hmnr(m,n,r)=dcmplx(re , im)
        enddo
      enddo
      !ws(:,r)=(/r1,r2,r3/)
     ! print *, (/r1,r2,r3/)
    enddo

!    print *, hmnr(:,:,1)
end subroutine readspin











