!---------------------------------------------------------
!
!calculate different Z invriant
!Zhang Zeying
!
!---------------------------------------------------------



subroutine parityz2
    !only for inversion sym
    ! not avilible yet
    use compara
    use constants,only:zzero,zero,one,zero3,half,pauy,pau0,zi,zone
    implicit none
    integer :: i,j,k
    integer :: nocc,nelectocc,hb
    real(acc) :: eig(nbands)
    complex(acc) :: wf(nbands)
    complex(acc)  :: pari
    complex(acc)  :: hk(nbands,nbands)
    complex(acc)   :: trs(nbands,nbands)
    real(acc)     :: trm(8,3)
    real(acc)     :: trmp(8)
    complex(acc),allocatable  :: eyehb(:,:)

    complex(acc),allocatable  :: treig(:)
    complex(acc),allocatable  :: trhk(:,:)
    complex(acc),allocatable   ::   trvec(:,:),trvecdg(:,:)
    
    hb=nbands/2
    allocate(eyehb(hb,hb))
    allocate(trhk(hb,hb))

    allocate(treig(nbands))
    allocate(trvec(nbands,nbands))
    allocate(trvecdg(nbands,nbands))
    eyehb=zzero
    trs=zzero
    trvec=zzero
    treig=zzero
    do i=1,hb
        eyehb(i,i)=zone
    enddo
    call kproduct(hb,2,nbands, eyehb, zi*pauy,trs)

!    print '(14F8.2)',real(trs)

    call eigensystem_non(nbands,trs,treig,trvec)
    call sortmirreigsys(nbands,treig,trvec)
    trvecdg=transpose(conjg(trvec))
!    print '(14F8.3)',imagpart(treig)
    

    nocc=0
    nelectocc=2
    hk=zzero
    call hamk(zero3,hk)
    call getnocc(nbands,hk,nocc)
    print 110,nocc
110 format ('Occ num for  ham ',I3)
    
    trm(1,:)=zero3
    trm(2,:)=(/half,zero,zero/)
    trm(3,:)=(/zero,half,zero/)
    trm(4,:)=(/zero,zero,half/)
    trm(5,:)=(/zero,half,half/)
    trm(6,:)=(/half,zero,half/)
    trm(7,:)=(/half,half,zero/)
    trm(8,:)=(/half,half,half/)

    do i=1,8
        hk=zero
        trmp(i)=one
        print *,trm(i,:)
        call hamk(trm(i,:),hk)
        hk=matmul(trs,hk)
        hk=matmul(transpose(conjg(hk)),trs)
        print '(14F8.3)',realpart(hk)
        do j=nelectocc+1,nocc,2
            wf=hk(:,j)
      !      wfv=hkv(:,j)
            !wfdg=conjg(wfdg)
         !   pari=dot_product(wfv,wf)
            pari=dot_product(wf,matmul(sym,wf))
            trmp(i)=trmp(i)*pari
          !  print *,real(pari)
          !  print *,eig(j)
        enddo
     !   print *,
     !   print *,trmp(i)
    enddo


end subroutine parityz2




























