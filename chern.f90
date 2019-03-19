!---------------------------------------------------------
!
!calculate different chern number
!Zhang Zeying
!
!---------------------------------------------------------


subroutine firstchern
    !Physical Review B 74, 195118
    !not avilable yet
    use compara
    use constants,only:zero
    implicit none
    integer :: i,j,k,m,n
    integer :: occ, nkpt1,nkpt2
    complex(acc)  :: vx(nbands,nbands),vxd(nbands,nbands)
    complex(acc)  :: vy(nbands,nbands),vyd(nbands,nbands)
    complex(acc)  :: vz(nbands,nbands)
    complex(acc)  :: hk(nbands,nbands)
    complex(acc)  :: vec(nbands,nbands),vecd(nbands,nbands)
    complex(acc)  :: invec(nbands,nbands)
    complex(acc)  :: tr
    real(acc) :: kv2(3)
    real(acc) :: kp(3)
    real(acc) :: eig(nbands),berry,chern
    real(acc) :: eigdiag(nbands,nbands)
    real(acc) :: kx,ky,dkx,dky
    nkpt1=100
    nkpt2=100
    vx=zero
    vy=zero
    vz=zero
    kx=zero
    ky=zero
    kp=zero
    tr=0
    eigdiag=zero
    occ=7
    chern=zero
    dkx=float((1))/nkpt1
    kp(3)=0.5_acc
    dky=float((1))/nkpt2
    do i=1,nkpt1
        kx= zero + float((i-1))/nkpt1
        do j=1,nkpt2
            kx= zero + float((j-1))/nkpt2
            kp(1)=kx
            kp(2)=ky
            call hamk(kp,hk)
            vec=hk
            call eigensystem('V',nbands,vec,eig)
            call dhamk(kp,2,vx)
            call dhamk(kp,3,vy)
            vx=matmul(vx,vec)
            vy=matmul(vy,vec)
            vecd=transpose(conjg(vec))
            vx=matmul(vecd,vx)
            vy=matmul(vecd,vy)
            vxd=vx
            vyd=vy
            do m=1,nbands
                do n=1,nbands
                    !if (m .ne. n) then
                    !if (eig(n) > 0d0 .and. eig(m)< 0d0) then
                     if (m> occ .and. n<= occ) then 
                        vx(n,m) = vx(n,m)/(eig(m)-eig(n))
                        vy(n,m) = vy(n,m)/(eig(m)-eig(n))
                    endif
                    !if (m .ne. n) then
                    !if (eig(m) > 0d0 .and. eig(n)< 0d0) then
                     if (n> occ .and. m<= occ) then 
                        vxd(n,m) = vx(n,m)/(eig(m)-eig(n))
                        vyd(n,m) = vy(n,m)/(eig(m)-eig(n))
                    endif
                    if (n.eq.m) then
                        vx(n,m) = zero
                        vy(n,m) = zero
                    endif
                enddo
            enddo
            vx=matmul(vec,matmul(vx,vecd))
            vxd=matmul(vec,matmul(vx,vecd))
            vy=matmul(vecd,matmul(vy,vec))
            vyd=matmul(vecd,matmul(vyd,vec))
            berry=zero
            vz=matmul(vx,vyd)

            call trace(nbands,vz,tr)
            berry=aimag(tr)
            chern=chern+berry
            print *,chern
        enddo
    enddo
    print *,chern
end subroutine firstchern


subroutine chernfirst
    use compara
    use constants,only:zzero,zero,one,zero3
    implicit none
    integer :: i,j,k
    integer :: nb,nocc,nw,m,n
    integer :: nkpt1,nkpt2
!    real(acc) :: kx,ky,dkx,dky,length
    real(acc) :: kp(3)
    real(acc) :: eig(nbands)
    complex(acc)  :: hk(nbands,nbands)
    complex(acc), allocatable :: wf(:, :, :)
    complex(acc), allocatable :: allwf(:,:,:, :)
    real(acc) :: berr_c,chern
!    real(acc) :: center(2)
    real(acc) :: space(2,3)
    

    space(1,:)=(/zero,one,one/)
    space(2,:)=(/one,zero,zero/)
    nb=nbands
    nw=5
    nkpt1=50
    nkpt2=50
    chern=zero
!    length=one
!    center(1)=zero
!    center(2)=zero

    nocc=0
    call hamk(zero3,hk)
    call getnocc(nbands,hk,nocc)
    print 110,nocc
110 format ('Occ num for  ham ',I3)

    allocate(wf(nw,nocc,nb))
    allocate(allwf(nkpt1,nkpt2,nocc,nb))
    allwf = zzero
    wf=zzero
!    dkx=float((1))/nkpt1
!    kp(3)=0.0_acc
!    dky=float((1))/nkpt2
    do i=1,nkpt1
        !kx= center(1)+length*(-.5_acc+float(i-1)/float(nkpt1-1))
        do j=1,nkpt2
            kp=((float(i-1)/float(nkpt1-1)))*space(1,:)&
                +((float(j-1)/float(nkpt2-1)))*space(2,:)
         !   ky= center(2)+length*(-.5_acc+float(j-1)/float(nkpt2-1))
           ! kp(1)=kx
           ! kp(2)=ky
            !print '(3F8.2)',kp
            call hamk(kp,hk)
            call eigensystem('V',nbands,hk,eig)
            allwf(i,j,:,:)=transpose(hk(:,:nocc))
        enddo
    enddo
    do i=1,nkpt1-1
        do j=1,nkpt2-1
            wf(1,:,:)=allwf(i,j,:,:)
            wf(2,:,:)=allwf(i+1,j,:,:)
            wf(3,:,:)=allwf(i+1,j+1,:,:)
            wf(4,:,:)=allwf(i,j+1,:,:)
            wf(5,:,:)=allwf(i,j,:,:)
            call ptberry(nb,nocc,nw,wf,berr_c)
            chern=chern+berr_c
!            print *,chern
        enddo
    enddo
    print 200,chern
200 format ('Chern = ', F12.6)
end subroutine chernfirst



subroutine chernmirr
    use compara
    use constants,only:eye3,zone,eps8,one,zzero,zero,eps2,zero3
    implicit none
    real(acc)   ::      nulspace(3,3)
    real(acc)   ::      invcenter(3,3)
    integer     ::      dimnul,nocc,nw
    integer     ::      nkpt1,nkpt2
    integer     ::      i,j,k
    integer     ::      mirrbands
    complex(acc)   ::   mirreig(nbands,nbands),mirreigdag(nbands,nbands)
    !complex(acc),allocatable   ::   eig(:)
    real(acc),allocatable   ::   eig(:)
    complex(acc)   ::   mirreval(nbands)
    complex(acc)  :: hk(nbands,nbands)
    complex(acc),allocatable  :: mirrhk(:,:)
    complex(acc), allocatable :: wf(:, :, :)
    complex(acc), allocatable :: allwf(:,:,:, :)
    real(acc) :: chern,berr_c
    real(acc) :: kp(3)
    integer        ::   eignocc
!    real(acc)   ::      sortwcc(3,nbands)
    mirreig=zzero
    mirreval=zzero
    nulspace=zero
    dimnul=0
    nw=5
    mirrbands=nbands/2
    invcenter=centering
    nkpt1=30
    nkpt2=30
    chern=zero

    call invr(3,invcenter)
    invcenter=matmul(matmul(centering,lacsym),invcenter)-eye3

!    print '(3F8.2)',invcenter
    call nullspace(3,invcenter,dimnul,nulspace)
    if (dimnul.ne.2) then
        print *,'not mirr symmetry'
        return 
    endif


    !mk pbc
    do i=1,2
        do j=1,3
            if ((dabs(nulspace(i,j)).gt.eps2) .and. ( nulspace(i,j).gt.eps2))then
!                print *,sign(nulspace(i,j),one)
                nulspace(i,j)=one
            else if  ((dabs(nulspace(i,j)).gt.eps2) .and. ( nulspace(i,j).lt.-eps2))then
                nulspace(i,j)=-one

            endif
        enddo
    enddo
!    print *,dimnul
    print 120,nulspace(1,:)
    print 121,nulspace(2,:)
120 format ('First  mirr plane vec',3F12.6)
121 format ('Second mirr plane vec',3F12.6)
!    print '(3F8.2)',nulspace(1,:)
!    print '(3F8.2)',nulspace(2,:)
    call eigensystem_non(nbands,sym,mirreval,mirreig)
    call sortmirreigsys(nbands,mirreval,mirreig)
    mirreigdag=transpose(conjg(mirreig))
   ! print '(14F8.3)',realpart(mirreig)



    hk=zzero
    eignocc=0
    call hamk(zero3,hk)
    call getnocc(nbands,hk,eignocc)
    nocc=eignocc/2
    print 110,nocc
110 format ('Occ num for mirr ham ',I3)
    allocate(eig(mirrbands))
    allocate(mirrhk(mirrbands,mirrbands))
    allocate(wf(nw,nocc,mirrbands))
    allocate(allwf(nkpt1,nkpt2,nocc,mirrbands))


    do i=1,nkpt1
        do j=1,nkpt2
            kp=((float(i-1)/float(nkpt1-1)))*nulspace(1,:)&
                +((float(j-1)/float(nkpt2-1)))*nulspace(2,:)
!            kp=1.0_acc*kp
!            print '(2I8)',i,j
!            print '(3F8.2)',kp
            call hamk(kp,hk)
            hk=matmul(hk,mirreig)
            hk=matmul(mirreigdag,hk)
            !mirrhk=hk(mirrbands+1:,mirrbands+1:)
            mirrhk=hk(:mirrbands,:mirrbands)
!            call checkhermit(mirrbands,mirrhk)
!            print '(7F8.2)',realpart(mirrhk)
!            print '(14F8.2)',realpart(hk)
            call eigensystem('V',mirrbands,mirrhk,eig)
!            print '(7F8.2)',(eig)
!            print '(7F8.2)',realpart(mirrhk)
!            print *,
            allwf(i,j,:,:)=transpose(mirrhk(:,:nocc))
!            print *,real(transpose(mirrhk(:,:nocc)))
        enddo
    enddo
    do i=1,nkpt1-1
        do j=1,nkpt2-1
            wf(1,:,:)=allwf(i,j,:,:)
            wf(2,:,:)=allwf(i+1,j,:,:)
            wf(3,:,:)=allwf(i+1,j+1,:,:)
            wf(4,:,:)=allwf(i,j+1,:,:)
            wf(5,:,:)=allwf(i,j,:,:)
            call ptberry(mirrbands,nocc,nw,wf,berr_c)
!            print '(1F12.6)',berr_c
            chern=chern+(berr_c)
!            if (i.eq.1 .or. j.eq.1) then
!                print *,berr_c
!            endif
        enddo
    enddo
!    print '('Mirr',1F12.6)',chern/(2*pi)
    print 100,chern/(2*pi)
100 format ('Mirror Chern = ', F12.6)




end subroutine chernmirr






















