!------------------------------------------------------------
!
!Calculate ham and check symmetry
!
!------------------------------------------------------------




subroutine hamk(k,hk)
    use compara
    use constants,only:zi,zero
    implicit none
    integer :: i,j,ir

    real(acc),intent(in) :: k(3)
    real(acc)  :: r(3)

    complex(acc),intent(out)  :: hk(nbands,nbands)
    real(acc) :: pij(3,nbands,nbands)

    hk=zero
    do ir=1,nrpts
        do i=1,nbands
            do j=1,nbands
                !r=wsmn(i,j,:,ir)
                !hk(i,j)=hk(i,j)+hmnr(i,j,ir)*exp(2_acc*pi*zi*dot_product(r,k))
                r=wsmn(j,i,:,ir)
                hk(j,i)=hk(j,i)+hmnr(j,i,ir)*exp(2_acc*pi*zi*dot_product(r,k))
            enddo
        enddo
    enddo

end subroutine hamk


subroutine dhamk(k,axis,dhk)
    use compara
    use constants,only:zi,zero
    implicit none
    integer :: i,j,ir

    real(acc),intent(in) :: k(3)
    integer,intent(in) :: axis

    real(acc)  :: r(3)
    complex(acc)  :: ra

    complex(acc),intent(out)  :: dhk(nbands,nbands)
    real(acc) :: pij(3,nbands,nbands)

    dhk=zero
    do ir=1,nrpts
        do i=1,nbands
            do j=1,nbands
                r=wsmn(i,j,:,ir)
                ra=zi*r(axis)
                dhk(i,j)=dhk(i,j)+hmnr(i,j,ir)*ra*exp(2_acc*pi*zi*dot_product(r,k))
            enddo
        enddo
    enddo
end subroutine dhamk


subroutine ptberry(nb,nocc,nw,wf,berr_c)
    use compara
    use constants,only:zzero,zi
    implicit none
    integer :: i,j,k
    integer :: m,n
    integer, intent(in) :: nb,nw,nocc
    complex(acc), intent(in) :: wf(nw,nocc,nb)
    complex(acc)  :: overlap(nocc,nocc),mover(nocc,nocc)
    complex(acc)  :: det
    real(acc),intent(out)  :: berr_c
    mover=zzero
    do j=1,nocc
        mover(j,j)=1.0_acc+0*zi
    enddo
    do i=1,nw-1
        !overlap=zzero
        do j=1,nocc
            do k=1,nocc
                overlap(j,k) = dot_product(wf(i,j,:),wf(i+1,k,:))
            enddo
        enddo
        mover=matmul(mover,overlap)
    enddo
    call detz(nocc,mover,det)
    berr_c = atan2(aimag(det),real(det))

    if (berr_c<-pi .or. berr_c > pi) then
        print *,berr_c
    endif

end subroutine ptberry










subroutine checksym
    use compara
    use constants,only:ten,zero
    implicit none
    integer :: i,j,k,m,n,num
    real(acc) :: x,y,z
    real(acc) :: kp1(3)
    real(acc) :: kp2(3)
    real(acc) :: mnorm
    real(acc) :: klist(ten,ten,ten,3)
    complex(acc)  :: hk1(nbands,nbands)
    complex(acc)  :: hk2(nbands,nbands)
    complex(acc)  :: tmp(nbands,nbands)
    complex(acc) :: inv_sym(nbands, nbands)
    complex(acc) :: commut(ten,ten,ten,nbands, nbands)

    commut=zero
    hk1=zero
    hk2=zero
    inv_sym=zero
    mnorm=zero
    num=0
    !print *,sym
    inv_sym=sym
    call inv(nbands,inv_sym)
   ! do i=1,nbands
   !    print '(14F8.2)',realpart(sym(i,:))
   ! enddo
  !  print *,''
    !do i=1,nbands
    !   print '(14F8.2)',realpart(inv_sym(:,i))
    !enddo
    
    do i=1,ten
        do j=1,ten
            do k=1,ten
                x=float(i)/(ten+1)
                y=float(j)/(ten+1)
                z=float(k)/(ten+1)
               ! kp1(1)=y+z
               ! kp1(2)=x+z
               ! kp1(3)=x+y
                kp1(1)=x
                kp1(2)=y
                kp1(3)=z
                kp2=matmul(lacsym,kp1)
                kp1=matmul(kp1,centering)
                kp2=matmul(kp2,centering)
              !  print '(3F8.2)',kp1
              !  print '(3F8.2)',kp2
              !  kp2(1)=-y+z
              !  kp2(2)=x+z
              !  kp2(3)=x-y
                call hamk(kp1,hk1)
                call hamk(kp2,hk2)
                !print *,hk1
                tmp=matmul(inv_sym,hk1)
          !      print *,tmp
                tmp=matmul(tmp,sym)
                tmp=tmp-hk2
                commut(i,j,k,:,:)=(tmp)
            enddo
        enddo
    enddo
    do i=1,ten
        do j=1,ten
            do k=1,ten
                do m=1,nbands
                    do n=1,nbands
                        mnorm=mnorm+abs(commut(i,j,k,m,n))
                        num=num+1
                    enddo
                enddo
            enddo
        enddo
    enddo
    !mnorm=mnorm/num



    write(*,'("Commutation check",f22.15," ")')  (mnorm)
    !write(*,'("Commutation check",f22.15')  (mnorm)
    if (mnorm > 0.00001_acc) then
        stop 'Symetry WORNG! Cannot Do Further Calculation'
    endif
   


end subroutine checksym



subroutine invarianthamk(k,hk)
    use compara
    use constants,only:zi,zero
    implicit none
    integer :: i,j,ir

    real(acc),intent(in) :: k(3)
    real(acc)  :: r(3)

    complex(acc),intent(out)  :: hk(nbands,nbands)
    real(acc) :: pij(3,nbands,nbands)

    hk=zero
    do ir=1,nrpts
        do i=1,nbands
            do j=1,nbands
                r=wsmn(i,j,:,ir)
                hk(i,j)=hk(i,j)+hmnr(i,j,ir)*exp(2_acc*pi*zi*dot_product(r,k))
            enddo
        enddo
    enddo

end subroutine invarianthamk



subroutine testham
    use compara
    implicit none
    integer :: i
    real(acc) :: k(3)
    real(acc) :: k2(3)
    complex(acc) :: inv_sym(nbands, nbands)
    complex(acc)  :: hk(nbands,nbands)
    complex(acc)  :: hk2(nbands,nbands)
    complex(acc)  :: tmp(nbands,nbands)

    k(1)=.1_acc
    k(2)=.2_acc
    k(3)=.4_acc
    k2(1)=.3_acc
    k2(2)=-.2_acc
    k2(3)=.4_acc
    call hamk(k,hk)
    call hamk(k2,hk2)
    inv_sym=sym
    call inv2(nbands,inv_sym)

    tmp=matmul(inv_sym,hk)
    tmp=matmul(tmp,sym)
    tmp=tmp-hk2
    print *,'re'
    do i=1,nbands
            print '(8F16.8)',realpart(sym(i,:))
    enddo
    print *,'re'
    do i=1,nbands
            print '(8F16.8)',realpart(hk(i,:))
    enddo
    print *,'im'
    do i=1,nbands
            print '(8F16.8)',imagpart(hk(i,:))
    enddo

end subroutine testham





















































