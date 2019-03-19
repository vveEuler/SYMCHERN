


!Zhang Zeying


module compara
    implicit none
    integer,parameter :: hrfile=11
    integer,parameter :: symfile=12
    integer,parameter :: stdout=8
    integer,parameter :: acc=kind(1.0d0)
    !integer,parameter :: acc=kind(1.0)
    real(acc),parameter :: pi= 3.14159265359_acc
!    real(acc),parameter :: zero= 0_acc
!    complex(acc),parameter    :: zzero=(0.0_acc, 0.0_acc)
!    complex(acc),parameter    :: zi=(0.0_acc, 1.0_acc)
    
    integer :: nbands
    integer :: nrpts
    real(acc) :: centering(3,3)
    real(acc) :: lac(3,3)
    real(acc) :: lacsym(3,3)
    real(acc) :: trans(3)
    real(acc), allocatable :: wcc(:, :)
    complex(acc), allocatable :: sym(:, :)
    complex(acc), allocatable :: hmnr(:,:,:)
    real(acc), allocatable     :: wsmn(:,:,:,:)
end module compara



module constants
    use compara,only: acc
    implicit none
    integer,parameter           :: ten= 5
    integer,parameter           :: five= 5
    real(acc),parameter         :: pi= 3.14159265359_acc
    real(acc),parameter         :: twopi= 2*pi
    real(acc),parameter         :: zero= 0.0_acc
    real(acc),parameter         :: one= 1.0_acc
    real(acc),parameter         :: half= 0.5_acc
    complex(acc),parameter      :: zzero=(0.0_acc, 0.0_acc)
    complex(acc),parameter      :: zone=(1.0_acc, 0.0_acc)
    complex(acc),parameter      :: zi=(0.0_acc, 1.0_acc)
    real(acc),parameter         :: zero3(3)=(/zero,zero,zero/)
    real(acc),parameter         :: eye3(3,3)=reshape((/one,zero,zero,zero,one,zero,zero,zero,one/),(/3,3/))
    real(acc),parameter         :: eye2(2,2)=reshape((/one,zero,zero,one/),(/2,2/))
    real(acc),parameter         :: eps8  = 1.0e-8_acc
    real(acc),parameter         :: eps5  = 1.0e-5_acc
    real(acc),parameter         :: eps2  = 1.0e-2_acc

    complex(acc),parameter      :: pau0(2,2)=reshape((/zone,zzero,zzero,zone/),(/2,2/))
    complex(acc),parameter      :: paux(2,2)=reshape((/zzero,zone,zone,zzero/),(/2,2/))
    complex(acc),parameter      :: pauy(2,2)=reshape((/zzero,zi,-zi,zzero/),(/2,2/))
    complex(acc),parameter      :: pauz(2,2)=reshape((/zone,zzero,zzero,-zone/),(/2,2/))
end module constants







module para
    use compara
    implicit none
!    integer,parameter :: ten= 5
    real(acc), allocatable     :: ws(:,:)
end module para

module spin
    use compara 
    integer,parameter :: spinfile=13
    complex(acc), allocatable :: smnr(:,:,:,:)
end module spin




