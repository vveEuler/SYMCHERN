





program main
    use compara
    implicit none
    integer :: i
    real(acc) :: start, finish
    call cpu_time(start)
    call readsym
    call readhr
!    call testham
!    call checksym
    call chernfirst
!    call chernmirr
!    call parityz2
    call cpu_time(finish)
    write(*,'("Time = ",f15.3," seconds.")')  (finish-start)

end program

