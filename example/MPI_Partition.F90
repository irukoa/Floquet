module MPI_Partition

  implicit none
  private

  public :: MPI_task_partition

contains

  subroutine MPI_task_partition(task_size, nodes, counts, displs)

    integer, intent(in) :: task_size, nodes
    integer, allocatable, intent(out) :: counts(:), &
                                         displs(:)

    integer :: ratio, remainder, i

    character(len=1024) :: errormsg
    integer :: istat

    if (task_size < 1) error stop "Floquet: Error #1: task_size must be a positive integer."
    if (nodes < 1) error stop "Floquet: Error #1: nodes must be a positive integer."

    allocate (counts(0:nodes - 1), displs(0:nodes - 1), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "Floquet: Error #3: failure allocating counts and displs. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    ratio = task_size/nodes
    remainder = modulo(task_size, nodes)

    do i = 0, nodes - 1
      if (i < remainder) then
        counts(i) = ratio + 1
        displs(i) = i*(ratio + 1)
      else
        counts(i) = ratio
        displs(i) = remainder*(ratio + 1) + (i - remainder)*ratio
      end if
    end do

  end subroutine MPI_task_partition

end module MPI_Partition
