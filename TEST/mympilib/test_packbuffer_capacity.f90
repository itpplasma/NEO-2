program test_packbuffer_capacity
  use mpi
  use packBuffer_module, only: packBuffer
  implicit none

  type(packBuffer) :: buffer
  double precision, allocatable :: matrix(:, :), received(:, :)
  integer :: ierr, rank, ranks, i, j, marker

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, ranks, ierr)
  if (ranks /= 2) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

  buffer = packBuffer()
  call buffer%allocateBuffers(1)
  allocate(matrix(-1:510, 2:129))
  do j = 2, 129
    do i = -1, 510
      matrix(i, j) = real(i + 512*j, kind=kind(1.0d0))
    end do
  end do

  if (rank == 0) then
    call buffer%add(37)
    ! Half a MiB fits the requested buffer but exceeds the old byte capacity.
    call buffer%add(matrix)
    call buffer%sendTo(1, 79)
  else
    call buffer%receiveFrom(0, 79)
    call buffer%get(marker)
    call buffer%get(received)
    if (marker /= 37) call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
    if (any(lbound(received) /= lbound(matrix))) &
      call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
    if (any(ubound(received) /= ubound(matrix))) &
      call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
    if (any(received /= matrix)) call MPI_Abort(MPI_COMM_WORLD, 5, ierr)
    print *, 'All tests passed!'
  end if

  call MPI_Finalize(ierr)
end program test_packbuffer_capacity
