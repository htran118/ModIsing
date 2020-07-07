*********************************************************************
*Analysis the data from lattice simulation
*********************************************************************

      program latticeAnal
	use constantClass, only : LDIM
	use crysProcClass, only : dataIsingCube,
     $	                          readChainIsing1d, writeChainIsing1d
	use dataProcClass, only : getLatticeFile, getDataFile
	implicit none

	integer*2, dimension(:), allocatable :: crys(:)
	character(len=64) :: arg
	integer :: lSze, lVol, evalStp, erg, mgnt,
     $	           dataID, crysID, ioS, i
	integer, dimension(:) :: lkErg(1:LDIM)
	real, dimension(:) :: kMgnt(1:LDIM)
	double precision :: temp

*********************************************************************
c	Get command-line argument
*********************************************************************
	i = iargc()
	if(i .ne. 3) then
	  call commandError()
	end if

	call getarg(1, arg)
	read(arg, *) lSze
	if(lSze .lt. 1) then
	  call commandError()
	end if

	call getarg(2, arg)
	read(arg, *) temp
	if(temp .le. 0.0) then
	  call commandError()
	end if

	call getarg(3, arg)
	read(arg, *) evalStp
	if(evalStp .lt. 1) then
	  call commandError()
	end if

*********************************************************************
c	Initialization
*********************************************************************
	crysID = 30
	dataID = 31
	lVol = lSze ** LDIM
	allocate(crys(1:lVol))

*********************************************************************
c	Data analysis
*********************************************************************
	call getLatticeFile(lSze, temp, crysID)
	call getDataFile(lSze, temp, dataID)

	do i = 1, evalStp
	  call readChainIsing1d(crys, crysID)
	  call dataIsingCube(crys, lSze, erg, mgnt, lkErg, kMgnt)
	  write(dataID, *), erg, mgnt, lkErg(1:LDIM), kMgnt(1:LDIM)
	end do

	deallocate(crys)
	close(crysID)
	close(dataID)

      contains

	subroutine commandError()
	  write(*,*) 'Invalid command-line arguments'
	  call exit(1)
	end subroutine

      end program
