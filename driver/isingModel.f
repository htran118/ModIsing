*********************************************************************
*Simulation of 2D Ising model using MC method
*********************************************************************

      program xy2dModel
	use randomSeed, only : getRandom
	use constantClass, only : LDIM, RANLSTSZE, CHKSTP, AVGEQSTP
	use crysProcClass, only : writeChainIsing1d
	use statProcClass, only : isingCluster
	use dataProcClass, only : getLatticeFile 
	implicit none

	integer*2, dimension(:), allocatable :: crys(:)
	double precision, dimension(:), allocatable :: rLst(:)
	double precision, dimension(:) :: expDif(0 : LDIM)
	integer, dimension(:) :: local(1 : (2 * LDIM))
	character(len=64) :: arg
	integer :: lSze, lVol,
     $	           stpNum, eqStp, calStp, 
     $	           pos, ergDif,
     $	           crysID,
     $	           rPos,
     $	           i, j, k,
     $	           success
	double precision :: temp, rVal

*********************************************************************
c	Get command-line argument
*********************************************************************
	i = iargc()
	if(i .ne. 4) then
	  call commandError()
	endif

	call getarg(1, arg)
	read(arg, *) lSze
	if(lSze .lt. 1) then
	  call commandError()
	endif

	call getarg(2, arg)
	read(arg, *) temp
	if(temp .le. 0.0) then
	  call commandError()
	endif

	call getarg(3, arg)
	read(arg, *) stpNum
	call getarg(4, arg)
	read(arg, *) calStp
	if((stpNum .lt. 1) .or. (calStp .gt. stpNum)) then
	  call commandError()
	endif

*********************************************************************
c	initialization
*********************************************************************
	crysID = 20
	lVol = lSze ** LDIM
	eqStp = stpNum / 2
	if(eqStp .lt. (AVGEQSTP * lVol)) then
	  eqStp = AVGEQSTP * lVol
	endif
	stpNum = stpNum + eqStp

	write(*,*), lSze, temp, stpNum, eqStp

	expDif(0) = 1 - exp(-2.0 / temp)
	do i = 1, LDIM
	  expDif(i) = exp(-4.0 * i / temp)
	enddo

	allocate(rLst(1 : lVol))
	rPos = lVol

	success = 0

*********************************************************************
c	Create lattice
*********************************************************************
	allocate(crys(1:lVol))

*********************************************************************
c	Initialize spin
*********************************************************************
c	Randomize spin
	do i = 1, lVol
	  call getRandom(rVal, rPos, rLst)
	  if(rVal .ge. 0.5) then
	    crys(i) = 1
	  else
	    crys(i) = -1
	  endif
	enddo

	deallocate(rLst)
	allocate(rLst(1 : RANLSTSZE))
	rPos = RANLSTSZE

*********************************************************************
c	Energy evaluation
*********************************************************************
c	call energyIsing2d(l, erg)

*********************************************************************
c	MC simulation
*********************************************************************
	call getLatticeFile(lSze, temp, crysID)

	do i = 1, stpNum
	  call isingCluster(crys, lSze, rLst, rPos, expDif(0))
c	  call isingMetropolis(crys, lSze, rLst, rPos, expDif(1:LDIM))
	  if(i .gt. eqStp) then
	    if(mod(i, calStp) .eq. 0) then
	      call writeChainIsing1d(crys, crysID)
	    endif
	  endif
	enddo 

	write(*,*), (real(success) / (stpNum * lVol)), lSze, lVol, temp

	close(crysID)
	deallocate(crys)
	deallocate(rLst)

      contains

	subroutine commandError()
	  write(*,*) 'Invalid command-line arguments'
	  write(*,*) 'Usage:'
	  write(*,*) './(programName) (latticeSize) (temperature) ',
     $	             '(stepNumber) (stepPerMeasure)'
	  write(*,*) '(latticeSize) must be even'
	  write(*,*) '(temperature) must be positive'
	  write(*,*) '(stepNumber) must be positive'
	  call exit(1)
	end subroutine
      end program
