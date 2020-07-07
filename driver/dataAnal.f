*********************************************************************
* statistical analysis of the raw data
*********************************************************************

      program dataAnal
	use constantClass, only : LDIM, TAU
	use dataProcClass, only : getDataFile, getStatisticsFile 
	implicit none

	integer :: lSze, lVol, evalStp, dataID, statID, ioSt, i
	double precision :: temp
	real :: ergAvg, ergVar, ergDens, heatDens, mgntAvg, mgntVar, 
     $	        mgntDens, suscept, mgntSqr, mgntQud, binder,
     $	        helix, helix4th, corrRate
	real, dimension(:) :: lkErgAvg(1:LDIM), lkErgSqr(1:LDIM),
     $	                      h(1:LDIM), h4(1:LDIM), corr(1:LDIM)
	real, dimension(:), allocatable :: erg(:), mgnt(:)
	real, dimension(:,:), allocatable :: lkErg(:,:), kMgnt(:,:)
	character(len=64) :: arg

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
	dataID = 31
	statID = 30
	lVol = lSze ** 2

	allocate(erg(1:evalStp))
	allocate(mgnt(1:evalStp))
	allocate(lkErg(1:evalStp,1:LDIM))
	allocate(kMgnt(1:evalStp,1:LDIM))
*********************************************************************
c	Get data from file
*********************************************************************
	call getDataFile(lSze, temp, dataID)
	call getStatisticsFile(statID)

	do i = 1, evalStp
	  read(dataID, *), erg(i), mgnt(i), lkErg(i,1:LDIM),
     $	                   kMgnt(i,1:LDIM) 
	end do

	ergAvg = real(sum(erg)) / evalStp
	ergVar = sum((erg - ergAvg) ** 2) / evalStp
	heatDens = ergVar / ((temp ** 2) * lVol)
	ergDens = ergAvg / lVol

	mgntAvg = real(sum(abs(mgnt))) / evalStp
	mgntVar = sum((abs(mgnt) - mgntAvg) ** 2) / evalStp 
	suscept = mgntVar / (temp * lVol)
	mgntDens = mgntAvg / lVol

c	mgntAbs = sum(abs(mgnt(1:evalStp))) / evalStp

	mgntSqr = real(sum(mgnt ** 2)) / evalStp
	mgntQud = real(sum(mgnt ** 4)) / evalStp
	binder = 1 - mgntQud / (3 * (mgntSqr ** 2))

	lkErgAvg = sum(lkErg, dim=1) / evalStp
	lkErgSqr = sum(lkErg ** 2, dim=1) / evalStp
	corr = sum(kMgnt, dim=1) / evalStp

	do i = 1, LDIM
	  h(i) = lkErgAvg(i) / lVol
          h4(i) = (3 * ((lkErgAvg(i) ** 2) - lkErgSqr(i)) / temp +
     $	           lkErgAvg(i)) / (lVol ** 2)
	  corr(i) = sqrt((mgntSqr / corr(i)) - 1) / 
     $	            (2 * sin(TAU / (2 * lSze)) * lSze)
	enddo
	helix = sum(h) / LDIM
	helix4th = sum(h4) / LDIM
	corrRate = sum(corr) / LDIM
	write(statID, *), lSze, temp, evalStp, ergDens, heatDens,
     $	                  mgntDens, suscept, binder, helix, helix4th,
     $	                  corrRate

	deallocate(erg)
	deallocate(mgnt)
	close(dataID)
	close(statID)

      contains
c	subroutine listFile(dir, lsFl)
c	  character(len=64) :: dir, lsFl
c	  lsFl = trim(dir) // 'ls.txt'
c	  call system('ls ' // dir // ' > ' // lsFl)
c	end subroutine
	subroutine commandError()
	  write(*,*) 'Invalid command-line arguments'
	  write(*,*) 'Usage:'
	  write(*,*) './(programName) (latticeSize) (temperature)',
     $	             '(stepNumber)'
	  write(*,*) '(latticeSize) must be even'
	  write(*,*) '(temperature) must be positive'
	  write(*,*) '(stepNumber) must be positive'
	  call exit(1)
	end subroutine
      end program

