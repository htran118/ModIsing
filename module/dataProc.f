*********************************************************************
* module for data-storing procedures
*********************************************************************

      module dataProcClass
	use constantClass, only : SUBDATNME, 
     $	                          RESPATH, RAWPATH, STATFILE, CORRFILE,
     $	                          WRITEFORM
	implicit none

      contains

	subroutine getLatticeFile(lSze, temp, flID)
	  character(len=64) :: flNm
	  integer, intent(in) :: flID, lSze
	  double precision, intent(in) :: temp
	  write(flNm, '(a, a, i0.3, a, f5.3, a)'), RAWPATH, '/sze',
     $	        lSze, '/lattice', temp, '.trj'
	  open(flID, file=flNm, form=WRITEFORM)
	end subroutine

	subroutine getDataFile(lSze, temp, flID)
	  character(len=64) :: flNm
	  integer, intent(in) :: flID, lSze
	  double precision, intent(in) :: temp
	  write(flNm, '(a, a, i0.3, a, f5.3, a)'), RESPATH, '/sze', 
     $	        lSze, '/data', temp, '.dat'
	  open(flID, file=flNm, form=WRITEFORM)
	end subroutine 

	subroutine getSubDataFile(lSze, temp, flID, ctrl)
	  character(len=64) :: flNm
	  integer, intent(in) :: lSze
	  integer, dimension(:), intent(in) :: flID
	  double precision, intent(in) :: temp
	  logical, dimension(:), intent(in) :: ctrl
	  integer :: i
	  do i = 1, size(flID)
	    if(ctrl(i)) then 
	      write(flNm, '(a, a, i0.3, a, a, f5.3, a)'),
     $	            RESPATH, '/sze', lSze, '/',
     $	            SUBDATNME(i), temp, '.dat'
	      open(flID(i), file=flNm, form=WRITEFORM)
	    endif
	  enddo
	end subroutine

	subroutine checkSubDataFile(lSze, temp, flID, ctrl, evalStp)
	  character(len=1) :: rdInp
	  integer, intent(in) :: lSze, evalStp
	  integer, dimension(:), intent(in) :: flID
	  double precision, intent(in) :: temp
	  logical, dimension(:), intent(inout) :: ctrl
	  integer :: flSze, i, n, ioSt
	  logical :: exst
	  do i = 1, size(flID)
	    if(ctrl(i)) then
	      ioSt = 0
	      n = 0
	      do while(ioSt .eq. 0)
	        read(flID(i)), rdInp
	        n = n + 1
	      enddo
	      if(n .eq. (evalStp + 1)) then
	        ctrl(i) = .false.
	        close(flID(i))
	      else
	        rewind(flID(i))
	      endif
	    endif
	  enddo
	end subroutine
c	subroutine getNewDataFile(lSze, temp, flID)
c	  character(len=64) :: flNm
c	  integer, intent(in) :: flID, lSze
c	  double precision, intent(in) :: temp
c	  write(flNm, '(a, a, i0.3, a, f5.3, a)'), RESPATH, '/sze', 
c     $	        lSze, '/data', temp, '.dat.new'
c	  open(flID, file=flNm, form=WRITEFORM)
c	end subroutine

c	subroutine clearDataFile(lSze, temp)
c	  character(len=64) :: oldNm, newNm
c	  character(len=256) :: command
c	  integer, intent(in) :: lSze
c	  double precision, intent(in) :: temp
c	  write(oldNm, '(a, a, i0.3, a, f5.3, a)'), RESPATH, '/sze', 
c     $	        lSze, '/data', temp, '.dat'
c	  write(newNm, '(a, a, i0.3, a, f5.3, a)'), RESPATH, '/sze', 
c     $	        lSze, '/data', temp, '.new'
c	  write(command, '(a,a,a,a)', 'mv ', newNm, ' ', oldNm 
c	  write(*,*), command
c	  call system(command)
c	  write(command, '(a,a)'), 'rm ', newNm
c	  write(*,*), command
c	  call system(command)
c	end subroutine

c	subroutine getExtDataFile(lSze, temp, flID)
c	  character(len=64) :: flNm
c	  integer, intent(in) :: flID, lSze
c	  double precision, intent(in) :: temp
c	  write(flNm, '(a, a, i0.3, a, f5.3, a)'), RESPATH, '/sze', 
c     $	        lSze, '/data', temp, '.ext'
c	  open(flID, file=flNm, form=WRITEFORM)
c	end subroutine

c	subroutine mergeDataFile(lSze, temp)
c	  character(len=64) :: oldNm, newNm
c	  character(len=256) :: command
c	  integer, intent(in) :: lSze
c	  double precision, intent(in) :: temp
c	  write(oldNm, '(a, a, i0.3, a, f5.3, a)'), RESPATH, '/sze', 
c     $	        lSze, '/data', temp, '.dat'
c	  write(newNm, '(a, a, i0.3, a, f5.3, a)'), RESPATH, '/sze', 
c     $	        lSze, '/data', temp, '.dat.new'
c	  write(command, 'a,a,a,a,a,a'), 'paste ', oldNm, ' ', newNm,
c     $	                                 ' > ', oldNm
c	  write(*,*), command
c	  call system(command)
c	  write(command, 'a,a'), 'rm ', newNm
c	  write(*,*), command
c	end subroutine

	subroutine writeStatisticsFile(flID)
	  character(len=64) :: flNm
	  integer, intent(in) :: flID
	  logical :: old
	  write(flNM, '(a, a)'), STATFILE, '.txt'
	  inquire(file=flNm, exist=old)
	  if(old) then
	    open(flID, file=flNm, status='old', position='append')
	  else
	    open(flID, file=flNm, status='new')
	  endif
	end subroutine

	subroutine writeBinaryStatisticsFile(flID)
	  character(len=64) :: flNm
	  integer, intent(in) :: flID
	  logical :: old
	  write(flNM, '(a, a)'), STATFILE, '.bin'
	  inquire(file=flNm, exist=old)
	  if(old) then
	    open(flId, file=flNm, status='old', position='append',
     $	         form=WRITEFORM)
	  else
	    open(flID, file=flNm, status='new', form=WRITEFORM)
	  endif
	end subroutine

	subroutine readBinaryStatisticsFile(flID)
	  character(len=64) :: flNm
	  integer, intent(in) :: flID
	  logical :: old
	  write(flNM, '(a, a)'), STATFILE, '.bin'
	  open(flId, file=flNm, form=WRITEFORM)
	end subroutine

	subroutine writeCorrelationFile(flID)
	  character(len=64) :: flNm
	  integer, intent(in) :: flID
	  logical :: old
	  write(flNM, '(a, a)'), CORRFILE, '.txt'
	  inquire(file=flNm, exist=old)
	  if(old) then
	    open(flID, file=flNm, status='old', position='append')
	  else
	    open(flID, file=flNm, status='new')
	  endif
	end subroutine

c	subroutine getEnergyFile(lSze, temp, flID)
c	  character(len=64) :: fl
c	  integer, intent(in) :: flID, lSze
c	  double precision, intent(in) :: temp
c	  write(fl, '(a, a, i0.3, a, f4.2, a)'), RESPATH, '/sze',
c    $	        lSze, '/energy', temp, '.txt'
c	  open(flID, file=fl, form=WRITEFORM)
c	end subroutine

c	subroutine getMagnetFile(lSze, temp, flID)
c	  character(len=64) :: flNm
c	  integer, intent(in) :: flID, lSze
c	  double precision, intent(in) :: temp
c	  write(flNm, '(a, a, i0.3, a, f4.2, a)'), RESPATH, '/sze',
c     $	        lSze, '/magnet', temp, '.txt'
c	  open(flID, file=flNm, form=WRITEFORM)
c	end subroutine

c	subroutine getLinkEnergyFile(lSze, temp, flID)
c	  character(len=64) :: flNm
c	  character :: lkNm
c	  integer, intent(in) :: flID, lSze
c	  double precision, intent(in) :: temp
c	  if(LKDIR .eq. 1) then
c	    lkNm = 'X'
c	  else if(LKDIR .eq. 2) then
c	    lkNm = 'Y'
c	  else if(LKDIR .eq. 3) then
c	    lkNm = 'Z'
c	  endif
c	  write(flNm, '(a, a, i0.3, a, a, a, f4.2, a)'), RESPATH,
c     $	        '/sze', lSze, '/link', lkNm, 'Energy', temp, '.txt'
c	  open(flID, file=flNm, form=WRITEFORM)
c	end subroutine

c	subroutine getLinkCurrentFile(lSze, temp, flID)
c	  character(len=64) :: flNm
c	  character :: lkNm
c	  integer, intent(in) :: flID, lSze
c	  double precision, intent(in) :: temp
c	  if(LKDIR .eq. 1) then
c	    lkNm = 'X'
c	  else if(LKDIR .eq. 2) then
c	    lkNm = 'Y'
c	  else if(LKDIR .eq. 3) then
c	    lkNm = 'Z'
c	  endif
c	  write(flNm, '(a, a, i0.3, a, a, a, f4.2, a)'), RESPATH,
c     $	        '/sze', lSze, '/link', lkNm, 'Current', temp, '.txt'
c	  open(flID, file=flNm, form=WRITEFORM)
c	end subroutine

	subroutine countLineNumber(flID, lnNum)
	  integer, intent(in) :: flID
	  integer, intent(out) :: lnNum
	  integer :: ioSt
	  character(len=1) :: rdInp
	  lnNum = 0
	  ioSt = 0
	  do while(ioSt .eq. 0)
	    read(flID, *, ioStat=ioSt) rdInp
	    lnNum = lnNum + 1
	  end do
	  rewind(flID)
	end subroutine
      end module

