*********************************************************************
* module for statistical procedures
*********************************************************************

      module statProcClass
	use randomSeed, only : getRandom
	use constantClass, only : ORRATE, ANISODIR, ANISOERG, 
     $	                          CLOCKDIR, CLOCKERG, TAU,
     $	                          MSKCUT, LDIM
	use crysProcClass, only : getPBCCube, localFieldXYCube,
     $	                          energyChangeXYCube,
     $	                          energyChangeIsingCube
	implicit none

      contains

	subroutine isingCluster(crys, lSze, rLst, rPos, expDif)
	  double precision, dimension(:), intent(inout) :: rLst(:)
	  integer*2, dimension(:), intent(inout) :: crys(:)
	  integer, intent(inout) :: rPos
	  double precision, intent(in) :: expDif
	  integer, intent(in) :: lSze
	  double precision :: rVal
	  integer, dimension(:) :: local(1:(2 * LDIM)),
*	  List of possible position for member of the cluster
     $	                           pLst(1:(2 * LDIM * size(crys)))
	  integer :: pPos, pPtr, pMax, pos,
     $	             seedSpin, countr,
     $	             i, j
	  logical, dimension(:) :: msk(1:size(crys))
	  do i = 1, size(crys)
	    msk(i) = .false. 
	  enddo
	  pMax = size(pLst)
	  pPos = pMax
	  pPtr = pMax
*	Seed for cluster
	  call getRandom(rVal, rPos, rLst)
	  pos = int(rVal * size(crys) + 1)
	  seedSpin = crys(pos)
	  msk(pos) = .true.
	  crys(pos) = -crys(pos)
	  call getPBCCube(pos, lSze, local)
	  do i = 1, (2 * LDIM)
	    pPos = merge(1, (pPos + 1), (pPos .eq. pMax))
	    pLst(pPos) = local(i)
	  enddo
	  countr = 2 * LDIM
c	Building cluster
	  do while(countr .gt. 0)
	    pPtr = merge(1, (pPtr + 1), (pPtr .eq. pMax))
	    countr = countr - 1
	    pos = pLst(pPtr)
	    if(.not. msk(pos)) then
	      if(crys(pos) .eq. seedSpin) then
	        call getRandom(rVal, rPos, rLst)
	        if(rVal .le. expDif) then
	          msk(pos) = .true.
	          crys(pos) = -crys(pos)
	          call getPBCCube(pos, lSze, local)
	          do i = 1, (2 * LDIM)
	            if(.not. msk(local(i))) then
	              pPos = merge(1, (pPos + 1), (pPos .eq. pMax))
	              pLst(pPos) = local(i)
	            endif
	          enddo
	          countr = countr + 2 * LDIM
	        endif 
	      endif
	    endif
	    if((countr .ge. pMax) .or. (countr .lt. 0)) then
	      write(*,*), 'Stack error'
	    endif
	  enddo
	end subroutine

	subroutine isingMetropolis(crys, lSze, rLst, rPos, expDif)
	  double precision, dimension(:), intent(inout) :: rLst(:)
	  integer*2, dimension(:), intent(inout) :: crys(:)
	  integer, intent(inout) :: rPos
	  double precision, dimension(:), intent(in) :: expDif(:)
	  integer, intent(in) :: lSze
	  integer :: i, j, pos, ergDif, stopCond
	  double precision :: rVal
	  do i = 1, size(crys)
	    call getRandom(rVal, rPos, rLst)
	    pos = int(rVal * size(crys) + 1)
	    call energyChangeIsingCube(crys, lSze, pos, ergDif)
	    if(ergDif .le. 0) then
	      crys(pos) = -crys(pos)
	    else
	      stopCond = 0
	      call getRandom(rVal, rPos, rLst)
	      do j = 1, LDIM
	        if(ergDif .eq. (j * 4)) then
	          stopCond = 1
	          if(expDif(j) .ge. rVal) then
	            crys(pos) = -crys(pos)
	          endif
	        endif
	      enddo
	      if(stopCond .eq. 0) then
	        write(*,*), 'Energy error'
	      endif
	    endif
	  enddo
	end subroutine

	subroutine XYOverrelax(crys, lSze, temp, rLst, rPos)
	  double precision, dimension(:), intent(inout) :: crys(:), 
     $	                                                   rLst(:)
	  integer, intent(inout) :: rPos
	  double precision, intent(in) :: temp
	  integer, intent(in) :: lSze
	  double precision :: localF, nwSpin, ergDif, rVal
	  integer :: i, j
	  do i = 1, nint(ORRATE * lSze)
	    do j = 1, size(crys)
	      call localFieldXYCube(crys, lSze, j, localF)
	      if(ANISOERG .eq. 0.0) then
	        crys(j) = 2 * localF - crys(j)
	      else
	        nwSpin = 2 * localF - crys(j)
	        call energyChangeXYCube(crys, lSze, j, nwSpin, ergDif)
	        call getRandom(rVal, rPos, rLst)
	        if(exp(-ergDif / temp) .ge. rVal) then
	          crys(j) = nwSpin
	        endif
	      endif
	    enddo
	  enddo
	end subroutine

	subroutine XYClusterOriginal(crys, lSze, temp, rLst, rPos)
	  double precision, dimension(:), intent(inout) :: crys(:),
     $	                                                   rLst(:)
	  integer, dimension(:) :: mLst(1:(size(crys) / 2)),
     $	                           local(1:(2 * LDIM)),
*	List of possible position for member of the cluster
     $	                           pLst(1:(2 * LDIM * size(crys)))
	  integer, intent(inout) :: rPos
	  double precision, intent(in) :: temp
	  integer, intent(in) :: lSze
	  double precision :: rVal, symPln, dirPln, seedProj, projProd,
     $	                      ergDif
	  integer :: pPos, pPtr, pMax, pos, mPos, mPtr, mMax,
     $	             countr, i, j
	  logical, dimension(:) :: msk(1:size(crys))
	  do i = 1, size(msk)
	    msk(i) = .false.
	  enddo
	  pMax = size(pLst) 
	  pPos = pMax
	  pPtr = pMax
	  mMax = size(mLst)
	  mPos = mMax
	  mPtr = mMax
	  ergDif = 0.0
*	Seed for reflection direction
	  call getRandom(rVal, rPos, rLst)
	  symPln = rVal * TAU
	  dirPln = symPln + TAU / 4
*	Seed for cluster
	  call getRandom(rVal, rPos, rLst)
	  pos = int(rVal * size(crys) + 1)
	  msk(pos) = .true.
	  mPos = merge(1, (mPos + 1), (mPos .eq. mMax))
	  mLst(mPos) = pos
	  call getPBCCube(pos, lSze, local)
	  do i = 1, (2 * LDIM)
	    pPos = merge(1, (pPos + 1), (pPos .eq. pMax))
	    pLst(pPos) = local(i)
	  enddo
	  countr = 2 * LDIM
*	Building cluster
	  do while(countr .gt. 0)
	    if(mod(pPtr, (2 * LDIM)) .eq. 0) then
	       mPtr = merge(1, (mPtr + 1), (mPtr .eq. mMax))
	       seedProj = dcos(crys(mLst(mPtr)) - dirPln)
	    endif
	    pPtr = merge(1, (pPtr + 1), (pPtr .eq. pMax))
	    pos = pLst(pPtr)
	    countr = countr - 1
	    if(.not. msk(pos)) then
	      projProd = dcos(crys(pos) - dirPln) * seedProj
	      if(projProd .gt. 0.0) then
	        call getRandom(rVal, rPos, rLst)
	        if((1 - exp(-2 * projProd / temp)) .ge. rVal) then
	          msk(pos) = .true.
	          mPos = merge(1, (mPos + 1), (mPos .eq. mMax))
	          mLst(mPos) = pos
	          call getPBCCube(pos, lSze, local)
	          do i = 1, (2 * LDIM)
	            pPos = merge(1, (pPos + 1), (pPos .eq. pMax))
	            pLst(pPos) = local(i)
	          enddo
	          countr = countr + 2 * LDIM 
	        endif 
	      endif
	    endif
	    if((countr .ge. pMax) .or. (countr .lt. 0)) then
	      write(*,*), 'Stack error'
	    endif
	  enddo
*	Flipping cluster
	  if(ANISOERG .ne. 0.0) then
	    do i = 1, size(crys)
	      if(msk(i)) then
	        ergDif = ergDif + (dcos(ANISODIR * crys(i)) -
     $	                           dcos(ANISODIR * (2 * symPln -
     $	                                            crys(i)))) *
     $	                           ANISOERG 
	      endif
	    enddo
	  endif
	  if(ergDif .le. 0.0) then
	    do i = 1, size(crys)
	      if(msk(i)) then
	        crys(i) = 2 * symPln - crys(i)
	      endif
	    enddo
	  else
	    call getRandom(rVal, rPos, rLst)
	    if(exp(-ergDif / temp) .ge. rVal) then
	      do i = 1, size(crys)
	        if(msk(i)) then
	          crys(i) = 2 * symPln - crys(i)
	        endif
	      enddo
	    endif
	  endif
	end subroutine

	subroutine XYClusterModified(crys, lSze, temp, rLst, rPos)
	  double precision, dimension(:), intent(inout) :: crys(:),
     $	                                                   rLst(:)
	  integer, dimension(:) :: mLst(1:(size(crys) / 2)),
     $	                           local(1:(2 * LDIM)),
*	List of possible position for member of the cluster
     $	                           pLst(1:(2 * LDIM * size(crys)))
	  integer, intent(inout) :: rPos
	  double precision, intent(in) :: temp
	  integer, intent(in) :: lSze
	  double precision :: rVal, symPln, dirPln, olSpin, nwSpin,
     $	                      pSpn, aniDif, ergDif
	  integer :: pPos, pPtr, pMax, pos, mPos, mPtr, mMax,
     $	             countr, i, j
	  logical, dimension(:) :: msk(1:size(crys))
	  do i = 1, size(msk)
	    msk(i) = .false.
	  enddo
	  pMax = size(pLst) 
	  pPos = pMax
	  pPtr = pMax
	  mMax = size(mLst)
	  mPos = mMax
	  mPtr = mMax
	  ergDif = 0.0
*	Seed for reflection direction
	  call getRandom(rVal, rPos, rLst)
	  symPln = rVal * TAU
*	Seed for cluster
	  call getRandom(rVal, rPos, rLst)
	  pos = int(rVal * size(crys) + 1)
	  msk(pos) = .true.
	  mPos = merge(1, (mPos + 1), (mPos .eq. mMax))
	  mLst(mPos) = pos
	  call getPBCCube(pos, lSze, local)
	  do i = 1, (2 * LDIM)
	    pPos = merge(1, (pPos + 1), (pPos .eq. pMax))
	    pLst(pPos) = local(i)
	  enddo
	  countr = 2 * LDIM
*	Building cluster
	  do while(countr .gt. 0)
	    if(mod(pPtr, (2 * LDIM)) .eq. 0) then
	       mPtr = merge(1, (mPtr + 1), (mPtr .eq. mMax))
	       olSpin = crys(mLst(mPtr))
	       nwSpin = 2 * symPln - olSpin
	    endif
	    pPtr = merge(1, (pPtr + 1), (pPtr .eq. pMax))
	    pos = pLst(pPtr)
	    countr = countr - 1
	    if(.not. msk(pos)) then
	      pSpn = crys(pos)
	      ergDif = dcos(olSpin - pSpn) - dcos(nwSpin - pSpn)
	      if(CLOCKERG .ne. 0.0) then
	      ergDif = CLOCKERG * (dcos(CLOCKDIR * (olSpin - pSpn)) -
     $	                           dcos(CLOCKDIR * (nwSpin - pSpn))) +
     $	               ergDif
	      endif
	      if(ergDif .gt. 0.0) then
	        call getRandom(rVal, rPos, rLst)
	        if((1 - exp(-ergDif / temp)) .ge. rVal) then
	          msk(pos) = .true.
	          mPos = merge(1, (mPos + 1), (mPos .eq. mMax))
	          mLst(mPos) = pos
	          call getPBCCube(pos, lSze, local)
	          do i = 1, (2 * LDIM)
	            pPos = merge(1, (pPos + 1), (pPos .eq. pMax))
	            pLst(pPos) = local(i)
	          enddo
	          countr = countr + 2 * LDIM 
	        endif 
	      endif
	    endif
	    if((countr .ge. pMax) .or. (countr .lt. 0)) then
	      write(*,*), 'Stack error'
	    endif
	  enddo
*       Flipping cluster
	  if(ANISOERG .ne. 0.0) then
	    aniDif = 0.0
	    do i = 1, size(crys)
	      if(msk(i)) then
	        aniDif = ANISOERG * (dcos(ANISODIR * crys(i)) -
     $	                             dcos(ANISODIR * (2 * symPln -
     $	                                              crys(i)))) + 
     $	                 aniDif
	      endif
	    enddo
	  endif
	  if(aniDif .le. 0.0) then
	    do i = 1, size(crys)
	      if(msk(i)) then
	        crys(i) = 2 * symPln - crys(i)
	      endif
	    enddo
	  else
	    call getRandom(rVal, rPos, rLst)
	    if(exp(-aniDif / temp) .ge. rVal) then
	      do i = 1, size(crys)
	        if(msk(i)) then
	          crys(i) = 2 * symPln - crys(i)
	        endif
	      enddo
	    endif
	  endif
	end subroutine

	subroutine XYMetropolis(crys, lSze, temp, rLst, rPos)
	  double precision, dimension(:), intent(inout) :: crys(:),
     $	                                                   rLst(:)
	  integer, intent(inout) :: rPos
	  double precision, intent(in) :: temp
	  integer, intent(in) :: lSze
	  integer :: i, pos
	  double precision :: nwSpin, ergDif, rVal
	  do i = 1, size(crys)
	    call getRandom(rVal, rPos, rLst)
	    pos = int(rVal * size(crys) + 1)
	    call getRandom(rVal, rPos, rLst)
	    nwSpin = rVal * TAU
	    call energyChangeXYCube(crys, lSze, pos, nwSpin, ergDif)
	    if(ergDif .le. 0.0) then
	      crys(pos) = nwSpin
	    else
	      call getRandom(rVal, rPos, rLst)
	      if(exp(-ergDif / temp) .ge. rVal) then
	        crys(pos) = nwSpin
	      endif
	    endif
	  enddo
	end subroutine

	subroutine XYHeatBath(crys, lSze, rLst, rPos)
	  double precision, dimension(:), intent(inout) :: crys(:),
     $	                                                   rLst(:)
	  integer, intent(inout) :: rPos
	  integer, intent(in) :: lSze
	end subroutine

	subroutine average(inp, avg)
	  real, dimension(:), intent(in) :: inp(:)
	  real, intent(out) :: avg
	  avg = sum(inp) / size(inp)
	end subroutine

	subroutine variance(inp, var)
	  real :: avg
	  real, dimension(:), intent(in) :: inp(:)
	  real, intent(out) :: var
	  avg = sum(inp) / size(inp)
	  var = sum((inp - avg) ** 2) / (size(inp) - 1)
	end subroutine

	subroutine rawMoment(inp, rwMmt, ord)
	  integer, intent(in) :: ord
	  real, dimension(:), intent(in) :: inp(:)
	  real, intent(out) :: rwMmt
	  rwMmt = sum(inp ** ord) / size(inp)
	end subroutine 

      end module 
