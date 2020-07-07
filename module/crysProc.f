*********************************************************************
* module for lattice procedures 
*********************************************************************

      module crysProcClass
	use constantClass, only : LDIM, SDIM, 
     $	                          SUBDATNME, SUBDATSZE, 
     $	                          XYFORM, ISINGFORM,
     $	                          TAU, RTAU,
     $	                          SPACE, ANISOERG, ANISODIR,
     $	                          POSITIVE, NEGATIVE
	implicit none

      contains

*********************************************************************
*	general procedure using LDIM
*********************************************************************
	subroutine getLinkPBCCube(pos, lSze, lnkLst)
	  integer, intent(in) :: pos, lSze
	  integer, dimension(:), intent(out) :: lnkLst(1:(2 * LDIM))
	  integer :: i, lVol, prvLnk, offset
	  lVol = lSze ** LDIM
	  offset = 1
	  do i = 1, LDIM
	    lnkLst(2 * i) = pos + (i - 1) * lVol
	    prvLnk = pos - offset
	    offset = offset * lSze
	    if(prvLnk .le. ((pos - 1) / offset * offset)) then
	      prvLnk = prvLnk + offset
	    endif
	    lnkLst((2 * i) - 1) = prvLnk + (i - 1) * lVol
	  enddo
	end subroutine

	subroutine getPBCCube(pos, lSze, local)
	  integer, intent(in) :: pos, lSze
	  integer, dimension(:) :: lVol(0:LDIM)
	  integer :: i, offset
	  integer, dimension(:), intent(out) :: local(1:(2 * LDIM))
	  offset = 1
	  do i = 1, LDIM
	    local((2 * i) - 1) = pos - offset
	    local(2 * i) = pos + offset
	    offset = offset * lSze
	    if(local((2 * i) - 1) .le. ((pos - 1) / offset * offset))
     $	    then
	      local((2 * i) - 1) = local((2 * i) - 1) + offset
	    endif
	    if(local(2 * i) .gt. (((pos - 1) / offset + 1) * offset)) 
     $	    then
	      local(2 * i) = local(2 * i) - offset
	    endif
	  enddo
	end subroutine

	subroutine getCoordinateCube(pos, lSze, coord)
	  integer, intent(in) :: pos, lSze
	  integer :: i, dummy
	  integer, dimension(:), intent(out) :: coord(0:LDIM)
	  dummy = pos
	  coord(0) = 0
	  do i = 1, LDIM
	    dummy = dummy - coord(i - 1) * (lSze ** (LDIM - i + 1))
	    coord(i) = dummy / (lSze ** (LDIM - i))
	  enddo
	end subroutine

	subroutine energyChangeXYCube(crys, lSze, pos, nwSpin, ergDif)
	  double precision, dimension(:), intent(in) :: crys
	  double precision, intent(in) :: nwSpin
	  integer, intent(in) :: lSze, pos
	  integer, dimension(:) :: local(1:(2 * LDIM))
	  integer :: i
	  double precision, intent(out) :: ergDif
	  ergDif = merge((ANISOERG * (dcos(ANISODIR * crys(pos)) -
     $	                              dcos(ANISODIR * nwSpin))),
     $	                 0.0D0, (ANISOERG .ne. 0.0))
	  call getPBCCube(pos, lSze, local)
	  do i = 1, (2 * LDIM)
	    ergDif = ergDif + dcos(crys(local(i)) - crys(pos)) - 
     $	                      dcos(crys(local(i)) - nwSpin)
	  enddo
	end subroutine

	subroutine localFieldXYCube(crys, lSze, pos, localF)
	  double precision, dimension(:), intent(in) :: crys
	  integer, intent(in) :: lSze, pos
	  double precision :: x, y
	  integer, dimension(:) :: local(1:(2 * LDIM))
	  integer :: i
	  double precision, intent(out) :: localF
	  x = 0.0
	  y = 0.0
	  call getPBCCube(pos, lSze, local)
	  do i = 1, (2 * LDIM)
	    x = x + dcos(crys(local(i)))
	    y = y + dsin(crys(local(i)))
	  enddo
	  localF = datan2(y, x)
	end subroutine

*********************************************************************
* Control is an integer, whose binary representation which data get
* analysed:
* 1 -- PRIME
*      Standard Package: Hamiltonian, magnetization, k-magnetization
* 2 -- HELIX
*      Helicity modulus: link current, link energy
* 3 -- CMGNT
*      Clock magnetization and its associated k-magnetization
*      (with q = ANISODIR)
* 4 -- ACMGN
*      Absolute clock magnetization and its associated
*      k-magnetization (with q = ANISODIR)
* 5 -- NMGNT
*      Nematic magnetization and its associated k-magnetization
*      (with q = ANISODIR)
* 6 -- 2MGNT
*      Half-vortex magnetization and its associated k-magnetization
*      (with q = 2)
* 7 -- A2MGN
* 8 -- N2MGN
*      Nematic half-vortex magnetization and its associated
*      k-magnetization (with q = 2)
* 9 -- DISTR
*      Distribution of spin direction
*********************************************************************

      subroutine dataXYCube(crys, lSze, ctrl, outData, distSze)
	  real, dimension(:), intent(in) :: crys(:)
	  integer, intent(in) :: lSze, distSze
	  logical, dimension(:), intent(in) :: ctrl(:)
	  real, dimension(:), intent(out) :: outData(:)
	  integer, dimension(:) :: coord(0:LDIM)
	  integer :: i, j, k, s, pos, adjPos, ptr
*	  hMgn is half-vortex magnetization. Fortran syntax does not
*	  allow the use of 2Mgn
	  real :: dif, mgnt, cMgn, nMgn, hMgn, n2Mg, erg, 
     $	          angle
	  real, dimension(:) :: kMgntFT(1:(LDIM * SDIM * 2)),
     $	                        kNMgnFT(1:(LDIM * SDIM * 2)),
     $	                        kN2MgFT(1:(LDIM * SDIM * 2)),
     $	                        trig(1:(lSze * 2)),
     $	                        spin(1:SDIM), spinTot(1:SDIM),
     $	                        nSpn(1:SDIM), nSpnTot(1:SDIM),
     $	                        n2Sp(1:SDIM), n2SpTot(1:SDIM),
     $	                        lkErg(1:LDIM), lkCur(1:LDIM),
     $	                        kMgnt(1:LDIM), kNMgn(1:LDIM),
     $	                        kN2Mg(1:LDIM), dist(1:distSze)
	  if(SDIM .ne. 2) then
	    write(*,*) 'Spin dimension error'
	    call exit(1)
	  endif
*	  initialization
	  coord(0) = 0
	  do s = 1, SDIM
	    spinTot(s) = 0.0
	    nSpnTot(s) = 0.0
	    n2SpTot(s) = 0.0
	  enddo
	  cMgn = 0.0
	  hMgn = 0.0
	  do i = 1, (LDIM * SDIM * 2)
	    kMgntFT(i) = 0.0
	    kNMgnFT(i) = 0.0
	    kN2MgFT(i) = 0.0
	  enddo
	  do i = 1, lSze
	    trig(2 * i - 1) = cos((i - 1) * TAU / lSze)
	    trig(2 * i) = sin((i - 1) * TAU / lSze)
	  enddo
	  do i = 1, LDIM
	    lkErg(i) = 0.0
	    lkCur(i) = 0.0
	  enddo
	  do i = 1, distSze
	    dist(i) = 0
	  enddo
c	  begin data analysis
	  do i = 1, size(crys)
	    if(ctrl(3)) then
	      cMgn = cMgn + cos(ANISODIR * crys(i))
	    endif
	    if(ctrl(6)) then
	      hMgn = hMgn + cos(2 * crys(i))
	    endif
	    do s = 1, SDIM
	      if(ctrl(1)) then
	        spin(s) = merge(cos(crys(i)), sin(crys(i)), 
     $	                        (s .eq. 1))
	        spinTot(s) = spinTot(s) + spin(s)
	      endif
	      if(ctrl(5)) then
	        nSpn(s) = merge(cos(ANISODIR * crys(i)),
     $	                        sin(ANISODIR * crys(i)),
     $	                        (s .eq. 1))
	        nSpnTot(s) = nSpnTot(s) + nSpn(s)
	      endif
	      if(ctrl(8)) then
	        n2Sp(s) = merge(cos(2 * crys(i)), sin(2 * crys(i)),
     $	                         (s .eq. 1))
	        n2SpTot(s) = n2SpTot(s) + n2Sp(s)
	      endif
	    enddo
	    call getCoordinateCube((i - 1), lSze, coord)
*	  component of the k-dependent magnetization
	    do j = 1, LDIM
	      do s = 1, SDIM
	        do k = 1, 2
	          if(ctrl(1)) then
	            kMgntFT(SDIM * (2 * j + s - 3) + k) =
     $	            kMgntFT(SDIM * (2 * j + s - 3) + k) +
     $	            spin(s) * trig(2 * coord(LDIM - j + 1) + k)
	          endif
	          if(ctrl(5)) then
	            kNMgnFT(SDIM * (2 * j + s - 3) + k) =
     $	            kNMgnFT(SDIM * (2 * j + s - 3) + k) +
     $	            nSpn(s) * trig(2 * coord(LDIM - j + 1) + k)
	          endif
	          if(ctrl(8)) then
	            kN2MgFT(SDIM * (2 * j + s - 3) + k) =
     $	            kN2MgFT(SDIM * (2 * j + s - 3) + k) +
     $	            n2Sp(s) * trig(2 * coord(LDIM - j + 1) + k)
	          endif
	        enddo
	      enddo
	      if(ctrl(1)) then
	        kMgnt(j) = sum(kMgntFT((SDIM * 2 * (j - 1) + 1):
     $	                             (SDIM * 2 * j)) ** 2)
	      endif
	      if(ctrl(5)) then
	        kNMgn(j) = sum(kNMgnFT((SDIM * 2 * (j - 1) + 1):
     $	                               (SDIM * 2 * j)) ** 2)
	      endif
	      if(ctrl(8)) then
	        kN2Mg(j) = sum(kN2MgFT((SDIM * 2 * (j - 1) + 1):
     $	                                 (SDIM * 2 * j)) ** 2)
	      endif
	    enddo
*	  component of the helicity modulus
*	  the sum of lkErg is the Hamiltonian due to the spin-spin
*	  interaction. So it's a good idea to compute it in both cases
	    if(ctrl(2) .or. ctrl(1)) then
	      do j = 1, LDIM
	        pos = coord(j)
	        coord(j) = mod((coord(j) + 1), lSze)
	        adjPos = 1
	        do k = 1, LDIM
	          adjPos = adjPos + coord(k) * (lSze ** (LDIM - k))
	        enddo
	        dif = crys(adjPos) - crys(i)
	        lkErg(j) = lkErg(j) + cos(dif)
	        if(ctrl(2)) then
	          lkCur(j) = lkCur(j) + sin(dif)
	        endif
	        coord(j) = pos 
	      enddo
	    endif
*	  spin distribution
*	  this one sucks. may need better approach
	    if(ctrl(9)) then
	      angle = modulo(crys(i), RTAU)
	      j = mod(ceiling(angle / TAU * distSze), distSze) + 1
	      dist(j) = dist(j) + 1
	    endif
	  enddo
c	  end data analysis
	  if(ctrl(1)) then
	    erg = sum(lkErg(1:LDIM))
	    if(ANISOERG .ne. 0.0) then
	      do i = 1, size(crys)
	        erg = erg + ANISOERG * cos(ANISODIR * crys(i))
	      enddo
	    endif
	    mgnt = sqrt(sum(spinTot(1:SDIM) ** 2))
	  endif
	  if(ctrl(5)) then
	    nMgn = sqrt(sum(nSpnTot(1:SDIM) ** 2))
	  endif
	  if(ctrl(8)) then
	    n2Mg = sqrt(sum(n2SpTot(1:SDIM) ** 2))
	  endif
*	  Write results to outData
	  ptr = 0
	  if(ctrl(1)) then
	    outData(ptr + 1) = erg
	    outData(ptr + 2) = mgnt
	    do i = 1, LDIM
	      outData(ptr + 2 + i) = kMgnt(i)
	    enddo
	    ptr = ptr + SUBDATSZE(1)
	  endif
	  if(ctrl(2)) then
	    do i = 1, LDIM
	      outData(ptr + i) = lkErg(i)
	      outData(ptr + LDIM + i) = lkCur(i)
	    enddo
	    ptr = ptr + SUBDATSZE(2)
	  endif
	  if(ctrl(3)) then
	    outData(ptr + 1) = cMgn
	    ptr = ptr + SUBDATSZE(3)
	  endif
	  if(ctrl(5)) then
	    outData(ptr + 1) = nMgn
	    do i = 1, LDIM
	      outData(ptr + 1 + i) = kNMgn(i)
	    enddo
	    ptr = ptr + SUBDATSZE(5)
	  endif
	  if(ctrl(6)) then
	    outData(ptr + 1) = hMgn
	    ptr = ptr + SUBDATSZE(6)
	  endif
	  if(ctrl(8)) then
	    outData(ptr + 1) = n2Mg
	    do i = 1, LDIM
	      outData(ptr + 1 + i) = kN2Mg(i)
	    enddo
	    ptr = ptr + SUBDATSZE(8)
	  endif
	  if(ctrl(9)) then
	    do i = 1, distSze
	      outData(ptr + i) = dist(i)
	    enddo
	    ptr = ptr + distSze
	  endif
	end subroutine

	subroutine energyChangeIsingCube(crys, lSze, pos, ergDif)
	  integer*2, dimension(:), intent(in) :: crys
	  integer, intent(in) :: lSze, pos
	  integer, dimension(:) :: local(1:(2 * LDIM))
	  integer :: i, localF
	  integer, intent(out) :: ergDif
	  ergDif = 0
	  localF = 0
	  call getPBCCube(pos, lSze, local)
	  do i = 1, (2 * LDIM)
	    localF = localF + crys(local(i))
	  enddo
	  ergDif = ergDif + 2 * crys(pos) * localF
	end subroutine

	subroutine dataIsingCube(crys, lSze, erg, mgnt, lkErg, kMgn)
	  integer*2, dimension(:), intent(in) :: crys(:)
	  integer, intent(in) :: lSze
	  integer, dimension(:) :: coord(0:LDIM)
	  integer :: i, j, k, pos, adjPos
	  real, dimension(:) :: kMgnFT(1:(LDIM * 2)),
     $	                        trig(1:(lSze * 2))
	  integer, intent(out) :: erg, mgnt
	  integer, intent(out), dimension(:) :: lkErg(1:LDIM)
	  real, intent(out), dimension(:) :: kMgn(1:LDIM)
	  if(SDIM .ne. 1) then
	    write(*,*), 'Spin dimension error'
	    call exit(1)
	  endif
	  mgnt = 0
	  erg = 0
	  coord(0) = 0
	  do i = 1, (LDIM * 2)
	    kMgnFT(i) = 0.0
	  enddo
	  do i = 1, lSze
	    trig(2 * i - 1) = cos((i - 1) * TAU / lSze)
	    trig(2 * i) = sin((i - 1) * TAU / lSze)
	  enddo
	  do i = 1, LDIM
	    lkErg(i) = 0
	  enddo
	  do i = 1, size(crys)
	    mgnt = mgnt + crys(i)
	    call getCoordinateCube((i - 1), lSze, coord)
*	  component of the k-dependent magnetization
	    do j = 1, LDIM
	      do k = 1, 2
	        kMgnFT(2 * (j - 1) + k) =
     $	        kMgnFT(2 * (j - 1) + k) +
     $	        crys(i) * trig(2 * coord(LDIM - j + 1) + k)
	      enddo
	      kMgn(j) = sum(kMgnFT((2 * j - 1):(2 * j)) ** 2)
	    enddo
*	  component of the helicity modulus
	    do j = 1, LDIM
	      pos = coord(j)
	      coord(j) = mod((coord(j) + 1), lSze)
	      adjPos = 1
	      do k = 1, LDIM
	        adjPos = adjPos + coord(k) * (lSze ** (LDIM - k))
	      enddo
	      lkErg(j) = lkErg(j) + crys(i) * crys(adjPos)
	      coord(j) = pos
	    enddo
	  enddo
	  erg = sum(lkErg(1:LDIM)) 
	end subroutine

*********************************************************************
*	routines for specific cases
*********************************************************************
	subroutine get2dPBC(pos, lSze, r, u, l, d)
	  integer, intent(in) :: pos, lSze
	  integer :: lVol, lwBound, offset
	  integer, intent(out) :: r, u, l, d
	  lVol = lSze ** 2
	  lwBound = lVol - lSze
	  offset = lSze - 1
	  r = pos + 1
	  u = pos - lSze
	  l = pos - 1
	  d = pos + lSze
*	  4 corners
	  if(pos .eq. 1) then
	    u = lwBound + 1
	    l = lSze
	  else if(pos .eq. lSze) then
	    r = 1
	    u = lVol
	  else if(pos .eq. (lwBound + 1)) then
	    l = lVol
	    d = 1
	  else if(pos .eq. lVol) then
	    r = lwBound + 1
	    d = lSze
*	  4 boundaries
c	  upper boundary
	  else if((pos .lt. lSze) .and. (pos .gt. 1)) then
	    u = lwBound + pos
c	  lower boundary
	  else if((pos .lt. lVol) .and. (pos .gt. lwBound)) then
	    d = pos - lwBound
c	  left boundary
	  else if((pos .lt. lwBound) .and. (pos .gt. 1) .and.
     $	          (mod(pos - 1, lSze) .eq. 0)) then
	    l = pos + offset
c	  right boundary
	  else if((pos .lt. lVol) .and. (pos .gt. lSze) .and.
     $	          (mod(pos, lSze) .eq. 0)) then
	    r = pos - offset
	  endif
	end subroutine

	subroutine get1dPBC(pos, lSze, r, l)
	  integer, intent(in) :: pos, lSze
	  integer, intent(out) :: r, l
	  r = pos + 1
	  if(pos .eq. lSze)  then
	    r = 1
	  endif
	  l = pos - 1
	  if(l .eq. 0)  then
	    l = lSze
	  endif
	end subroutine

	subroutine get2dPBCWithMod(pos, lSze, r, u, l, d)
	  integer, intent(in) :: pos, lSze
	  integer :: lVol
	  integer, intent(out) :: r, u, l, d
	  lVol = lSze ** 2
	  u = pos - lSze
	  if(u .le. 0) then
	    u = lVol + u
	  endif
	  d = pos + lSze
	  if(d .gt. lVol) then
	    d = d - lVol
	  endif
	  r = pos + 1
	  if(mod(pos, lSze) .eq. 0) then
	    r = r - lSze
	  endif
	  l = pos - 1
	  if(mod(l, lSze) .eq. 0) then
	    l = l + lSze
	  endif
	end subroutine

	subroutine writeChainXY1d(crys, flID)
	  double precision, dimension(:), intent(in) :: crys(:)
	  integer, intent(in) :: flID
	  integer :: i
	  do i = 1, (size(crys) - 1)
	    write(flID, XYFORM, advance='no'), crys(i), SPACE
	  enddo
	  write(flID, XYFORM), crys(size(crys)), SPACE
	end subroutine

	subroutine readChainXY1d(crys, flID)
	  double precision, dimension(:), intent(in) :: crys(:)
	  integer, intent(in) :: flID
	end subroutine

	subroutine writeChainXY2d(crys, flID)
	  double precision, dimension(:), intent(in) :: crys(:)
	  integer, intent(in) :: flID
	  integer :: lSze, i
	  lSze = nint(sqrt(real(size(crys))))
	  if((lSze ** 2) .ne. size(crys)) then
	    write(*,*), 'Lattice size error'
	    call exit(1)
	  endif
	  do i = 1, lSze
	    call writeChainXY1d(crys(((i - 1) * lSze + 1) :
     $	                             (i * lSze)), flID)
	  enddo
	end subroutine

	subroutine readChainXY2d(crys, flID)
	  double precision, dimension(:), intent(in) :: crys(:)
	  integer, intent(in) :: flID
	end subroutine

	subroutine writeChainIsing1d(crys, flID)
	  integer*2, dimension(:), intent(in) :: crys(:)
	  integer, intent(in) :: flID
	  integer :: i
	  do i = 1, (size(crys) - 1)
	    if(crys(i) .eq. 1) then
	      write(flID, ISINGFORM, advance='no'), POSITIVE
	    else if(crys(i) .eq. -1) then
	      write(flID, ISINGFORM, advance='no'), NEGATIVE
	    else
	      write(*,*) 'Write error'
	      call exit(1)
	    endif
	  enddo
	  if(crys(i) .eq. 1) then
	    write(flID, ISINGFORM), POSITIVE
	  else if(crys(i) .eq. -1) then
	    write(flID, ISINGFORM), NEGATIVE
	  else
	    write(*,*) 'Write error'
	    call exit(1)
	  endif
	end subroutine

	subroutine readChainIsing1d(crys, flID)
	  integer*2, dimension(:), intent(out) :: crys
	  integer, intent(in) :: flID
	  integer :: i
	  character :: spn
	  do i = 1, (size(crys) - 1)
	    read(flId, ISINGFORM, advance='no'), spn
	    if(spn .eq. POSITIVE) then
	      crys(i) = 1
	    else if(spn .eq. NEGATIVE) then
	      crys(i) = -1
	    else
	      write(*,*) 'Read error'
	      call exit(1)
	    endif
	  enddo
	  read(flID, ISINGFORM), spn
	  if(spn .eq. POSITIVE) then
	    crys(size(crys)) = 1
	  else if (spn .eq. NEGATIVE) then
	    crys(size(crys)) = -1
	  else
	    write(*,*) 'Read error'
	    call exit(1)
	  endif
	end subroutine

	subroutine writeChainIsing2d(crys, flID)
	  integer*2, dimension(:), intent(in) :: crys
	  integer, intent(in) :: flID
	  integer :: lSze, i
	  lSze = nint(sqrt(real(size(crys))))
	  if((lSze ** 2) .ne. size(crys)) then
	    write(*,*), 'Lattice size error'
	    call exit(1)
	  endif
	  do i = 1, lSze
	    call writeChainIsing1d(crys(((i - 1) * lSze + 1) :
     $	                                (i * lSze)), flID)
	  enddo
	end subroutine

	subroutine readChainIsing2d(crys, flID)
	  integer*2, dimension(:), intent(in) :: crys
	  integer, intent(in) :: flID
	end subroutine

	subroutine writeFormatXY2d(crys, flID)
	  double precision, dimension(:), intent(in) :: crys(:)
	  integer, intent(in) :: flID
	  integer :: lSze, i, j
	  lSze = nint(sqrt(real(size(crys))))
	  if((lSze ** 2) .ne. size(crys)) then
	    write(*,*), 'Lattice size error'
	    call exit(1)
	  endif
	  do i = 1, lSze
	    do j = 1, (lSze - 1)
	      write(flID, XYFORM, advance='no'),
     $	           (crys((i - 1) * lSze + j) / TAU), SPACE
	    enddo 
	    write(flID, XYFORM), (crys(i * lSze) / TAU), SPACE
	  enddo  
	  write(flID, '(A1)'), SPACE
	end subroutine

	subroutine readFormatXY2d(crys, flID, ioSt)
	  double precision, dimension(:), intent(out) :: crys(:)
	  integer, intent(in) :: flID
	  integer :: lSze, i, j
	  character :: spn
	  integer, intent(out) :: ioSt
	  lSze = nint(sqrt(real(size(crys))))
	  if((lSze ** 2) .ne. size(crys)) then
	    write(*,*) 'Lattice size error'
	    call exit(1)
	  endif
	  do i = 1, lSze
	    do j = 1, (lSze - 1)
	      read(flID, XYFORM, advance='no', ioStat=ioSt),
     $	          crys((i - 1) * lSze + j), spn
	    enddo
	    read(flID, XYFORM, ioStat=ioSt), crys(i * lSze)
	  enddo
	  do i = 1, lSze ** 2
	    crys(i) = crys(i) * TAU
	  enddo
	  read(flID, '(A1)'), spn
	end subroutine

	subroutine writeUnformatXY2d(crys, flID)
	  double precision, dimension(:), intent(in) :: crys(:)
	  integer, intent(in) :: flID
	  write(flID), crys(1:size(crys))
	end subroutine

	subroutine readUnformatXY2d(crys, flID, ioSt)
	  double precision, dimension(:), intent(out) :: crys(:)
	  integer, intent(in) :: flID
	  integer :: ioSt
	  read(flID, ioStat=ioSt), crys(1:size(crys))
	end subroutine

	subroutine energyXY2d(crys, erg)
	  real, dimension(:), intent(in) :: crys(:)
	  integer :: lSze, i, l, r, u, d
	  double precision :: dif1, dif2
	  double precision, intent(out) :: erg
	  erg = 0.0
	  lSze = nint(sqrt(real(size(crys))))
	  if((lSze ** 2) .ne. size(crys)) then
	    write(*,*) 'Lattice size error'
	    call exit(1)
	  endif
	  do i = 1, size(crys)
	    call get2dPBCWithMod(i, lSze, r, u, l, d)
	    dif1 = crys(r) - crys(i)
	    dif2 = crys(d) - crys(i)
	    erg = erg - cos(dif1) - cos(dif2)
	  enddo
	end subroutine

      end module
