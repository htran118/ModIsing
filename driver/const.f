*********************************************************************
* module of global variables
*********************************************************************

      module constantClass
c	Dimension of lattice and of each spin
	integer, parameter :: LDIM = 2
	integer, parameter :: SDIM = 2
c	Size of random list, number of equilibrium step per spin
	integer, parameter :: RANLSTSZE = 10000
	integer, parameter :: AVGEQSTP = 50
c	Number of data part, size and name of each part 
c	(see xyAnal.f for further explanation)
	integer, parameter :: MAXDAT = 9
	integer, parameter, dimension(1:MAXDAT) :: 
     $	SUBDATSZE = (/ (3 + LDIM), (2 * LDIM), 1, 1, (1 + LDIM),
     $	                                       1, 1, (1 + LDIM),
     $	                                       0 /)
	character(len=5), parameter, dimension(1:MAXDAT) ::
     $	SUBDATNME = (/ 'prime', 'helix', 'cMgnt', 'acMgn', 'nMgnt',
     $	               '2Mgnt', 'a2Mgn', 'n2Mgn', 'distr' /)
*	Mask for unavailable data part
	integer, parameter, dimension(1:MAXDAT) :: 
     $	SUBDISABL = (/ 0, 0, 0, 1, 0, 0, 1, 0, 0 /)
c	2 pies
	double precision, parameter :: TAU = 8.0D0 * datan(1.D0)
	real, parameter :: RTAU = real(TAU)
c	Number of overrelaxation steps per spin, cutoff threshold
c	for Wolff algorithm (as portion of total number of spins)
	double precision, parameter :: ORRATE = 5.0D0 / 2
	double precision, parameter :: MSKCUT = 0.5D0
c	Number of anisotropic direction and strength of the
c	anisotropic Hamiltonian (take the standard energy
c	coupling to be 1.0)
	integer, parameter :: ANISODIR = 4
	double precision, parameter :: ANISOERG = 0.0D0 
c	Number of clock direction and strength of the
c	clock Hamiltonian (take the standard energy
c	coupling to be 1.0)
	integer, parameter :: CLOCKDIR = 2
	double precision, parameter :: CLOCKERG = 1.0D0
c	Misc.
	character, parameter :: POSITIVE = '+'
	character, parameter :: NEGATIVE = 'o'
	character, parameter :: SPACE = ' '
	character(len=64), parameter :: XYFORM = '(F9.6, A1)'
	character(len=64), parameter :: ISINGFORM = '(L1)'
	character(len=64), parameter :: WRITEFORM = 'unformatted'
	character(len=8), parameter :: STATFILE = 'statData'
	character(len=8), parameter :: CORRFILE = 'corrData'
	character(len=6), parameter :: RESPATH = 'result'
	character(len=3), parameter :: RAWPATH = 'raw'
      end module
