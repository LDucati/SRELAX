! Copyright (C) 2024-2025 L. N. Vidal
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

! ******************************************************************************
! Module containing shared variables
! ******************************************************************************
MODULE shared_vars

! Constants (SI units: CODATA)
real(8), parameter :: hbar=1.054571817D-34, hartree=4.3597447222060D-18, a0=5.29177210544D-11

! Static allocation
real(8)            :: timestep, splstep, limtime, Ispin, Qmom, B0ext, gammaI, &
                      sigmaEFG(5), GlmRe(5), GlmENR(5), SDmin, SDmax
integer            :: nQNUC, nframes, splframes, natom, nacfs
integer, parameter :: DYNPYU=10, acfflen=160
logical            :: ORCA, DYNPY, ACF_ENSAV, LTIME, RACF, ACF_AVG, PrintSD
character          :: trajdir*50, symbol*3, asymbol*2

! vector to label componentes of the 2nd-rank 
! spherical tensor R_2m: R_20 R_2+1 R_2-1 R_2+2 R_2-2
integer, dimension(5) :: mvec=(/0,1,-1,2,-2/)

! Dynamic allocation
complex(8), dimension(:,:), allocatable :: T2m
real(8),    dimension(:),   allocatable :: Xspl, Fspl, newX, newF
real(8),    dimension(:,:), allocatable :: Txyz, ACFre, ACFim
integer,    dimension(:),   allocatable :: atmlist
character(len=acfflen),  dimension(:),   allocatable :: acffiles

END MODULE

! ******************************************************************************
! Main code: Computes T1 and T2 relaxation times of the quadrupolar mechanism
! ******************************************************************************
PROGRAM Qrelax
Use shared_vars
Implicit None

integer :: ios, natm, alist
logical :: EFGexist

!$      integer :: nthreads, tid
!$      integer omp_get_num_threads,omp_get_thread_num
!$      external omp_get_num_threads,omp_get_thread_num
!$omp parallel private(tid,nthreads)
!$      nthreads = omp_get_num_threads()
!$      tid = omp_get_thread_num()
!$      if (tid == 1) then
!$         print '("QRELAX - Quadrupolar Relaxation program")'
!$         print '(" Parallel execution with number of threads: ",i3)', nthreads
!$      end if
!$omp end parallel

write(*,'(A/A/A//A/A/)')  &
" ************************* ", &
" OUTPUT OF PROGRAM QRELAX  ", &
" ************************* ", &
" This program computes relaxation times and related ", &
" quantities for the quadrupolar mechanism."

! Reading input file "qrelax.inp"
call read_input()

! Reading ACFs from a file
if (RACF) then
   call AVG_ACF()
else

   ! Extracting EFG from electronic structure outputs and saving them in EFG.xyz file
   inquire(file='EFG.xyz', exist=EFGexist )
   if ( EFGexist ) then
      write(*,'(/A/)') " => EFG.xyz file exists. Reading EFGs from it"
   else
      open(20,file='EFG.xyz',status='new')
        if (ORCA) then
           write(20,'(A)') " EFGs computed using ORCA code "
        else
           write(20,'(A)') " EFGs computed using ADF code "
        endif
        call WRITE_EFG()
      close(20)
   endif

   ! Timestep used for interpolation of EFGs
   if (splstep > 0.d0 ) then
   ! Add EFGs starting and ending at the same times of those EFGs 
   ! computed by the electronic structure code, but with timestep equal to 'splstep'
      splframes=NINT(timestep*DBLE(nframes-1)/splstep)+1
   else
      splframes=nframes
      splstep=timestep
   endif

   ! allocating memory for the EFG tensors (Txyz and T2m)
   if (splframes < nframes) then
      allocate(Txyz(nQNUC*6,nframes), STAT=ios)
      allocate(T2m(nQNUC*5,nframes), STAT=ios)
   else
      allocate(Txyz(nQNUC*6,splframes), STAT=ios)
      allocate(T2m(nQNUC*5,splframes), STAT=ios)
   endif
   call AllocStat(ios,'Txyz','Qrelax','A')
   call AllocStat(ios,'T2m','Qrelax','A')

   ! Creating DYNPY input file (CSV format)
   if (DYNPY) then
      write(*,'(/A/)') " Creating DYNPY 'ADF-efg.csv' input file"
      open(DYNPYU,file='ADF-efg.csv',status='unknown')
        write(DYNPYU,'(A)') &
"system,traj,frame,time,label,symbol,Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,V11,V22,V33,Eta,ComputeTime"
        call READ_EFG()
      close(DYNPYU)
      Stop " ADF-efg.csv successfully created!"
   endif

   ! read EFGs from EFG.xyz
   call READ_EFG()

   ! Using spline to increase the number of Cartesian EFG tensors 
   ! (thus the number of frames) through cubic spline interpolation
   call EFG_spliner()
   if ( splframes /= nframes ) then
       write(*,'(/A/A/A//A,I4,A,I6/)') &
   " --------------------------------", &
   " Using Cubic Spline Interpolation", &
   " --------------------------------", &
   " Number of frames changed from ",nframes, &
   " to ", splframes
   ! new timestep and nframes after interpolation 
      timestep=splstep
      nframes=splframes
   endif
   
   ! deallocating Txyz
   if (allocated(Txyz)) deallocate(Txyz, STAT=ios)
   call AllocStat(ios,'Txyz','Qrelax','D')
   
   ! EFG-autocorrelation functions are averaged over 'nQNUC' equivalent nuclei
   ! by ACF_R2m and WK_ACF_R2m routines
   if (nQNUC > 1) then
      write(*,'(/A/A,I4,A,A,A/A)') &
   " ---------------------------------------------------------------------------",     &
   " Averaged properties over ",nQNUC," equivalent ",symbol(:2)," nuclei from the following list:", &
   " ---------------------------------------------------------------------------"
      natm=1
      do while ( natm <= nQNUC)
         if (natm+10 < nQNUC) then
            write(*,'(100I4)') (atmlist(alist),alist=natm,natm+10)
         else
            write(*,'(100I4)') (atmlist(alist),alist=natm,nQNUC)
         end if
         natm = natm +10
      end do
   else
      write(*,'(/A/A,I4,A,A,A/A)') &
   " --------------------------------------------------------- ", &
   " Relaxation properties evaluated for atom number ", &
   atmlist(1)," (",symbol(:2),") ", &
   " --------------------------------------------------------- "
   endif
   write(*,'(/A,F15.5)') " Trajectory duration (ps) : ",(nframes-1)*timestep
   write(*,'(A,F15.5)') " Trajectory step (ps)     : ",timestep
   
   ! Allocating memory to store the Real and Imag 
   ! parts of the ACF_2m(EFG) functions
   allocate(ACFre(nQNUC*5,nframes), STAT=ios)
   call AllocStat(ios,'ACFre','Qrelax','A')
   allocate(ACFim(nQNUC*5,nframes), STAT=ios)
   call AllocStat(ios,'ACFim','Qrelax','A')
   
   ! Ensemble Average 
   if (ACF_ENSAV) then
       write(*,"(/A)") " Autocorrelation functions computed from Ensemble Averages! "
      call ACF_R2m()
   else
   ! FFT/Wiener-Khinchin Theorem
       write(*,"(/A)") " Autocorrelation functions computed using the Wiener-Khinchin Theorem "
      call WK_ACF_R2m()
   endif
   
   ! deallocating T2m 
   if (allocated(T2m)) deallocate(T2m, STAT=ios)
   call AllocStat(ios,'T2m','Qrelax','D')

endif

! Spectral Densities of the EFG-ACFs 
call SpectDens()

! deallocating ACF_2m
if (allocated(ACFre)) deallocate(ACFre, STAT=ios)
call AllocStat(ios,'ACFre','Qrelax','D')
if (allocated(ACFim)) deallocate(ACFim, STAT=ios)
call AllocStat(ios,'ACFim','Qrelax','D')

! Quadrupolar relaxation properties
call QrelaxProp()

! Date and Time
call LOCALINF
write(*,'(A/A/A/)')  &
" *********************** ",&
" End of all calculations ",&
" *********************** "
END PROGRAM Qrelax

! ******************************************************************************
! Routine to read the input file 'qrelax.inp'
! ******************************************************************************
SUBROUTINE read_input()
Use shared_vars
Implicit None

integer   :: ios, loop, leniso, Amass, strgsz
character :: line*80, emsg*80
real*8    :: ProtonFq=0.d0, Hspin, gammaH
logical   :: Isotp=.FALSE.

! Initializing some input vars
nQNUC=0; splframes=0; nframes=0
limtime=0.d0; Ispin=0.d0; Qmom=0.d0; timestep=0.d0; splstep=0.d0; B0ext=0.d0; gammaI=0.d0
ORCA=.FALSE.; DYNPY=.FALSE.; ACF_ENSAV=.FALSE.; LTIME=.FALSE.; ACF_AVG=.FALSE.; PrintSD=.FALSE.

! Opening input file 'qrelax.inp'
open(01,file='qrelax.inp',status='old',iostat=ios)

! Checking if file exists
 if (ios /= 0) call ERRMSG(01,"Could not open input file 'qrelax.inp")

! Endless loop to read input options
 do 
    read(01,'(A)',end=01) line

! Each component of the spherical EFG-ACF will be replaced by the mean
! value taken from its five components
    if ( line(1:13) == 'ForceIsotropy' ) ACF_AVG=.TRUE.
    
! relaxation properties will be calculated for those nuclei
! contained in the .EquivNuclei list
    if ( line(1:12) == '.EquivNuclei' ) then
       call COMMENT(01)
       read(01,*) nQNUC
       if (allocated(atmlist)) deallocate(atmlist, STAT=ios)
       call AllocStat(ios,'atmlist','read_input','D')
       allocate(atmlist(nQNUC), STAT=ios)
       call AllocStat(ios,'atmlist','read_input','A')
       call COMMENT(01)
       read(01,*) (atmlist(loop),loop=1,nQNUC)
    endif

! Reading EFGs from ORCA outputs
    if ( line(1:4) == 'ORCA' ) ORCA=.TRUE.

! DYNPY
! Will create an input file for dynpy code
    if ( line(1:5) == 'DYNPY' ) DYNPY=.TRUE.

! ACFs from ensemble averages (for testing purposes)
    if ( line(1:9) == 'ACF_ENSAV' ) ACF_ENSAV=.TRUE.

! Read ACFs to average them 
    if ( line(1:9) == '.READ_ACF' ) then
       RACF=.TRUE.
       call COMMENT(01)
       read(01,*) nacfs
       if (allocated(acffiles)) deallocate(acffiles, STAT=ios)
       call AllocStat(ios,'acffiles','read_input','D')
       allocate(acffiles(nacfs), STAT=ios)
       call AllocStat(ios,'acffiles','read_input','A')
       do loop=1,nacfs
          call COMMENT(01)
          read(01,'(A)') acffiles(loop)
       enddo
       write(*,'(/A,I4,A/A//A/A)') &
" This job will read ",nacfs," EGF-ACFs from the files listed below ", &
" and average them to compute relaxation properties ", &
" EFG_ACFs to be read:", &
" -------------------"
       do loop=1,nacfs
          call KBLANK(acffiles(loop),acfflen,strgsz)
          write(*,'(A)') acffiles(loop)(:strgsz)
       enddo
       write(*,*)
    endif

! Number of MD frames 
    if ( line(1:7) == '.frames' ) then
       call COMMENT(01)
       read(01,*) nframes
       if ( MOD(nframes,2).NE.0 ) then
          write(*,"(/A,I4/)") " @ Warning: 'nframes' MUST be even! Setting it to ",nframes-1
          nframes=nframes-1
       endif
    endif

! timestep between consecutive frames
    if ( line(1:9) == '.timestep' ) then
       call COMMENT(01)
       read(01,*) timestep
    endif

! Path to ORCA/ADF outputs files
    if ( line(1:11) == '.EFGoutputs' ) then
       call COMMENT(01)
       read(01,'(A)') trajdir
    endif

! Limit Time for the Correlation Function of EFG (ps)
    if ( line(1:10) == '.LimitTime' ) then
       call COMMENT(01)
       LTIME=.TRUE.
       read(01,*) limtime
    end if

! Get nuclear spin, magnetogyric ratio and eletric quadrupole moment for a given isotope
    if ( line(1:8) == '.Isotope' ) then
       Isotp=.TRUE.
       call COMMENT(01)
       read(01,'(A)') line
       call KBLANK(line,80,leniso)
       call getIsotope(line,leniso,Amass,asymbol)
       call SpinGammaI(Amass,asymbol,Ispin,gammaI)
       if ( gammaI == 0.d0 ) then
          emsg='Could not get nuclear spin and magnetogyric ratio for the isotope: '//line(:leniso)
          call ERRMSG(01,emsg)
       endif
       call EQuadrupol(Amass,asymbol,Qmom)
       if ( Qmom == 0.d0 ) then
          emsg='Could not get nuclear eletric quadrupole for the isotope: '//line(:leniso)
          call ERRMSG(01,emsg)
       endif
    endif

! Larmor frequency of the proton (MHz)
    if ( line(1:17) == '.Proton_Frequency' ) then
       call COMMENT(01)
       read(01,*) ProtonFq
    endif

! Printing spectral densities
    if ( line(1:16) == '.SpectralDensity' ) then
       PrintSD = .TRUE.
       call COMMENT(01)
       read(01,*) SDmin, SDmax
       write(*,"(/A/A,F10.2,A,F8.2,A/)") & 
" Spectral Densities will be calculated for magnetic fields (in T)",&
" within the interval: [", SDmin,",", SDmax,"]"
    endif

 end do
 01 continue

close(01)
if ( ProtonFq > 0.d0 ) then
   call SpinGammaI(1,'H ',Hspin,gammaH)
   B0ext = ProtonFq/gammaH
endif
if ( B0ext == 0.d0 ) then
   write(*,'(/A/A/)') &
" => Proton Larmor frequency not found in input file: setting it to zero.", &
"    As the Larmor frequency will be zero, this calculation implies the Extreme Narrowing Regime."
endif

if ( .NOT.Isotp )      call ERRMSG(0," Missing input option: '.Isotope' ")
if ( RACF ) return

if (nQNUC == 0 )       call ERRMSG(0," Missing input option: '.EquivNuclei' ")
if (nframes == 0 )     call ERRMSG(0," Missing input option: '.frames' ")
if (timestep == 0.d0 ) call ERRMSG(0," Missing input option: '.timestep' ")
END SUBROUTINE read_input

! ******************************************************************************
! Read EFG's from the EFG.xyz file
! ******************************************************************************
SUBROUTINE READ_EFG()
Use shared_vars
Implicit None

real(8)   :: EFGxyz(3,3), time
integer   :: frame, lin, col, nuc, ax, ixyz, ios
character :: string*80, emsg*125
logical   :: rewd=.FALSE., found

open(20,file='EFG.xyz',status='old',iostat=ios)
 if (ios /= 0) call ERRMSG(20,'File EFG.xyz does not exist!')

! reading EFGs of each frame
 do frame=1,nframes
    
! Frame time
    time = DBLE(frame-1)*timestep

! loop over equivalent nuclei
    do nuc=1,nQNUC
       write(string,'(A,D20.10,A,I4)') ' Time (ps): ',time,' ; EFG of atom number:',atmlist(nuc)
! Time (ps):     0.6250000000D-01 ; EFG of atom number:   1   ;   type: C
       call FINDSTRING (string,58,20,rewd,found)
       if (found) then
          read(20,'(A)') string
          symbol=string(72:73)//','
! Checking if read EFG belongs to the correct atom type
          if ( symbol(:2) /= asymbol ) then
             write(emsg,'(A,I3,A,D20.10,5A)') &
' Atom number ',atmlist(nuc),' from timestep ',time,' (type=',symbol(:2), &
') differs from that specified in qrelax.inp (type=',asymbol,')'
             call ERRMSG(20,emsg)
          endif
          read(20,'(A)') string
          do lin=1,3
             read(20,'(3X,3D30.15)') (EFGxyz(lin,col),col=1,3)
          end do

! --------------
! DYNPY CSV FILE
! --------------
!system, traj, frame, time, label, symbol, Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,V11,V22,V33,Eta,ComputeTime
! <time between frames in (ps)>
          if (DYNPY) then
             time = DBLE(frame-1)*timestep
             write(DYNPYU,'(A,I4,A,E14.6,A,I3,A,A,9(E14.6,A),A)') &
"ADF,1,", frame,',', time,',', atmlist(nuc),',', symbol, &
EFGxyz(1,1),',',EFGxyz(1,2),',',EFGxyz(1,3),',',EFGxyz(2,1),',',EFGxyz(2,2),',', &
EFGxyz(2,3),',',EFGxyz(3,1),',',EFGxyz(3,2),',',EFGxyz(3,3),',', &
" 0.,0.,0.,0.,0."
          endif

! Loading the Cartesian EFG (symmetric tensor) of atom "nuc" in Txyz(ax,frame), where:
! ax = 1, 2, 3, 4, 5, 6 for xx, xy, xz, yy, yz, zz of atom 1
!      7, 8, 9,10,11,12 for xx, xy, xz, yy, yz, zz of atom 2 ... 
          ixyz=0
          do lin=1,3
             do col=lin,3
                ixyz = ixyz+1
                ax   = (nuc-1)*6 +ixyz
                Txyz(ax,frame) = EFGxyz(lin,col)
             enddo
          enddo
       else
          write(emsg,'(3A)') "The following string is missing in EFG.xyz: ",string
          call ERRMSG(20,emsg)
       endif
    end do

 end do
close(20)
END SUBROUTINE READ_EFG

! ******************************************************************************
! READ EFG's from ADF/ORCA outputs and WRITE then into the EFG.xyz file
! ******************************************************************************
SUBROUTINE WRITE_EFG()
Use shared_vars
Implicit None

integer                 :: frame, lchar
integer, parameter      :: lEFGfile=80
character               :: fname*10
character(len=lEFGfile) :: EFGfile

! Path to the EFG file 
call KBLANK(trajdir,50,lchar)
do frame=1,nframes
   write(20,'(/A/A,I5/A,D20.10/A/)') &
" -------------------------------", &
" FRAME:   ",frame, &
" TIME (ps): ",DBLE(frame-1)*timestep, &
" -------------------------------"
   if (frame < 10                     ) write(fname,'(A,I1)') '/frame000',frame
   if (frame < 100   .and. frame >   9) write(fname,'(A,I2)') '/frame00',frame
   if (frame < 1000  .and. frame >  99) write(fname,'(A,I3)') '/frame0',frame
   if (frame < 10000 .and. frame > 999) write(fname,'(A,I4)') '/frame',frame
   EFGfile=trajdir(:lchar)//fname

! Reading EFG's
   if (ORCA) then
      call EFG_ORCA(EFGfile,lEFGfile,frame)
   else
      call EFG_ADF(EFGfile,lEFGfile,frame)
   end if
end do
END SUBROUTINE WRITE_EFG

! ******************************************************************************
! Read EFGs from ORCA output files
! ******************************************************************************
SUBROUTINE EFG_ORCA(EFGfile,lEFGfile,frame)
Use shared_vars
Implicit None

integer                 :: frame, lEFGfile, lline
character(len=lEFGfile) :: EFGfile
real(8)                 :: time
real(8), dimension(3,3) :: EFGxyz
character               :: strg*120, charXYZ(3)*1=(/'X','Y','Z'/), emsg*160
integer                 :: ios, lin, col, lchar, nuc
logical                 :: found, rewd, ORCA6

! Removing possible blank spaces from the output file/directory name
call KBLANK(EFGfile,lEFGfile,lchar)

open(01,file=EFGfile(:lchar),status='old',iostat=ios)

! Check if output file exists
 if (ios /= 0) then
    write(emsg,'(A,A)') "Failed to open the following ORCA output file: ",EFGfile(:lchar)
    call system('rm -f ./EFG.xyz')
    call ERRMSG(01,emsg)
 endif

! Finding ORCA's version (5 or 6)
 rewd=.True.
 call FINDSTRING ('                         Program Version',40,01,rewd,found)
 emsg="in function EFG_ORCA: Could not find ORCA's version within the output file"//EFGfile(:lchar)
 if (.not.found) then 
    call system('rm -f ./EFG.xyz')
    call ERRMSG(01,emsg)
 endif
 read(01,'(A)') strg 
 ORCA6=.FALSE.
 if ( strg(42:42) == '6') ORCA6=.TRUE.

! Find the number of atoms
 rewd=.True.
 call FINDSTRING ('Number of atoms                             ...',47,01,rewd,found)
 emsg="in function EFG_ORCA: Could not find the number of atoms in file "//EFGfile(:lchar)
 if (.not.found) then 
    call system('rm -f ./EFG.xyz')
    call ERRMSG(01,emsg)
 endif
 read(01,'(47X,I7)') natom

! Loop over those atoms from the '.EquivNuclei' list
 do nuc=1,natom

! Searching for the string: " Nucleus   5N : A:ISTP=   14 I=  1.0 P= 38.5677 MHz/au**3"
    if      (                nuc <   10 ) then
            write(strg,'(A,I1)') ' Nucleus   ',nuc-1 
    else if ( nuc >  9 .and. nuc <  100 ) then
            write(strg,'(A,I2)') ' Nucleus  ',nuc-1 
    else if ( nuc > 99 .and. nuc < 1000 ) then
            write(strg,'(A,I3)') ' Nucleus ',nuc-1 
    else
        stop ' @ Error in EFG_ORCA: system has more than 999 atoms!'
    end if
    rewd=.True.
    call FINDSTRING (strg(1:12),12,01,rewd,found)
    if (.not.found) then
       call system('rm -f ./EFG.xyz')
       write(emsg,'(A,A,A,A)') &
"could not find the EFG string '",strg(1:12),"' within ORCA output file: ",EFGfile(:lchar)
       call system('rm -f EFG.xyz')
       call ERRMSG(01,emsg)
    endif
    
    if (found) then

! Reading the EFG:

! Nucleus   5N : A:ISTP=   14 I=  1.0 P= 38.5677 MHz/au**3
!                Q:ISTP=   14 I=  1.0 Q=  0.0204 barn     
! -----------------------------------------------------------
! Raw EFG matrix (all values in a.u.**-3):
! -----------------------------------
!            -0.8172      -0.5297      -0.1931
!            -0.5297       0.4576       0.0818
!            -0.1931       0.0818       0.3596
       time = DBLE(frame-1)*timestep
       lline=5
       IF (ORCA6) lline=8
       do lin=1,lline
          read(01,'(A)') strg
          if (lin == 1) then
             symbol=strg(13:14)//','
             if (strg(14:14) == ' ') symbol=strg(13:13)//', '
             write(*,'(A,I3,4A,A,L,A)') &
"Extracting EFG for atom=(",nuc,", symbol ",symbol(1:index(symbol,',')-1),") from file: ",EFGfile(:lchar),' (ORCA6=',ORCA6,')'
          endif
       end do
! Writing data in EFG.xyz file
       write(20,'(A,D20.10,A,I4,2A)') &
" Time (ps): ",time," ; EFG of atom number:",nuc,"   ;   type: ",symbol(1:index(symbol,',')-1)
       write(20,'(A)') "            X                             Y                             Z"
       do lin=1,3
          IF (ORCA6)      read(01, '(7X,3F16.7)') (EFGxyz(lin,col),col=1,3)
          IF (.not.ORCA6) read(01, '(6X,3F13.4)') (EFGxyz(lin,col),col=1,3)
          write(20,'(1X,A,1X,3D30.15)') charXYZ(lin),(EFGxyz(lin,col),col=1,3)
       end do
    endif
 enddo
close(01)
END SUBROUTINE EFG_ORCA

! ******************************************************************************
! Read EFGs from ADF output files
! ******************************************************************************
SUBROUTINE EFG_ADF(EFGfile,lEFGfile,frame)
Use shared_vars
Implicit None

integer                 :: frame, lEFGfile
character(len=lEFGfile) :: EFGfile

real*8                  :: time
real(8), dimension(3,3) :: EFGxyz
character               :: strg*120, charXYZ(3)*1=(/'X','Y','Z'/), emsg*80
integer                 :: ios, lin, col, lchar, nuc
logical                 :: found, rewd

call KBLANK(EFGfile,lEFGfile,lchar)
open(01,file=EFGfile(:lchar),status='old',iostat=ios)

! Check if output file exists
 if (ios /= 0) then
    write(emsg,'(A,A)') "Failed to open the following ADF output file: ",EFGfile(:lchar)
    call system('rm -f ./EFG.xyz')
    call ERRMSG(01,emsg)
 endif

! Find the number of atoms
 rewd=.True.
 call FINDSTRING ('x',47,01,rewd,found)
 if (.not.found) then
    call system('rm -f ./EFG.xyz')
    call ERRMSG(01,"Function EFG_ADF: Could not find the number of atoms in file frame0001!")
 endif
 read(01,'(47X,I7)') natom

! Loop over those atoms to read EFG's
 do nuc=1,natom

! Searching for the string " Atom    1         N Input number    1  xyznuc(Angstrom) ="
    strg = ' Atom    '
    if ( atmlist(nuc) > 9) then
       write(strg(9:10), '(I2)') atmlist(nuc)
    else
       write(strg(10:10),'(I1)') atmlist(nuc)
    endif
    strg(11:11) = ' '
    rewd=.True.
    call FINDSTRING (strg(1:11),11,01,rewd,found)
    if (.not.found) then
       call system('rm -f ./EFG.xyz')
       write(emsg,'(A,A,A,A)') &
"Could not find string '",strg(1:11),"' within ADF output file: ",EFGfile(:lchar)
       call ERRMSG(01,emsg)
    endif

! Reading the EFG:
! Atom    1         N Input number    1  xyznuc(Angstrom) =     5.5953    7.9404    7.0538
!
! ==== EFG-tensor (a.u.) N
!            X             Y             Z
!     X   -0.164832E+00  0.685588E-01  0.398570E-01
!     Y    0.685588E-01  0.181437E+00  0.191942E+00
!     Z    0.398570E-01  0.191942E+00 -0.166043E-01
    time = DBLE(frame-1)*timestep
    do lin=1,4
       read(01,'(A)') strg
       if (lin == 1) then
          symbol=strg(13:14)//','
          if (strg(14:14) == ' ') symbol=strg(13:13)//', '
          write(*,'(A,I3,4A)') &
"Extracting EFG for atom=(",atmlist(nuc),", symbol ",symbol,") from file: ",EFGfile
       endif
    end do
    do lin=1,3
       read(01,'(8X,3E14.6)') (EFGxyz(lin,col),col=1,3)
!      write(*,'(8X,3E14.6)') (EFGxyz(lin,col),col=1,3)
    end do

! Writing data in EFG.xyz file
    write(20,'(A,D20.10,A,I4,2A)') &
" Time (ps): ",time," ; EFG of atom number:",nuc,"   ;   type: ",symbol(1:index(symbol,',')-1)
    write(20,'(A)') "            X                             Y                             Z"
    do lin=1,3
       write(20,'(1X,A,1X,3D30.15)') charXYZ(lin),(EFGxyz(lin,col),col=1,3)
    end do
 enddo
close(01)
END SUBROUTINE EFG_ADF

! ******************************************************************************
! Computing EFG spherical tensor R_lm (complex)
! from the Cartesian representation R_ab (real)
! ******************************************************************************
SUBROUTINE RxRm(frm,EFG,atm)
Use shared_vars
Implicit None

integer                 :: frm, atm
real(8), dimension(3,3) :: EFG

! T2m(a,b): 'a' index:  1  2  3  4  5 (atom 1) => R20 R21 R2_1 R22 R2_2
! T2m(a,b): 'a' index:  6  7  8  9 10 (atom 2) => R20 R21 R2_1 R22 R2_2
! T2m(a,b): 'a' index: 11 12 13 14 15 (atom 3) ...
! T2m(a,b): 'b' index: frame nmber
T2m(5*atm-4,frm) = 3.D0*(1/DSQRT(6.D0))*DCMPLX( EFG(3,3),                      0.D0)
T2m(5*atm-3,frm) =                      DCMPLX(-EFG(1,3),                 -EFG(2,3))
T2m(5*atm-2,frm) =                      DCMPLX( EFG(1,3),                 -EFG(2,3))
T2m(5*atm-1,frm) =                      DCMPLX( 0.5D0*(EFG(1,1)-EFG(2,2)), EFG(1,2))
T2m(5*atm,  frm) =                      DCMPLX( 0.5D0*(EFG(1,1)-EFG(2,2)),-EFG(1,2))
END SUBROUTINE RxRm

! ******************************************************************************
! ACF from ensemble averages: this routine computes ACFs from Eq.21 of ref:
!
! dx.doi.org/10.1021/ct400419s 
! J. Chem. Theory Comput. 2013, 9, 4074−4086
! Quadrupolar NMR Spin Relaxation Calculated Using Ab Initio Molecular Dynamics: 
! Group 1 and Group 17 Ions in Aqueous Solution
! Shyam Badu,† Lionel Truflandier,‡ and Jochen Autschbach
! ******************************************************************************
SUBROUTINE ACF_R2m() 
Use shared_vars
Implicit None

integer    :: frm1, frm2, nuc, am, tm, mcom, norm(nframes)
real(8)    :: meanRe, meanIm
complex(8) :: z

! ACF2m(m,n):
! n = 1 (tau =   0)
! n = 2 (tau =  dt)
! n = 3 (tau = 2dt)
! n = 4 (tau = 3dt) ...

! Initializing ACF with zeros
do frm1=1,nframes
   norm(frm1) = 0
   do nuc=1,5*nQNUC
      ACFre(nuc,frm1) = 0.d0
      ACFim(nuc,frm1) = 0.d0
   end do
end do

! Double Looping over frames to get time differences (tau)
do frm1=1,nframes
   do frm2=frm1,nframes
      tm = frm2-frm1+1
      norm(tm) = norm(tm) +1
! Loop over equivalent nuclei
      do nuc=1,nQNUC
!ACF_lm: lm = 20 21 2-1 22 2-2
         do mcom=1,5
            am = 5*(nuc-1) +mcom
            z = DCONJG(T2m(am,frm1)) * T2m(am,frm2)
            ACFre(am,tm) = ACFre(am,tm) +DBLE(z)
            ACFim(am,tm) = ACFim(am,tm) +DIMAG(z)
         end do
      end do
   end do
end do

! Averaging normalized ACFs over equivalent atoms
! loop over the 'm'-components of R2m
do mcom=1,5
! loop over trajectory frames
   do frm1=1,nframes
      meanRe=0.d0
      meanIm=0.d0
! loop over equivalent nuclei
      do nuc=1,nQNUC
          am = 5*(nuc-1) +mcom
          meanRe = meanRe +ACFre(am,frm1) / DBLE(norm(frm1))
          meanIm = meanIm +ACFim(am,frm1) / DBLE(norm(frm1))
      end do
! Overwriting the first atom's ACF with the mean value
      ACFre(mcom,frm1) = meanRe / DBLE(nQNUC)
      ACFim(mcom,frm1) = meanIm / DBLE(nQNUC)
   end do
end do

! Saving averaged over nuclei ACFs 
open(20,file='acf_efg-ensembleav.dat',status='unknown')
 write(20,'(A,1X,I5,1X,E20.10)') "# ", nframes, timestep
 write(20,'(2A)') &
"#            Time/ps     Re(f_{20})/a.u.    Re(f_{2+1})/a.u.    Re(f_{2-1})/a.u.    Re(f_{2+2})/a.u.    Re(f_{2-2})/a.u.",&
                     "    Im(f_{20})/a.u.    Im(f_{2+1})/a.u.    Im(f_{2-1})/a.u.    Im(f_{2+2})/a.u.    Im(f_{2-2})/a.u."
 do frm1=1,nframes
    write(20,'(11E20.10)') DBLE(frm1-1)*timestep, (ACFre(mcom,frm1),mcom=1,5), (ACFim(mcom,frm1),mcom=1,5)
 end do
close(20)
END SUBROUTINE ACF_R2m

! ******************************************************************************
! Autocorrelation Function of R_2m from Wiener-Khinchin Theorem
! Corr(g, g) ⇐⇒ |G(f)|
! ******************************************************************************
SUBROUTINE WK_ACF_R2m() 
Use shared_vars
Implicit None

integer                            :: nuc, frm, odd, mcom, am, ios, ISIGN, zpad
real(8), dimension(:), allocatable :: Rdata
real(8)                            :: meanRe, meanIm, pi

pi = 4.0d0*DATAN(1.d0)

zpad=2*nframes

! Array used in FFT-routine 'FOUR1'
allocate(Rdata(2*zpad), STAT=ios)
call AllocStat(ios,'Rdata','WK_ACF_R2m','A')

! Loop over equivalent nuclei
do nuc=1,nQNUC

! Loop over the "m" component of R2m: 
   do mcom=1,5

! Atom Number plus 'm'
! am : atom : m
! 1    1      0
! 2    1      1
! 3    1     -1
! 4    1      2
! 5    1     -2
! 6    2      0
      am = 5*(nuc-1)+mcom

! Computing the FT of R2m tensor
      do frm=1,nframes

! NR: Real      Part -> odd  numbers
!     Imaginary Part -> even numbers
         odd = 2*frm -1
         Rdata(odd)   =  DBLE(T2m(am,frm))   
         Rdata(odd+1) = DIMAG(T2m(am,frm))
      end do
      do frm=nframes+1,zpad
         odd = 2*frm -1
         Rdata(odd)   = 0.d0
         Rdata(odd+1) = 0.d0
      end do
      ISIGN=1
!      call FOUR1(Rdata,zpad,ISIGN)
      call DFourierT(Rdata,zpad,ISIGN)
   
! ACF = inverseFFT( |FFT(R2m)|^2 ) / N*(N-1)  : WK
      do frm=1,zpad
         odd = 2*frm -1
         Rdata(odd)   = Rdata(odd)**2 +Rdata(odd+1)**2
         Rdata(odd+1) = 0.d0
      enddo
      ISIGN=-1
!      call FOUR1(Rdata,zpad,ISIGN)
      call DFourierT(Rdata,zpad,ISIGN)

! Re[ ACF2m(am,n) ]:
! n = 0 (tau =    0)
! n = 1 (tau =   dt)
! n = 2 (tau =  2dt)
! nframes-1 = tf - ti
      do frm=1,nframes
         odd = 2*frm-1
         ACFre(am,frm) = Rdata(odd)  /DBLE(nframes-1)
         ACFim(am,frm) = Rdata(odd+1)/DBLE(nframes-1)
      end do
   end do
end do

! deallocating Rdata array
if (allocated(Rdata)) deallocate(Rdata, STAT=ios)
call AllocStat(ios,'Rdata','WK_ACF_R2m','D')

! Averaging ACFs over equivalent atoms
do mcom=1,5
   do frm=1,nframes
      meanRe=0.d0
      meanIm=0.d0
      do nuc=1,nQNUC
          am = 5*(nuc-1) +mcom
          meanRe = meanRe +ACFre(am,frm) 
          meanIm = meanIm +ACFim(am,frm) 
      end do
! Overwriting the first atom's ACF with the mean value
      ACFre(mcom,frm) = meanRe / DBLE(nQNUC)
      ACFim(mcom,frm) = meanIm / DBLE(nQNUC)
   end do
end do

! Saving ACFs within 'acf_r2m_wk.dat' file
open(20,file='acf_efg_wk.dat',status='unknown')
 write(20,'(A,1X,I5,1X,E20.10)') "# ", nframes, timestep
 write(20,'(2A)') &
"#            Time/ps     Re(f_{20})/a.u.    Re(f_{2+1})/a.u.    Re(f_{2-1})/a.u.    Re(f_{2+2})/a.u.    Re(f_{2-2})/a.u.",&
                     "    Im(f_{20})/a.u.    Im(f_{2+1})/a.u.    Im(f_{2-1})/a.u.    Im(f_{2+2})/a.u.    Im(f_{2-2})/a.u."
 do frm=1,nframes
    write(20,'(11E20.10)') DBLE(frm-1)*timestep, (ACFre(mcom,frm),mcom=1,5), (ACFim(mcom,frm),mcom=1,5)
 end do
close(20)
END SUBROUTINE WK_ACF_R2m

! ******************************************************************************
! Read ACFs from user given files and average them 
! ******************************************************************************
SUBROUTINE AVG_ACF()
Use shared_vars
Implicit None

real*8  :: time, tslist(nacfs), ReSum, ImSum, thres
integer :: acf, flen, ios, acfloc, mcom, frm, frmlist(nacfs)
character ::  cstring*160, cmt*2, io_emsg*160

do acf=1,nacfs

! acfloc: 
! 1, 2, 3, 4, 5 for m=0,1,-1,2,-2 of file 1
! 6, 7, 8, 9,10 for m=0,1,-1,2,-2 of file 2 and so on
   acfloc=(acf-1)*5

   call KBLANK(acffiles(acf),acfflen,flen)
   open(20,file=acffiles(acf)(:flen),status='old',iostat=ios,iomsg=io_emsg)
     cstring="Could not open file '"//acffiles(acf)(:flen)//"'"
     if (ios /= 0) write(*,*) io_emsg
     if (ios /= 0) call ERRMSG(20,cstring)

     ! reading ACFs
     read(20,*) cmt, nframes, timestep
     frmlist(acf) = nframes
     tslist(acf)  = timestep

     if ( acf == 1 ) then
        ! Allocating memory to store the Real and Imag parts of the ACF_2m(EFG) functions
        allocate(ACFre(nacfs*5,nframes), STAT=ios)
        call AllocStat(ios,'ACFre','AVG_ACF','A')
        allocate(ACFim(nacfs*5,nframes), STAT=ios)
        call AllocStat(ios,'ACFim','AVG_ACF','A')
     endif

     write(*,'(3A,I5,F15.5)') &
" Number of frames and timestep (ps) read from '",acffiles(acf)(:flen),"' file: ",nframes, timestep
     read(20,'(A)') cstring
     do frm=1,nframes
       read(20,*) time, (ACFre(acfloc+mcom,frm),mcom=1,5), (ACFim(acfloc+mcom,frm),mcom=1,5)
     end do
   close(20)
enddo

! Checking if nframes and timestep are the same in all files
thres=timestep/1.D4
do acf=2,nacfs
   if (     frmlist(1) /= frmlist(acf)           .or. &
        DABS(tslist(1)    -tslist(acf)) > thres )     &
call ERRMSG(00," Number of frames and timestep (first line) must be equal in all ACFs files!")
enddo

write(*,'(/A,F15.5)') " Trajectory duration (ps) : ",(nframes-1)*timestep
write(*,'(A,F15.5)') " Trajectory step (ps)     : ",timestep

! Averaging ACFs
nQNUC=1
do mcom=1,5
   do frm=1,nframes
      ReSum=0.d0
      ImSum=0.d0
      do acf=1,nacfs
         acfloc = (acf-1)*5 +mcom
         ReSum = ReSum +ACFre(acfloc,frm)
         ImSum = ImSum +ACFim(acfloc,frm)
      enddo
      ACFre(mcom,frm) = ReSum / DBLE(nacfs)
      ACFim(mcom,frm) = ImSum / DBLE(nacfs)
   enddo
enddo

! Saving ACFs within 'acf_r2m_wk.dat' file
open(20,file='acf_efg_averaged.dat',status='unknown')
 write(20,'(A,1X,I5,1X,E20.10)') "# ", nframes, timestep
 write(20,'(2A)') &
"#            Time/ps     Re(f_{20})/a.u.    Re(f_{2+1})/a.u.    Re(f_{2-1})/a.u.    Re(f_{2+2})/a.u.    Re(f_{2-2})/a.u.",&
                     "    Im(f_{20})/a.u.    Im(f_{2+1})/a.u.    Im(f_{2-1})/a.u.    Im(f_{2+2})/a.u.    Im(f_{2-2})/a.u."
 do frm=1,nframes
    write(20,'(11E20.10)') DBLE(frm-1)*timestep, (ACFre(mcom,frm),mcom=1,5), (ACFim(mcom,frm),mcom=1,5)
 end do
close(20)

END SUBROUTINE AVG_ACF

! ******************************************************************************
! Frequency Dependent Spectral Densities from the EFG-ACFs
! ******************************************************************************
SUBROUTINE SpectDens()
Use shared_vars
Implicit None

integer :: ios, frm, mcom, iw0
real(8) :: time, cs, sn, w0, wstep
real(8) :: AvgRe, AvgIm

real(8)           :: SINT, y2(nframes), pi
real(8), external :: ACFfunc

integer, parameter :: NptsSD=1000
real(8), dimension(:,:), allocatable :: SDen

pi = 4.d0*DATAN(1.d0)

if (LTIME) write(*,'(/A,F10.2/A/)') &
" Limit Time (in ps): t_lim = ",limtime, &
" ACF(EFG)(tau > t_lim) = Zero "

! Allocating memory to store spline arrays:
allocate(Xspl(nframes), STAT=ios)
call AllocStat(ios,'Xspl','SpectDens','A')
allocate(Fspl(nframes), STAT=ios)
call AllocStat(ios,'Fspl','SpectDens','A')

! Larmor Frequency (rad/s)
w0 = B0ext * gammaI * 1.d6 * 2.d0*pi

! Isotropic average of ACFs
if ( ACF_AVG ) then
   write(*,'(/A/)')" Option Force Isotropy is on! "
   do frm=1,nframes
      AvgRe=0.d0
      AvgIm=0.d0
      do mcom=1,5
         AvgRe=AvgRe+ACFre(mcom,frm)
         AvgIm=AvgIm+ACFim(mcom,frm)
      enddo
      do mcom=1,5
         ACFre(mcom,frm) = AvgRe / 5.d0
         ACFim(mcom,frm) = AvgIm / 5.d0
      enddo
   enddo
   open(20,file='acf_efg.iso',status='unknown')
    write(20,'(A,1X,I5,1X,E20.10)') "# ", nframes, timestep
    write(20,'(2A)') &
   "#            Time/ps     Re(f_{20})/a.u.    Re(f_{2+1})/a.u.    Re(f_{2-1})/a.u.    Re(f_{2+2})/a.u.    Re(f_{2-2})/a.u.",&
                        "    Im(f_{20})/a.u.    Im(f_{2+1})/a.u.    Im(f_{2-1})/a.u.    Im(f_{2+2})/a.u.    Im(f_{2-2})/a.u."
    do frm=1,nframes
       write(20,'(11E20.10)') DBLE(frm-1)*timestep, (ACFre(mcom,frm),mcom=1,5), (ACFim(mcom,frm),mcom=1,5)
    end do
   close(20)
endif

! Loop over 'm' components of ACF_2m: 0, 1, -1, 2, -2
do mcom=1,5

! EFG variance
   sigmaEFG(mcom) = ACFre(mcom,1)

   do frm=1,nframes
       time = DBLE(frm-1)*timestep
! Time (ps)
       Xspl(frm) = time
! for m =  0 -> w = 0
! for m =  1 -> w =   w0
! for m = -1 -> w =  -w0
! for m =  2 -> w =  2w0
! for m = -2 -> w = -2w0
       cs = DCOS( DBLE(mvec(mcom)) * w0 * time * 1.D-12 )
       sn = DSIN( DBLE(mvec(mcom)) * w0 * time * 1.D-12 )
       Fspl(frm) = ACFre(mcom,frm) * cs -ACFim(mcom,frm) * sn
! Limit Time option
       if ( LTIME .and. time > limtime ) Fspl(frm) = 0.d0
   end do

! Using Simpson's rule from Numerical Recipes in Fortran 77
! to integrate the ACFs x exp(i m w0 t)
!   call spline(Xspl,Fspl,nframes,yp1,ypn,y2)
   call y2spline(Xspl,Fspl,nframes,y2)
!   call qsimp(ACFfunc,nframes,Xspl,Fspl,y2,Xspl(1),Xspl(nframes),SINT)
   call simpson(ACFfunc,nframes,Xspl,Fspl,y2,Xspl(1),Xspl(nframes),SINT)
   GlmRe(mcom) = SINT

! Integrating pure ACFs (without the exponential factor)  
! to get correlation times
   do frm=1,nframes
       Fspl(frm) = ACFre(mcom,frm) 
   end do
!   call spline(Xspl,Fspl,nframes,yp1,ypn,y2)
   call y2spline(Xspl,Fspl,nframes,y2)
!   call qsimp(ACFfunc,nframes,Xspl,Fspl,y2,Xspl(1),Xspl(nframes),SINT)
   call simpson(ACFfunc,nframes,Xspl,Fspl,y2,Xspl(1),Xspl(nframes),SINT)
   GlmENR(mcom) = SINT

end do

! Calculating Spectral Densities for a selected range of Larmor frequencies
if (PrintSD) then

! Allocating memory to store spectral densities
   allocate(SDen(5,NptsSD), STAT=ios)
   call AllocStat(ios,'SDen','SpectDens','A')

! Loop over frequencies (in rad/s)
   SDmin = SDmin * gammaI * 1.d6 * 2.d0*pi
   SDmax = SDmax * gammaI * 1.d6 * 2.d0*pi
   wstep = (SDmax-SDmin)/dfloat(NptsSD-1)
   do mcom=1,5
      do iw0=1,NptsSD
         w0 = (SDmin + wstep*dfloat(iw0-1))
         do frm=1,nframes
             time = DBLE(frm-1)*timestep
             Xspl(frm) = time
             cs = DCOS( DBLE(mvec(mcom)) * w0 * time * 1.D-12 )
             sn = DSIN( DBLE(mvec(mcom)) * w0 * time * 1.D-12 )
             Fspl(frm) = ACFre(mcom,frm) * cs -ACFim(mcom,frm) * sn
         end do
         call y2spline(Xspl,Fspl,nframes,y2)
         call simpson(ACFfunc,nframes,Xspl,Fspl,y2,Xspl(1),Xspl(nframes),SINT)
         SDen(mcom,iw0) = SINT
      end do
   end do
   open(20,file='spectden.dat',status='unknown')
       SDmin = SDmin / 1.d6 / 2.d0 / pi
       SDmax = SDmax / 1.d6 / 2.d0 / pi
       wstep = wstep / 1.d6 / 2.d0 / pi
       write(20,"(/A/A,F10.2,A,F10.2,A/)") &
"# Spectral Densities, in units of (Eh x e^-1 x a0^-2)^2 x ps, calculated for", &
"# Larmor frequencies (in MHz) within the range: [", SDmin,",", SDmax,"]"
       write(20,'(2A)') &
"# Larmor Freq./MHz      Re(g_{20})/a.u.     Re(g_{2+1})/a.u.    Re(g_{2-1})/a.u.    Re(g_{2+2})/a.u.    Re(g_{2-2})/a.u."
       do iw0=1,NptsSD
          w0 = (SDmin + wstep*dfloat(iw0-1))  
          write(20,'(F20.10,5E20.10)') w0, (SDen(mcom,iw0),mcom=1,5)
       end do
   close(20)
   if (allocated(SDen)) deallocate(SDen, STAT=ios)
   call AllocStat(ios,'SDen','SpectDens','D')
endif

! deallocating Spline arrays
if (allocated(Xspl)) deallocate(Xspl, STAT=ios)
call AllocStat(ios,'Xspl','SpectDens','D')
if (allocated(Fspl)) deallocate(Fspl, STAT=ios)
call AllocStat(ios,'Fspl','SpectDens','D')

END SUBROUTINE SpectDens

! ******************************************************************************
! Quadrupolar relaxation properties
! ******************************************************************************
SUBROUTINE QrelaxProp()
Use shared_vars
Implicit None

integer :: mcom
real(8) :: tauEFG(5)
real(8) :: G20, G21, au2SI, Rfac, w0, pi

pi = 4.0d0*DATAN(1.d0)

write(*,'(/A/A/A//A)') & 
" --------------------------------- ", &
" EFG variance and Correlation Time ", &
" --------------------------------- ", &
"  m         <V(0)^2>/a.u.          tau/ps"

! Larmor Frequency (MHz)
w0 = B0ext*gammaI

! Loop over 'm' components of ACF_2m: 0, 1, -1, 2, -2
do mcom=1,5
! correlation time
   tauEFG(mcom) = GlmENR(mcom)/sigmaEFG(mcom)
   write(*,'(I3,2F20.7)') mvec(mcom), sigmaEFG(mcom), tauEFG(mcom)
end do

! Linear combinations of spectral densities
! R20 R21 R2_1 R22 R2_2
! 1   2   3    4   5
!       4*g_{2,2}          +g_{2,1}       +g_{2,-1} +4*g_{2,-2}
G20   = 4.d0*GlmRe(4)      +GlmRe(2)      +GlmRe(3) +4.d0*GlmRe(5)
!       2*g_{2,-2}    +3*g_{2,-1}    +2*g_{2,1}     +3*g_{2,0}
G21   = 2.d0*GlmRe(5) +3.d0*GlmRe(3) +2.d0*GlmRe(2) +3.d0*GlmRe(1)

! Conversion of Q (mbarn) and EFG (Eh/e/a0^2) to SI units
au2SI = 1.D-31 * hartree / a0**2 / hbar
! Time is given in ps
au2SI = 1.d-12 * au2SI**2
Rfac = au2SI * Qmom**2 * (2.d0*Ispin+3.d0)/40.d0/Ispin**2/(2.d0*Ispin-1.d0)

write(*,"(/A/A/A/)") &
" ----------------------- ", &
" Nuclear spin properties ", &
" ------------------------ "
write(*,'(A,I3//A,F12.4//A,F12.4//A,F12.4//A,F12.4//A,E12.4/)')  &
" Nuclear spin quantum number: ", nint(Ispin), &
" Electric quadrupole moment (mbarn): ", Qmom, &
" External magnetic field (T):        ", B0ext, &
" Magnetogyric ratio (MHz/T):         ", gammaI, &
" Larmor frequency (MHz):             ", w0, &
" Larmor angular frequency (rad/s):   ", w0*1.d6*2.d0*pi

write(*,"(/A/A/A)") &
" -------------------------- ", &
" Relaxation Times and Rates ", &
" -------------------------- "
write(*,'(/A,F13.3/A,F13.3//A,F13.5/A,F13.5/)') &
" T1 (ms) : ", 1.d3/Rfac/G20, &
" T2 (ms) : ", 1.d3/Rfac/G21, &
" R1 (Hz) : ", Rfac*G20, &
" R2 (Hz) : ", Rfac*G21
END SUBROUTINE QrelaxProp

! ******************************************************************************
! This 'external' function uses spline interpolation to estimate ACF_R2m at 'arg'
! It's used in simpson's qsimp routine
! ******************************************************************************
DOUBLE PRECISION FUNCTION ACFfunc(arg,ndim,Xspl,Fspl,y2)
Implicit None

real*8  :: arg, Xspl(*), Fspl(*), y2(*), SPL
integer :: ndim

!call splint(Xspl,Fspl,y2,ndim,arg,SPL)
call ncspline(arg,Xspl,Fspl,ndim,y2,SPL)
ACFfunc = SPL
END FUNCTION ACFfunc

! ******************************************************************************
! Splining the Cartesian EFGs computed by eletronic structure codes
! ******************************************************************************
SUBROUTINE EFG_spliner()
Use shared_vars
Implicit None

real*8    :: EFGxyz(3,3), time
integer   :: nuc, xcom, frm, lin, col, ax, ios
character :: routine*11='EFG_spliner'

! Allocating memory to store spline arrays:
allocate(Xspl(nframes), STAT=ios)
call AllocStat(ios,'Xspl',routine,'A')
allocate(Fspl(nframes), STAT=ios)
call AllocStat(ios,'Fspl',routine,'A')
allocate(newX(splframes), STAT=ios)
call AllocStat(ios,'newX',routine,'A')
allocate(newF(splframes), STAT=ios)
call AllocStat(ios,'newF',routine,'A')

! Loop over the '6' unequal Cartesian componentes of the EFG 
! (actually there is five independent components)
do nuc=1,nQNUC
   do xcom=1,6
      ax = (nuc-1)*6 +xcom

! loading input function
      do frm=1,nframes
         time      = DBLE(frm-1)*timestep
         Xspl(frm) = time
         Fspl(frm) = Txyz(ax,frm)
!write(30+xcom,*) time, Fspl(frm)
      end do
! calling the spliner for the EFG of a given nucleus
      call SpliningF()

! loading spline-generated EFGs
      do frm=1,splframes
         Txyz(ax,frm) = newF(frm)
      end do
   end do
   do frm=1,splframes
      ax = (nuc-1)*6 
      do lin=1,3
         do col=lin,3
            ax = ax+1
            EFGxyz(lin,col) = Txyz(ax,frm)
         end do
      end do
      EFGxyz(2,1) = EFGxyz(1,2)
      EFGxyz(3,1) = EFGxyz(1,3)
      EFGxyz(3,2) = EFGxyz(2,3)
!time = dble(frm-1)*splstep
!write(41,*) time,Txyz((nuc-1)*6+1,frm)
!write(42,*) time,Txyz((nuc-1)*6+2,frm)
!write(43,*) time,Txyz((nuc-1)*6+3,frm)
!write(44,*) time,Txyz((nuc-1)*6+4,frm)
!write(45,*) time,Txyz((nuc-1)*6+5,frm)
!write(46,*) time,Txyz((nuc-1)*6+6,frm)
! computing the spherical tensors 
      call RxRm(frm,EFGxyz,nuc)
   end do
end do

! Deallocating Spline arrays
if (allocated(Xspl)) deallocate(Xspl, STAT=ios)
call AllocStat(ios,'Xspl',routine,'D')
if (allocated(Fspl)) deallocate(Fspl, STAT=ios)
call AllocStat(ios,'Fspl',routine,'D')
if (allocated(newX)) deallocate(newX, STAT=ios)
call AllocStat(ios,'newX',routine,'D')
if (allocated(newF)) deallocate(newF, STAT=ios)
call AllocStat(ios,'newF',routine,'D')
END SUBROUTINE EFG_spliner

! ******************************************************************************
! Splining a function
! ******************************************************************************
SUBROUTINE SpliningF()
Use shared_vars
Implicit None

real*8 :: time

! Numerical Recipes in Fortran77 spline (natural spline)
integer :: frm
real*8  :: SPL, y2(nframes)

! spline from Numerical Recipes in Fortran77
!call spline(Xspl,Fspl,nframes,yp1,ypn,y2)
call y2spline(Xspl,Fspl,nframes,y2)
do frm=1,splframes
   time = splstep*DBLE(frm-1)
   newX(frm) = time
!   call splint(Xspl,Fspl,y2,nframes,time,SPL)
   call ncspline(time,Xspl,Fspl,nframes,y2,SPL)
   newF(frm) = SPL
end do
END SUBROUTINE SpliningF

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Check if dynamic (de)allocation succeeded
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE AllocStat(ios,var,routine,dalloc)
Implicit None

integer   :: ios
character :: var*(*), routine*(*), dalloc*1

if (ios /= 0) then
   if (dalloc == 'A') then
      write(*,'(/4A/)') " @ Error: Allocation failed for '",var,"' array in routine '",routine,"'!"
   else
      write(*,'(/4A/)') " @ Error: Deallocation failed for '",var,"' array in routine '",routine,"'!"
   endif
   Stop
endif
END SUBROUTINE AllocStat

! ******************************************************************************
! Print an error message and stop the code
! ******************************************************************************
SUBROUTINE ERRMSG(UNIT,emessage)
Implicit None
integer   :: UNIT
character :: emessage*(*)

write(*,'(/A,A/)') " @ Error: ", emessage
if (UNIT /= 0) close(UNIT)
Stop
END SUBROUTINE ERRMSG

! ******************************************************************************
! Discrete Fourier Transform from Numerical Recipes in Fortran 77 book
! ******************************************************************************
! Real      Part -> odd  numbers
! Imaginary Part -> even numbers
! ******************************************************************************
SUBROUTINE DFourierT(Rdata,npts,ISIGN)
Implicit None
real*8     :: Rdata(*), HnRe(npts), HnIm(npts)
real*8     :: Theta, c, s, Nfac, Rodd, Reven, sRe, sIm, Ofac
integer    :: n, k, odd, npts, ISIGN

Ofac = 2.d0 * 4.0d0*DATAN(1.d0) * DBLE(ISIGN)/DBLE(npts)

! For inverse transform the results must be divided by 1/npts
Nfac=1.d0
if ( ISIGN == -1 ) Nfac = DBLE(npts)

! Loop over discrete frequencies
!$OMP PARALLEL SHARED(npts,Rdata,HnRe,HnIm,Nfac) PRIVATE(sRe,sIm,n,odd,Theta,Rodd,Reven,c,s)
!$OMP DO
do n=0,npts-1
   sRe = 0.d0
   sIm = 0.d0
! Loop over function points
   do k=0,npts-1
      odd = 2*(k+1) -1
      Rodd  = Rdata(odd)
      Reven = Rdata(odd+1)
! Frequency here is not circular, i.e., it has unit of Hz only (without radians)
      Theta = Ofac * DBLE(k * n)
      c = DCOS(Theta)
      s = DSIN(Theta)
      sRe = sRe + c * Rodd  -s * Reven
      sIm = sIm + c * Reven +s * Rodd
   end do
   HnRe(n+1) = sRe / Nfac
   HnIm(n+1) = sIm / Nfac
enddo
!$OMP END DO
!$OMP END PARALLEL
do k=0,npts-1
   odd = 2*(k+1) -1
   Rdata(odd)   = HnRe(k+1)
   Rdata(odd+1) = HnIm(k+1)
enddo
END SUBROUTINE DFourierT

! ******************************************************************************
! Natural Cubic Spline Interpolation
! 
! Based on sections "2.3 LU Decomposition and Its Applications"
! and "3.3 Cubic Spline Interpolation" of Numerical Recipes in Fortran 77
! ******************************************************************************
SUBROUTINE ncspline(arg,x,y,ndim,y2,fx)
Implicit None

Integer :: ndim, j, jhi, jlo
Real*8  :: x(*), y(*), y2(ndim), arg, fx, A,B, xx

! Finding shortest interval where x_j \le arg \le x_{j+1},
! by means of bisection.
jlo=1
jhi=ndim
do while (jhi-jlo.gt.1)
   j=(jhi+jlo)/2
   if (x(j).gt.arg) then
      jhi=j
   else
      jlo=j
   endif
end do
 
xx = x(jhi) -x(jlo)
if ( xx < 1.d-10 ) call ERRMSG(0,"x(jhi) -x(jlo) < 1e-20 in ncspline")
A  = (x(jhi) -arg)/xx
B  = 1.d0 -A
fx = A*y(jlo) +B*y(jhi) +( (A**3-A) * xx**2 * y2(jlo)            &
                          +(B**3-B) * xx**2 * y2(jhi) ) / 6.d0
END SUBROUTINE ncspline

! ******************************************************************************
! Generating second derivatives of y(x_i) by
! solving the tridiagonal system of equations (Eq. 3.3.7 of NR)
! using LU decomposition 
! ******************************************************************************
SUBROUTINE y2spline(x,y,ndim,y2)
Implicit None

Integer :: ndim, i
Real*8  :: x(*), y(*), y2(ndim) 
Real*8  :: xi, yi, xp, yp, xm, ym, z(ndim), alpha

! Natural spline has y"_1 = y"_N = 0
y2(1)    = 0.d0
y2(ndim) = 0.d0

! forward substitution
! "beta_ii" is stored in "y2" 
z(2)  = (y(3)-y(2))/(x(3)-x(2)) &
       -(y(2)-y(1))/(x(2)-x(1)) 
y2(2) = (x(3)-x(1))/3.d0
do i=3,ndim-1
 yi = y(i)
 yp = y(i+1)
 ym = y(i-1)
 xi = x(i)
 xp = x(i+1)
 xm = x(i-1)
 alpha = (xi-xm)/6.d0/y2(i-1)
 y2(i) = (xp-xm)/3.d0 -alpha * (xi-xm)/6.d0
 z(i)  = (yp-yi)/(xp-xi) -(yi-ym)/(xi-xm) -alpha * z(i-1)
end do

!backsubstitution
do i=ndim-1,2,-1
 y2(i) = ( z(i) -y2(i+1) * (x(i+1)-x(i))/6.d0 ) / y2(i)
end do
END SUBROUTINE y2spline

!*******************************************************************************
! Extended Simpson’s rule, Eq.(4.1.13) and extended Trapezoidal rule, 
! Eq.(4.1.11), of Numerical Recipes in Fortran 77
! Simpson's rule is implemented through Eq.(4.2.4): S = 4/3 * S2N -1/3 * SN,
! where S2N and SN are integrals computed using Eq.(4.1.11) with 2N and N points
!*******************************************************************************
SUBROUTINE simpson(func,npts,Xspl,Fspl,y2,x1,xN,Sint)
Implicit None

Integer   :: Nmax, npts, ngrid, i, j
Real*8    :: x, step, x1, xN, pSint, sum, f1N, SN, S2N, func, Thresh, Sint
Real*8    :: Xspl(*), Fspl(*), y2(*)
Parameter (Thresh=1.d-6, Nmax=20)
External func

!----------------------------------------------------
!Extended Trapezoidal rule for 9 points evenly spaced
!----------------------------------------------------
! Contribution from f(x_1) and f(x_N) 
f1N = ( func(x1,npts,Xspl,Fspl,y2) &
       +func(xN,npts,Xspl,Fspl,y2) ) * 0.5d0
! Middle points now: f(x_2) ... f(x_8)
ngrid = 8
step = (xN-x1)/Dble(ngrid)
sum=0.d0
do i=2,ngrid
   x = Dble(i-1) * step
   sum = sum +func(x,npts,Xspl,Fspl,y2)
end do
SN = (f1N+sum) * step
!------------------------------------------------------------
!Increasing (doubling) the number of points until convergence
!------------------------------------------------------------
pSint=-1.d80
do i = 1,Nmax
   ngrid = ngrid * 2
   step = step / 2.d0
! Only even x_i are considered here to avoid calculating 
! f(x) many times at the same point
   do j=2,ngrid+1,2
      x = Dble(j-1) * step
      sum = sum +func(x,npts,Xspl,Fspl,y2)
   end do
   S2N = (f1N+sum) * step
   Sint = (4.d0 * S2N -SN) / 3.d0
   if ( Dabs(Sint-pSint).LT.Thresh ) return
   SN=S2N
   pSint=Sint
end do
call ERRMSG(0,"Simpson integration did not converge!")
END SUBROUTINE simpson

!*******************************************************************************
! Input: atomic mass and atomic symbol 
! Output: nuclear spin and magnetogyric ratio (MHz/T)
! Data from:
! EasySpin Nuclear isotope database :
! The associated file is easyspin/private/isotopedata.txt.  
!
! Stefan Stoll, Arthur Schweiger
! EasySpin, a comprehensive software package for spectral simulation and analysis in EPR
! J. Magn. Reson. 178(1), 42-55 (2006) 
!*******************************************************************************
SUBROUTINE SpinGammaI(mass,symbol,Inuc,gnuc)
Implicit None

integer, parameter :: lisotope=351
integer            :: Amass(lisotope), mass, isotope
character(len=2)   :: Asymbol(lisotope), symbol
real(8)            :: gnI(lisotope), spin(lisotope), Inuc, gnuc, eCharge, pMass, Pi

! Atomic Mass & Symbol & Spin & Gyromagnetic ratio
data (Amass(isotope),isotope=1,lisotope) / &
  1,  2,  3,  3,  4,  6,  7,  9, 10, 11, 12, 13, 14, 14, 15, 16, 17, 18, 19, 20, & 
 21, 22, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 36, 35, 36, 37, 36, &
 38, 39, 40, 39, 40, 41, 40, 41, 42, 43, 44, 46, 48, 45, 46, 47, 48, 49, 50, 50, &
 51, 50, 52, 53, 54, 53, 55, 54, 56, 57, 58, 59, 60, 58, 60, 61, 62, 64, 63, 65, &
 64, 66, 67, 68, 70, 69, 71, 70, 72, 73, 74, 76, 75, 74, 76, 77, 78, 79, 80, 82, &
 79, 81, 78, 80, 82, 83, 84, 85, 86, 85, 87, 84, 86, 87, 88, 89, 90, 91, 92, 94, &
 96, 93, 92, 94, 95, 96, 97, 98,100, 99, 96, 98, 99,100,101,102,104,103,102,104, &
105,106,108,110,107,109,106,108,110,111,112,113,114,116,113,115,112,114,115,116, &
117,118,119,120,122,124,121,123,125,120,122,123,124,125,126,128,130,127,129,124, &
126,128,129,130,131,132,134,136,133,134,135,137,130,132,133,134,135,136,137,138, &
137,138,139,136,138,140,142,141,142,143,144,145,146,148,150,147,144,147,148,149, &
150,151,152,154,151,152,153,154,155,152,154,155,156,157,158,160,157,159,160,156, &
158,160,161,162,163,164,165,162,164,166,167,168,170,169,171,168,170,171,172,173, &
174,176,173,174,175,176,174,176,177,178,179,180,180,181,180,182,183,184,186,185, &
187,184,186,187,188,189,190,192,191,193,190,192,194,195,196,198,197,196,198,199, &
200,201,202,204,203,204,205,204,206,207,208,207,209,209,  0,  0,  0,  0,227,229, &
232,231,234,235,238,237,239,243,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
data (Asymbol(isotope),isotope=1,lisotope) / &
"H ", "H ", "H ", "He", "He", "Li", "Li", "Be", "B ", "B ", "C ", "C ", "C ", & 
"N ", "N ", "O ", "O ", "O ", "F ", "Ne", "Ne", "Ne", "Na", "Na", "Mg", "Mg", &
"Mg", "Al", "Si", "Si", "Si", "P ", "S ", "S ", "S ", "S ", "Cl", "Cl", "Cl", &
"Ar", "Ar", "Ar", "Ar", "K ", "K ", "K ", "Ca", "Ca", "Ca", "Ca", "Ca", "Ca", &
"Ca", "Sc", "Ti", "Ti", "Ti", "Ti", "Ti", "V ", "V ", "Cr", "Cr", "Cr", "Cr", &
"Mn", "Mn", "Fe", "Fe", "Fe", "Fe", "Co", "Co", "Ni", "Ni", "Ni", "Ni", "Ni", &
"Cu", "Cu", "Zn", "Zn", "Zn", "Zn", "Zn", "Ga", "Ga", "Ge", "Ge", "Ge", "Ge", &
"Ge", "As", "Se", "Se", "Se", "Se", "Se", "Se", "Se", "Br", "Br", "Kr", "Kr", &
"Kr", "Kr", "Kr", "Kr", "Kr", "Rb", "Rb", "Sr", "Sr", "Sr", "Sr", "Y ", "Zr", &
"Zr", "Zr", "Zr", "Zr", "Nb", "Mo", "Mo", "Mo", "Mo", "Mo", "Mo", "Mo", "Tc", &
"Ru", "Ru", "Ru", "Ru", "Ru", "Ru", "Ru", "Rh", "Pd", "Pd", "Pd", "Pd", "Pd", &
"Pd", "Ag", "Ag", "Cd", "Cd", "Cd", "Cd", "Cd", "Cd", "Cd", "Cd", "In", "In", &
"Sn", "Sn", "Sn", "Sn", "Sn", "Sn", "Sn", "Sn", "Sn", "Sn", "Sb", "Sb", "Sb", &
"Te", "Te", "Te", "Te", "Te", "Te", "Te", "Te", "I ", "I ", "Xe", "Xe", "Xe", &
"Xe", "Xe", "Xe", "Xe", "Xe", "Xe", "Cs", "Cs", "Cs", "Cs", "Ba", "Ba", "Ba", &
"Ba", "Ba", "Ba", "Ba", "Ba", "La", "La", "La", "Ce", "Ce", "Ce", "Ce", "Pr", &
"Nd", "Nd", "Nd", "Nd", "Nd", "Nd", "Nd", "Pm", "Sm", "Sm", "Sm", "Sm", "Sm", &
"Sm", "Sm", "Sm", "Eu", "Eu", "Eu", "Eu", "Eu", "Gd", "Gd", "Gd", "Gd", "Gd", & 
"Gd", "Gd", "Tb", "Tb", "Tb", "Dy", "Dy", "Dy", "Dy", "Dy", "Dy", "Dy", "Ho", &
"Er", "Er", "Er", "Er", "Er", "Er", "Tm", "Tm", "Yb", "Yb", "Yb", "Yb", "Yb", &
"Yb", "Yb", "Lu", "Lu", "Lu", "Lu", "Hf", "Hf", "Hf", "Hf", "Hf", "Hf", "Ta", &
"Ta", "W ", "W ", "W ", "W ", "W ", "Re", "Re", "Os", "Os", "Os", "Os", "Os", &
"Os", "Os", "Ir", "Ir", "Pt", "Pt", "Pt", "Pt", "Pt", "Pt", "Au", "Hg", "Hg", &
"Hg", "Hg", "Hg", "Hg", "Hg", "Tl", "Tl", "Tl", "Pb", "Pb", "Pb", "Pb", "Bi", &
"Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Th", "Pa", "U ", "U ", "U ", &
"Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", &
"Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og" /
data (spin(isotope),isotope=1,lisotope) / &
0.5d0, 1.0d0, 0.5d0, 0.5d0, 0.0d0, 1.0d0, 1.5d0, 1.5d0, 3.0d0, 1.5d0, 0.0d0, & 
0.5d0, 0.0d0, 1.0d0, 0.5d0, 0.0d0, 2.5d0, 0.0d0, 0.5d0, 0.0d0, 1.5d0, 0.0d0, &
3.0d0, 1.5d0, 0.0d0, 2.5d0, 0.0d0, 2.5d0, 0.0d0, 0.5d0, 0.0d0, 0.5d0, 0.0d0, & 
1.5d0, 0.0d0, 0.0d0, 1.5d0, 2.0d0, 1.5d0, 0.0d0, 0.0d0, 3.5d0, 0.0d0, 1.5d0, & 
4.0d0, 1.5d0, 0.0d0, 3.5d0, 0.0d0, 3.5d0, 0.0d0, 0.0d0, 0.0d0, 3.5d0, 0.0d0, & 
2.5d0, 0.0d0, 3.5d0, 0.0d0, 6.0d0, 3.5d0, 0.0d0, 0.0d0, 1.5d0, 0.0d0, 3.5d0, & 
2.5d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 3.5d0, 5.0d0, 0.0d0, 0.0d0, 1.5d0, 0.0d0, & 
0.0d0, 1.5d0, 1.5d0, 0.0d0, 0.0d0, 2.5d0, 0.0d0, 0.0d0, 1.5d0, 1.5d0, 0.0d0, & 
0.0d0, 4.5d0, 0.0d0, 0.0d0, 1.5d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 3.5d0, 0.0d0, & 
0.0d0, 1.5d0, 1.5d0, 0.0d0, 0.0d0, 0.0d0, 4.5d0, 0.0d0, 4.5d0, 0.0d0, 2.5d0, &
1.5d0, 0.0d0, 0.0d0, 4.5d0, 0.0d0, 0.5d0, 0.0d0, 2.5d0, 0.0d0, 0.0d0, 0.0d0, &
4.5d0, 0.0d0, 0.0d0, 2.5d0, 0.0d0, 2.5d0, 0.0d0, 0.0d0, 4.5d0, 0.0d0, 0.0d0, &
2.5d0, 0.0d0, 2.5d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0, 2.5d0, 0.0d0, 0.0d0, &
0.0d0, 0.5d0, 0.5d0, 0.0d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0, &
4.5d0, 4.5d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 0.5d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0, &
0.0d0, 2.5d0, 3.5d0, 3.5d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0, &
0.0d0, 2.5d0, 3.5d0, 0.0d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 1.5d0, 0.0d0, 0.0d0, &
0.0d0, 3.5d0, 4.0d0, 3.5d0, 3.5d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 1.5d0, 0.0d0, &
1.5d0, 0.0d0, 3.5d0, 5.0d0, 3.5d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 2.5d0, 0.0d0, &
3.5d0, 0.0d0, 3.5d0, 0.0d0, 0.0d0, 0.0d0, 3.5d0, 0.0d0, 3.5d0, 0.0d0, 3.5d0, &
0.0d0, 2.5d0, 0.0d0, 0.0d0, 2.5d0, 3.0d0, 2.5d0, 3.0d0, 2.5d0, 0.0d0, 0.0d0, &
1.5d0, 0.0d0, 1.5d0, 0.0d0, 0.0d0, 1.5d0, 1.5d0, 3.0d0, 0.0d0, 0.0d0, 0.0d0, &
2.5d0, 0.0d0, 2.5d0, 0.0d0, 3.5d0, 0.0d0, 0.0d0, 0.0d0, 3.5d0, 0.0d0, 0.0d0, &
0.5d0, 0.5d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 2.5d0, 0.0d0, 0.0d0, 3.5d0, 1.0d0, &
3.5d0, 7.0d0, 0.0d0, 0.0d0, 3.5d0, 0.0d0, 4.5d0, 0.0d0, 9.0d0, 3.5d0, 0.0d0, &
0.0d0, 0.5d0, 0.0d0, 0.0d0, 2.5d0, 2.5d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 1.5d0, &
0.0d0, 0.0d0, 1.5d0, 1.5d0, 0.0d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0, 1.5d0, &
0.0d0, 0.0d0, 0.5d0, 0.0d0, 1.5d0, 0.0d0, 0.0d0, 0.5d0, 2.0d0, 0.5d0, 0.0d0, &
0.0d0, 0.5d0, 0.0d0, 4.5d0, 4.5d0, 0.5d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0,    &
1.5d0, 2.5d0, 0.0d0, 1.5d0, 0.0d0, 3.5d0, 0.0d0, 2.5d0, 0.5d0, 2.5d0, -1.0d0,   &
-1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
-1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, &
-1.0d0, -1.0d0 /
data (gnI(isotope),isotope=1,lisotope) / &
  5.58569468d0,  0.85743820d0,  5.95799369d0, -4.25499544d0,  0.00000000d0,  & 
  0.82204730d0,  2.17095100d0, -0.78495000d0,  0.60021500d0,  1.79243260d0,  & 
  0.00000000d0,  1.40482360d0,  0.00000000d0,  0.40376100d0, -0.56637768d0,  & 
  0.00000000d0, -0.75751600d0,  0.00000000d0,  5.25773600d0,  0.00000000d0,  & 
 -0.44119800d0,  0.00000000d0,  0.58200000d0,  1.47834800d0,  0.00000000d0,  & 
 -0.34218000d0,  0.00000000d0,  1.45660280d0,  0.00000000d0, -1.11058000d0,  & 
  0.00000000d0,  2.26320000d0,  0.00000000d0,  0.42921400d0,  0.00000000d0,  & 
  0.00000000d0,  0.54791620d0,  0.64273500d0,  0.45608240d0,  0.00000000d0,  & 
  0.00000000d0, -0.45370000d0,  0.00000000d0,  0.26098000d0, -0.32452500d0,  & 
  0.14324670d0,  0.00000000d0, -0.45565170d0,  0.00000000d0, -0.37637000d0,  & 
  0.00000000d0,  0.00000000d0,  0.00000000d0,  1.35899000d0,  0.00000000d0,  & 
 -0.31539000d0,  0.00000000d0, -0.31547700d0,  0.00000000d0,  0.55761480d0,  & 
  1.47106000d0,  0.00000000d0,  0.00000000d0, -0.31636000d0,  0.00000000d0,  & 
  1.43900000d0,  1.38130000d0,  0.00000000d0,  0.00000000d0,  0.18090000d0,  & 
  0.00000000d0,  1.32200000d0,  0.75980000d0,  0.00000000d0,  0.00000000d0,  & 
 -0.50001000d0,  0.00000000d0,  0.00000000d0,  1.48240000d0,  1.58780000d0,  & 
  0.00000000d0,  0.00000000d0,  0.35019200d0,  0.00000000d0,  0.00000000d0,  & 
  1.34439000d0,  1.70818000d0,  0.00000000d0,  0.00000000d0, -0.19543730d0,  & 
  0.00000000d0,  0.00000000d0,  0.95965000d0,  0.00000000d0,  0.00000000d0,  & 
  1.07008000d0,  0.00000000d0, -0.29000000d0,  0.00000000d0,  0.00000000d0,  & 
  1.40426700d0,  1.51370800d0,  0.00000000d0,  0.00000000d0,  0.00000000d0,  & 
 -0.21570400d0,  0.00000000d0, -0.22330000d0,  0.00000000d0,  0.54119200d0,  & 
  1.83421000d0,  0.00000000d0,  0.00000000d0, -0.24284000d0,  0.00000000d0,  & 
 -0.27483080d0,  0.00000000d0, -0.52144800d0,  0.00000000d0,  0.00000000d0,  & 
  0.00000000d0,  1.37120000d0,  0.00000000d0,  0.00000000d0, -0.36570000d0,  & 
  0.00000000d0, -0.37340000d0,  0.00000000d0,  0.00000000d0,  1.26320000d0,  & 
  0.00000000d0,  0.00000000d0, -0.25600000d0,  0.00000000d0, -0.28800000d0,  & 
  0.00000000d0,  0.00000000d0, -0.17680000d0,  0.00000000d0,  0.00000000d0,  & 
 -0.25700000d0,  0.00000000d0,  0.00000000d0,  0.00000000d0, -0.22714000d0,  & 
 -0.26112000d0,  0.00000000d0,  0.00000000d0,  0.00000000d0, -1.18977000d0,  & 
  0.00000000d0, -1.24460200d0,  0.00000000d0,  0.00000000d0,  1.22860000d0,  & 
  1.23130000d0,  0.00000000d0,  0.00000000d0, -1.83770000d0,  0.00000000d0,  & 
 -2.00208000d0,  0.00000000d0, -2.09456000d0,  0.00000000d0,  0.00000000d0,  & 
  0.00000000d0,  1.34540000d0,  0.72851000d0,  0.75100000d0,  0.00000000d0,  & 
  0.00000000d0, -1.47389600d0,  0.00000000d0, -1.77701020d0,  0.00000000d0,  & 
  0.00000000d0,  0.00000000d0,  1.12531000d0,  0.74886000d0,  0.00000000d0,  & 
  0.00000000d0,  0.00000000d0, -1.55595000d0,  0.00000000d0,  0.46100000d0,  & 
  0.00000000d0,  0.00000000d0,  0.00000000d0,  0.73772140d0,  0.74843000d0,  & 
  0.78069000d0,  0.81466000d0,  0.00000000d0,  0.00000000d0, -1.54330000d0,  & 
  0.00000000d0,  0.55863000d0,  0.00000000d0,  0.62491000d0,  0.00000000d0,  & 
  0.77140000d0,  0.74272900d0,  0.79515600d0,  0.00000000d0,  0.00000000d0,  & 
  0.00000000d0,  0.00000000d0,  1.71020000d0,  0.00000000d0, -0.30430000d0,  & 
  0.00000000d0, -0.18700000d0,  0.00000000d0,  0.00000000d0,  0.00000000d0,  & 
  0.73700000d0,  0.00000000d0, -0.23200000d0,  0.00000000d0, -0.19080000d0,  & 
  0.00000000d0,  0.14440000d0,  0.00000000d0,  0.00000000d0,  1.38870000d0,  & 
 -0.64670000d0,  0.61340000d0, -0.66830000d0,  0.60800000d0,  0.00000000d0,  & 
  0.00000000d0, -0.17150000d0,  0.00000000d0, -0.22650000d0,  0.00000000d0,  & 
  0.00000000d0,  1.34000000d0,  1.34300000d0,  0.59670000d0,  0.00000000d0,  & 
  0.00000000d0,  0.00000000d0, -0.19200000d0,  0.00000000d0,  0.26900000d0,  & 
  0.00000000d0,  1.66800000d0,  0.00000000d0,  0.00000000d0,  0.00000000d0,  & 
 -0.16110000d0,  0.00000000d0,  0.00000000d0, -0.46200000d0, -0.45600000d0,  & 
  0.00000000d0,  0.00000000d0,  0.98734000d0,  0.00000000d0, -0.25920000d0,  & 
  0.00000000d0,  0.00000000d0,  0.65170000d0,  1.98800000d0,  0.63780000d0,  & 
  0.45170000d0,  0.00000000d0,  0.00000000d0,  0.22670000d0,  0.00000000d0,  & 
 -0.14240000d0,  0.00000000d0,  0.53610000d0,  0.67729000d0,  0.00000000d0,  & 
  0.00000000d0,  0.23556950d0,  0.00000000d0,  0.00000000d0,  1.27480000d0,  & 
  1.28790000d0,  0.00000000d0,  0.00000000d0,  0.12930380d0,  0.00000000d0,  & 
  0.43995600d0,  0.00000000d0,  0.00000000d0,  0.10050000d0,  0.10910000d0,  & 
  0.00000000d0,  0.00000000d0,  0.00000000d0,  1.21900000d0,  0.00000000d0,  & 
  0.00000000d0,  0.09716400d0,  0.00000000d0,  0.00000000d0,  1.01177100d0,  & 
  0.00000000d0, -0.37348400d0,  0.00000000d0,  0.00000000d0,  3.24451574d0,  & 
  0.04500000d0,  3.27642920d0,  0.00000000d0,  0.00000000d0,  1.18512000d0,  & 
  0.00000000d0,  0.90920000d0,  0.91340000d0,  1.50000000d0,  0.00000000d0,  & 
  0.00000000d0,  0.00000000d0,  0.00000000d0,  0.73000000d0,  0.18000000d0,  & 
  0.00000000d0,  1.34000000d0,  0.00000000d0, -0.10900000d0,  0.00000000d0,  & 
  1.25600000d0,  0.40600000d0,  0.60000000d0,  0.00000000d0,  0.00000000d0,  & 
  0.00000000d0,  0.00000000d0,  0.00000000d0,  0.00000000d0,  0.00000000d0,  & 
  0.00000000d0,  0.00000000d0,  0.00000000d0,  0.00000000d0,  0.00000000d0,  & 
  0.00000000d0,  0.00000000d0,  0.00000000d0,  0.00000000d0,  0.00000000d0,  & 
  0.00000000d0,  0.00000000d0,  0.00000000d0,  0.00000000d0,  0.00000000d0,  & 
  0.00000000d0 /

! Souce: CODATA 2022 (https://physics.nist.gov/) SI units
eCharge=1.602176634d-19
pMass=1.67262192595d-27

Pi=4.d0*ATAN(1.d0)

gnuc=1.d40
Inuc=1.d40
do isotope=1,lisotope
   if ( mass == Amass(isotope) .and. symbol == Asymbol(isotope) ) then
      gnuc=gnI(isotope) * eCharge / 4.d0 / Pi / pMass / 1.D6
      Inuc=spin(isotope)
      return
   end if
enddo
if (Inuc == 1.d40) call ERRMSG(01," SUBROUTINE SpinGammaI: unknown isotope!")
END SUBROUTINE SpinGammaI

!*******************************************************************************
! Input: atomic mass and atomic symbol 
! Output: Nuclear eletric quadrupole (mbarn, 10^-31 m^2)
! Data from:
! Pekka Pyykkö (2018) Year-2017 nuclear quadrupole moments, Molecular Physics, 
! 116:10, 1328-1338, DOI: 10.1080/00268976.2018.1426131
!*******************************************************************************
SUBROUTINE EQuadrupol(mass,symbol,Qnuc)
Implicit None

integer, parameter :: nEQuad=111
integer            :: Amass(nEQuad), mass, isotope
character(len=2)   :: Asymbol(nEQuad), symbol
real(8)            :: EQuad(nEQuad), Qnuc

data (Amass(isotope),isotope=1,nEQuad) / &
  2,   6,   7,   9,  10,  11,  11,  14,  17,  19,  21,  23,  25,  27,  33,  35,  35,  37,  39,  40, &
 41,  41,  43,  45,  47,  49,  50,  51,  53,  55,  57,  59,  61,  63,  65,  67,  69,  71,  73,  75, &
 77,  79,  81,  83,  85,  87,  87,  90,  91,  93,  95,  97,  99,  99, 101, 100, 105, 113, 115, 119, &
121, 123, 127, 131, 133, 135, 137, 138, 139, 141, 143, 145, 147, 147, 149, 151, 153, 155, 157, 159, &
161, 163, 165, 167, 173, 175, 176, 177, 179, 181, 185, 187, 189, 191, 193, 197, 201, 209, 209, 209, &
211, 227, 229, 231, 233, 235, 237, 241, 241, 243, 253 /
data (Asymbol(isotope),isotope=1,nEQuad) / &
"H ", "Li", "Li", "Be", "B ", "B ", "C ", "N ", "O ", "F ", "Ne", "Na", "Mg", "Al", "S ", &
"S ", "Cl", "Cl", "K ", "K ", "K ", "Ca", "Ca", "Sc", "Ti", "Ti", "V ", "V ", "Cr", "Mn", &
"Fe", "Co", "Ni", "Cu", "Cu", "Zn", "Ga", "Ga", "Ge", "As", "Se", "Br", "Br", "Kr", "Rb", &
"Rb", "Sr", "Y ", "Zr", "Nb", "Mo", "Mo", "Tc", "Ru", "Ru", "Rh", "Pd", "In", "In", "Sn", &
"Sb", "Sb", "I ", "Xe", "Cs", "Ba", "Ba", "La", "La", "Pr", "Nd", "Nd", "Pm", "Sm", "Sm", &
"Eu", "Eu", "Gd", "Gd", "Tb", "Dy", "Dy", "Ho", "Er", "Yb", "Lu", "Lu", "Hf", "Hf", "Ta", &
"Re", "Re", "Os", "Ir", "Ir", "Au", "Hg", "Pb", "Bi", "Rn", "Fr", "Ac", "Th", "Pa", "U ", &
"U ", "Np", "Pu", "Am", "Am", "Es" / 
data (EQuad(isotope),isotope=1,nEQuad) / &
2.85783d0, -0.808d0, -40.1d0, 52.88d0, 84.59d0, 40.59d0, 33.27d0, 20.44d0, -25.58d0, -94.2d0, &
101.55d0, 104d0, 199.4d0, 146.6d0, -69.4d0, 48.3d0, -81.12d0, -63.93d0, 60.3d0, -75.0d0, &
73.4d0, -66.5d0, -40.8d0, -220d0, 302d0, 247d0, 210d0, -52d0, -150d0, 330d0, &
160d0, 420d0, 162d0, -220d0, -204d0, 122d0, 171d0, 107d0, -196d0, 311d0, &
760d0, 308.7d0, 257.9d0, +259d0, 276d0, 133.5d0, 305d0, -125d0, -176d0, -320d0, &
-22d0, 255d0, -129d0, 79d0, 457d0, 153d0, 660d0, 761d0, 772d0, -132d0, &
-543d0, -692d0, -688.22d0, -114.6d0, -3.43d0, 153d0, 236d0, 450d0, 206d0, -58.9d0, &
-630d0, -330d0, +740d0, -259d0, +75d0, 903d0, 2412d0, 1270d0, 1350d0, 1432d0, &
2507d0, 2648d0, 3580d0,  3565d0, 2800d0, 3490d0, 4970d0, 3365d0, 3793d0, 3170d0, & 
2180d0, 2070d0, 856d0, 816d0, 751d0, 547d0, 387d0, -269d0, -516d0, 311d0, &
-210d0, 1700d0, 3110d0, -1720d0, 3663d0,  4936d0, 3886d0, 5600d0, 4200d0, 4200d0, &
6700d0 /
Qnuc=-1.d0
do isotope=1,nEQuad
   if ( mass == Amass(isotope) .and. symbol == Asymbol(isotope) ) then
      Qnuc=EQuad(isotope)
      return
   end if
enddo
if (Qnuc == -1.d0) call ERRMSG(01," SUBROUTINE EQuadrupol: unknown isotope!")
END SUBROUTINE EQuadrupol

!*******************************************************************************
SUBROUTINE getIsotope(line,linelen,mass,asymbol)
!*******************************************************************************
Implicit None

character :: line*(*), asymbol*2, number*1
integer   :: char, Masslen, num, linelen, mass

Masslen=0
do char=1,linelen
  do num=0,9
     write(number,'(I1)') num
     if ( line(char:char) == number ) Masslen=Masslen+1
  end do
enddo
if ( linelen-Masslen == 1) then
   asymbol=line(Masslen+1:Masslen+1)//' '
else
   asymbol=line(Masslen+1:Masslen+2)//' '
endif
mass=ICHAR(line(:Masslen))
if (Masslen==1) read(line(:Masslen),'(I1)') mass
if (Masslen==2) read(line(:Masslen),'(I2)') mass
if (Masslen==3) read(line(:Masslen),'(I3)') mass
END SUBROUTINE getIsotope

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! F77 Codes extracted from PLACZEK
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ******************************
! REMOVE BLANK SPACES IN STRINGS
! ******************************
      SUBROUTINE KBLANK(TNAM,ILEN,LEN)
      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (MAXLEN=160)
      CHARACTER(len=MAXLEN) NAM
      CHARACTER(len=*) TNAM

      NAM  = ' '
      IDEL = 1
      LEN  = 0

      IF (ILEN.GT.MAXLEN) STOP "Error in KBLANK: Too long string!"
      DO JDEL = 1,ILEN
         IF( TNAM(JDEL:JDEL).NE.' ' ) THEN
             LEN = LEN + 1
             NAM(IDEL:IDEL) = TNAM(JDEL:JDEL)
             IDEL = IDEL + 1
         END IF
      END DO

      TNAM = NAM(1:LEN)

      RETURN
      END SUBROUTINE KBLANK

! *********************
! SUBROUTINE FINDSTRING
! *********************
      SUBROUTINE FINDSTRING (STRING,LEN,NOP,NBACK,FIND)
      IMPLICIT REAL*8 (A-H,O-Z)

      CHARACTER(len=*) STRING
      CHARACTER LINE*80
      LOGICAL   FIND, NBACK

      IF ( LEN > 80 ) THEN
         WRITE(*,'(/A/)') " @ ERROR: Too long STRING in FINDSTRING!"
         STOP
      END IF
      IF ( NBACK ) REWIND(NOP)
 01   READ(NOP,'(A)',END=02,ERR=02) LINE
      IF (LINE(1:LEN) .EQ. STRING(1:LEN)) THEN
          FIND = .TRUE.
          BACKSPACE(NOP)
          RETURN
      ELSE
          GOTO 01
      END IF
 02   FIND = .FALSE.

      RETURN
      END SUBROUTINE FINDSTRING

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Search for comments in a file
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE COMMENT(UNIT)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER UNIT

      CHARACTER LINE*80

 01   READ(UNIT,'(A)',END=02) LINE
!      IF (LINE(1:1) .EQ. '#') GO TO 01
      IF (LINE(1:1) .EQ. '!') GO TO 01
      BACKSPACE(UNIT)
      RETURN
 02   CONTINUE
      END SUBROUTINE COMMENT

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C   Date and Time / Hostname
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE LOCALINF

CHARACTER*(24) FFDATE
CHARACTER*40   HSTNAM

CALL FDATE(FFDATE)
CALL HOSTNM(HSTNAM)
WRITE(*,'(//1X,2A/1X,2A/)') &
'Host name     : ',HSTNAM,  &
'Date and time : ',FFDATE
END SUBROUTINE LOCALINF
