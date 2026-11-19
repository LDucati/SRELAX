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
real(8), parameter :: hbar=1.054571817D-34, hartree=4.3597447222060D-18, a0=5.29177210544D-11, &
                      kBtz=1.380649d-23, Navog=6.02214076d23, E_Hartree=4.3597447222060d-18
!vacuum magnetic permeability
real*8, parameter :: mu0 =  1.25663706127D-6 ! N/A^2 (CODATA)

! Static allocation
real(8)            :: timestep, splstep, limtime, Ispin, Qmom, B0ext, gammaI, &
                      sigmaEFG(5), GlmRe(5), GlmENR(5), SDmin, SDmax, MDdt, Rmax, Lbox, DEcoef, Temper, etaCoef, &
                      SaveR1(3), SaveR2(3), SaveTau, TauInf, SRtensor(3)
integer            :: nQNUC, nframes, splframes, natom, nacfs, MolAtm, nspin, maxmcom=5, &
                      MDstartProd, MDendProd, MDsampleFreq, inttype, fft_size, MolFrame(4), SR_nmol
integer, parameter :: DYNPYU=10, acfflen=160, lMDfile=120
logical            :: ORCA, DYNPY, ACF_ENSAV, LTIME, RACF, ACF_AVG, PrintSD, GMX, EQRelax, & 
                      MSD_COM, CP2K, DRELAX, NoInterACF, NucProp, UNWRAPMD, HwangFreed, SaveR1R2, &
                      WK_divbyN, RDDACF, TINKER, SR_RELAX, Read_Vel, AllInterPairs
character          :: trajdir*50, symbol*3, Asymbol*2, acfisofn*20='acf.iso'
character(len=lMDfile) :: MDtrajfile, MDvelfile

! vector to label componentes of the 2nd-rank 
! spherical tensor R_2m: R_20 R_2+1 R_2-1 R_2+2 R_2-2
integer, dimension(5) :: mvec=(/0,1,-1,2,-2/)

! Dynamic allocation
complex(8), dimension(:,:), allocatable :: T2m
real(8),    dimension(:),   allocatable :: Xspl, Fspl, newX, newF, AtomMass
real(8),    dimension(:,:), allocatable :: Txyz, AtmPos, ACFre, ACFim, AtmVel
integer,    dimension(:),   allocatable :: atmlist
character(len=acfflen),  dimension(:),   allocatable :: acffiles

END MODULE

! ***************************************************************************
! Main code: Computes T1 and T2 relaxation times of the following mechanisms: 
!            (i) Quadrupolar, 
!            (ii) Dipolar of identical spins
! ***************************************************************************
PROGRAM Srelax
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
!$         print '(" SRELAX - Spin Relaxation program")'
!$         print '(" Parallel execution with number of threads: ",i3)', nthreads
!$      end if
!$omp end parallel

write(*,'(A/A/A//A/A/A/A/A/)')  &
" ***************************** ", &
" OUTPUT OF PROGRAM SPIN_RELAX  ", &
" ***************************** ", &
" This program computes relaxation times and ", &
" related quantities for the following mechanisms:", &
" (i)   Dipolar (equal nuclei), ", &
" (ii)  Quadrupolar, ", &
" (iii) Spin-Rotation."

! Reading input file "srelax.inp"
call read_input()

! ************************
! Spin-Rotation Relaxation
! ************************
if (SR_RELAX) then

   ! SR relaxation involves a rank-1 tensor, i.e., 3 components (here, x,y,z)
   maxmcom = 3

   ! Reading MD trajectory file
   if      (GMX) then
      call ReadGRO()
   else if (CP2K) then
      call ReadCP2K()
   else if (TINKER) then
      call ReadTINKER()
   else
      call ERRMSG(0,"The trajectory file was not defined in 'srelax.inp' via the '.CP2K'/'.GROMACS'/'.Tinker' options.")
   endif
   write(*,'(A/A/A,20I4)') &
   " ------------------------------------------------------------------------------------- ", &
   " Spin-Rotation relaxation will be calculated for the following atoms in each molecule: ", &
   " atom type = "//Asymbol(:2)//", atom number = ",(atmlist(alist),alist=1,nspin)
   write(*,'(A)') &
   " ------------------------------------------------------------------------------------- "
   write(*,'(A,F15.5,A,I7,A)') " Start sampling trajectories at time (ps) : ", &
        (MDstartProd-1)*MDdt," (",MDstartProd,"-th frame)"
   write(*,'(A,F15.5,A,I7,A)') " Stop  sampling trajectories at time (ps) : ", &
        (MDendProd-1)*MDdt," (",MDendProd,"-th frame)"
   write(*,'(A,F15.5)') " Trajectory duration (ps)                 : ",(nframes-1)*timestep
   write(*,'(A,F15.5)') " Trajectory step (ps)                     : ",timestep
   write(*,'(A,I15)')   " Number of frames of the current job      : ",nframes
   write(*,'(A,I4,A)')  " Sampling frequency                       : (every ",MDsampleFreq," frames)"

! -------------------------------------------------------------
! Moving the CoM of each molecule, in each frame, to the origin
! -------------------------------------------------------------
   write(*,"(/A/)") " Moving the CoM of each molecule, in each frame, to the origin."
   call molLab2CoM()

! Reading the linear velocity of atoms 
   if (Read_Vel) then
!      if (GMX)    call ReadGROvel()
!      if (CP2K)   call ReadCP2Kvel()
      if (TINKER) then
         call ReadTINKERvel()
         write(*,"(/A/)") " Removing the CoM velocity from the atomic velocities "
         call molLab2CoMvel()
      else 
         call ERRMSG(0,"The current implementation only allows for reading atomic velocities generated by Tinker")
      endif
   else
      ! Estimates atomic linear velocities in the molecular CoM frame 
      ! previously established by molLab2CoM routine
      call velCoMframe()
   endif

! ----------------------------------------------------------------------------
! Determining the molecular angular momentum from the molecular inertia tensor 
! and angular velocity. The last one is obtained using equation 13 from the 
! reference: https://doi.org/10.1021/acs.jpcb.2c06411
!            J. Phys. Chem. B 2023, 127, 1167−1177
! 
! Computes also the corresponding ACF, the spectral densities e relaxation 
! properties
! ----------------------------------------------------------------------------
   call J_MolFrame()

! Date and Time
   call LOCALINF
   write(*,'(A/A/A/)')  &
" *********************** ",&
" End of all calculations ",&
" *********************** "
   stop
endif

! ******************
! Dipolar Relaxation 
! ******************
if (DRELAX) then
! Dipolar relaxation involves a spherical function with m=0,1,2,-1,-2
! Here, only the first 3 of them are calculated
   mvec    = (/0,1,2,0,0/)
   maxmcom = 3

! Reading MD trajectory file 
   if      (GMX) then
      call ReadGRO()
   else if (CP2K) then
      call ReadCP2K()
   else if (TINKER) then
      call ReadTINKER()
   else
      call ERRMSG(0,"The trajectory file was not defined in 'srelax.inp' via the '.CP2K'/'.GROMACS'/'.Tinker' options.")
   endif

! Computes the Mean Squared Displacement for the center of mass (COM)
   if (MSD_COM) then
      call MSDcom()
      write(*,"(/A/)") " MSD(t) of the center of mass written to file 'msd_com.dat' "
      stop
   endif

! Interacting nuclei and time information
   write(*,'(A/A/A,20I4)') &
   " ------------------------------------------------------------------------------- ", &
   " Dipolar relaxation will be calculated for the following atoms in each molecule: ", &
   " atom type = "//Asymbol(:2)//", atom number = ",(atmlist(alist),alist=1,nspin)
   write(*,'(A)') &
   " ------------------------------------------------------------------------------- "
   write(*,'(A,F15.5,A,I7,A)') " Start sampling trajectories at time (ps) : ", &
        (MDstartProd-1)*MDdt," (",MDstartProd,"-th frame)"
   write(*,'(A,F15.5,A,I7,A)') " Stop  sampling trajectories at time (ps) : ", &
        (MDendProd-1)*MDdt," (",MDendProd,"-th frame)"
   write(*,'(A,F15.5)') " Trajectory duration (ps)                 : ",(nframes-1)*timestep
   write(*,'(A,F15.5)') " Trajectory step (ps)                     : ",timestep
   write(*,'(A,I15)')   " Number of frames of the current job      : ",nframes
   write(*,'(A,I4,A)')  " Sampling frequency                       : (every ",MDsampleFreq," frames)"

! Autocorrelation functions for the F^IJ_m auxiliary functions
! from reference Phys. Chem. Chem. Phys., 2019, 21, 26621 (DOI: 10.1039/c9cp04976b)
   write(*,'(/A/)') &
   " Computing F^ij_m(t) of the dipolar interaction and the corresponding autocorrelation functions"
   call FIJmPCCP()

! Freeing up memory
   if (allocated(AtmPos)) then
      deallocate(AtmPos, STAT=ios)
      call AllocStat(ios,'AtmPos','SRelax','D')
   endif

! Intermolecular dipolar relaxation according to J. Chem. Phys. 160, 074102 (2024); doi: 10.1063/5.0191052
   if (.not.NoInterACF .and. HwangFreed) then
      call HF_ACF()
   endif

! Date and Time
   call LOCALINF
   write(*,'(A/A/A/)')  &
" *********************** ",&
" End of all calculations ",&
" *********************** "
   stop
endif

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
   call AllocStat(ios,'Txyz','Srelax','A')
   call AllocStat(ios,'T2m','Srelax','A')

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
   if (allocated(Txyz)) then
      deallocate(Txyz, STAT=ios)
      call AllocStat(ios,'Txyz','Srelax','D')
   endif
   
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
   call AllocStat(ios,'ACFre','Srelax','A')
   allocate(ACFim(nQNUC*5,nframes), STAT=ios)
   call AllocStat(ios,'ACFim','Srelax','A')
   
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
   if (allocated(T2m)) then
      deallocate(T2m, STAT=ios)
      call AllocStat(ios,'T2m','Srelax','D')
   endif

endif

! Spectral Densities of the EFG-ACFs 
call SpectDens()

! deallocating ACF_2m
if (allocated(ACFre)) then
   deallocate(ACFre, STAT=ios)
   call AllocStat(ios,'ACFre','Srelax','D')
endif
if (allocated(ACFim)) then
   deallocate(ACFim, STAT=ios)
   call AllocStat(ios,'ACFim','Srelax','D')
endif

! Quadrupolar relaxation properties
call QrelaxProp()

! Date and Time
call LOCALINF
write(*,'(A/A/A/)')  &
" *********************** ",&
" End of all calculations ",&
" *********************** "
END PROGRAM Srelax

! ******************************************************************************
! Routine to read the input file 'srelax.inp'
! ******************************************************************************
SUBROUTINE read_input()
Use shared_vars
Implicit None

integer   :: ios, loop, leniso, Amass, strgsz, Zatm, curr_len, char_pos
integer, parameter      :: line_len=512
character(len=line_len) :: line, emsg
character :: IspinSymbol*2, atomsb*2
logical   :: Isotp=.FALSE., in_word=.false.
real(8)    :: ProtonFq=0.d0, Hspin, gammaH

! Initializing some input vars
nQNUC=0; splframes=0; nframes=0; MDstartProd=0; MDendProd=0; MDsampleFreq=0
limtime=0.d0; Ispin=0.d0; Qmom=0.d0; timestep=0.d0; splstep=0.d0; B0ext=0.d0; gammaI=0.d0; 
MDdt=0.d0; MDtrajfile=""; MDvelfile=""; MolAtm=0; Rmax=1.d0; Lbox=0.d0
ORCA=.FALSE.; DYNPY=.FALSE.; ACF_ENSAV=.FALSE.; LTIME=.FALSE.; ACF_AVG=.FALSE.; PrintSD=.FALSE.; GMX=.FALSE.; EQRelax=.FALSE.
MSD_COM=.FALSE.; CP2K=.FALSE.; DRELAX=.FALSE.; NoInterACF=.FALSE.; UNWRAPMD=.FALSE.; DEcoef=0.d0; HwangFreed=.FALSE.;
etaCoef=0.d0; Temper=0.d0; TauInf=0.d0; fft_size=0; SaveR1R2=.FALSE.; WK_divbyN=.FALSE.; RDDACF=.FALSE.;
NucProp=.TRUE.; TINKER=.FALSE.; SR_RELAX=.FALSE.; SRtensor=0.d0; MolFrame=0; SR_nmol=0; Read_Vel=.FALSE.
AllInterPairs=.FALSE.; RACF=.FALSE.; nacfs=0; Asymbol="  "; IspinSymbol="  " 

! Opening input file 'srelax.inp'
open(01,file='srelax.inp',status='old',iostat=ios)

! Checking if file exists
 if (ios /= 0) call ERRMSG(01,"Could not open input file 'srelax.inp")

! Endless loop to read input options
 do 
    read(01,'(A)',end=01) line

! Reading GROMACS trajectory file (.gro)
!
! MDstartProd:  SRELAX starts sampling snapshots at the "MDstartProd"-th frame contained in the trajectory file
! MDendProd:    SRELAX stops  sampling snapshots at the "MDendProd"-th   frame contained in the trajectory file
! MDsampleFreq: Number of frames that elapse between reading coordinates from the trajectory file (after reading a frame, SRELAX 
!               will skip the next "MDsampleFreq-1" frames before reading a frame again
! MDdt:         Timestep in picoseconds between consecutive frames in the trajectory file
! After ".GROMACS", you must provide the following information:
! MDstartProd, MDendProd, MDsampleFreq, MDdt, MDtrajfile
    if ( line(1:8) == '.GROMACS' ) then 
       MDtrajfile = ''
       GMX=.TRUE.
       call COMMENT(01)
! dt nstxout start_prod end_prod sample_freq filename
       read(01,*) MDstartProd, MDendProd, MDsampleFreq, MDdt, MDtrajfile
! If the MD filename has slash bar '/'
       if ( MDtrajfile.eq.'' ) then
          backspace(01)
          read(01,'(A)') line 
          MDtrajfile = line(index(line,'/'):)
       endif 
       nframes   = (MDendProd-MDstartProd)/MDsampleFreq+1
       timestep  = MDdt*MDsampleFreq
    endif

! Reading atomic velocities from file 'MDvelfile'
    if ( line(1:11) == '.Velocities' ) then 
       Read_Vel=.TRUE.
       call COMMENT(01)
       read(01,*) MDvelfile
       if ( MDvelfile.eq.'' ) then
          backspace(01)
          read(01,'(A)') line
          MDvelfile = line(index(line,'/'):)
       endif
    endif

! Reading TINKER trajectory file (.arc)
! MDstartProd:  SRELAX starts sampling snapshots at the "MDstartProd"-th frame contained in the trajectory file
! MDendProd:    SRELAX stops  sampling snapshots at the "MDendProd"-th   frame contained in the trajectory file
! MDsampleFreq: Number of frames that elapse between reading coordinates from the trajectory file (after reading a frame, SRELAX
!               will skip the next "MDsampleFreq-1" frames before reading a frame again
! MDdt:         Timestep in picoseconds between consecutive frames
! After ".TINKER", you must provide the following information:
! MDdt, MDstartProd, MDendProd, MDsampleFreq, MD_trajectory_file_name
    if ( line(1:8) == '.Tinker' ) then
       MDtrajfile = ''
       TINKER=.TRUE.
       call COMMENT(01)
       read(01,*) MDstartProd, MDendProd, MDsampleFreq, MDdt, MDtrajfile
! If the MD filename has slash bar '/'
       if ( MDtrajfile.eq.'' ) then
          backspace(01)
          read(01,'(A)') line
          MDtrajfile = line(index(line,'/'):)
       endif
       nframes   = (MDendProd-MDstartProd)/MDsampleFreq+1
       timestep  = MDdt*MDsampleFreq
    endif

! Number of molecules to sample when computing angular velocities, correlation functions, and SR rate. 
! Can be up to total number of molecules in the system but may be reduced for performance.
    if (line(:14).eq.'.SR_nmol') then
       call COMMENT(01)
       read(01,*) SR_nmol
       write(*,'(A/A,I5/)') &
" Number of molecules to sample when computing angular velocities, ", &
" correlation functions, and SR rate: ",SR_nmol
    endif

! Reading the Spin-Rotation Tensor (principal components) in kHz
    if (line(:14).eq.'.SpinRotTensor') then
       SR_RELAX=.TRUE.
       call COMMENT(01)
       read(01,*) (SRtensor(loop),loop=1,3)
       write(*,'(A,3F12.5/)') " Principal components of the SR tensor (kHz): ",(SRtensor(loop),loop=1,3)
    endif

! The SR mechanism requires the definition of a coordinate system fixed to the molecule. 
!
! The first line after '.MolFrame' specifies which atoms are to be used to define the z-axis. 
! In this example (CH4 molecule), the z-axis will be parallel to the vector: 
! vec{k} = vec{H1} - vec{C}, where vec{H1} and vec{C} represent the positions of 
! hydrogen 1 and carbon, respectively.
! The second line defines the vector used to establish the y-axis; in this case, 
! vec{j} = vec{k} x (vec{H2} -vec{C}). 
! Finally, the x-axis is defined as: vec{i} = vec{k} x vec{j}.
    if (line(:9).eq.'.MolFrame') then
       call COMMENT(01)
       read(01,*) (MolFrame(loop),loop=1,2)
       call COMMENT(01)
       read(01,*) (MolFrame(loop),loop=3,4)
       write(*,'(A/A/A,I3,A,I3,A/A,I3,A,I3,A/A/)') &
" Definition of the coordinate system fixed to the molecule:", &
" --------------------------------------------------------- ", &
" z-axis: vec{k} propto            vec{atom ",MolFrame(1),"} -vec{atom ",MolFrame(2),"}", &
" y-axis: vec{j} propto  vec{k} x (vec{atom ",MolFrame(3),"} -vec{atom ",MolFrame(4),"})", &
" x-axis: vec{i} =       vec{k} x  vec{j}"
    endif


! The autocorrelation functions will be computed by dividing all terms by N, 
! the number of frames.
    if (line(:9).eq.'WK_divbyN') then
       write(*,"(A/)") &
" Note: the autocorrelation functions will be computed by dividing all terms by N, the number of frames."
       WK_divbyN=.TRUE.
    endif

! Unwrap the MD simulation box
    if ( line(1:8) == 'UNWRAPMD' ) then
       write(*,'(A/)') " Unwrapping the molecular dynamics trajectory=.T."
       UNWRAPMD=.TRUE.
    endif

! Correlation time (ps) of the infinite-size simulation
! obtained via the empirical function: tau_i = tau_inf + b*L^-1
    if ( line(1:7) == '.TauInf' ) then 
       call COMMENT(01)
       read(01,*) TauInf
    endif

! Read Dipolar ACFs
    if (line(:11).eq.'READ_DD_ACF') then
       RDDACF=.TRUE.
       write(*,"(A/)") " Dipolar ACFs will be read from avg-acf-intra_wk.dat and avg-acf-inter_wk.dat files. "
    endif

! Reading Self-Diffusion Coeficient (angstrom^2/ps)
    if ( line(1:14) == '.DiffusionCoef' ) then 
       HwangFreed=.TRUE.
       ACF_AVG=.TRUE.
       SaveR1R2=.TRUE.
       call COMMENT(01)
       read(01,*) DEcoef
! Converting DEcoef to nm^2/ps
       DEcoef = DEcoef/100.d0
    endif

! shear viscosity η (in mPa.s) and its Temperature  (K)
    if ( line(1:15) == '.ShearViscosity' ) then
       call COMMENT(01)
       read(01,*) etaCoef, Temper
    endif
! 
! Reading CP2K trajectory file (*-pos-1.xyz)
!
! MDstartProd:  SRELAX starts sampling snapshots at the "MDstartProd"-th frame contained in the trajectory file
! MDendProd:    SRELAX stops  sampling snapshots at the "MDendProd"-th   frame contained in the trajectory file
! MDsampleFreq: Number of frames that elapse between reading coordinates from the trajectory file (after reading a frame, SRELAX
!               will skip the next "MDsampleFreq-1" frames before reading a frame again
! MDdt:         Timestep in picoseconds between consecutive frames in the trajectory file
!
! After ".CP2K", provide the following information:
! MDdt, MDstartProd, MDendProd, MDsampleFreq, MD_trajectory_file_name
    if ( line(1:5) == '.CP2K' ) then 
       MDtrajfile = ''
       CP2K=.TRUE.
       call COMMENT(01)
       read(01,*) MDstartProd, MDendProd, MDsampleFreq, MDdt, MDtrajfile
! If the MD filename has slash bar '/'
       if ( MDtrajfile.eq.'' ) then
          backspace(01)
          read(01,'(A)') line
          MDtrajfile = line(index(line,'/'):)
       endif
       nframes  = (MDendProd-MDstartProd)/MDsampleFreq+1
       timestep = MDdt*MDsampleFreq
    endif

! Enabling the calculation of dipolar relaxation properties
    if ( line(1:9) == 'DD_RELAX' ) then
       DRELAX=.TRUE.
       write(*,'(A/)') &
" This job will compute properties related to the (direct) dipolar interaction mechanism."
    endif

! ***************************************
! Defining the molecule: (.Molecule flag)
!
! A single line witj the atomic symbols in the order they occur within the trajectory file.
! An asterisk after the element symbol indicates that it will contribute to 
! the dipolar/SR relaxation mechanism. The isotope corresponding to this/these 
! atoms is defined using the ".Isotope" option below.
! Notice: Only dipolar relaxation of like nuclei (same isotope) is available 
! in the current implementation.
    if ( line(1:9) == '.Molecule' ) then
       call COMMENT(01)
       read(01,"(A)") line
       
       ! Getting the number of atoms and interacting spins 
       MolAtm   = 0
       nspin    = 0
       char_pos = 0
       in_word  = .false.
       curr_len = len_trim(line)

       ! Count atoms and interacting nuclei
       do loop = 1, curr_len
         if (line(loop:loop) == '*') then
           nspin = nspin + 1
         else if (line(loop:loop) /= ' ') then
           if (.not. in_word) then
             MolAtm = MolAtm + 1
             in_word = .true.
           end if
         else
           in_word = .false.
         end if
       end do

       if (MolAtm == 0) call ERRMSG(01, "No atoms found in '.Molecule' option.")

       ! Safe allocation for AtomMass and atmlist
       if (allocated(AtomMass)) then
          deallocate(AtomMass, STAT=ios)
          call AllocStat(ios,'AtomMass','read_input','D')
       endif
       allocate(AtomMass(MolAtm),STAT=ios)
       call AllocStat(ios,'AtomMass','read_input','A')

       if (allocated(atmlist)) then
          deallocate(atmlist, STAT=ios)
          call AllocStat(ios,'atmlist','read_input','D')
       endif
       allocate(atmlist(nspin),STAT=ios)
       call AllocStat(ios,'atmlist','read_input','A')

       ! Extract and store atom masses 
       MolAtm      = 0
       nspin       = 0
       in_word     = .false.
       atomsb      = "  "
       IspinSymbol = "  " 

       do loop = 1, curr_len
         if (line(loop:loop) == '*') then
            nspin = nspin + 1
            atmlist(nspin) = MolAtm
            
            ! Store IspinSymbol on the first spin, or check consistency on subsequent ones
            if (nspin == 1) then
               IspinSymbol = atomsb
            else if (atomsb /= IspinSymbol) then
               emsg = "Interacting spins (*) in '.Molecule' must belong to the same element: " &
                      // trim(IspinSymbol) // " vs " // trim(atomsb)
               call ERRMSG(01, emsg)
            end if

            ! Check if the asterisk is the very last character of the line
            if (loop == curr_len) then
               call GetAtomZ(atomsb,Zatm)
               call GetAtomM(Zatm,AtomMass(MolAtm))
            end if
            cycle
         endif

         if (line(loop:loop) /= ' ') then
           if (.not. in_word) then
             MolAtm   = MolAtm + 1
             in_word  = .true.
             char_pos = 1
             atomsb   = "  " 
           else
             char_pos = char_pos + 1
           end if

           if (char_pos <= 2) then
              atomsb(char_pos:char_pos) = line(loop:loop)
           end if

           ! Compute Z and Mass when reaching the end of the line
           if (loop == curr_len) then
              call GetAtomZ(atomsb,Zatm)
              call GetAtomM(Zatm,AtomMass(MolAtm))
           end if

         else
           ! Compute Z and Mass once the full atom symbol has been read
           if (in_word) then
              call GetAtomZ(atomsb,Zatm)
              call GetAtomM(Zatm,AtomMass(MolAtm))
           end if
           in_word = .false.
         end if
       end do
    end if

!If spins i and j are at a distance greater than Rmax (in Angstrom), 
! then Fij_m=0
    if ( line(1:8) == '.Rij_max' ) then
       call COMMENT(01)
       read(01,*) Rmax
! conventing Rmax to nm
       Rmax = Rmax/10.d0
    endif

! Simulation box size (Angstrom)
    if ( line(1:5) == '.LBOX' ) then
       call COMMENT(01)
       read(01,*) Lbox
! conventing Lbox to nm
       Lbox = Lbox/10.d0
       write(*,'(A,F10.5,A/)') " Edge length of the (cubic) simulation box: ",Lbox*10.d0," angstrom"
    endif

! Dipolar ACFs for intermolecular interactions will not be calculated
    if ( line(1:10) == 'NoInterACF' ) NoInterACF=.TRUE.

! All itermolecular dipolar interactions will be calculated
    if ( line(1:13) == 'AllInterPairs' ) AllInterPairs=.TRUE.

! Calculates the root-mean-square displacement of the molecule's center of mass, 
! which is saved in the file "msd_com.dat"
    if ( line(1:7) == 'MSD_COM' ) MSD_COM=.TRUE.

! Each component of the spherical EFG-ACF will be replaced by the mean
! value taken from its five components
    if ( line(1:13) == 'ForceIsotropy' ) ACF_AVG=.TRUE.
    
! relaxation properties will be calculated for those nuclei
! contained in the .EquivNuclei list
    if ( line(1:12) == '.EquivNuclei' ) then
       EQRelax=.TRUE.
       call COMMENT(01)
       read(01,*) nQNUC
       if (allocated(atmlist)) then
          deallocate(atmlist, STAT=ios)
          call AllocStat(ios,'atmlist','read_input','D')
       endif
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
       if (allocated(acffiles)) then
          deallocate(acffiles, STAT=ios)
          call AllocStat(ios,'acffiles','read_input','D')
       endif
       allocate(acffiles(nacfs), STAT=ios)
       call AllocStat(ios,'acffiles','read_input','A')
       do loop=1,nacfs
          call COMMENT(01)
          read(01,'(A)') acffiles(loop)
       enddo
       write(*,'(A,I4,A/A//A/A)') &
" This job will read ",nacfs," EFG-ACFs from the files listed below ", &
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
       call getIsotope(line,leniso,Amass,Asymbol)
       call SpinGammaI(Amass,Asymbol,Ispin,gammaI)
       if ( gammaI == 0.d0 ) then
          emsg='Could not get nuclear spin and magnetogyric ratio for the isotope: '//line(:leniso)
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
       write(*,"(A/A,F10.2,A,F8.2,A/)") & 
" Spectral Densities will be calculated for magnetic fields (in T)",&
" within the interval: [", SDmin,",", SDmax,"]"
    endif

 end do
 01 continue

close(01)

if (EQRelax .and. Isotp) then
   call EQuadrupol(Amass,Asymbol,Qmom)
   if ( Qmom == 0.d0 ) then
      emsg='Could not get nuclear eletric quadrupole for the isotope: '//line(:leniso)
      call ERRMSG(01,emsg)
   endif
endif

if ( ProtonFq > 0.d0 ) then
   call SpinGammaI(1,'H ',Hspin,gammaH)
   B0ext = ProtonFq/gammaH
endif

if (SR_RELAX) then
   write(*,'(A/A/)') &
" The Spin-Rotation mechanism was implemented assuming the extreme-narrowing condition ", &
" (rapid motion regime); therefore, B_0 is set to zero."
endif

if ( B0ext == 0.d0 ) then
   write(*,'(A/A/)') &
" => Proton Larmor frequency not found in input file: setting it to zero.", &
"    As the Larmor frequency will be zero, this calculation implies the Extreme Narrowing Regime."
endif

! FFT array size
fft_size=2
do while (fft_size.lt.2*nframes) 
   fft_size=fft_size*2
enddo

if ( .NOT.Isotp ) call ERRMSG(0," Missing input option: '.Isotope' ")
if ( RACF ) return

if (EQRelax .and. nQNUC == 0 ) call ERRMSG(0," Missing input option: '.EquivNuclei' ")
if (nframes == 0 )     call ERRMSG(0," Missing input option: '.frames' ")
if (timestep == 0.d0 ) call ERRMSG(0," Missing input option: '.timestep' ")

! Consistency check between the atoms selected with the "*" flag and the specified isotope.
if ((DRELAX .or. SR_RELAX) .and. IspinSymbol.ne.Asymbol) then
   write(*,'(A)') " @ The selected isotope (option .Isotope) and the atomic symbol marked ", &
                  "   with a '*' (option .Molecule) are incompatible!"
   call ERRMSG(0," SUBROUTINE read_input")
endif
if (DRELAX) then
   write(*,'(A,F10.4,A)') &
   " Cutoff distance for the intermolecular interaction of a pair of spins:",Rmax*10.d0," angstrom"
   if (Lbox.eq.0.d0) call ERRMSG(0," Missing input option: '.LBOX' ")
endif
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
          if ( symbol(:2) /= Asymbol ) then
             write(emsg,'(A,I3,A,D20.10,5A)') &
' Atom number ',atmlist(nuc),' from timestep ',time,' (type=',symbol(:2), &
') differs from that specified in srelax.inp (type=',Asymbol,')'
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
real(8)                            :: meanRe, meanIm, pi, denon

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

! Normalization factors associated with number of time origins
! ACF normalization
! for ACF(1) : divides by N
! for ACF(2) : divides by N-1
! for ACF(3) : divides by N-2
! ...
! for ACF(N) : divides by 1
      if (WK_divbyN) then
         denon = DBLE(nframes)
         do frm=1,nframes
            odd = 2*frm-1
            ACFre(mcom,frm) = Rdata(odd)  /denon
            ACFim(mcom,frm) = Rdata(odd+1)/denon
         end do
      else
         do frm=1,nframes
            odd = 2*frm-1
            denon = DBLE(nframes-frm+1)
            ACFre(mcom,frm) = Rdata(odd)  /denon
            ACFim(mcom,frm) = Rdata(odd+1)/denon
         end do
      endif
   end do
end do

! deallocating Rdata array
if (allocated(Rdata)) then
   deallocate(Rdata, STAT=ios)
   call AllocStat(ios,'Rdata','WK_ACF_R2m','D')
endif

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
! Frequency Dependent Spectral Densities
! ******************************************************************************
SUBROUTINE SpectDens()
 USE shared_vars
 IMPLICIT NONE

 integer :: ios, frm, mcom, iw0, lstrg, npts
 real(8) :: time, cs, w0, wstep, w_min, w_max
 real(8) :: AvgRe, AvgIm, inv_maxmcom
 real(8) :: SINT, y2(nframes), Fspl_loc(nframes), pi
 real(8), external :: SplineExtfunc

 integer, parameter :: NptsSD = 1000
 real(8), dimension(:,:), allocatable :: SDen

pi = 4.0d0 * atan(1.0d0)

if (LTIME) write(*,'(/A,F10.2/A/)') &
     " Limit Time (in ps): t_lim = ", limtime, &
     " ACF(EFG)(tau > t_lim) = Zero "

! Allocate memory for spline arrays
allocate(Xspl(nframes), STAT=ios)
call AllocStat(ios, 'Xspl', 'SpectDens', 'A')
allocate(Fspl(nframes), STAT=ios)
call AllocStat(ios, 'Fspl', 'SpectDens', 'A')

! Larmor frequency (rad s^-1) - gammaI in MHz/T
w0 = B0ext * gammaI * 1.0d6 * 2.0d0 * pi

! Isotropic averaging of the ACFs
if (ACF_AVG) then
   write(*,'(/A)') " Option Force Isotropy is on! "
   inv_maxmcom = 1.0d0 / real(maxmcom, kind=8)

   do frm = 1, nframes
      AvgRe = sum(ACFre(1:maxmcom, frm)) * inv_maxmcom
      AvgIm = sum(ACFim(1:maxmcom, frm)) * inv_maxmcom
      ACFre(1:maxmcom, frm) = AvgRe
      ACFim(1:maxmcom, frm) = AvgIm
   end do

   call KBLANK(acfisofn, 20, lstrg)
   open(20, file=acfisofn(:lstrg), status='unknown')
   write(20,'(A,1X,I5,1X,E20.10)') "# ", nframes, timestep

! Quadrupolar ACFs
   if (maxmcom == 5) write(20,'(2A)') &
        "#            Time/ps     Re(f_{20})/a.u.    Re(f_{2+1})/a.u.    Re(f_{2-1})/a.u.    Re(f_{2+2})/a.u.    Re(f_{2-2})/a.u.",&
        "    Im(f_{20})/a.u.    Im(f_{2+1})/a.u.    Im(f_{2-1})/a.u.    Im(f_{2+2})/a.u.    Im(f_{2-2})/a.u."
! Dipolar ACFs
   if (DRELAX) write(20,'(2A)') &
        "#            Time/ps     Re(F^ij_0)/a.u.    Re(F^ij_1)/a.u.     Re(F^ij_2)/a.u.     ", &
        "Im(F^ij_0)/a.u.    Im(F^ij_1)/a.u.     Im(F^ij_2)/a.u."
! Spin-Rotation ACFs
   if (SR_RELAX) write(20,'(2A)') &
        "#            Time/ps     Re(J_x)/a.u.       Re(J_y)/a.u.        Re(J_z)/a.u.        ", &
        "Im(J_x)/a.u.       Im(J_y)/a.u.        Im(J_z)/a.u."

   do frm = 1, nframes
      write(20,'(11E20.10)') real(frm-1, kind=8)*timestep, (ACFre(mcom,frm), mcom=1,maxmcom), (ACFim(mcom,frm), mcom=1,maxmcom)
   end do
   close(20)
end if

! Populate time axis (Xspl) and determine npts 
npts = nframes
do frm = 1, nframes
   time = real(frm-1, kind=8) * timestep
   Xspl(frm) = time
   if (LTIME .and. time > limtime) then
      npts = frm
      exit
   end if
end do

! Loop over 'm' components of ACF_2m
do mcom = 1, maxmcom
   sigmaEFG(mcom) = ACFre(mcom, 1)

   ! 1. Pure ACF integration (correlation time)
   Fspl(1:npts) = ACFre(mcom, 1:npts)
   call y2spline(Xspl, Fspl, npts, y2)
   call simpson(SplineExtfunc, npts, Xspl, Fspl, y2, Xspl(1), Xspl(npts), SINT)
   GlmENR(mcom) = SINT

   ! 2. ACF * cos(m * w0 * t) integration
   if (mvec(mcom) == 0) then
      ! For m = 0, cos(0) = 1, so the result equals the correlation time
      GlmRe(mcom) = GlmENR(mcom)
   else
      do frm = 1, npts
         cs = cos(real(mvec(mcom), kind=8) * w0 * Xspl(frm) * 1.0d-12)
         Fspl(frm) = ACFre(mcom, frm) * cs
      end do
      call y2spline(Xspl, Fspl, npts, y2)
      call simpson(SplineExtfunc, npts, Xspl, Fspl, y2, Xspl(1), Xspl(npts), SINT)
      GlmRe(mcom) = SINT
   end if
end do

! Compute Spectral Densities over a frequency range
if (PrintSD) then
   allocate(SDen(maxmcom, NptsSD), STAT=ios)
   call AllocStat(ios, 'SDen', 'SpectDens', 'A')

   ! Use local variables to avoid corrupting global SDmin and SDmax
   w_min = SDmin * gammaI * 1.0d6 * 2.0d0 * pi
   w_max = SDmax * gammaI * 1.0d6 * 2.0d0 * pi
   wstep = (w_max - w_min) / real(NptsSD - 1, kind=8)

   do mcom = 1, maxmcom
      if (mvec(mcom) == 0) then
         ! Avoid 1000 identical integrations when m = 0
         SDen(mcom, 1:NptsSD) = GlmENR(mcom)
      else
         ! Parallel region: each thread gets private copies of Fspl_loc, y2, and SINT
         !$omp parallel do default(none) &
         !$omp shared(w_min, wstep, npts, mvec, mcom, Xspl, ACFre, SDen) &
         !$omp private(iw0, w0, frm, cs, Fspl_loc, y2, SINT) &
         !$omp schedule(static)
         do iw0 = 1, NptsSD
            w0 = w_min + wstep * real(iw0 - 1, kind=8)

            do frm = 1, npts
               cs = cos(real(mvec(mcom), kind=8) * w0 * Xspl(frm) * 1.0d-12)
               Fspl_loc(frm) = ACFre(mcom, frm) * cs
            end do

            call y2spline(Xspl, Fspl_loc, npts, y2)
            call simpson(SplineExtfunc, npts, Xspl, Fspl_loc, y2, Xspl(1), Xspl(npts), SINT)
            SDen(mcom, iw0) = SINT
         end do
         !$omp end parallel do
      end if
   end do

   open(20, file='spectden.dat', status='unknown')
   if (maxmcom == 5) write(20,'(2A)') &
        "# Larmor Freq./MHz      Re(g_{20})/a.u.     Re(g_{2+1})/a.u.    Re(g_{2-1})/a.u.    Re(g_{2+2})/a.u.    Re(g_{2-2})/a.u."

   wstep = (SDmax - SDmin) / real(NptsSD - 1, kind=8)
   do iw0 = 1, NptsSD
      w0 = SDmin + wstep * real(iw0 - 1, kind=8)
      write(20,'(F20.10,5E20.10)') w0, (SDen(mcom, iw0), mcom=1,maxmcom)
   end do
   close(20)

   deallocate(SDen, STAT=ios)
   call AllocStat(ios, 'SDen', 'SpectDens', 'D')
end if

! Deallocate spline arrays
deallocate(Xspl, STAT=ios)
call AllocStat(ios, 'Xspl', 'SpectDens', 'D')
deallocate(Fspl, STAT=ios)
call AllocStat(ios, 'Fspl', 'SpectDens', 'D')
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
! gammaI in MHz/T
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
if (allocated(Xspl)) then
   deallocate(Xspl, STAT=ios)
   call AllocStat(ios,'Xspl',routine,'D')
endif
if (allocated(Fspl)) then
   deallocate(Fspl, STAT=ios)
   call AllocStat(ios,'Fspl',routine,'D')
endif
if (allocated(newX)) then
   deallocate(newX, STAT=ios)
   call AllocStat(ios,'newX',routine,'D')
endif
if (allocated(newF)) then
   deallocate(newF, STAT=ios)
   call AllocStat(ios,'newF',routine,'D')
endif
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
SUBROUTINE ncspline(arg, x, y, ndim, y2, fx)
 IMPLICIT NONE

 integer, intent(in)  :: ndim
 real(8), intent(in)  :: x(*), y(*), y2(ndim), arg
 real(8), intent(out) :: fx

 integer :: j, jhi, jlo
 real(8) :: A, B, xx, xx2, dx, x_span
 real(8), parameter :: inv6 = 1.0d0 / 6.0d0

! Fast O(1) index lookup for uniform grids (like Xspl)
dx = x(2) - x(1)
x_span = x(ndim) - x(1)

if (abs(x_span - dx * real(ndim - 1, kind=8)) <= 1.0d-10 * abs(x_span)) then
   jlo = int((arg - x(1)) / dx) + 1
   if (jlo < 1) jlo = 1
   if (jlo >= ndim) jlo = ndim - 1
   jhi = jlo + 1
else
   ! Fallback to bisection for non-uniform grids
   jlo = 1
   jhi = ndim
   do while (jhi - jlo > 1)
      j = (jhi + jlo) / 2
      if (x(j) > arg) then
         jhi = j
      else
         jlo = j
      end if
   end do
end if

xx = x(jhi) - x(jlo)
if (xx < 1.0d-10) call ERRMSG(0, "x(jhi) - x(jlo) < 1e-10 in ncspline")

A   = (x(jhi) - arg) / xx
B   = 1.0d0 - A
xx2 = xx * xx

fx = A * y(jlo) + B * y(jhi) + &
     ((A * (A * A - 1.0d0)) * y2(jlo) + &
      (B * (B * B - 1.0d0)) * y2(jhi)) * xx2 * inv6
END SUBROUTINE ncspline

! ******************************************************************************
! Generating second derivatives of y(x_i) by
! solving the tridiagonal system of equations (Eq. 3.3.7 of NR)
! using LU decomposition 
! ******************************************************************************
SUBROUTINE y2spline(x, y, ndim, y2)
 IMPLICIT NONE

 integer, intent(in)  :: ndim
 real(8), intent(in)  :: x(*), y(*)
 real(8), intent(out) :: y2(*)

 integer :: i
 real(8) :: dx_prev, dx_curr, slope_prev, slope_curr
 real(8) :: h6, alpha, z(ndim)
 real(8), parameter :: inv3 = 1.0d0 / 3.0d0
 real(8), parameter :: inv6 = 1.0d0 / 6.0d0

! Natural spline boundary conditions: y''_1 = y''_N = 0
y2(1) = 0.0d0
if (ndim < 2) return
y2(ndim) = 0.0d0
if (ndim == 2) return

! Initialization for i = 2
dx_prev    = x(2) - x(1)
dx_curr    = x(3) - x(2)
slope_prev = (y(2) - y(1)) / dx_prev
slope_curr = (y(3) - y(2)) / dx_curr

z(2)  = slope_curr - slope_prev
y2(2) = (x(3) - x(1)) * inv3

! Forward substitution reusing previous step intervals and slopes
do i = 3, ndim - 1
   dx_prev    = dx_curr
   slope_prev = slope_curr

   dx_curr    = x(i+1) - x(i)
   slope_curr = (y(i+1) - y(i)) / dx_curr

   h6    = dx_prev * inv6
   alpha = h6 / y2(i-1)

   y2(i) = (x(i+1) - x(i-1)) * inv3 - alpha * h6
   z(i)  = (slope_curr - slope_prev) - alpha * z(i-1)
end do

! Backsubstitution
do i = ndim - 1, 2, -1
   y2(i) = (z(i) - y2(i+1) * (x(i+1) - x(i)) * inv6) / y2(i)
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
Parameter (Thresh=1.d-8, Nmax=20)
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
   if ( Dabs((Sint-pSint)/pSint).LT.Thresh ) return
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
if (Inuc == 1.d40) then
   write(*,*) mass, symbol
   call ERRMSG(01," SUBROUTINE SpinGammaI: unknown isotope!")
endif
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
SUBROUTINE getIsotope(line,linelen,mass,Asymbol)
!*******************************************************************************
Implicit None

character :: line*(*), Asymbol*2, number*1
integer   :: char, Masslen, num, linelen, mass

Masslen=0
do char=1,linelen
  do num=0,9
     write(number,'(I1)') num
     if ( line(char:char) == number ) Masslen=Masslen+1
  end do
enddo
if ( linelen-Masslen == 1) then
   Asymbol=line(Masslen+1:Masslen+1)//' '
else
   Asymbol=line(Masslen+1:Masslen+2)
endif
mass=ICHAR(line(:Masslen))
if (Masslen==1) read(line(:Masslen),'(I1)') mass
if (Masslen==2) read(line(:Masslen),'(I2)') mass
if (Masslen==3) read(line(:Masslen),'(I3)') mass
END SUBROUTINE getIsotope

! ******************************************************************************
! Reading Gromacs output file (.gro)
! ******************************************************************************
SUBROUTINE ReadGRO()
Use shared_vars
Implicit None

integer   :: lfname, ios, iatm, iframe, jframe, line, unit_id
integer   :: frame_base, skip_lines
character :: emsg*80

! Remove trailing blank spaces from the filename
call KBLANK(MDtrajfile, lMDfile, lfname)

open(newunit=unit_id, file=MDtrajfile(:lfname), status='old', action='read', iostat=ios)

! Check if GROMACS trajectory file exists
if (ios /= 0) then
   write(emsg,'(A,A)') "Failed to open GROMACS trajectory file named: ", MDtrajfile(:lfname)
   call ERRMSG(unit_id, emsg)
endif
write(*,"(A,A)") " Reading GROMACS trajectory file ", MDtrajfile(:lfname)

! Read the number of atoms from the second line of the first frame
read(unit_id, '(A)', iostat=ios)
call read_ios(unit_id, ios, MDtrajfile(:lfname))
read(unit_id, *, iostat=ios) natom
call read_ios(unit_id, ios, MDtrajfile(:lfname))
backspace(unit_id)
backspace(unit_id)

if (RDDACF .and. .not. MSD_COM) then
   close(unit_id)
   return
endif

! Allocate XYZ coordinates array with contiguous layout (3, total_atoms)
allocate(AtmPos(3, nframes*natom), STAT=ios)
call AllocStat(ios, 'AtmPos', 'ReadGRO', 'A')

! Skip frames before MDstartProd using '(A)' to avoid list-directed parsing overhead
skip_lines = (MDstartProd - 1) * (natom + 3)
do line = 1, skip_lines
   read(unit_id, '(A)', iostat=ios)
   call read_ios(unit_id, ios, MDtrajfile(:lfname))
enddo

! Pre-compute number of lines to skip when MDsampleFreq > 1
skip_lines = (MDsampleFreq - 1) * (natom + 3)

! Read frames from MDstartProd to MDendProd
do iframe = 0, nframes - 1
   frame_base = iframe * natom

   ! Skip title and atom count header lines of the GROMACS frame
   read(unit_id, '(A)', iostat=ios)
   call read_ios(unit_id, ios, MDtrajfile(:lfname))
   read(unit_id, '(A)', iostat=ios)
   call read_ios(unit_id, ios, MDtrajfile(:lfname))

   ! Read atomic coordinates directly into contiguous slice (already in nm)
   do iatm = 1, natom
      read(unit_id, '(20X,3F8.3)', iostat=ios) AtmPos(1:3, frame_base + iatm)
      call read_ios(unit_id, ios, MDtrajfile(:lfname))
   end do

   ! Read box length
   read(unit_id, *, iostat=ios) Lbox
   call read_ios(unit_id, ios, MDtrajfile(:lfname))

   ! Skip intermediate frames if MDsampleFreq > 1 (except after the last frame)
   if (MDsampleFreq > 1 .and. iframe /= nframes - 1) then
      do jframe = 1, skip_lines
         read(unit_id, '(A)', iostat=ios)
         call read_ios(unit_id, ios, MDtrajfile(:lfname))
      enddo
   endif
end do
close(unit_id)
END SUBROUTINE ReadGRO

! ******************************************************************
! Print an error message and close the unit if ios not equal to zero
! ******************************************************************
SUBROUTINE read_ios(UNIT, ios, filename)
Implicit None
integer,          intent(in) :: UNIT, ios
character(len=*), intent(in) :: filename

! Fast path: return immediately when read succeeds (99.99% of calls)
if (ios == 0) return

! Check if end-of-file was reached (negative value)
if (ios < 0) then
   write(*,'(A,A)') " ERROR: End of file reached while reading ", trim(filename)
   close(UNIT)
   error stop 1
else
   write(*,'(A,A)') " ERROR: An error occurred during reading of file ", trim(filename)
   close(UNIT)
   error stop 1
end if
END SUBROUTINE read_ios

! ******************************************************************************
! Reading Tinker output file .arc
! ******************************************************************************
SUBROUTINE ReadTINKER()
Use shared_vars
Implicit None

integer   :: lfname, ios, iatm, iframe, jframe, line, unit_id
integer   :: frame_base, skip_lines
character :: emsg*80

! Remove trailing blank spaces from the filename
call KBLANK(MDtrajfile, lMDfile, lfname)

open(newunit=unit_id, file=MDtrajfile(:lfname), status='old', action='read', iostat=ios)

! Check if trajectory file exists
if (ios /= 0) then
   write(emsg,'(A,A)') "Failed to open Tinker trajectory file named: ", MDtrajfile(:lfname)
   call ERRMSG(unit_id, emsg)
endif
write(*,"(A,A/)") " Reading Tinker trajectory file: ", MDtrajfile(:lfname)

! Read the number of atoms
read(unit_id, *, iostat=ios) natom
call read_ios(unit_id, ios, MDtrajfile(:lfname))
backspace(unit_id)

if (RDDACF) then
   close(unit_id)
   return
endif

! Allocate XYZ coordinates array with contiguous layout (3, total_atoms)
allocate(AtmPos(3, nframes*natom), STAT=ios)
call AllocStat(ios, 'AtmPos', 'ReadTINKER', 'A')

! Skip frames before MDstartProd using '(A)' to avoid list-directed parsing overhead
skip_lines = (MDstartProd - 1) * (natom + 2)
do line = 1, skip_lines
   read(unit_id, '(A)', iostat=ios)
   call read_ios(unit_id, ios, MDtrajfile(:lfname))
enddo

! Pre-compute number of lines to skip when MDsampleFreq > 1
skip_lines = (MDsampleFreq - 1) * (natom + 2)

! Read frames from MDstartProd to MDendProd
do iframe = 0, nframes - 1
   frame_base = iframe * natom

   ! Skip the two header lines of each Tinker frame
   read(unit_id, '(A)', iostat=ios)
   call read_ios(unit_id, ios, MDtrajfile(:lfname))
   read(unit_id, '(A)', iostat=ios)
   call read_ios(unit_id, ios, MDtrajfile(:lfname))

   ! Read atomic coordinates directly into contiguous slice
   do iatm = 1, natom
      read(unit_id, '(11X,3F12.6)', iostat=ios) AtmPos(1:3, frame_base + iatm)
      call read_ios(unit_id, ios, MDtrajfile(:lfname))
   end do

   ! Skip intermediate frames if MDsampleFreq > 1
   if (MDsampleFreq > 1) then
      do jframe = 1, skip_lines
         read(unit_id, '(A)', iostat=ios)
         call read_ios(unit_id, ios, MDtrajfile(:lfname))
      enddo
   endif
end do
close(unit_id)

! Convert coordinates from Angstrom to nm in a single vectorized pass
AtmPos(:,:) = AtmPos(:,:) * 0.1d0
END SUBROUTINE ReadTINKER

! ******************************************************************************
! Reading Tinker velocity file *.vel
! ******************************************************************************
SUBROUTINE ReadTINKERvel()
Use shared_vars
Implicit None

integer   :: lfname, ios, iatm, iframe, jframe, line, unit_id
integer   :: frame_base, skip_lines
character :: emsg*80

! Remove trailing blank spaces from the filename
call KBLANK(MDvelfile, lMDfile, lfname)

open(newunit=unit_id, file=MDvelfile(:lfname), status='old', action='read', iostat=ios)

! Check if velocity file exists
if (ios /= 0) then
   write(emsg,'(A,A)') "Failed to open Tinker trajectory file named: ", MDvelfile(:lfname)
   call ERRMSG(unit_id, emsg)
endif
write(*,"(/A,A/)") " Reading Tinker velocity file: ", MDvelfile(:lfname)

! Read the number of atoms
read(unit_id, *, iostat=ios) natom
call read_ios(unit_id, ios, MDvelfile(:lfname))
backspace(unit_id)

! Return if ACF will be read from results of a previous calculation
if (RDDACF) then
   close(unit_id)
   return
endif

! Allocate velocity array with contiguous layout 
allocate(AtmVel(3, nframes*natom), STAT=ios)
call AllocStat(ios, 'AtmVel', 'ReadTINKERvel', 'A')

! Skip frames before MDstartProd using '(A)' to avoid list-directed parsing overhead
skip_lines = (MDstartProd - 1) * (natom + 2)
do line = 1, skip_lines
   read(unit_id, '(A)', iostat=ios)
   call read_ios(unit_id, ios, MDvelfile(:lfname))
enddo

! Pre-compute number of lines to skip when MDsampleFreq > 1
skip_lines = (MDsampleFreq - 1) * (natom + 2)

! Read frames from MDstartProd to MDendProd
do iframe = 0, nframes - 1
   frame_base = iframe * natom

   ! Skip the two header lines of each Tinker frame
   read(unit_id, '(A)', iostat=ios)
   call read_ios(unit_id, ios, MDvelfile(:lfname))
   read(unit_id, '(A)', iostat=ios)
   call read_ios(unit_id, ios, MDvelfile(:lfname))

   ! Read atomic velocities directly into contiguous slice
   do iatm = 1, natom
      read(unit_id, '(11X,3F12.6)', iostat=ios) AtmVel(1:3, frame_base + iatm)
      call read_ios(unit_id, ios, MDvelfile(:lfname))
   end do

   ! Skip intermediate frames if MDsampleFreq > 1
   if (MDsampleFreq > 1) then
      do jframe = 1, skip_lines
         read(unit_id, '(A)', iostat=ios)
         call read_ios(unit_id, ios, MDvelfile(:lfname))
      enddo
   endif
end do
close(unit_id)

! Convert units from A/ps to nm/ps in a single vectorized pass
AtmVel(:,:) = AtmVel(:,:) * 0.1d0
END SUBROUTINE ReadTINKERvel

! ******************************************************************************
! Reading CP2K output file (*-pos-1.xyz)
! ******************************************************************************
SUBROUTINE ReadCP2K()
Use shared_vars
Implicit None

integer   :: lfname, ios, iatm, iframe, jframe, line, unit_id
integer   :: frame_base, skip_lines
character :: emsg*80

! Remove trailing blank spaces from the filename
call KBLANK(MDtrajfile, lMDfile, lfname)

open(newunit=unit_id, file=MDtrajfile(:lfname), status='old', action='read', iostat=ios)

! Check if trajectory file exists
if (ios /= 0) then
   write(emsg,'(A,A)') "Failed to open CP2K trajectory file named: ", MDtrajfile(:lfname)
   call ERRMSG(unit_id, emsg)
endif
write(*,"(A,A/)") " Reading CP2K trajectory file: ", MDtrajfile(:lfname)

! Read the number of atoms
read(unit_id, *, iostat=ios) natom
call read_ios(unit_id, ios, MDtrajfile(:lfname))
backspace(unit_id)

if (RDDACF) then
   close(unit_id)
   return
endif

! Allocate XYZ coordinates array with contiguous layout (3, total_atoms)
allocate(AtmPos(3, nframes*natom), STAT=ios)
call AllocStat(ios, 'AtmPos', 'ReadCP2K', 'A')

! Skip frames before MDstartProd using '(A)' to avoid list-directed parsing overhead
skip_lines = (MDstartProd - 1) * (natom + 2)
do line = 1, skip_lines
   read(unit_id, '(A)', iostat=ios)
   call read_ios(unit_id, ios, MDtrajfile(:lfname))
enddo

! Pre-compute number of lines to skip when MDsampleFreq > 1
skip_lines = (MDsampleFreq - 1) * (natom + 2)

! Read frames from MDstartProd to MDendProd
do iframe = 0, nframes - 1
   frame_base = iframe * natom

   ! Skip the two header lines of each XYZ frame
   read(unit_id, '(A)', iostat=ios)
   call read_ios(unit_id, ios, MDtrajfile(:lfname))
   read(unit_id, '(A)', iostat=ios)
   call read_ios(unit_id, ios, MDtrajfile(:lfname))

   ! Read atomic coordinates directly into contiguous slice
   do iatm = 1, natom
      read(unit_id, '(4X,3F20.10)', iostat=ios) AtmPos(1:3, frame_base + iatm)
      call read_ios(unit_id, ios, MDtrajfile(:lfname))
   end do

   ! Skip intermediate frames if MDsampleFreq > 1
   if (MDsampleFreq > 1) then
      do jframe = 1, skip_lines
         read(unit_id, '(A)', iostat=ios)
         call read_ios(unit_id, ios, MDtrajfile(:lfname))
      enddo
   endif
end do
close(unit_id)

! Convert coordinates from Angstrom to nm in a single vectorized pass
AtmPos(:,:) = AtmPos(:,:) * 0.1d0
END SUBROUTINE ReadCP2K

! ******************************************************************************
! Moving the CoM of each molecule, in each frame, to the origin
! ******************************************************************************
SUBROUTINE molLab2CoM()
Use shared_vars
Use omp_lib
Implicit None

integer :: Nmol, imol, iMolAtm, iframe, pos, frame_base, mol_base
real(8) :: Tmass, inv_Tmass, mass, com_x, com_y, com_z

! Number of molecules in each frame
Nmol = natom / MolAtm

! Total mass of a single molecule and its inverse
Tmass = SUM(AtomMass(1:MolAtm))
inv_Tmass = 1.d0 / Tmass

! Parallelize over frames 
!$omp parallel do default(none) &
!$omp shared(nframes, Nmol, natom, MolAtm, AtomMass, AtmPos, inv_Tmass) &
!$omp private(iframe, frame_base, imol, mol_base, com_x, com_y, com_z, iMolAtm, pos, mass) &
!$omp schedule(static)
do iframe = 0, nframes - 1
   frame_base = iframe * natom

   do imol = 1, Nmol
      mol_base = frame_base + (imol - 1) * MolAtm

      com_x = 0.d0
      com_y = 0.d0
      com_z = 0.d0

      ! 1. Accumulate mass-weighted atomic positions
      do iMolAtm = 1, MolAtm
         pos  = mol_base + iMolAtm
         mass = AtomMass(iMolAtm)
         com_x = com_x + AtmPos(1, pos) * mass
         com_y = com_y + AtmPos(2, pos) * mass
         com_z = com_z + AtmPos(3, pos) * mass
      end do

      ! Compute center of mass coordinates
      com_x = com_x * inv_Tmass
      com_y = com_y * inv_Tmass
      com_z = com_z * inv_Tmass

      ! 2. Shift atomic positions to the molecular center-of-mass frame
      do iMolAtm = 1, MolAtm
         pos = mol_base + iMolAtm
         AtmPos(1, pos) = AtmPos(1, pos) - com_x
         AtmPos(2, pos) = AtmPos(2, pos) - com_y
         AtmPos(3, pos) = AtmPos(3, pos) - com_z
      end do

   end do
end do
!$omp end parallel do
END SUBROUTINE molLab2CoM

! ******************************************************************************
! Removing the CoM velocity from the atomic velocities 
! ******************************************************************************
SUBROUTINE molLab2CoMvel()
Use shared_vars
Use omp_lib
Implicit None

integer :: Nmol, imol, iMolAtm, iframe, pos, frame_base, mol_base
real(8) :: Tmass, inv_Tmass, mass, vcom_x, vcom_y, vcom_z

! Number of molecules in each frame
Nmol = natom / MolAtm

! Total mass of a single molecule and its inverse (computed once)
Tmass = SUM(AtomMass(1:MolAtm))
inv_Tmass = 1.d0 / Tmass

! Parallelize over frames (each frame is completely independent)
!$omp parallel do default(none) &
!$omp shared(nframes, Nmol, natom, MolAtm, AtomMass, AtmVel, inv_Tmass) &
!$omp private(iframe, frame_base, imol, mol_base, vcom_x, vcom_y, vcom_z, iMolAtm, pos, mass) &
!$omp schedule(static)
do iframe = 0, nframes - 1
   frame_base = iframe * natom

   do imol = 1, Nmol
      mol_base = frame_base + (imol - 1) * MolAtm

      vcom_x = 0.d0
      vcom_y = 0.d0
      vcom_z = 0.d0

      ! 1. Accumulate mass-weighted atomic velocities
      do iMolAtm = 1, MolAtm
         pos  = mol_base + iMolAtm
         mass = AtomMass(iMolAtm)
         vcom_x = vcom_x + AtmVel(1, pos) * mass
         vcom_y = vcom_y + AtmVel(2, pos) * mass
         vcom_z = vcom_z + AtmVel(3, pos) * mass
      end do

      ! Compute center-of-mass velocity
      vcom_x = vcom_x * inv_Tmass
      vcom_y = vcom_y * inv_Tmass
      vcom_z = vcom_z * inv_Tmass

      ! 2. Shift atomic velocities to the molecular center-of-mass frame
      do iMolAtm = 1, MolAtm
         pos = mol_base + iMolAtm
         AtmVel(1, pos) = AtmVel(1, pos) - vcom_x
         AtmVel(2, pos) = AtmVel(2, pos) - vcom_y
         AtmVel(3, pos) = AtmVel(3, pos) - vcom_z
      end do

   end do
end do
!$omp end parallel do
END SUBROUTINE molLab2CoMvel

! ****************************************************************************
! Determining the molecular angular momentum from the molecular inertia tensor 
! and angular velocity. The last one is obtained using equation 13 from the 
! reference: https://doi.org/10.1021/acs.jpcb.2c06411
!            J. Phys. Chem. B 2023, 127, 1167−1177
! ****************************************************************************
SUBROUTINE J_MolFrame()
 Use shared_vars
 Use omp_lib
 Use, intrinsic :: iso_c_binding
 Implicit None
 Include 'fftw3.f03'

 real(8) :: vecM(3,4), ivec(3), jvec(3), kvec(3), norm, &
            rLab(3), vLab(3), rMol(3), vMol(3), vec(3), I_ts(3,3), mass, conv
 real(8) :: prod11, prod22, prod33, prod12, prod13, prod23, &
            prod_diag1, prod_diag2, prod_diag3
 real(8) :: Mxyz(3,3), b_vec(3), omega(3)
 real(8) :: det, inv_det, c11, c12, c13, c22, c23, c33
 real(8), parameter :: ZeroNorm = 1.d-6
 real(8), dimension(:,:), allocatable :: meanRe, meanIm, ACFre_loc, ACFim_loc
 complex(8), dimension(:,:), allocatable :: T2m_loc
 complex(C_DOUBLE_COMPLEX), dimension(:), allocatable :: in_arr, out_arr, dummy_in, dummy_out
 type(C_PTR) :: plan_fwd, plan_bwd
 integer :: i, iframe, imol, Nmol, pos, iMolAtm, ios
 integer :: mol_base, frame_base, atom0(4)
 logical :: is_zero(4)

allocate(ACFre(maxmcom,nframes), STAT=ios)
call AllocStat(ios,'ACFre','J_MolFrame','A')
allocate(ACFim(maxmcom,nframes), STAT=ios)
call AllocStat(ios,'ACFim','J_MolFrame','A')
allocate(meanRe(maxmcom,nframes), STAT=ios)
call AllocStat(ios,'meanRe','J_MolFrame','A')
allocate(meanIm(maxmcom,nframes), STAT=ios)
call AllocStat(ios,'meanIm','J_MolFrame','A')

Nmol = natom / MolAtm
if (SR_nmol .eq. 0) SR_nmol = Nmol
SR_nmol = min(Nmol, SR_nmol)

meanRe(:,:) = 0.d0
meanIm(:,:) = 0.d0
is_zero(:)  = (MolFrame(:) .eq. 0)

! SINGLE creation of FFTW plans before the parallel region
allocate(dummy_in(fft_size), dummy_out(fft_size))
plan_fwd = fftw_plan_dft_1d(fft_size, dummy_in, dummy_out, FFTW_FORWARD, FFTW_ESTIMATE)
plan_bwd = fftw_plan_dft_1d(fft_size, dummy_out, dummy_in, FFTW_BACKWARD, FFTW_ESTIMATE)
deallocate(dummy_in, dummy_out)

!$omp parallel default(none) &
!$omp shared(SR_nmol, MolAtm, MolFrame, is_zero, nframes, natom, AtmPos, AtmVel, AtomMass, maxmcom, fft_size, plan_fwd, plan_bwd) &
!$omp private(imol, mol_base, atom0, iframe, frame_base, i, vecM, kvec, norm, vec, jvec, ivec, &
!$omp         Mxyz, I_ts, b_vec, iMolAtm, mass, pos, rLab, vLab, rMol, vMol, &
!$omp         prod11, prod22, prod33, prod12, prod13, prod23, prod_diag1, prod_diag2, prod_diag3, &
!$omp         c11, c12, c13, det, inv_det, c22, c23, c33, omega, T2m_loc, ACFre_loc, ACFim_loc, in_arr, out_arr) &
!$omp reduction(+:meanRe, meanIm)

! Thread-private buffers
allocate(T2m_loc(maxmcom,nframes))
allocate(ACFre_loc(maxmcom,nframes))
allocate(ACFim_loc(maxmcom,nframes))
allocate(in_arr(fft_size))
allocate(out_arr(fft_size))

!$omp do schedule(static)
do imol = 1, SR_nmol
   mol_base = (imol -1) *MolAtm
   atom0(:) = mol_base +MolFrame(:)

   do iframe = 0, nframes-1
      frame_base = iframe * natom + mol_base

! Gathering the coordinates of the atoms that will define the molecular coordinate system
      do i = 1,4
! Use CoM position 
         if (is_zero(i)) then
            vecM(:,i) = 0.d0
! Use Atom position
         else
            vecM(:,i) = AtmPos(:, iframe *natom +atom0(i))
         end if
      end do

! z-axis: vec{k} propto vec{atom MolFrame(1)} -vec{atom MolFrame(2)}
      kvec(:) = vecM(:,1) -vecM(:,2)
      norm    = sqrt(kvec(1)*kvec(1) + kvec(2)*kvec(2) + kvec(3)*kvec(3))
      if (norm .lt. ZeroNorm) then
         call ERRMSG(0, "vec{k} has zero norm!")
      else
         kvec(:) = kvec(:) / norm
      end if

! y-axis: vec{j} propto  vec{k} x (vec{atom MolFrame(3)} -vec{atom MolFrame(4)})
      vec(:)  = vecM(:,3) - vecM(:,4)
      !Cross product
      jvec(1) = kvec(2)*vec(3) - kvec(3)*vec(2)
      jvec(2) = kvec(3)*vec(1) - kvec(1)*vec(3)
      jvec(3) = kvec(1)*vec(2) - kvec(2)*vec(1)
      norm    = sqrt(jvec(1)*jvec(1) + jvec(2)*jvec(2) + jvec(3)*jvec(3))
      if (norm .lt. ZeroNorm) then
         call ERRMSG(0, "vec{j} has zero norm!")
      else
         jvec(:) = jvec(:) / norm
      end if

! x-axis: vec{i} = vec{k} x  vec{j}
      ! Cross product
      ivec(1) = kvec(2)*jvec(3) - kvec(3)*jvec(2)
      ivec(2) = kvec(3)*jvec(1) - kvec(1)*jvec(3)
      ivec(3) = kvec(1)*jvec(2) - kvec(2)*jvec(1)

      Mxyz(:,:) = 0.d0
      I_ts(:,:) = 0.d0
      b_vec(:)  = 0.d0

      do iMolAtm = 1, MolAtm
         mass    = AtomMass(iMolAtm)
         pos     = frame_base + iMolAtm
         rLab(:) = AtmPos(:,pos)
         vLab(:) = AtmVel(:,pos)

         rMol(1) = ivec(1)*rLab(1) + ivec(2)*rLab(2) + ivec(3)*rLab(3)
         rMol(2) = jvec(1)*rLab(1) + jvec(2)*rLab(2) + jvec(3)*rLab(3)
         rMol(3) = kvec(1)*rLab(1) + kvec(2)*rLab(2) + kvec(3)*rLab(3)

         vMol(1) = ivec(1)*vLab(1) + ivec(2)*vLab(2) + ivec(3)*vLab(3)
         vMol(2) = jvec(1)*vLab(1) + jvec(2)*vLab(2) + jvec(3)*vLab(3)
         vMol(3) = kvec(1)*vLab(1) + kvec(2)*vLab(2) + kvec(3)*vLab(3)

         prod11 = rMol(1) * rMol(1)
         prod22 = rMol(2) * rMol(2)
         prod33 = rMol(3) * rMol(3)
         prod12 = rMol(1) * rMol(2)
         prod13 = rMol(1) * rMol(3)
         prod23 = rMol(2) * rMol(3)

         prod_diag1 = prod22 + prod33
         prod_diag2 = prod11 + prod33
         prod_diag3 = prod11 + prod22

         Mxyz(1,1) = Mxyz(1,1) + prod_diag1
         I_ts(1,1) = I_ts(1,1) + prod_diag1 * mass

         Mxyz(2,2) = Mxyz(2,2) + prod_diag2
         I_ts(2,2) = I_ts(2,2) + prod_diag2 * mass

         Mxyz(3,3) = Mxyz(3,3) + prod_diag3
         I_ts(3,3) = I_ts(3,3) + prod_diag3 * mass

         Mxyz(1,2) = Mxyz(1,2) - prod12
         I_ts(1,2) = I_ts(1,2) - prod12 * mass

         Mxyz(1,3) = Mxyz(1,3) - prod13
         I_ts(1,3) = I_ts(1,3) - prod13 * mass

         Mxyz(2,3) = Mxyz(2,3) - prod23
         I_ts(2,3) = I_ts(2,3) - prod23 * mass

         b_vec(1) = b_vec(1) + (rMol(2)*vMol(3) - rMol(3)*vMol(2))
         b_vec(2) = b_vec(2) + (rMol(3)*vMol(1) - rMol(1)*vMol(3))
         b_vec(3) = b_vec(3) + (rMol(1)*vMol(2) - rMol(2)*vMol(1))
      end do

! Solving 3x3 linear system of equations using Cramer's rule

! Compute the main determinant 
      c11 = Mxyz(2,2)*Mxyz(3,3) - Mxyz(2,3)*Mxyz(2,3)
      c12 = Mxyz(1,3)*Mxyz(2,3) - Mxyz(1,2)*Mxyz(3,3)
      c13 = Mxyz(1,2)*Mxyz(2,3) - Mxyz(1,3)*Mxyz(2,2)
      det = Mxyz(1,1)*c11 + Mxyz(1,2)*c12 + Mxyz(1,3)*c13

      if (det .le. 1.d-20) then
         PRINT *, "Error: Mxyz matrix is singular or not positive-definite at mol/frame:", imol, iframe
         omega(:) = 0.d0
      else
         inv_det = 1.d0 / det
!Cofactors of the Symmetric Matrix 
         c22 = Mxyz(1,1)*Mxyz(3,3) - Mxyz(1,3)*Mxyz(1,3)
         c23 = Mxyz(1,2)*Mxyz(1,3) - Mxyz(1,1)*Mxyz(2,3)
         c33 = Mxyz(1,1)*Mxyz(2,2) - Mxyz(1,2)*Mxyz(1,2)
! Angular velocity : rad/ps
         omega(1) = (c11*b_vec(1) + c12*b_vec(2) + c13*b_vec(3)) * inv_det
         omega(2) = (c12*b_vec(1) + c22*b_vec(2) + c23*b_vec(3)) * inv_det
         omega(3) = (c13*b_vec(1) + c23*b_vec(2) + c33*b_vec(3)) * inv_det
      end if

! Angular momentum (Eq.3) in u.nm^2/ps
      vec(1) = I_ts(1,1)*omega(1) + I_ts(1,2)*omega(2) + I_ts(1,3)*omega(3)
      vec(2) = I_ts(1,2)*omega(1) + I_ts(2,2)*omega(2) + I_ts(2,3)*omega(3)
      vec(3) = I_ts(1,3)*omega(1) + I_ts(2,3)*omega(2) + I_ts(3,3)*omega(3)

      T2m_loc(:,iframe+1) = dcmplx(vec(:))
   end do

! ACF of the current molecul
   call FFT_OMP(T2m_loc, ACFre_loc, ACFim_loc, plan_fwd, plan_bwd, in_arr, out_arr)

! Storing them to compute ensemble averages
   meanRe(:,:) = meanRe(:,:) + ACFre_loc(:,:)
   meanIm(:,:) = meanIm(:,:) + ACFim_loc(:,:)
end do
!$omp end do

deallocate(T2m_loc, ACFre_loc, ACFim_loc, in_arr, out_arr)
!$omp end parallel

! Destroy the plans once after the parallel loop
call fftw_destroy_plan(plan_fwd)
call fftw_destroy_plan(plan_bwd)

! converting ACFs to a.u., i.e. hbar^2
conv  = ((1.d-18 / Navog / 1.d3 / 1.d-12) / hbar)**2 / dble(SR_nmol)
ACFre = meanRe * conv
ACFim = meanIm * conv

! Saving Angular Momentum ACFs
open(20,file='acf-J.dat',status='unknown')
write(20,'(A,1X,I5,1X,E20.10)') "# ", nframes, timestep
write(20,'(A,A)') &
      "#            Time/ps     Re(J_x)/a.u.        Re(J_y)/a.u.        Re(J_z)/a.u.        ", &
                               "Im(J_x)/a.u.        Im(J_y)/a.u.        Im(J_z)/a.u."
do iframe=1,nframes
   write(20,'(7E20.10)')  DBLE(iframe-1)*timestep, (ACFre(i,iframe),i=1,maxmcom), (ACFim(i,iframe),i=1,maxmcom)
end do
close(20)

! Spectral Densities
acfisofn='acf-J.iso'
call SpectDens()

! Relaxation Properties
call SRrelaxProp()

! Releasing memory
deallocate(meanRe, STAT=ios)
call AllocStat(ios,'meanRe','J_MolFrame','D')
deallocate(meanIm, STAT=ios)
call AllocStat(ios,'meanIm','J_MolFrame','D')
deallocate(ACFre, STAT=ios)
call AllocStat(ios,'ACFre','J_MolFrame','D')
deallocate(ACFim, STAT=ios)
call AllocStat(ios,'ACFim','J_MolFrame','D')
END SUBROUTINE J_MolFrame

! ******************************************************************************
! Cross product
! ******************************************************************************
SUBROUTINE crossp(u, v, w) 
REAL(8), DIMENSION(3), INTENT(IN) :: u, v
REAL(8), DIMENSION(3) :: w
w(1) = u(2)*v(3) - u(3)*v(2)
w(2) = u(3)*v(1) - u(1)*v(3)
w(3) = u(1)*v(2) - u(2)*v(1)
END SUBROUTINE crossp

! ******************************************************************************
! Routine to estimate the linear velocity of atoms between consecutive frames
! velocities in nm/ps
! ******************************************************************************
SUBROUTINE velCoMframe()
Use shared_vars
Use omp_lib
Implicit None

integer :: ios, iatom, iframe, frame_base1, frame_base2
integer :: src_start, src_end, dst_start, dst_end
real(8) :: inv_dt

! Allocate atomic velocity array
allocate(AtmVel(3, nframes*natom), STAT=ios)
call AllocStat(ios, 'AtmVel', 'velCoMframe', 'A')

inv_dt = 1.d0 / timestep

! 1. Compute atomic velocities via forward finite differences (parallelized over frames)
!$omp parallel do default(none) &
!$omp shared(nframes, natom, AtmPos, AtmVel, inv_dt) &
!$omp private(iframe, frame_base1, frame_base2, iatom) &
!$omp schedule(static)
do iframe = 0, nframes - 2
   frame_base1 =  iframe      * natom
   frame_base2 = (iframe + 1) * natom

   do iatom = 1, natom
      AtmVel(1, frame_base1 + iatom) = (AtmPos(1, frame_base2 + iatom) - AtmPos(1, frame_base1 + iatom)) * inv_dt
      AtmVel(2, frame_base1 + iatom) = (AtmPos(2, frame_base2 + iatom) - AtmPos(2, frame_base1 + iatom)) * inv_dt
      AtmVel(3, frame_base1 + iatom) = (AtmPos(3, frame_base2 + iatom) - AtmPos(3, frame_base1 + iatom)) * inv_dt
   end do
end do
!$omp end parallel do

! 2. Shift trajectory sequentially frame-by-frame to avoid temporary array allocation on overlapping slices
do iframe = 1, nframes - 1
   dst_start = (iframe - 1) * natom + 1
   dst_end   =  iframe      * natom
   src_start =  iframe      * natom + 1
   src_end   = (iframe + 1) * natom

   AtmPos(:, dst_start:dst_end) = AtmPos(:, src_start:src_end)
end do

! 3. Update the frame counter
nframes = nframes - 1

write(*,'(A/A/)') " Note: The number of frames has been reduced by one (the first one was discarded)", &
                  "       so that atomic velocities could be estimated."
END SUBROUTINE velCoMframe

! ******************************************************************************
! Computes the Mean Square Displacement (MSD) for the Center of Mass (com)
! ******************************************************************************
SUBROUTINE MSDcom()
Use shared_vars
Use omp_lib
Implicit None

integer, parameter :: MAX_LAG_CAP = 25000
integer :: ios, iatom, Nmol, imol, iMolAtm, iframe, pos1, pos2, frm1, frm2, tm, unit_id, max_lag
integer :: frame_base, prev_base, mol_base, com_base, com_base1, com_base2
real(8) :: Wmass, inv_Wmass, inv_Lbox, deltar, dx, dy, dz, cx, cy, cz, mass
real(8), dimension(:,:), allocatable :: COMXYZ
real(8), dimension(:),   allocatable :: MSDarray, amass

! 1. Unfolding the periodic box (sequential in time, parallel over atoms)
if (UNWRAPMD) then
   write(*,"(A)") " Unwrapping the simulation box "
   inv_Lbox = 1.d0 / Lbox

   do iframe = 1, nframes - 1
      frame_base =  iframe      * natom
      prev_base  = (iframe - 1) * natom

      !$omp parallel do default(none) &
      !$omp shared(natom, AtmPos, frame_base, prev_base, Lbox, inv_Lbox) &
      !$omp private(iatom, pos1, pos2, dx, dy, dz) &
      !$omp schedule(static)
      do iatom = 1, natom
         pos1 = frame_base + iatom
         pos2 = prev_base  + iatom

         dx = AtmPos(1, pos1) - AtmPos(1, pos2)
         dy = AtmPos(2, pos1) - AtmPos(2, pos2)
         dz = AtmPos(3, pos1) - AtmPos(3, pos2)

         AtmPos(1, pos1) = AtmPos(1, pos2) + (dx - ANINT(dx * inv_Lbox) * Lbox)
         AtmPos(2, pos1) = AtmPos(2, pos2) + (dy - ANINT(dy * inv_Lbox) * Lbox)
         AtmPos(3, pos1) = AtmPos(3, pos2) + (dz - ANINT(dz * inv_Lbox) * Lbox)
      end do
      !$omp end parallel do
   end do
endif

! Number of molecules in each frame
Nmol = natom / MolAtm

! Normalized atomic mass weights
allocate(amass(MolAtm), STAT=ios)
call AllocStat(ios, 'amass', 'MSDcom', 'A')

Wmass = SUM(AtomMass(1:MolAtm))
inv_Wmass = 1.d0 / Wmass
amass(1:MolAtm) = AtomMass(1:MolAtm) * inv_Wmass

! Allocate COM coordinates array
allocate(COMXYZ(3, nframes*Nmol), STAT=ios)
call AllocStat(ios, 'COMXYZ', 'MSDcom', 'A')

! 2. Computing the Center of Mass for each molecule in each frame (Parallelized)
write(*,"(/A/)") " Computing the Center of Mass for each molecule in each frame "

!$omp parallel do default(none) &
!$omp shared(nframes, Nmol, natom, MolAtm, amass, AtmPos, COMXYZ) &
!$omp private(iframe, frame_base, com_base, imol, mol_base, cx, cy, cz, iMolAtm, pos2, mass) &
!$omp schedule(static)
do iframe = 0, nframes - 1
   frame_base = iframe * natom
   com_base   = iframe * Nmol

   do imol = 1, Nmol
      mol_base = frame_base + (imol - 1) * MolAtm
      cx = 0.d0
      cy = 0.d0
      cz = 0.d0

      do iMolAtm = 1, MolAtm
         pos2 = mol_base + iMolAtm
         mass = amass(iMolAtm)
         cx = cx + AtmPos(1, pos2) * mass
         cy = cy + AtmPos(2, pos2) * mass
         cz = cz + AtmPos(3, pos2) * mass
      end do

      COMXYZ(1, com_base + imol) = cx
      COMXYZ(2, com_base + imol) = cy
      COMXYZ(3, com_base + imol) = cz
   end do
end do
!$omp end parallel do

max_lag = nframes
! max_lag capped at MAX_LAG_CAP for 8GB RAM / CPU safety)
if (max_lag > MAX_LAG_CAP) then
   write(*,'(A,I0,A,I0,A)') " Warning: max_lag (", max_lag, ") exceeds recommended cap. Capping at ", MAX_LAG_CAP, " frames."
   max_lag = MAX_LAG_CAP
end if

allocate(MSDarray(max_lag), STAT=ios)
call AllocStat(ios, 'MSDarray', 'MSDcom', 'A')
MSDarray(:) = 0.d0

write(*,"(A,I0,A)") " Computing the Mean Square Displacement up to ", max_lag, " frames "

!$omp parallel do default(none) &
!$omp shared(nframes, max_lag, Nmol, COMXYZ) &
!$omp private(frm1, frm2, tm, com_base1, com_base2, deltar, imol, dx, dy, dz) &
!$omp reduction(+:MSDarray) &
!$omp schedule(dynamic, 1)
do frm1 = 1, nframes
   com_base1 = (frm1 - 1) * Nmol

   do frm2 = frm1, min(nframes, frm1 + max_lag - 1)
      tm = frm2 - frm1 + 1
      com_base2 = (frm2 - 1) * Nmol
      deltar = 0.d0

      ! Vectorizable loop over molecules
      do imol = 1, Nmol
         dx = COMXYZ(1, com_base2 + imol) - COMXYZ(1, com_base1 + imol)
         dy = COMXYZ(2, com_base2 + imol) - COMXYZ(2, com_base1 + imol)
         dz = COMXYZ(3, com_base2 + imol) - COMXYZ(3, com_base1 + imol)
         deltar = deltar + dx*dx + dy*dy + dz*dz
      end do

      MSDarray(tm) = MSDarray(tm) + deltar
   end do
end do
!$omp end parallel do

! 4. Saving MSD into msd_com.dat file (using analytical pair count: nframes - tm + 1)
open(newunit=unit_id, file="msd_com.dat", status="unknown")
write(unit_id, "(A)") "#    Time/ps     MSD/nm^2"
do tm = 1, max_lag
   MSDarray(tm) = MSDarray(tm) / (DBLE(nframes - tm + 1) * DBLE(Nmol))
   write(unit_id, "(F12.3,1X,F12.4)") DBLE(tm - 1) * timestep, MSDarray(tm)
end do
close(unit_id)

! Releasing Memory
if (allocated(COMXYZ)) then
   deallocate(COMXYZ, STAT=ios)
   call AllocStat(ios, 'COMXYZ', 'MSDcom', 'D')
endif
if (allocated(MSDarray)) then
   deallocate(MSDarray, STAT=ios)
   call AllocStat(ios, 'MSDarray', 'MSDcom', 'D')
endif
if (allocated(amass)) then
   deallocate(amass, STAT=ios)
   call AllocStat(ios, 'amass', 'MSDcom', 'D')
endif
END SUBROUTINE MSDcom

! ********************************************************************************
! Dipolar Interaction Function from Phys. Chem. Chem. Phys., 2019, 21, 26621-26629
! ********************************************************************************
SUBROUTINE FIJmPCCP()
Use shared_vars
Use omp_lib
Use, intrinsic :: iso_c_binding
Implicit None
Include 'fftw3.f03'

integer :: frm, mcom, ios, unit_id
integer :: inuc, jnuc, posi, posj, iframe, Nmol, atom1, atom2, imol, jmol, ncycle
integer :: frame_base
real(8), dimension(:,:), allocatable :: meanRe, meanIm, IntraRe, IntraIm, ACFre_loc, ACFim_loc
complex(8), dimension(:,:), allocatable :: T2m_loc
complex(C_DOUBLE_COMPLEX), dimension(:), allocatable :: in_arr, out_arr, dummy_in, dummy_out
integer, dimension(:,:,:), allocatable :: rand_pairs
type(C_PTR) :: plan_fwd, plan_bwd
real(8) :: xij, yij, zij, r1, r2, inv_r5, fac0, fac1, fac2, Rmax2, inv_Lbox, time, randnb(2), conv_au

write(*,"(A)") " Autocorrelation functions computed using the Wiener-Khinchin Theorem "

! Constants necessary to compute F^ij_m
fac0     =  0.5d0
fac1     = -dsqrt(3.d0 / 2.d0)
fac2     = -fac1 * 0.5d0
Rmax2    =  Rmax * Rmax
inv_Lbox =  1.d0 / Lbox
conv_au  = (a0 * 1.d9)**6

! Number of molecules
Nmol = natom / MolAtm

! Allocate global output arrays
allocate(ACFre(3, nframes), STAT=ios)
call AllocStat(ios, 'ACFre', 'FIJmPCCP', 'A')
allocate(ACFim(3, nframes), STAT=ios)
call AllocStat(ios, 'ACFim', 'FIJmPCCP', 'A')
allocate(meanRe(3, nframes), STAT=ios)
call AllocStat(ios, 'meanRe', 'FIJmPCCP', 'A')
allocate(meanIm(3, nframes), STAT=ios)
call AllocStat(ios, 'meanIm', 'FIJmPCCP', 'A')
allocate(IntraRe(3, nframes), STAT=ios)
call AllocStat(ios, 'IntraRe', 'FIJmPCCP', 'A')
allocate(IntraIm(3, nframes), STAT=ios)
call AllocStat(ios, 'IntraIm', 'FIJmPCCP', 'A')

! Create FFTW plans ONCE before parallel regions 
allocate(dummy_in(fft_size), dummy_out(fft_size))
plan_fwd = fftw_plan_dft_1d(fft_size, dummy_in, dummy_out, FFTW_FORWARD, FFTW_ESTIMATE)
plan_bwd = fftw_plan_dft_1d(fft_size, dummy_out, dummy_in, FFTW_BACKWARD, FFTW_ESTIMATE)
deallocate(dummy_in, dummy_out)

! *******************
! Intramolecular ACFs
! *******************
IntraRe(:,:) = 0.d0
IntraIm(:,:) = 0.d0
meanRe(:,:)  = 0.d0
meanIm(:,:)  = 0.d0
ACFre(:,:)   = 0.d0
ACFim(:,:)   = 0.d0

if (nspin .eq. 1) then
   write(*,'(/A/)') " @ This system has no Intramolecular Dipolar ACFs"
else if (nspin .gt. 1) then
   write(*,"(/A/A/A)") &
      " ---------------------------------------------------- ", &
      " Dipolar relaxation from INTRA-molecular interactions ", &
      " ---------------------------------------------------- "
   write(*,'(/A)'  ) " @ Intramolecular Dipolar ACFs"
   write(*,"(A,I5)") "   Number of molecules: ", Nmol
   write(*,"(A,I8)") "   Total number of interacting pairs: ", Nmol * nspin * (nspin - 1) / 2

! Read previously calculated ACFs
   if (RDDACF) then
      open(newunit=unit_id, file='avg-acf-intra_wk.dat', status='old', action='read', iostat=ios)
      if (ios /= 0) call ERRMSG(unit_id, "Could not open file 'avg-acf-intra_wk.dat'")
      read(unit_id, '(A)')
      read(unit_id, '(A)')
      do frm = 1, nframes
         read(unit_id, *) time, (ACFre(mcom, frm), mcom=1, 3), (ACFim(mcom, frm), mcom=1, 3)
      end do
      close(unit_id)
! Calculate ACFs
   else
      !$omp parallel default(none) &
      !$omp shared(Nmol, nspin, MolAtm, atmlist, nframes, natom, maxmcom, AtmPos, fac0, fac1, fac2, fft_size, plan_fwd, plan_bwd) &
      !$omp private(imol, inuc, jnuc, atom1, atom2, iframe, frame_base, posi, posj, xij, yij, zij, r2, r1, inv_r5, &
      !$omp         T2m_loc, ACFre_loc, ACFim_loc, in_arr, out_arr) &
      !$omp reduction(+:meanRe, meanIm)

      allocate(T2m_loc(maxmcom, nframes), ACFre_loc(maxmcom, nframes), ACFim_loc(maxmcom, nframes))
      allocate(in_arr(fft_size), out_arr(fft_size))

      !$omp do schedule(static)
      !Loop over molecules
      do imol = 1, Nmol
         !Loop over spin-i
         do inuc = 1, nspin
            atom1 = MolAtm * (imol - 1) + atmlist(inuc)
            !Loop over spin-j
            do jnuc = inuc + 1, nspin
               atom2 = MolAtm * (imol - 1) + atmlist(jnuc)
               !Loop over frames
               do iframe = 0, nframes - 1
                  frame_base = iframe * natom
                  posi = frame_base + atom1
                  posj = frame_base + atom2

                  xij = AtmPos(1, posj) - AtmPos(1, posi)
                  yij = AtmPos(2, posj) - AtmPos(2, posi)
                  zij = AtmPos(3, posj) - AtmPos(3, posi)

                  r2     = xij*xij + yij*yij + zij*zij
                  r1     = dsqrt(r2)
                  inv_r5 = 1.d0 / (r2 * r2 * r1)

                  ! F^ij_m(t) functions:
                  ! T2m(1,...) = F^ij_0
                  ! T2m(2,...) = F^ij_1
                  ! T2m(3,...) = F^ij_2
                  ! Spherical harmonics Y2_m in Cartesian coordinates
                  T2m_loc(1, iframe + 1) = dcmplx(fac0 * (3.d0 * zij*zij - r2) * inv_r5, 0.d0)
                  T2m_loc(2, iframe + 1) = dcmplx(fac1 * zij * xij * inv_r5, fac1 * zij * yij * inv_r5)
                  T2m_loc(3, iframe + 1) = dcmplx(fac2 * (xij*xij - yij*yij) * inv_r5, fac2 * 2.d0 * xij * yij * inv_r5)
               end do

               ! ACFs of the atom pair
               call FFT_OMP(T2m_loc, ACFre_loc, ACFim_loc, plan_fwd, plan_bwd, in_arr, out_arr)

               ! Storing them to compute averages
               meanRe(:,:) = meanRe(:,:) + ACFre_loc(:,:)
               meanIm(:,:) = meanIm(:,:) + ACFim_loc(:,:)
            end do
         end do
      end do
      !$omp end do

      deallocate(T2m_loc, ACFre_loc, ACFim_loc, in_arr, out_arr)
      !$omp end parallel

      ! Average over equivalent pairs and different molecules 
      ! Convert to a.u. (a0^-6) in a single pass
      !
      ! The total number of intramolecular interaction pairs is (N/2)*(N-1),
      ! where N is the number of interacting spins within a given molecule;
      ! However, each nuclei will interact only with (N-1) neighbors;
      ! In a molecule, the average ACF is obtained by summing the ACFs computed
      ! for all "(N/2)*(N-1)" pairs and multiplying the resulting ACF by 2/N.
      !
      ! For a set of M molecules, the number of ACFs of this type is simply M.
      ! Therefore, the sum of all possible intramolecular ACFs, considering all molecules,
      ! must be multiplied by 2/(M*N) to obtain an average.
      ACFre(:,:) = meanRe(:,:) * (2.d0 / DBLE(Nmol * nspin)) * conv_au
      ACFim(:,:) = meanIm(:,:) * (2.d0 / DBLE(Nmol * nspin)) * conv_au

      ! Saving INTRAMOLECULAR ACFs within 'avg-acf-intra_wk.dat' file
      open(newunit=unit_id, file='avg-acf-intra_wk.dat', status='unknown')
      write(unit_id, '(A,1X,I5,1X,E20.10)') "# ", nframes, timestep
      write(unit_id, '(A,A)') &
         "#            Time/ps     Re(Fij_0)/a.u.      Re(Fij_1)/a.u.      Re(Fij_2)/a.u.      ", &
                                  "Im(Fij_0)/a.u.      Im(Fij_1)/a.u.      Im(Fij_2)/a.u."
      do frm = 1, nframes
         write(unit_id, '(7E20.10)') DBLE(frm - 1) * timestep, (ACFre(mcom, frm), mcom=1, 3), (ACFim(mcom, frm), mcom=1, 3)
      end do
      close(unit_id)
   endif

   IntraRe(:,:) = ACFre(:,:)
   IntraIm(:,:) = ACFim(:,:)
endif

! Spectral Densities & Relaxation Properties (Intramolecular)
acfisofn = 'avg-acf-intra_wk.iso'
call SpectDens()

NucProp = .FALSE.
if (NoInterACF) NucProp = .TRUE.
inttype = 1
call DDrelaxProp()

! Intermolecular relaxation will not be calculated
if (NoInterACF) then
   call fftw_destroy_plan(plan_fwd)
   call fftw_destroy_plan(plan_bwd)
   deallocate(meanRe, meanIm, IntraRe, IntraIm, ACFre, ACFim)
   return
endif

! *************************
! INTERmolecular Relaxation
! *************************
ACFre(:,:)  = 0.d0
ACFim(:,:)  = 0.d0
meanRe(:,:) = 0.d0
meanIm(:,:) = 0.d0

write(*,"(/A/A/A)") &
   " ---------------------------------------------------- ", &
   " Dipolar relaxation from INTER-molecular interactions ", &
   " ---------------------------------------------------- "
write(*,'(/A)'  ) " @ Intermolecular Dipolar ACFs"
write(*,"(A,I5)") "   Total number of molecules: ", Nmol
write(*,"(A,I8)") "   Total number of interacting pairs: ", (Nmol**2 - Nmol) * nspin**2

! Read previously calculated ACFs
if (RDDACF) then
   open(newunit=unit_id, file='avg-acf-inter_wk.dat', status='old', action='read', iostat=ios)
   if (ios /= 0) call ERRMSG(unit_id, "Could not open file 'avg-acf-inter_wk.dat'")
   read(unit_id, '(A)')
   read(unit_id, '(A)')
   do frm = 1, nframes
      read(unit_id, *) time, (ACFre(mcom, frm), mcom=1, 3), (ACFim(mcom, frm), mcom=1, 3)
   end do
   close(unit_id)
! Calculate ACFs
else
   if (AllInterPairs) then
      write(*,"(A,I10,A)") "   ACFs will be computed for ALL ", (Nmol * (Nmol - 1) / 2) * nspin**2, " intermolecular pairs"

      !$omp parallel default(none) &
      !$omp shared(Nmol, nspin, MolAtm, atmlist, nframes, natom, AtmPos, Lbox, inv_Lbox, Rmax2, &
      !$omp        fac0, fac1, fac2, maxmcom, fft_size, plan_fwd, plan_bwd) &
      !$omp private(imol, jmol, inuc, jnuc, atom1, atom2, ncycle, iframe, frame_base, posi, posj, &
      !$omp         xij, yij, zij, r2, r1, inv_r5, T2m_loc, ACFre_loc, ACFim_loc, in_arr, out_arr) &
      !$omp reduction(+:meanRe, meanIm)

      allocate(T2m_loc(maxmcom, nframes), ACFre_loc(maxmcom, nframes), ACFim_loc(maxmcom, nframes))
      allocate(in_arr(fft_size), out_arr(fft_size))

      !$omp do schedule(dynamic, 1)
      do imol = 1, Nmol - 1
         do jmol = imol + 1, Nmol
            ! Loop over all nspin x nspin intermolecular atom pairs
            do inuc = 1, nspin
               atom1 = MolAtm * (imol - 1) + atmlist(inuc)
               do jnuc = 1, nspin
                  atom2 = MolAtm * (jmol - 1) + atmlist(jnuc)

                  ncycle = 0
                  do iframe = 0, nframes - 1
                     frame_base = iframe * natom
                     posi = frame_base + atom1
                     posj = frame_base + atom2

                     ! Minimum Image Convention
                     xij = AtmPos(1, posj) - AtmPos(1, posi)
                     yij = AtmPos(2, posj) - AtmPos(2, posi)
                     zij = AtmPos(3, posj) - AtmPos(3, posi)

                     xij = xij - ANINT(xij * inv_Lbox) * Lbox
                     yij = yij - ANINT(yij * inv_Lbox) * Lbox
                     zij = zij - ANINT(zij * inv_Lbox) * Lbox

                     r2 = xij*xij + yij*yij + zij*zij

                     ! F^ij_m(t) is not calculated for spin pairs at distances greater than Rmax
                     if (r2 .gt. Rmax2) then
                        ncycle = ncycle + 1
                        T2m_loc(:, iframe + 1) = (0.d0, 0.d0)
                        cycle
                     end if

                     r1     = dsqrt(r2)
                     inv_r5 = 1.d0 / (r2 * r2 * r1)

                     T2m_loc(1, iframe + 1) = dcmplx(fac0 * (3.d0 * zij*zij - r2) * inv_r5, 0.d0)
                     T2m_loc(2, iframe + 1) = dcmplx(fac1 * zij * xij * inv_r5, fac1 * zij * yij * inv_r5)
                     T2m_loc(3, iframe + 1) = dcmplx(fac2 * (xij*xij - yij*yij) * inv_r5, fac2 * 2.d0 * xij * yij * inv_r5)
                  end do

                  if (ncycle .eq. nframes) cycle

                  ! ACFs of the atom pair
                  call FFT_OMP(T2m_loc, ACFre_loc, ACFim_loc, plan_fwd, plan_bwd, in_arr, out_arr)

                  ! Storing them to computes averages
                  meanRe(:,:) = meanRe(:,:) + ACFre_loc(:,:)
                  meanIm(:,:) = meanIm(:,:) + ACFim_loc(:,:)
               end do
            end do
         end do
      end do
      !$omp end do

      deallocate(T2m_loc, ACFre_loc, ACFim_loc, in_arr, out_arr)
      !$omp end parallel

      ! Averaging ACFs over equivalent pairs
      !
      ! The intermolecular interaction is calculated for ALL the M molecules in the box -that is,
      ! there are M ACFs of this type, which are accumulated so that an average can subsequently be calculated.
      ! If the molecule has more than one -say, N-spins involved in dipolar relaxation, then one has a total of
      ! M × N different ACFs in the sum of all calculated ACFs.
      ! The factor "2.d0" below appears because the loop over molecules-j starts at "molecule-i+1".
      !
      ! Exact normalization when all nspin^2 pairs are explicitly summed
      ACFre(:,:) = meanRe(:,:) * (2.d0 / DBLE(Nmol * nspin)) * conv_au
      ACFim(:,:) = meanIm(:,:) * (2.d0 / DBLE(Nmol * nspin)) * conv_au

   else
      write(*,"(A,I8,A)") "   ACFs will be computed for ", (Nmol * (Nmol - 1) / 2), " intermolecular pairs taken at random"

      ! Pre-generate independent random spin pairs (inuc, jnuc) per (imol, jmol) for deterministic OpenMP execution
      call init_random_seed()
      allocate(rand_pairs(2, Nmol, Nmol))
      do imol = 1, Nmol - 1
         do jmol = imol + 1, Nmol
            call random_number(randnb)
            rand_pairs(1, imol, jmol) = 1 + NINT((nspin - 1) * randnb(1))
            rand_pairs(2, imol, jmol) = 1 + NINT((nspin - 1) * randnb(2))
         end do
      end do

      !$omp parallel default(none) &
      !$omp shared(Nmol, MolAtm, atmlist, rand_pairs, nframes, natom, AtmPos, Lbox, inv_Lbox, Rmax2, &
      !$omp        fac0, fac1, fac2, maxmcom, fft_size, plan_fwd, plan_bwd) &
      !$omp private(imol, jmol, inuc, jnuc, atom1, atom2, ncycle, iframe, frame_base, posi, posj, &
      !$omp         xij, yij, zij, r2, r1, inv_r5, T2m_loc, ACFre_loc, ACFim_loc, in_arr, out_arr) &
      !$omp reduction(+:meanRe, meanIm)

      allocate(T2m_loc(maxmcom, nframes), ACFre_loc(maxmcom, nframes), ACFim_loc(maxmcom, nframes))
      allocate(in_arr(fft_size), out_arr(fft_size))

      !$omp do schedule(dynamic, 1)
      do imol = 1, Nmol - 1
         do jmol = imol + 1, Nmol
            inuc  = rand_pairs(1, imol, jmol)
            jnuc  = rand_pairs(2, imol, jmol)
            atom1 = MolAtm * (imol - 1) + atmlist(inuc)
            atom2 = MolAtm * (jmol - 1) + atmlist(jnuc)

            ncycle = 0
            do iframe = 0, nframes - 1
               frame_base = iframe * natom
               posi = frame_base + atom1
               posj = frame_base + atom2

               ! Minimum Image Convention
               xij = AtmPos(1, posj) - AtmPos(1, posi)
               yij = AtmPos(2, posj) - AtmPos(2, posi)
               zij = AtmPos(3, posj) - AtmPos(3, posi)

               xij = xij - ANINT(xij * inv_Lbox) * Lbox
               yij = yij - ANINT(yij * inv_Lbox) * Lbox
               zij = zij - ANINT(zij * inv_Lbox) * Lbox

               r2 = xij*xij + yij*yij + zij*zij

               if (r2 .gt. Rmax2) then
                  ncycle = ncycle + 1
                  T2m_loc(:, iframe + 1) = (0.d0, 0.d0)
                  cycle
               end if

               r1     = dsqrt(r2)
               inv_r5 = 1.d0 / (r2 * r2 * r1)

               T2m_loc(1, iframe + 1) = dcmplx(fac0 * (3.d0 * zij*zij - r2) * inv_r5, 0.d0)
               T2m_loc(2, iframe + 1) = dcmplx(fac1 * zij * xij * inv_r5, fac1 * zij * yij * inv_r5)
               T2m_loc(3, iframe + 1) = dcmplx(fac2 * (xij*xij - yij*yij) * inv_r5, fac2 * 2.d0 * xij * yij * inv_r5)
            end do

            if (ncycle .eq. nframes) cycle

            call FFT_OMP(T2m_loc, ACFre_loc, ACFim_loc, plan_fwd, plan_bwd, in_arr, out_arr)

            meanRe(:,:) = meanRe(:,:) + ACFre_loc(:,:)
            meanIm(:,:) = meanIm(:,:) + ACFim_loc(:,:)
         end do
      end do
      !$omp end do

      deallocate(T2m_loc, ACFre_loc, ACFim_loc, in_arr, out_arr)
      !$omp end parallel

      deallocate(rand_pairs)

      ! Statistical normalization when 1 random pair is sampled per molecule pair
      ACFre(:,:) = meanRe(:,:) * (2.d0 * DBLE(nspin) / DBLE(Nmol)) * conv_au
      ACFim(:,:) = meanIm(:,:) * (2.d0 * DBLE(nspin) / DBLE(Nmol)) * conv_au
   end if

   ! Saving INTER-MOLECULAR ACFs in a.u., i.e., a0^-6
   open(newunit=unit_id, file='avg-acf-inter_wk.dat', status='unknown')
   write(unit_id, '(A,1X,I5,1X,E20.10)') "# ", nframes, timestep
   write(unit_id, '(A,A)') &
      "#            Time/ps     Re(Fij_0)/a.u.      Re(Fij_1)/a.u.      Re(Fij_2)/a.u.      ", &
                               "Im(Fij_0)/a.u.      Im(Fij_1)/a.u.      Im(Fij_2)/a.u."
   do frm = 1, nframes
      write(unit_id, '(7E20.10)') DBLE(frm - 1) * timestep, (ACFre(mcom, frm), mcom=1, 3), (ACFim(mcom, frm), mcom=1, 3)
   end do
   close(unit_id)
endif

! Clean up FFTW plans
call fftw_destroy_plan(plan_fwd)
call fftw_destroy_plan(plan_bwd)

! Spectral Densities & Relaxation Properties (Intermolecular)
acfisofn = 'avg-acf-inter_wk.iso'
call SpectDens()

inttype = 2
SaveTau = GlmENR(1) / sigmaEFG(1)
call DDrelaxProp()

! ****************************************************************
! ALL (inter + intra)-molecular interactions
! ****************************************************************
write(*,"(/A/A/A)") &
   " ---------------------------------------------------------------- ", &
   " Dipolar relaxation from ALL (inter+intra)-molecular interactions ", &
   " ---------------------------------------------------------------- "

ACFre(:,:) = ACFre(:,:) + IntraRe(:,:)
ACFim(:,:) = ACFim(:,:) + IntraIm(:,:)

! Note: ACFre and IntraRe are already in a.u. (a0^-6)
open(newunit=unit_id, file='avg-acf_wk.dat', status='unknown')
write(unit_id, '(A,1X,I5,1X,E20.10)') "# ", nframes, timestep
write(unit_id, '(A,A)') &
   "#            Time/ps     Re(Fij_0)/a.u.      Re(Fij_1)/a.u.      Re(Fij_2)/a.u.      ", &
                            "Im(Fij_0)/a.u.      Im(Fij_1)/a.u.      Im(Fij_2)/a.u."
do frm = 1, nframes
   write(unit_id, '(7E20.10)') DBLE(frm - 1) * timestep, (ACFre(mcom, frm), mcom=1, 3), (ACFim(mcom, frm), mcom=1, 3)
end do
close(unit_id)

acfisofn = 'avg-acf_wk.iso'
call SpectDens()

NucProp = .TRUE.
inttype = 3
call DDrelaxProp()

! Freeing up memory space
deallocate(meanRe, STAT=ios)
call AllocStat(ios, 'meanRe', 'FIJmPCCP', 'D')
deallocate(meanIm, STAT=ios)
call AllocStat(ios, 'meanIm', 'FIJmPCCP', 'D')
deallocate(ACFre, STAT=ios)
call AllocStat(ios, 'ACFre', 'FIJmPCCP', 'D')
deallocate(ACFim, STAT=ios)
call AllocStat(ios, 'ACFim', 'FIJmPCCP', 'D')
deallocate(IntraRe, STAT=ios)
call AllocStat(ios, 'IntraRe', 'FIJmPCCP', 'D')
deallocate(IntraIm, STAT=ios)
call AllocStat(ios, 'IntraIm', 'FIJmPCCP', 'D')

END SUBROUTINE FIJmPCCP

! ******************************************************************************
! Dipolar Autocorrelation Function from Wiener-Khinchin Theorem
! ******************************************************************************
SUBROUTINE WK_ACF_DD()
Use shared_vars
Implicit None

integer :: frm, odd, mcom, ios, ISIGN, zpad
real(8), dimension(:), allocatable :: Rdata
real(8) :: denon

zpad=2*nframes

! Array used in FFT-routine 'FOUR1'
allocate(Rdata(2*zpad), STAT=ios)
call AllocStat(ios,'Rdata','WK_ACF_DD','A')

! Loop over the "m" component of F^{ij}_m(t): 
do mcom=1,maxmcom

! Computing the FT of R2m tensor
   do frm=1,nframes

! NR: Real      Part -> odd  numbers
!     Imaginary Part -> even numbers
      odd = 2*frm -1
      Rdata(odd)   =  DBLE(T2m(mcom,frm))   
      Rdata(odd+1) = DIMAG(T2m(mcom,frm))
   end do
   do frm=nframes+1,zpad
      odd = 2*frm -1
      Rdata(odd)   = 0.d0
      Rdata(odd+1) = 0.d0
   end do
   ISIGN=1
   call DFourierT(Rdata,zpad,ISIGN)
   
! ACF = inverseFFT( |FFT(R2m)|^2 ) / N*(N-1)  : WK
   do frm=1,zpad
      odd = 2*frm -1
      Rdata(odd)   = Rdata(odd)**2 +Rdata(odd+1)**2
      Rdata(odd+1) = 0.d0
   enddo
   ISIGN=-1
   call DFourierT(Rdata,zpad,ISIGN)

! Re[ ACF2m(mcom,n) ]:
! n = 0 (tau =    0)
! n = 1 (tau =   dt)
! n = 2 (tau =  2dt)
! nframes-1 = tf - ti

! Normalization factors associated with number of time origins
! ACF normalization
! for ACF(1) : divides by N
! for ACF(2) : divides by N-1
! for ACF(3) : divides by N-2
! ...
! for ACF(N) : divides by 1
   if (WK_divbyN) then
      denon = DBLE(nframes)
      do frm=1,nframes
         odd = 2*frm-1
         ACFre(mcom,frm) = Rdata(odd)  /denon
         ACFim(mcom,frm) = Rdata(odd+1)/denon
      end do
   else
      do frm=1,nframes
         odd = 2*frm-1
         denon = DBLE(nframes-frm+1)
         ACFre(mcom,frm) = Rdata(odd)  /denon
         ACFim(mcom,frm) = Rdata(odd+1)/denon
      end do
   endif
end do

! deallocating Rdata array
if (allocated(Rdata)) then
   deallocate(Rdata, STAT=ios)
   call AllocStat(ios,'Rdata','WK_ACF_DD','D')
endif
END SUBROUTINE WK_ACF_DD

! ******************************************************************************
! Dipolar Autocorrelation Function from Wiener-Khinchin Theorem
! FFT algorithm using the FFTW3 library
! Based on "https://github.com/Allen-Tildesley/examples/blob/master/corfun.f90"
! ******************************************************************************
SUBROUTINE FFT_ACF_DD()
Use shared_vars
Use, intrinsic :: iso_c_binding
Implicit none
Include 'fftw3.f03'  ! Contains FFTW constants and types

real(8) :: denon
complex(C_DOUBLE_COMPLEX), dimension(:), allocatable :: in_arr, out_arr
integer :: frm, mcom, ios
type(C_PTR) :: plan

! Array used in FFT-routine 'FOUR1'
allocate(in_arr(fft_size), STAT=ios)
call AllocStat(ios,'in_arr','FFT_ACF_DD','A')
allocate(out_arr(fft_size), STAT=ios)
call AllocStat(ios,'out_arr','FFT_ACF_DD','A')

! Loop over the "m" component of F^{ij}_m(t):
do mcom=1,maxmcom

! Computing the FT of R2m tensor
   in_arr = DCMPLX(0.0d0)
   do frm=1,nframes
      in_arr(frm) = cmplx(DBLE(T2m(mcom,frm)),DIMAG(T2m(mcom,frm)), kind=C_DOUBLE_COMPLEX)
   enddo

! Create a plan (FFTW_FORWARD = -1)
   plan = fftw_plan_dft_1d(fft_size, in_arr, out_arr, FFTW_FORWARD, FFTW_ESTIMATE) ! Set up plan
   CALL fftw_execute_dft ( plan, in_arr, out_arr )                       ! Execute FFT
   CALL fftw_destroy_plan ( plan )                                        ! Release plan

   out_arr = out_arr * DCONJG ( out_arr ) ! Square modulus

  ! Reverse FFT
   plan = fftw_plan_dft_1d ( fft_size, out_arr, in_arr, FFTW_BACKWARD, FFTW_ESTIMATE) ! Set up plan
   CALL fftw_execute_dft ( plan, out_arr, in_arr )                                   ! Execute FFT
   CALL fftw_destroy_plan ( plan )                                                    ! Release plan

   in_arr = in_arr / REAL ( fft_size ) ! Normalization factor associated with FFT itself

! Normalization factors associated with number of time origins 
! ACF normalization
! for ACF(1) : divides by N
! for ACF(2) : divides by N-1
! for ACF(3) : divides by N-2
! ...
! for ACF(N) : divides by 1
   if (WK_divbyN) then
      denon = DBLE(nframes)
      do frm=1,nframes
         ACFre(mcom,frm) =  DBLE(in_arr(frm))/denon
         ACFim(mcom,frm) = DIMAG(in_arr(frm))/denon
      end do
   else
      do frm=1,nframes
         denon = DBLE(nframes-frm+1)
         ACFre(mcom,frm) =  DBLE(in_arr(frm))/denon
         ACFim(mcom,frm) = DIMAG(in_arr(frm))/denon
      end do
   endif
end do

! deallocating Rdata array
if (allocated(in_arr)) then
   deallocate(in_arr, STAT=ios)
   call AllocStat(ios,'in_arr','FFT_ACF_DD','D')
endif
if (allocated(out_arr)) then
   deallocate(out_arr, STAT=ios)
   call AllocStat(ios,'out_arr','FFT_ACF_DD','D')
endif
END SUBROUTINE FFT_ACF_DD
! ******************************************************************************
! Dipolar Autocorrelation Function from Wiener-Khinchin Theorem
! FFT algorithm using the FFTW3 library
! Based on "https://github.com/Allen-Tildesley/examples/blob/master/corfun.f90"
! ******************************************************************************
SUBROUTINE FFT_OMP(T2m_loc, ACFre_loc, ACFim_loc, plan_fwd, plan_bwd, in_arr, out_arr)
 Use shared_vars
 Use, intrinsic :: iso_c_binding
 Implicit none
 Include 'fftw3.f03'
 
 complex(8), intent(in)  :: T2m_loc(maxmcom, nframes)
 real(8),    intent(out) :: ACFre_loc(maxmcom, nframes)
 real(8),    intent(out) :: ACFim_loc(maxmcom, nframes)
 type(C_PTR), intent(in) :: plan_fwd, plan_bwd
 complex(C_DOUBLE_COMPLEX), intent(inout) :: in_arr(fft_size), out_arr(fft_size)
 
 real(8) :: denon, inv_fft_size
 integer :: frm, mcom

inv_fft_size = 1.d0 / DBLE(fft_size)

do mcom = 1, maxmcom
   ! "Fill in the input vector and complete the zero-padding
   in_arr(1:nframes) = T2m_loc(mcom, 1:nframes)
   in_arr(nframes+1:fft_size) = (0.d0, 0.d0)

   ! Execute the forward FFT
   call fftw_execute_dft(plan_fwd, in_arr, out_arr)

   ! Magnitude squared
   out_arr(:) = out_arr(:) * DCONJG(out_arr(:))

   ! Execute the inverse FFT
   call fftw_execute_dft(plan_bwd, out_arr, in_arr)

   ! Double-precision normalization
   in_arr(:) = in_arr(:) * inv_fft_size

   ! Time-origin normalization
   if (WK_divbyN) then
      denon = 1.d0 / DBLE(nframes)
      ACFre_loc(mcom, 1:nframes) = DBLE(in_arr(1:nframes)) * denon
      ACFim_loc(mcom, 1:nframes) = DIMAG(in_arr(1:nframes)) * denon
   else
      do frm = 1, nframes
         denon = 1.d0 / DBLE(nframes -frm +1)
         ACFre_loc(mcom, frm) = DBLE(in_arr(frm)) * denon
         ACFim_loc(mcom, frm) = DIMAG(in_arr(frm)) * denon
      end do
   end if
end do
END SUBROUTINE FFT_OMP

! ******************************************************************************
! Dipolar of like spins relaxation properties
! ******************************************************************************
SUBROUTINE DDrelaxProp()
 USE shared_vars
 IMPLICIT NONE

 integer :: mcom
 real(8) :: tauDD(maxmcom), J2SI, Rfac, w0, w0_rad, gamma_rad, pi
 real(8) :: G12, G012, R1, R2, T1_ms, T2_ms

pi = 4.0d0 * atan(1.0d0)

write(*,'(/A/A/A/A)') &
     " ------------------------------------- ", &
     " Dipolar Variance and Correlation Time ", &
     " --------------------------------------", &
     "  m      <Fij^2>/a0^6            tau/ps"

! Larmor Frequency (MHz) - gammaI in MHz/T
w0 = B0ext * gammaI

! Correlation time calculation
do mcom = 1, maxmcom
   if (abs(sigmaEFG(mcom)) > 1.0d-30) then
      tauDD(mcom) = GlmENR(mcom) / sigmaEFG(mcom)
   else
      tauDD(mcom) = 0.0d0
   end if
   write(*,'(I3,E18.5,F18.5)') mvec(mcom), sigmaEFG(mcom), tauDD(mcom)
end do

! Spectral density combinations for Dipolar Relaxation (m = 0, 1, 2 -> indices 1, 2, 3)
G12  = 2.0d0 * (                           GlmRe(2) + 4.0d0 * GlmRe(3))
! For transversal relaxation, adapted from J. Phys. C: Solid State Phys., 14 (1981) 447-464
G012 = 2.0d0 * (1.5d0 * GlmRe(1) + 2.5d0 * GlmRe(2) +         GlmRe(3))

! Conversion of DD spectral density to SI units (ps/a0^6 -> s/m^6)
! gammaI in MHz/T -> rad / (T * s)
gamma_rad = gammaI * 1.0d6 * 2.0d0 * pi
J2SI      = 1.0d-12 / (a0**6)
Rfac      = J2SI * ((mu0 / (4.0d0 * pi))**2) * (gamma_rad**4) * (hbar**2) * (Ispin * (Ispin + 1.0d0))

R1 = Rfac * G12
R2 = Rfac * G012

if (R1 > 1.0d-30) then
   T1_ms = 1.0d3 / R1
else
   T1_ms = 0.0d0
end if

if (R2 > 1.0d-30) then
   T2_ms = 1.0d3 / R2
else
   T2_ms = 0.0d0
end if

write(*,"(/A/A/A)") &
     " -------------------------- ", &
     " Relaxation Times and Rates ", &
     " -------------------------- "
write(*,'(A,F13.3/A,F13.3/A,F13.5/A,F13.5)') &
     " T1 (ms) : ", T1_ms, &
     " T2 (ms) : ", T2_ms, &
     " R1 (Hz) : ", R1,    &
     " R2 (Hz) : ", R2

if (SaveR1R2) then
   SaveR1(inttype) = R1
   SaveR2(inttype) = R2
end if

if (NucProp) then
   w0_rad = w0 * 1.0d6 * 2.0d0 * pi
   write(*,"(/A/A/A)") &
        " ----------------------- ", &
        " Nuclear spin properties ", &
        " ------------------------ "
   write(*,'(A,14X,F5.1/A,F12.4/A,F12.4/A,E12.4/A,F12.4/A,E12.4)')  &
        " Nuclear spin quantum number: ", Ispin,     &
        " External magnetic field (T):        ", B0ext,     &
        " Magnetogyric ratio (MHz/T):         ", gammaI,    &
        " Magnetogyric ratio (rad/(T*s)):     ", gamma_rad, &
        " Larmor frequency (MHz):             ", w0,        &
        " Larmor angular frequency (rad/s):   ", w0_rad
end if
END SUBROUTINE DDrelaxProp

! ******************************************************************************
! Spin-Rotation relaxation properties
! ******************************************************************************
SUBROUTINE SRrelaxProp()
 USE shared_vars
 IMPLICIT NONE

 integer :: mcom
 real(8) :: tauSR(maxmcom), w0, w0_rad, gamma_rad, pi
 real(8) :: SR_rad, J_SI, R1, R2, T1_ms, T2_ms
 character(len=1), parameter :: xyz(3) = (/ 'x', 'y', 'z' /)

pi = 4.0d0 * atan(1.0d0)

! Larmor Frequency (MHz) - gammaI in MHz/T
w0 = B0ext * gammaI

write(*,'(/A/A/A/A)') &
     " ------------------------------------------- ", &
     " Spin-Rotation Variance and Correlation Time ", &
     " --------------------------------------------", &
     "  i      <J_i^2>/a.u.          tau_i/ps"

! Correlation time calculation
do mcom = 1, maxmcom
   if (abs(sigmaEFG(mcom)) > 1.0d-30) then
      tauSR(mcom) = GlmENR(mcom) / sigmaEFG(mcom)
   else
      tauSR(mcom) = 0.0d0
   end if
   write(*,'(A3,E18.5,F18.5)') xyz(mcom), sigmaEFG(mcom), tauSR(mcom)
end do

! Spin-Rotation relaxation rate calculation (preserving global arrays)
! SRtensor converted from kHz to rad/s; GlmRe converted from ps to s
R1 = 0.0d0
do mcom = 1, maxmcom
   SR_rad = SRtensor(mcom) * 1.0d3 * 2.0d0 * pi
   J_SI   = GlmRe(mcom) * 1.0d-12
   R1     = R1 + J_SI * (SR_rad**2)
end do
R1 = (2.0d0 / 3.0d0) * R1
R2 = R1

if (R1 > 1.0d-30) then
   T1_ms = 1.0d3 / R1
   T2_ms = 1.0d3 / R2
else
   T1_ms = 0.0d0
   T2_ms = 0.0d0
end if

write(*,"(/A/A/A)") &
     " -------------------------- ", &
     " Relaxation Times and Rates ", &
     " -------------------------- "
write(*,'(A,F13.3/A,F13.3/A,F13.5/A,F13.5)') &
     " T1 (ms) : ", T1_ms, &
     " T2 (ms) : ", T2_ms, &
     " R1 (Hz) : ", R1,    &
     " R2 (Hz) : ", R2

!if (SaveR1R2) then
!   SaveR1(inttype) = R1
!   SaveR2(inttype) = R2
!end if

if (NucProp) then
   gamma_rad = gammaI * 1.0d6 * 2.0d0 * pi
   w0_rad    = w0 * 1.0d6 * 2.0d0 * pi
   write(*,"(/A/A/A)") &
        " ----------------------- ", &
        " Nuclear spin properties ", &
        " ------------------------ "
   write(*,'(A,14X,F5.1/A,F12.4/A,F12.4/A,E12.4/A,F12.4/A,E12.4)')  &
        " Nuclear spin quantum number: ", Ispin,     &
        " External magnetic field (T):        ", B0ext,     &
        " Magnetogyric ratio (MHz/T):         ", gammaI,    &
        " Magnetogyric ratio (rad/(T*s)):     ", gamma_rad, &
        " Larmor frequency (MHz):             ", w0,        &
        " Larmor angular frequency (rad/s):   ", w0_rad
end if
END SUBROUTINE SRrelaxProp

! ******************************************************************************
! Intermolecular dipolar relaxation will be calculated using the
! Hwang-Freed diffusive model, as described in the reference
! J. Chem. Phys. 160, 074102 (2024); doi: 10.1063/5.0191052
! ******************************************************************************
SUBROUTINE HF_ACF()
 USE shared_vars
 IMPLICIT NONE

 integer :: Nmol, ipts, ios, npts_tr, ipts_tr, iframe, npts, mcom
 real(8) :: dij, pi, rhoJ, Dinf, w0, wu, x, x2, x4, x6, denom
 real(8) :: SINT, JnHFw(3), trtime, utime, fac, G0Re, inv_G0Re, c54_pi
 real(8) :: J2SI, Rfac, GT1, GT2, tauG, DeltaJn, tau_inf_eff, gamma_rad
 real(8), parameter :: zeta = 2.837297d0
 real(8), external :: SplineExtfunc
 real(8), allocatable :: y2(:), Xspl_loc(:), Fspl_loc(:), y2_loc(:)

pi     = 4.0d0 * atan(1.0d0)
c54_pi = 54.0d0 / pi

! Number of points used in numerical integration
npts = max(1000, nframes)

! Allocating memory to store spline arrays
allocate(Xspl(npts), STAT=ios)
call AllocStat(ios, 'Xspl', 'HF_ACF', 'A')
allocate(Fspl(npts), STAT=ios)
call AllocStat(ios, 'Fspl', 'HF_ACF', 'A')
allocate(y2(npts), STAT=ios)
call AllocStat(ios, 'y2', 'HF_ACF', 'A')
allocate(ACFre(3, nframes), STAT=ios)
call AllocStat(ios, 'ACFre', 'HF_ACF', 'A')

write(*,"(/A/A/A/A)") &
     " ---------------------------------------------------------------------------- ", &
     " Intermolecular dipolar relaxation will be calculated by combining MD results ", &
     " with the Hwang-Freed (HF) diffusive model (see: DOI 10.1063/5.0191052)       ", &
     " ---------------------------------------------------------------------------- "

! Number of molecules and numerical density of spins "i" (in nm^-3)
Nmol = natom / MolAtm
rhoJ = real(nspin * Nmol, kind=8) / (Lbox**3)

! Reading MD intermolecular isotropic ACFs from file 'avg-acf-inter_wk.iso'
open(20, file='avg-acf-inter_wk.iso', status='old')
read(20,*)
read(20,*)
do iframe = 1, nframes
   read(20,*) x, (ACFre(mcom, iframe), mcom = 1, 3)
end do
close(20)

! Converting ACF from a0^-6*ps to nm^-6*ps
ACFre(:,:) = ACFre(:,:) * (1.0d0 / ((a0 * 1.0d9)**6))

! Variance of Fij_m in nm^-6
G0Re     = ACFre(1, 1)
inv_G0Re = 1.0d0 / G0Re

! Distance of closest approach (nm)
dij = (4.0d0 * pi * rhoJ / (15.0d0 * G0Re))**(1.0d0 / 3.0d0)
write(*,"(/A,F10.5/)") " Distance of closest approach (Angs): ", dij * 10.0d0

! Size-independent inter-diffusion coefficient (Yeh-Hummer) (in nm^2/ps)
if (etaCoef == 0.0d0) then
   Dinf = 2.0d0 * DEcoef
   write(*,"(/A/A//A/A,F10.5,A/)") &
        " Warning: size independent inter-diffusion coefficient, D', cannot be estimated!", &
        " (missing shear viscosity, option '.ShearViscosity')", &
        " SPIN-RELAX will take D' = 2 D_E, where D_E is the self-diffusion coefficient from the COM MSD(t) ", &
        " D' = ", Dinf, " nm^2/ps"
else
   Dinf = 2.0d0 * (DEcoef + kBtz * Temper * zeta * 1.0d15 / (6.0d0 * pi * 1.0d-3 * etaCoef * Lbox))
   write(*,"(A,F10.5/)") " Size independent inter-diffusion coefficient (nm^2/ps): ", Dinf
end if

! ------------------------------------------------------------------------------------
! Normalized Hwang-Freed spectral densities J(0), J(w), J(2w) (in ps) for t in [0,inf]
! ------------------------------------------------------------------------------------
gamma_rad = gammaI * 1.0d6 * 2.0d0 * pi
w0        = B0ext * gamma_rad
fac       = 20.0d0 / real(npts, kind=8)

do mcom = 1,maxmcom
   wu = real(mcom - 1, kind=8) * w0 * (dij**2) * 1.0d-12 / Dinf
   do ipts = 1, npts
      x          = real(ipts, kind=8) * fac
      x2         = x * x
      x4         = x2 * x2
      x6         = x4 * x2
      Xspl(ipts) = x
      Fspl(ipts) = 1.0d0 / ((81.0d0 + 9.0d0 * x2 - 2.0d0 * x4 + x6) * (1.0d0 + (wu * wu) / x4))
   end do
   call y2spline(Xspl, Fspl, npts, y2)
   call simpson(SplineExtfunc, npts, Xspl, Fspl, y2, Xspl(1), Xspl(npts), SINT)
   JnHFw(mcom) = c54_pi * ((dij**2) / Dinf) * SINT
end do

! Trust time: tr (ps)
trtime = (Lbox**2) / (16.0d0 * pi * (2.0d0 * DEcoef))
write(*,"(A,F10.5,A/)") &
     " Time limit t_tr (ps): ", trtime, " (for t > t_tr, G^AB_MD(t) = G^AB_HF(t)) "

! Pre-determine number of frames up to trtime so the loop can be parallelized
npts_tr = nframes
do ipts_tr = 1, nframes
   if (timestep * real(ipts_tr, kind=8) > trtime) then
      npts_tr = ipts_tr
      exit
   end if
end do

! -----------------------------------------------------------------------------
! Parallel computation of Delta Gn = GnMD - GnHF over [0, t_tr]
! -----------------------------------------------------------------------------
!$omp parallel default(none) &
!$omp shared(npts_tr, npts, timestep, Dinf, dij, ACFre, inv_G0Re, c54_pi) &
!$omp private(ipts_tr, ipts, utime, fac, x, x2, x4, x6, denom, SINT, Xspl_loc, Fspl_loc, y2_loc)
allocate(Xspl_loc(npts), Fspl_loc(npts), y2_loc(npts))

! First, one needs GnHF for times t < t_tr, where GnMD was calculated
! loop over time
!$omp do schedule(static)
do ipts_tr = 1, npts_tr
   ! dimensionless time "u"
   utime = timestep * real(ipts_tr - 1, kind=8) * Dinf / (dij**2)

   ! controling the increment in "x" according to "u", where x = dij / rij
   if (ipts_tr > 1) then
      fac = sqrt(log(1.0d5) / utime) / real(npts, kind=8)
   else
      fac = 20.0d0 / real(npts, kind=8)
   end if

   ! building Gn_HF : Loop over x, with x in [0,inf]
   do ipts = 1, npts
      x              = real(ipts, kind=8) * fac
      x2             = x * x
      x4             = x2 * x2
      x6             = x4 * x2
      denom          = 81.0d0 + 9.0d0 * x2 - 2.0d0 * x4 + x6
      Xspl_loc(ipts) = x
      Fspl_loc(ipts) = x2 * exp(-utime * x2) / denom
   end do
   call y2spline(Xspl_loc, Fspl_loc, npts, y2_loc)
   call simpson(SplineExtfunc, npts, Xspl_loc, Fspl_loc, y2_loc, Xspl_loc(1), Xspl_loc(npts), SINT)
  !Delta Gn          = GnMD                         - GnHF
   ACFre(1, ipts_tr) = ACFre(1, ipts_tr) * inv_G0Re - SINT * c54_pi
end do
!$omp end do

deallocate(Xspl_loc, Fspl_loc, y2_loc)
!$omp end parallel

! -----------------------------------------------------------------------------
! Integration of Delta Gn to get Delta Jn and Delta tauG (in ps)
! -----------------------------------------------------------------------------
! Only the real portion of that integral is considered:
! "For random stationary motions, the correlation times flm (r) are real" (Spiess)
do ipts = 1, npts_tr
   Xspl(ipts) = timestep * real(ipts - 1, kind=8)
end do

! First integrate for m = 0 (mcom = 1): DeltaJn(0) and Delta_tauG are identical
Fspl(1:npts_tr) = ACFre(1, 1:npts_tr)
call y2spline(Xspl, Fspl, npts_tr, y2)
call simpson(SplineExtfunc, npts_tr, Xspl, Fspl, y2, Xspl(1), Xspl(npts_tr), SINT)

! Computing the full correlation time
!      tauHF                               + Delta tau
tauG = (4.0d0 * (dij**2)) / (9.0d0 * Dinf) + SINT

if (TauInf == 0.0d0) then
   tau_inf_eff = tauG
   write(*,"(/A/A/)") &
        " Note: correlation time for the infinity-size system was not given. ", &
        "       You may use the input option '.TauInf' to provide it.    "
else
   tau_inf_eff = TauInf
   write(*,"(A,F13.5/)") " Correlation time (ps) for the infinity-size system: ", tau_inf_eff
end if

! Apply finite-size correction for m = 0 : spectral density in ps/a0^6
JnHFw(1) = G0Re * (JnHFw(1) + (tau_inf_eff / tauG) * SINT)

! Now integrate for m = 1 and m = 2 (mcom = 2, 3)
do mcom = 2, maxmcom
   do ipts = 1, npts_tr
      Fspl(ipts) = ACFre(1, ipts) * cos(real(mcom - 1, kind=8) * w0 * Xspl(ipts) * 1.0d-12)
   end do
   call y2spline(Xspl, Fspl, npts_tr, y2)
   call simpson(SplineExtfunc, npts_tr, Xspl, Fspl, y2, Xspl(1), Xspl(npts_tr), DeltaJn)
   ! Apply finite-size correction : spectral density in ps/a0^6
   JnHFw(mcom) = G0Re * (JnHFw(mcom) + (tau_inf_eff / tauG) * DeltaJn)
end do

! -----------------------------------------------------------------------------
! Computing relaxation rates (in Hz)
! -----------------------------------------------------------------------------
J2SI = 1.0d-12 / ((1.0d-9)**6)
Rfac = J2SI * ((mu0 / (4.0d0 * pi))**2) * (gamma_rad**4) * (hbar**2) * (Ispin * (Ispin + 1.0d0))

GT1 = 2.0d0 * (                           JnHFw(2) + 4.0d0 * JnHFw(3))
GT2 = 2.0d0 * (1.5d0 * JnHFw(1) + 2.5d0 * JnHFw(2) +         JnHFw(3))

write(*,"(A/A/A)") &
     " -------------------------------------- ", &
     " Correlation Times and Relaxation Rates ", &
     " -------------------------------------- "
write(*,'(2(A,F13.5/)/5(A,F13.5/)/5(A,F13.5/))') &
     " tauG (inter,   ps)     : ", SaveTau, &
     " tauG (inter_HF,ps)     : ", tauG,    &
     " R1 (intra,         Hz) : ", SaveR1(1), &
     " R1 (inter,         Hz) : ", SaveR1(2), &
     " R1 (inter_HF,      Hz) : ", Rfac * GT1,  &
     " R1 (inter+intra,   Hz) : ", SaveR1(3), &
     " R1 (inter_HF+intra,Hz) : ", Rfac * GT1 + SaveR1(1), &
     " R2 (intra,         Hz) : ", SaveR2(1), &
     " R2 (inter,         Hz) : ", SaveR2(2), &
     " R2 (inter_HF,      Hz) : ", Rfac * GT2,  &
     " R2 (inter+intra,   Hz) : ", SaveR2(3), &
     " R2 (inter_HF+intra,Hz) : ", Rfac * GT2 + SaveR2(1)

! Freeing up memory space
deallocate(Xspl, STAT=ios)
call AllocStat(ios, 'Xspl', 'HF_ACF', 'D')
deallocate(Fspl, STAT=ios)
call AllocStat(ios, 'Fspl', 'HF_ACF', 'D')
deallocate(y2, STAT=ios)
call AllocStat(ios, 'y2', 'HF_ACF', 'D')
deallocate(ACFre, STAT=ios)
call AllocStat(ios, 'ACFre', 'HF_ACF', 'D')

END SUBROUTINE HF_ACF

! ******************************************************************************
! This external function uses spline interpolation to estimate f(arg) at 'arg'
! It's used in simpson's qsimp routine
! ******************************************************************************
FUNCTION SplineExtfunc(arg, ndim, Xspl, Fspl, y2) RESULT(SPL)
 IMPLICIT NONE
 integer, intent(in) :: ndim
 real(8), intent(in) :: arg, Xspl(*), Fspl(*), y2(*)
 real(8)             :: SPL

 call ncspline(arg, Xspl, Fspl, ndim, y2, SPL)
END FUNCTION SplineExtfunc

! --------------------------------------------------
! Seeding the random number generator (random_number)
! --------------------------------------------------
subroutine init_random_seed()
  implicit none
  integer :: i, n, clock, rate, cmax
  integer, allocatable :: seed(:)
  call random_seed(size = n)
  allocate(seed(n))
  call system_clock(count=clock, count_rate=rate, count_max=cmax)
  do i = 1, n
     seed(i) = clock + 37 * (i - 1)
  end do
  call random_seed(put = seed)
  write(*,*) 
  write(*,*) "Random number generator seed(s) : ",seed
  write(*,*) 
  deallocate(seed)
end subroutine init_random_seed

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
WRITE(*,'(/1X,2A/1X,2A/)') &
'Host name     : ',HSTNAM,  &
'Date and time : ',FFDATE
END SUBROUTINE LOCALINF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C GET ATOM CHARGE FROM ATOM NAME
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE GetAtomZ(Asymbol,Zatm)
     
INTEGER            :: iatm, Zatm
INTEGER, PARAMETER :: Natom=118
CHARACTER          :: Asymbol*2
     
CHARACTER(LEN=2), DIMENSION(Natom) :: Aname = (/ &
"H ","He", &
"Li","Be","B ","C ","N ","O ","F ","Ne", &
"Na","Mg","Al","Si","P ","S ","Cl","Ar", &
"K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn", &
    "Ga","Ge","As","Se", "Br","Kr", &
"Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd", &
     "In","Sn","Sb","Te","I ","Xe", &
"Cs","Ba", &
     "La","Ce","Pr","Nd","Pm","Sm","Eu", &
                    "Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", &
          "Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg", &
     "Tl","Pb","Bi","Po","At","Rn", &
"Fr","Ra", &
"Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",  &
"Rf","Db","Sg","Bh","Hs", &
"Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"/)

CHARACTER(LEN=2), DIMENSION(Natom) :: AnameUC= (/ &
"H ","HE", &
"LI","BE","B ","C ","N ","O ","F ","NE", &
"NA","MG","AL","SI","P ","S ","CL","AR", &
"K ","CA","SC","TI","V ","CR","MN","FE","CO","NI","CU","ZN", &
    "GA","GE","AS","SE", "BR","KR", &
"RB","SR","Y ","ZR","NB","MO","TC","RU","RH","PD","AG","CD", &
     "IN","SN","SB","TE","I ","XE", &
"CS","BA", &
     "LA","CE","PR","ND","PM","SM","EU", &
                    "GD","TB","DY","HO","ER","TM","YB","LU", &
          "HF","TA","W ","RE","OS","IR","PT","AU","HG", &
     "TL","PB","BI","PO","AT","RN", &
"FR","RA","AC","TH","PA","U ","NP","PU","AM","CM","BK","CF","ES","FM","MD","NO","LR",  &
"RF","DB","SG","BH","HS", &
"MT","DS","RG","CN","NH","FL","MC","LV","TS","OG" /)

Zatm=0
do iatm=1,Natom
   if ( Aname(iatm).EQ.Asymbol ) then
      Zatm = iatm
      return
   endif
end do
do iatm=1,Natom
   if ( AnameUC(iatm).EQ.Asymbol ) then
      Zatm = iatm
      return
   endif
end do
END SUBROUTINE GetAtomZ

!******************************************
!  Atomic masses and abundances taken from:
!
!  CRC Handbook of Chemistry and Physics
!   
!  Chief Editor:      David R. Lide  
!  Associate Editor:  H. P. R. Frederikse
!
!  78th edition; 1997-1998; CRC Press; USA
!******************************************
! ............................................................
! Computes atomic weights from isotopic masses and abundances
! ............................................................
SUBROUTINE GetAtomM(Zatm,AtomMass)

INTEGER            :: Zatm, niso, loop, Z
INTEGER, PARAMETER :: NMASS=906
character          :: emsg*80
real*8 :: AtomMass, fpZatm
REAL(8), DIMENSION(NMASS) :: AMASS = (/ &
!  Z           A. MASS           ABUNDANCE      
   1.D0,       1.007825D0,     99.985000D0,&   
   1.D0,       2.014102D0,      0.015000D0,&   
   1.D0,       3.016049D0,     -1.000000D0,&   
   2.D0,       4.002603D0,     99.999863D0,&   
   2.D0,       3.016029D0,      0.000137D0,&   
   3.D0,       7.016003D0,     92.500000D0,&   
   3.D0,       6.015121D0,      7.500000D0,&   
   4.D0,       9.012182D0,    100.000000D0,&   
   5.D0,      11.009305D0,     80.100000D0,&   
   5.D0,      10.012937D0,     19.900000D0,&   
   6.D0,      12.000000D0,     98.900000D0,&   
   6.D0,      13.003355D0,      1.100000D0,&   
   6.D0,      14.003242D0,     -1.000000D0,&   
   7.D0,      14.003074D0,     99.634000D0,&   
   7.D0,      15.000109D0,      0.366000D0,&   
   8.D0,      15.994915D0,     99.762000D0,&   
   8.D0,      17.999160D0,      0.200000D0,&   
   8.D0,      16.999131D0,      0.038000D0,&   
   9.D0,      18.998403D0,    100.000000D0,&   
  10.D0,      19.992436D0,     90.480000D0,&   
  10.D0,      21.991383D0,      9.250000D0,&   
  10.D0,      20.993843D0,      0.270000D0,&   
  11.D0,      22.989770D0,    100.000000D0,&   
  12.D0,      23.985042D0,     78.990000D0,&   
  12.D0,      25.982594D0,     11.010000D0,&   
  12.D0,      24.985837D0,     10.000000D0,&   
  13.D0,      26.981539D0,    100.000000D0,&   
  14.D0,      27.976927D0,     92.230000D0,&   
  14.D0,      28.976495D0,      4.670000D0,&   
  14.D0,      29.973771D0,      3.100000D0,&   
  15.D0,      30.973762D0,    100.000000D0,&   
  16.D0,      31.972071D0,     95.020000D0,&   
  16.D0,      33.967867D0,      4.210000D0,&   
  16.D0,      32.971459D0,      0.750000D0,&   
  16.D0,      35.967081D0,      0.020000D0,&   
  17.D0,      34.968853D0,     75.770000D0,&   
  17.D0,      36.965903D0,     24.230000D0,&   
  18.D0,      39.962384D0,     99.600000D0,&   
  18.D0,      35.967546D0,      0.337000D0,&   
  18.D0,      37.962733D0,      0.063000D0,&   
  19.D0,      38.963707D0,     93.258100D0,&   
  19.D0,      40.961825D0,      6.730200D0,&   
  19.D0,      39.963999D0,      0.011700D0,&   
  20.D0,      39.962591D0,     96.941000D0,&   
  20.D0,      43.955481D0,      2.086000D0,&   
  20.D0,      41.958618D0,      0.647000D0,&   
  20.D0,      47.952533D0,      0.187000D0,&   
  20.D0,      42.958767D0,      0.135000D0,&   
  20.D0,      45.953693D0,      0.004000D0,&   
  21.D0,      44.955910D0,    100.000000D0,&   
  22.D0,      47.947947D0,     73.800000D0,&   
  22.D0,      45.952629D0,      8.000000D0,&   
  22.D0,      46.951764D0,      7.300000D0,&   
  22.D0,      48.947871D0,      5.500000D0,&   
  22.D0,      49.944792D0,      5.400000D0,&   
  23.D0,      50.943962D0,     99.750000D0,&   
  23.D0,      49.947161D0,      0.250000D0,&   
  24.D0,      51.940510D0,     83.789000D0,&   
  24.D0,      52.940651D0,      9.501000D0,&   
  24.D0,      49.946046D0,      4.345000D0,&   
  24.D0,      53.938883D0,      2.365000D0,&   
  25.D0,      54.938047D0,    100.000000D0,&   
  26.D0,      55.934939D0,     91.720000D0,&   
  26.D0,      53.939613D0,      5.800000D0,&   
  26.D0,      56.935396D0,      2.100000D0,&   
  26.D0,      57.933277D0,      0.280000D0,&   
  27.D0,      58.933198D0,    100.000000D0,&   
  28.D0,      57.935346D0,     68.077000D0,&   
  28.D0,      59.930788D0,     26.223000D0,&   
  28.D0,      61.928346D0,      3.634000D0,&   
  28.D0,      60.931058D0,      1.140000D0,&   
  28.D0,      63.927968D0,      0.926000D0,&   
  29.D0,      62.929599D0,     69.170000D0,&   
  29.D0,      64.927793D0,     30.830000D0,&   
  30.D0,      63.929145D0,     48.600000D0,&   
  30.D0,      65.926035D0,     27.900000D0,&   
  30.D0,      67.924846D0,     18.800000D0,&   
  30.D0,      66.927129D0,      4.100000D0,&   
  30.D0,      69.925325D0,      0.600000D0,&   
  31.D0,      68.925580D0,     60.108000D0,&   
  31.D0,      70.924701D0,     39.892000D0,&   
  32.D0,      73.921177D0,     35.940000D0,&   
  32.D0,      71.922080D0,     27.660000D0,&   
  32.D0,      69.924250D0,     21.230000D0,&   
  32.D0,      72.923463D0,      7.730000D0,&   
  32.D0,      75.921402D0,      7.440000D0,&   
  33.D0,      74.921594D0,    100.000000D0,&   
  34.D0,      79.916520D0,     49.610000D0,&   
  34.D0,      77.917308D0,     23.780000D0,&   
  34.D0,      75.919212D0,      9.360000D0,&   
  34.D0,      81.916698D0,      8.730000D0,&   
  34.D0,      76.919913D0,      7.630000D0,&   
  34.D0,      73.922475D0,      0.890000D0,&   
  35.D0,      78.918336D0,     50.690000D0,&   
  35.D0,      80.916289D0,     49.310000D0,&   
  36.D0,      83.911507D0,     57.000000D0,&   
  36.D0,      85.910616D0,     17.300000D0,&   
  36.D0,      81.913482D0,     11.600000D0,&   
  36.D0,      82.914135D0,     11.500000D0,&   
  36.D0,      79.916380D0,      2.250000D0,&   
  36.D0,      77.920396D0,      0.350000D0,&   
  37.D0,      84.911794D0,     72.165000D0,&   
  37.D0,      86.909187D0,     27.835000D0,&   
  38.D0,      87.905619D0,     82.580000D0,&   
  38.D0,      85.909267D0,      9.860000D0,&   
  38.D0,      86.908884D0,      7.000000D0,&   
  38.D0,      83.913430D0,      0.560000D0,&   
  39.D0,      88.905849D0,    100.000000D0,&   
  40.D0,      89.904703D0,     51.450000D0,&   
  40.D0,      93.906315D0,     17.380000D0,&   
  40.D0,      91.905039D0,     17.150000D0,&   
  40.D0,      90.905644D0,     11.220000D0,&   
  40.D0,      95.908275D0,      2.800000D0,&   
  41.D0,      92.906377D0,    100.000000D0,&   
  42.D0,      97.905407D0,     24.130000D0,&   
  42.D0,      95.904679D0,     16.680000D0,&   
  42.D0,      94.905841D0,     15.920000D0,&   
  42.D0,      91.906808D0,     14.840000D0,&   
  42.D0,      99.907477D0,      9.630000D0,&   
  42.D0,      96.906021D0,      9.550000D0,&   
  42.D0,      93.905085D0,      9.250000D0,&   
  43.D0,      97.907215D0,     -1.000000D0,&   
  44.D0,     101.904349D0,     31.600000D0,&   
  44.D0,     103.905424D0,     18.700000D0,&   
  44.D0,     100.905582D0,     17.000000D0,&   
  44.D0,      98.905939D0,     12.700000D0,&   
  44.D0,      99.904219D0,     12.600000D0,&   
  44.D0,      95.907599D0,      5.520000D0,&   
  44.D0,      97.905287D0,      1.880000D0,&   
  45.D0,     102.905500D0,    100.000000D0,&   
  46.D0,     105.903478D0,     27.330000D0,&   
  46.D0,     107.903895D0,     26.460000D0,&   
  46.D0,     104.905079D0,     22.330000D0,&   
  46.D0,     109.905167D0,     11.720000D0,&   
  46.D0,     103.904029D0,     11.140000D0,&   
  46.D0,     101.905634D0,      1.020000D0,&   
  47.D0,     106.905092D0,     51.839000D0,&   
  47.D0,     108.904757D0,     48.161000D0,&   
  48.D0,     113.903357D0,     28.730000D0,&   
  48.D0,     111.902758D0,     24.130000D0,&   
  48.D0,     110.904182D0,     12.800000D0,&   
  48.D0,     109.903005D0,     12.490000D0,&   
  48.D0,     112.904400D0,     12.220000D0,&   
  48.D0,     115.904754D0,      7.490000D0,&   
  48.D0,     105.906461D0,      1.250000D0,&   
  48.D0,     107.904176D0,      0.890000D0,&   
  49.D0,     114.903880D0,     95.700000D0,&   
  49.D0,     112.904061D0,      4.300000D0,&   
  50.D0,     119.902199D0,     32.590000D0,&   
  50.D0,     117.901609D0,     24.230000D0,&   
  50.D0,     115.901747D0,     14.530000D0,&   
  50.D0,     118.903310D0,      8.590000D0,&   
  50.D0,     116.902956D0,      7.680000D0,&   
  50.D0,     123.905274D0,      5.790000D0,&   
  50.D0,     121.903440D0,      4.630000D0,&   
  50.D0,     111.904826D0,      0.970000D0,&   
  50.D0,     113.902784D0,      0.650000D0,&   
  50.D0,     114.903348D0,      0.340000D0,&   
  51.D0,     120.903821D0,     57.360000D0,&   
  51.D0,     122.904216D0,     42.640000D0,&   
  52.D0,     129.906229D0,     33.800000D0,&   
  52.D0,     127.904463D0,     31.690000D0,&   
  52.D0,     125.903314D0,     18.950000D0,&   
  52.D0,     124.904433D0,      7.139000D0,&   
  52.D0,     123.902823D0,      4.816000D0,&   
  52.D0,     121.903054D0,      2.603000D0,&   
  52.D0,     122.904271D0,      0.908000D0,&   
  52.D0,     119.904048D0,      0.096000D0,&   
  53.D0,     126.904473D0,    100.000000D0,&   
  54.D0,     131.904144D0,     26.900000D0,&   
  54.D0,     128.904780D0,     26.400000D0,&   
  54.D0,     130.905072D0,     21.200000D0,&   
  54.D0,     133.905395D0,     10.400000D0,&   
  54.D0,     135.907219D0,      8.900000D0,&   
  54.D0,     129.903509D0,      4.100000D0,&   
  54.D0,     127.903531D0,      1.910000D0,&   
  54.D0,     123.905894D0,      0.100000D0,&   
  54.D0,     125.904281D0,      0.090000D0,&   
  55.D0,     132.905429D0,    100.000000D0,&   
  56.D0,     137.905232D0,     71.700000D0,&   
  56.D0,     136.905812D0,     11.230000D0,&   
  56.D0,     135.904553D0,      7.854000D0,&   
  56.D0,     134.905665D0,      6.592000D0,&   
  56.D0,     133.904486D0,      2.417000D0,&   
  56.D0,     131.905042D0,      0.101000D0,&   
  56.D0,     129.906282D0,     -2.000000D0,&   
  57.D0,     138.906347D0,     99.909800D0,&   
  57.D0,     137.907105D0,      0.090200D0,&   
  58.D0,     139.905433D0,     88.480000D0,&   
  58.D0,     141.909241D0,     11.080000D0,&   
  58.D0,     137.905985D0,      0.250000D0,&   
  58.D0,     135.907140D0,      0.190000D0,&   
  59.D0,     140.907647D0,    100.000000D0,&   
  60.D0,     141.907719D0,     27.130000D0,&   
  60.D0,     143.910083D0,     23.800000D0,&   
  60.D0,     145.913113D0,     17.190000D0,&   
  60.D0,     142.909810D0,     12.180000D0,&   
  60.D0,     144.912570D0,      8.300000D0,&   
  60.D0,     147.916889D0,      5.760000D0,&   
  60.D0,     149.920887D0,      5.640000D0,&   
  61.D0,     144.912743D0,     -1.000000D0,&   
  62.D0,     151.919729D0,     26.700000D0,&   
  62.D0,     153.922206D0,     22.700000D0,&   
  62.D0,     146.914895D0,     15.000000D0,&   
  62.D0,     148.917181D0,     13.800000D0,&   
  62.D0,     147.914819D0,     11.300000D0,&   
  62.D0,     149.917273D0,      7.400000D0,&   
  62.D0,     143.911998D0,      3.100000D0,&   
  63.D0,     152.921225D0,     52.200000D0,&   
  63.D0,     150.919847D0,     47.800000D0,&   
  64.D0,     157.924019D0,     24.840000D0,&   
  64.D0,     159.927049D0,     21.860000D0,&   
  64.D0,     155.922118D0,     20.470000D0,&   
  64.D0,     156.923956D0,     15.650000D0,&   
  64.D0,     154.922618D0,     14.800000D0,&   
  64.D0,     153.920861D0,      2.180000D0,&   
  64.D0,     151.919786D0,      0.200000D0,&   
  65.D0,     158.925342D0,    100.000000D0,&   
  66.D0,     163.929171D0,     28.200000D0,&   
  66.D0,     161.926795D0,     25.500000D0,&   
  66.D0,     162.928728D0,     24.900000D0,&   
  66.D0,     160.926930D0,     18.900000D0,&   
  66.D0,     159.925193D0,      2.340000D0,&   
  66.D0,     157.924403D0,      0.100000D0,&   
  66.D0,     155.924277D0,      0.060000D0,&   
  67.D0,     164.930319D0,    100.000000D0,&   
  68.D0,     165.930290D0,     33.600000D0,&   
  68.D0,     167.932368D0,     26.800000D0,&   
  68.D0,     166.932046D0,     22.950000D0,&   
  68.D0,     169.935461D0,     14.900000D0,&   
  68.D0,     163.929198D0,      1.610000D0,&   
  68.D0,     161.928775D0,      0.140000D0,&   
  69.D0,     168.934212D0,    100.000000D0,&   
  70.D0,     173.938859D0,     31.800000D0,&   
  70.D0,     171.936378D0,     21.900000D0,&   
  70.D0,     172.938208D0,     16.120000D0,&   
  70.D0,     170.936323D0,     14.300000D0,&   
  70.D0,     175.942564D0,     12.700000D0,&   
  70.D0,     169.934759D0,      3.050000D0,&   
  70.D0,     167.933894D0,      0.130000D0,&   
  71.D0,     174.940770D0,     97.410000D0,&   
  71.D0,     175.942679D0,      2.590000D0,&   
  72.D0,     179.946546D0,     35.100000D0,&   
  72.D0,     177.943696D0,     27.297000D0,&   
  72.D0,     176.943217D0,     18.606000D0,&   
  72.D0,     178.945812D0,     13.629000D0,&   
  72.D0,     175.941406D0,      5.206000D0,&   
  72.D0,     173.940044D0,      0.162000D0,&   
  73.D0,     180.947992D0,     99.988000D0,&   
  73.D0,     179.947462D0,      0.012000D0,&   
  74.D0,     183.950928D0,     30.670000D0,&   
  74.D0,     185.954357D0,     28.600000D0,&   
  74.D0,     181.948202D0,     26.300000D0,&   
  74.D0,     182.950220D0,     14.300000D0,&   
  74.D0,     179.946701D0,      0.130000D0,&   
  75.D0,     186.955744D0,     62.600000D0,&   
  75.D0,     184.952951D0,     37.400000D0,&   
  76.D0,     191.961467D0,     41.000000D0,&   
  76.D0,     189.958436D0,     26.400000D0,&   
  76.D0,     188.958137D0,     16.100000D0,&   
  76.D0,     187.955860D0,     13.300000D0,&   
  76.D0,     186.955741D0,      1.600000D0,&   
  76.D0,     185.953830D0,      1.580000D0,&   
  76.D0,     183.952488D0,      0.020000D0,&   
  77.D0,     192.962917D0,     62.700000D0,&   
  77.D0,     190.960584D0,     37.300000D0,&   
  78.D0,     194.964766D0,     33.800000D0,&   
  78.D0,     193.962655D0,     32.900000D0,&   
  78.D0,     195.964926D0,     25.300000D0,&   
  78.D0,     197.967869D0,      7.200000D0,&   
  78.D0,     191.961019D0,      0.790000D0,&   
  78.D0,     189.959917D0,      0.010000D0,&   
  79.D0,     196.966543D0,    100.000000D0,&   
  80.D0,     201.970617D0,     29.860000D0,&   
  80.D0,     199.968300D0,     23.100000D0,&   
  80.D0,     198.968254D0,     16.870000D0,&   
  80.D0,     200.970277D0,     13.180000D0,&   
  80.D0,     197.966743D0,      9.970000D0,&   
  80.D0,     203.973467D0,      6.870000D0,&   
  80.D0,     195.965807D0,      0.150000D0,&   
  81.D0,     204.974401D0,     70.476000D0,&   
  81.D0,     202.972320D0,     29.524000D0,&   
  82.D0,     207.976627D0,     52.400000D0,&   
  82.D0,     205.974440D0,     24.100000D0,&   
  82.D0,     206.975872D0,     22.100000D0,&   
  82.D0,     203.973020D0,      1.400000D0,&   
  83.D0,     208.980374D0,    100.000000D0,&   
  84.D0,     208.982404D0,     -1.000000D0,&   
  85.D0,     209.987126D0,     -1.000000D0,&   
  86.D0,     222.017570D0,     -1.000000D0,&   
  87.D0,     223.019733D0,     -1.000000D0,&   
  88.D0,     226.025402D0,     -1.000000D0,&   
  89.D0,     227.027750D0,     -1.000000D0,&   
  90.D0,     232.038054D0,    100.000000D0,&   
  91.D0,     231.035880D0,     -1.000000D0,&   
  92.D0,     238.050785D0,     99.274500D0,&   
  92.D0,     235.043924D0,      0.720000D0,&   
  92.D0,     234.040947D0,      0.005500D0,&   
  93.D0,     237.048168D0,     -1.000000D0,&   
  94.D0,     239.052157D0,     -1.000000D0,&   
  94.D0,     244.064199D0,     -1.000000D0,&   
 119.D0,       1.D+06,        100.000000D0 /)

AtomMass=0.d0
fpZatm = DBLE(Zatm)

niso=0
do loop=1,302
   Z = 3*loop-2
   if (fpZatm.EQ.AMASS(Z) .and. AMASS(Z+2).gt.0.d0) then
      niso = niso +1
      AtomMass = AtomMass +AMASS(Z+1)*AMASS(Z+2)/100.d0
    endif
enddo
if (niso.eq.0) then
    write(emsg,'(A,I3)') "Failed to get Atomic Mass for the atomic number Z=",Zatm
    call ERRMSG(01,emsg)
endif
END SUBROUTINE GetAtomM
