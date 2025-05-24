
The QRELAX code is composed of a single file: 

i)   qrelax.f90: program written in Fortran 90 language

To generate the executable (qrelax.x), edit the BASH script compile.sh, 
inserting the Fortran compiler and compilation flags of choice, then type: 

% ./compile.sh

The QRELAX input file is called "qrelax.inp," and its available options are 
contained in the subroutine "read_input()." 

Here, we provide a sample input file with several available options commented out. 
We also included the "ORCA" directory containing 400 frames from trajectory 2 
of the article, used to generate the last line of Table 2. 
Running the program will initially extract ALL the EFGs available in the "framexxxx" 
files and insert them into the EFG.xyz file. 
Thus, if it is necessary to run the program again, the EFGs will be read from the 
EFG.xyz file, making the calculation of the relaxation properties much faster. 

First type

% ./qrelax.x 

to perform the commented extraction operation. By doing so, the R1 and R2 rates 
shown in Table 2 are calculated. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   R1      R2      R1^EN  s^2_iso  tau_iso
% $\omega$B97MV/SCNL/ma-def2-TZVP & 3.105 & 3.190 & 3.159 & 2.521 & 5.631 \\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

To generate R1^EN, s^2_iso, and tau_iso, uncomment the "ForceIsotropy" option 
and replace the Larmor frequency of 500 MHz with zero. 
Then, rerun the program. 

Please let me know if you experience any problems:

(Luciano N. Vidal, lnvidal@utfpr.edu.br).

