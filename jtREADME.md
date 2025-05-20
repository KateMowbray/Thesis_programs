# 2016-01-27 party_orb code (JT)

I have modified the non-relativistic code of Paolo Guiliani and others (documented in Guiliani et al. ApJ 635:pp636-646, 2005) to be relativistic, use a range of input environments and use simple Makefile.

See manual (./manual/) for full details of source code and IDL widgets.

Careful on the Makefile flags:
On my machine (oldjock.mcs.st-and.ac.uk) code compiles against older version of gfortran (gfortran44). 
Several debugging options have also been commented out - including these can severly slow down run-time.

To prepare the code for use, type "make" in the root directory.
This creates several additional directories:
./Data/	::	location of output data
./mod/	:: 	location of *.mod files
./obj/	::	location of *.o files
./bin/	::	binary file location

(the contents of these folders are ignored by git as some will contain binary info or large datafiles)

To run the code, simply type "./bin/rparty" on the command line (for the relativistic version).

"make clean" 	removes all these extra directories ready for recompilation.
"make datatidy"	deletes all the files in the ./Data/ directory.

The manual also details how to switch between setups. Sorry, you'll have to compile the pdf file yourself though.

To load lots of useful IDL routines to process the data, load the Start.pro file on starting IDL, with the command
idl Start.pro

Oh, the remaining subfolders contain stuff carried over from previous experiments, following the work of P. Giuliani, K. Grady, S. Oskoui, T. Neukirch, and perhaps others.
./text/	contains snippets of other peoples work, which I'm just holding onto for now in case they become important later.
Neither folder is affected by the Makefile.

2016 saw the addition of several plugins for different datafiles from Lare, in either 2D (but pretending to be 3D for the particles) or 3D data snapshots.
These features are still being tested, and should NOT be taken for granted.

Good luck - any problems, email me at jwt9@st-andrews.ac.uk.. Ciao!
