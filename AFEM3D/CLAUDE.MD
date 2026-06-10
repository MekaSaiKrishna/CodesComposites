# Instructions 

You are a software engineer proficient in FORTRAN and ABAQUS. I need you to go through the '.for' files and have them commented with brevity and separate out the code chunks to make them easily comprehensible by any non-coder.

The codes in this folder for an ABAQUS subroutine to implement 'incremental' progressive damage using the smeared crack approach for composite materials. 

Currently the codes are written for an isotropic material. 

Once it is tested on various ABAQUS simulations of a test block under tension, compression, and shear it should be deemed ready to be used for general purpose simulations.

The element type in ABAQUS for which this is compatible is C3D8.

Always update the 'README.MD' file which should be a summary of what all the fortran files do and which one is the latest version that is currently in use.

The test-cases in ABAQUS for which these fortran files will be tested are written in 'CHECKLIST.MD'

For each of the test cases things that we will be checking is, the amount of time it takes for the simulation to run, does the peak stress match with the max-stress that the material property is defined as, the failure strain, and the area under the stress vs strain curve. Best if all of these are tabulated in the same order as the CHECKLIST.MD file



