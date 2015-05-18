# Jacobi Solver

<img src="http://appgallery.maxeler.com/v0.1/app/Jacobi%20Solver/icon" alt="Jacobi Solver">

## Description

This App implements a solver for equations of the type Ax=b, where A is constant but where we have a set of possible b. Both the matrix A and the vectors b are randomly generated. This App always runs both the DFE and the CPU implementation and compares respective times taken.

####This App requires [maxpower](https://github.com/maxeler/maxpower) Maxeler standard library.####

## Content

The repo root directory contains the following items:

- APP
- LICENCE.txt

### APP

Directory containing project sources.
  
### LICENSE.txt

License of the project.

## Information to compile

Ensure the environment variables below are correctly set:
  * `MAXELEROSDIR`
  * `MAXCOMPILERDIR`

To compile the application, run:

    make RUNRULE="<ProfileName>"

If would like to remove the distributed maxfiles before recompiling the application run the following command before compilation:

    make RUNRULE="<ProfileName>" distclean

## Makefile targets

### build  

Compiles the application

### clean  

Removes results of compilation from build directories  

### distclean  

Removes all results of comakempilation from build directories, including all maxfiles

Jacobi Solver on [AppGallery](http://appgallery.maxeler.com/)   

