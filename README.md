# Non-Central Chi Sqare

<img src="http://appgallery.maxeler.com/v0.1/app/Non-Central%20Chi%20Square/icon" alt="Non-Central Chi Sqare">

## Description

This App generates lots of sample from complex distribution.

####This App requires maxjlibs Maxeler library.####

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

Non-Central Chi Sqare on [AppGallery](http://appgallery.maxeler.com/)   

