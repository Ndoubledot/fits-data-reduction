# fits-data-reduction
make is a python script written to clean astronomical data that is fits file, fits is flexible image transport system which means briefly that fits file is image taken from telescopes using a CCD. CCD device can be thought of as an matrix array made up of small pixels, in similar sense fits file is data collected from those pixels which can be thought of as an matrix itself.

Make can be a useful tool to subtract biase, dark substraction, flat fielding, cosmic ray correction and stacking of fits files.

make has various functions defined as masterbias, masterflat, masterdark and so on which take lists of ccd data and return a single cleaned ccd data

A thorough example is given which can walk the reader through the make functions

the raw data for the example has been provided by Devanand Ullas, which is taken using Devastal optical Telescope

Thank you so much devanand for the data
