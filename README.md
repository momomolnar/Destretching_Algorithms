# Destretching_Algorithms
Compilation of routines for removing (seeing induced) optical distortions, based on local correlation tracking (LCT).

To-do list (8/17/21):
  1) Make the xcorrelation box size to be independent parameter from the grid spacing;
  7) Make a function to calculate the destretching on a user-defined (arbitraty geometry) grid; 
    - Find a function to interpolate the arbitrary destretching geometry;
  9)  Break down the definition *doref* to be able to take a list of _arbitrary_ grid points;
  10) Make a control point index that contains all the subfield Control Point locations and remove;
  the extra dimension in all the fft routines;
  12) [X] Change the interpolation function with a modern implementation (faster);


  4) Insert a white noise estimator to see when the box size is small enough; 
  5) A box size estimator based on destretching vector correlation coefficient; 
  6) Add a function for (limb) observations, to zero destretching vector, as there is no signal in the box;
 
  11) [X] Make a jupyter notebook with example of the code;
  8) Add napari built-in support to be able to show the result(ing movies) in real time;
  13) Add control point visualization suite. 
  14) Add tests for the accuracy of the code:
    i. [X] Monolithic shift of an image; 
    ii. [X] Rotation of an image; 
    
    
  2) Test which approach to reg_loop-ing is better for restoring the original image: 
    - finer grid with large boxes and oversampled image; 
    - nested reg_looping; 
  3) Test for the atmospheric smearing by comparison with MHD simulation + atmospheric OTF; 

