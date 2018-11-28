# SLURP

This repository contains a pre-release version of SLURP, a Speech and Language Ultrasound Research Package.  At the moment, SLURP primarily consists of software to semi-automatically delineate and track the contour of the tongue in ultrasound recordings of speech.  Among other things, it contains an implementation of the particle filtering algorithm described in the following paper:

Laporte C, Ménard L. (2018). Multi-hypothesis tracking of the tongue surface in ultrasound video recordings of normal and impaired speech. Medical Image Analysis. 44: 98-114.

Please cite this paper if you decide to use SLURP in your research.

The software was built on top of an old version of GetContours (the newest version is available at https://github.com/mktiede/GetContours) and shares many of its user interface features.  It may (no, it DOES) also contain a few bugs as a result of our building on top of it for research purposes, which we are planning to fix prior to a more official release.  There is also no user manual at the moment other than this ReadMe file, so here are some very basic instructions to get our particle filtering method running:

	1) Start MATLAB and change the path to the directory where you installed SLURP
	2) Start the SLURP software by calling SLURP
	3) Open a video (.avi file) using the SLURP menu on the right
	4) Select your initial frame using the scroll bar at the bottom of the figure (or just use the first one)
	5) Choose about 6 anchor points just below the tongue surface, as you would have done using Edgetrak or GetContour
	6) Choose the Fit Snake option in the SLURP menu (this fits a snake to only this initial frame)
	7) To track the tongue contour using particle filtering, choose the “Track snake using particle filter” item from the SLURP menu

This should track the snake backwards to frame #1, then forwards to the end of the video (or perhaps the second to last frame), and show you a weird greyscale image afterwards (the Energy map), which shows a measure of contour quality as a function of frame # and position along the contour, white meaning potentially poor quality and black meaning probably good quality.  Browsing this map with the mouse will show you the corresponding frame in the main figure.  There are also SLURP menu items that can be used to semi-automatically fix tracking errors, but admittedly, they are not as easy to use as We would like at the moment and they are also a bit unstable.  You can always give them a try (no instructions provided for now).  These tools are described in the following paper:

Ghrenassia S, Ménard L, Laporte C. (2014). Interactive segmentation of tongue contours in ultrasound video sequences using quality maps. SPIE Medical Imaging, San Diego, United States (90344O-1 – 90344O-7)
	
	8) To save the contours, choose Save State from the SLURP menu. 

This will produce a .mat file containing a structure called “state”.  The contour information will be in the XY field (state.XY), a 3D array.  The three dimensions are: position along the contour (39 points per contour), coordinate (two coordinates: x and y) and frame number.
 
The software comes with default tongue shape and motion models (ShapeModel.mat and MotionModel.mat), which were built from Canadian French VCV tongue contour sequences.  We don’t expect too many problems with tracking tongue motion in other languages, but there may be some language dependence here.  If needed, new models more suited to your data can be built (the code is also included in this repository (without instructions for the time being, however).  The default parameters of the snake and particle filters can be changed using items in the SLURP menu; you might want to try this first if you find the software doesn’t work well on your data.  
