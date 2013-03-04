gaborjets
=========
Java code for computing "Gabor Jet" features from images (Stacked gabor filter responses at multiple orientations/scales). The included example shows how to track features across two included images. 
Based on Matlab code from [Xiaomin Yue](http://nmr.mgh.harvard.edu/~xiaomin/). Depends on Piotr Wendykier's [JTransforms](https://sites.google.com/site/piotrwendykier/software/jtransforms) (included, but outdated).


Usage:
======
Compile:  
> javac GWTGrid.java

Run: 
> java GWTGrid.java

> maxSim=0.9109851412448251, maxSimIdx=23,26
> elapsed time=666 ms

The default code prints the coordinates in corner2.jpg that are maximally similar to corner1.jpg at x,y = (64,58).
