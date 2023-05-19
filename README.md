# Weiszfeld approximation algorithm for approximating the Fermat point for arbitrary dimension of points.

This is a modification and an implementation of Weiszfeld approximation algorithm from the paper: https://cs.nyu.edu/exact/doc/fermat-esa21.pdf. There are 2 version of the implementation––– one using standard cpp numerical precision and the other using "exact" precision (CORE). CORE can be installed from https://cs.nyu.edu/exact/core_pages/intro.html. 

 The dataset and the corresponding dimension can be changed in the main.cpp file.
 
 This internship involved extending the 2-D epsilon Weiszfeld approximation algorithm for approximating the Fermat point to arbitrary dimensions and using the CORE library for reliable exact computation. The approximation algorithm was introduced by the paper: “Certified Approximation Algorithms for the Fermat Point and n-Ellipses”, which only dealt with 2-D points. The first task for this project was to implement the 2D epsilon approximation algorithm which formed a backbone for the extension. Significant portion of the time was spent on formulating and implementing the Newton operator for interval matrices, which involved experimenting with various Newton operators such as (Moore’s & Nickel’s) and Krawzcyk’s and making approximations/changes to the formulations which were specifically intended to handle 2-dim points. Newton operators such as (Moore’s & Nickel’s) consisted of a term which computed the inverse of an interval matrix, which is non-viable for arbitrary dimension matrices so Krawzcyk’s inverse was used. This required the computation of gradient of the Fermat distance function and the hessian of the Fermat distance function, which were also modified.
The algorithm was then implemented in C++ using the CORE Library. Many experiments were performed using different points sets of different sizes. The algorithm struggled to converge for points of dimension 26 and above, when using exact computation. Hence for points of dimension 26 and above, certain computations were done exactly until 120 bits of accuracy.
 
## Run
```
make
./main
```

## using standard cpp numerical precision
![output1](/img/out1.png?raw=true "out1")

##using 120 bits of accuracy ~ 30 digits (CORE)
![output2](/img/out2.png?raw=true "out2")
