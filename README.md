# FracSimOpt
**Estimating Fractured rock properties via ML**

This is a code I developed in MATLAB environment for my masters thesis in summer 2012. The aim was to forecast the properties of the fractured material based on numerical simulations coupled with Machine Learning algorithms. 

The error function would be calculated and via iterative numerical optimization algorithms, machine would learn to move down the "valley" of error and converge. I developed the numerical optimization algorithms A-Z and from scratch and didn't use embedded functions of MATLAB. Such algorithms including gradient descent, Levenberg Marquardts, and Newton method are widely used nowadays in machine learning cores to allow model convergence over epochs. 


**Error surface:**
![image](https://user-images.githubusercontent.com/121983512/223020779-1152af14-febf-49ff-a0d1-187919d2ed72.png)


**Gradient descent (GD):**
![image](https://user-images.githubusercontent.com/121983512/223021199-2daa5a1d-3f62-4466-a18b-aaa5608fc0d9.png)

**Newton:**
![image](https://user-images.githubusercontent.com/121983512/223021257-4e08b358-f934-40a2-b939-a4ba677e9e90.png)


**Levenberg-Marquardt (also known as Damped Least Squares):**
![image](https://user-images.githubusercontent.com/121983512/223021765-b1c32d87-efc1-4992-b436-9600e66b4076.png)

