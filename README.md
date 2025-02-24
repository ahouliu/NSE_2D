# NSE_2D

用 Fourier-Galerkin计算二维Navier-Stokes方程，streamfunction-velocity formulation.
时间外力f(x,y)g(t), 周期T。若T=Inf，视为时不变外力。时间演化方法选择的ETDRK4。
空间区域：[0,Lx]*[0,Ly],但只用Lx=Ly=2\pi。
BC：周期性边界条件.

各个文件夹说明：

DNS——二维Naver-Stokes方程的数值计算，外力F = (F1*g(t),F2*g(t)), Fj= ej(1-cos(wx))(1-cos(wy))/4+CONST, j=1,2. 但因为选择的streamfunction-velocity formulation，实际上只需要f:=F2_x-F1_y，然后再乘上时变g(t)即可。
主文件DNS_1_MAIN。
certain_forces——同DNS，差别是外力f.选择较多，其中一种见f(x,y).jpg.




注意：vorticity-stream --> only stream (4阶方程), may unconditionally unstable!!
参考文献： David Gottlieb, Steven A. Orszag: Numerical Analysis of Spectral Methods - Theory and Applications
