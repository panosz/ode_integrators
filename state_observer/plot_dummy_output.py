
import scipy as sp
import matplotlib.pyplot as plt


a = sp.loadtxt('../cmake-build-release/dummy.txt')
plt.plot(a[:,1],a[:,2])
plt.show()