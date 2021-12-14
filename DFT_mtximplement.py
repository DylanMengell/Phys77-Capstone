"""This was one of our preliminary tests to see if
our implementation of the Discrete Fourier Transform 
and Inverse Discrete Fourier Transforms was working."""

from cmath import exp
import numpy as np

def DFT(N):
    output = np.empty([N,N],dtype=complex)
    omega = exp(2*np.pi*1j/N)
    for i in range(N):
        for j in range(N):
            output[i][j]= omega**(-i*j)
    output *= N**-0.5
    return output

def IDFT(N):
    output = np.empty([N,N],dtype=complex)
    omega = exp(2*np.pi*1j/N)
    for i in range(N):
        for j in range(N):
            output[i][j]= omega**(i*j)
    output *= N**-0.5
    return output

print(DFT(4))
print("")
print(IDFT(4))
print("")
print(np.dot(DFT(4),IDFT(4))) #should print get the identity matrix
