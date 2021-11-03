import numpy as np

#Hadamard Gate
H = 0.5**0.5 * np.array([[1, 1],\
                         [1, -1]]) 

#NOT Gate
NOT = np.array([[0,1],\
                [1,0]])

#Controlled NOT Gate
CNOT = np.array([[1,0,0,0],\
                 [0,1,0,0],\
                 [0,0,0,1],\
                 [0,0,1,0]])
  
