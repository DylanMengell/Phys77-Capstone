"""This was a piece of code that we thought may be necessary
to implement Shor's Algorithm on a classical computer. It takes in
a sequence of quantum gates and the wires they should be applied to.
It then returns the result of a quantum computation with those gates
applied to the specified wires. However, not only is this code quite 
slow, it turns out to be overkill for our version of Shor's algorithm. 
Nevertheless, I think the code is very insightful from a pedagogical 
standpoint for those who are interested in learning about quantum 
computing. Feel free to peruse the code at your leasure."""

import numpy as np
from cmath import polar
from functools import reduce

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

#Hadamard Transform
def HFT(n):
    if n == 1:
        return H
    return np.kron(HFT(n-1),H)

#Measurement Functions
def amp_to_prob(cx_number):
    """takes a complex number a+bi and returns its norm squared a^2+b^2
    
    this is how we get from an amplitude (the numbers stored in the column vector) 
    to an actual probability"""
    return polar(cx_number)[0]**2
def measure(state):
    """takes a state (column vector) and returns a bit-string, simulating a measurement of the state
    
    remember that a column vector [a_1, a_2, ..., a_N] has probabilty amp_to_prob(a_i) of resulting in 
    a measurement of |bit(i)> where bit(i) is the bitstring associated with the number i"""
    
    #this is fairly standard notation 
    #n --> number of qubits; N --> dimension of vectorspace we're dealing with
    N = len(state)
    n = int(np.log2(N))
    #generate random number from 0 to 1
    random = np.random.random()
    #generate pdf from quantum state
    pdf = list(map(amp_to_prob,state))
    #run through each element to determine measurement
    for i in range(len(pdf)):
        if random<pdf[i]:
            rawoutput = format(i,"b")
            return '0'*(n-len(rawoutput))+rawoutput
        random -= pdf[i]

#Other Miscellaneous, Useful Functioncs
def binaryToDecimal(bitstring):
    """takes bitstring 
    returns decimal equivalent"""
    binary = int(bitstring)
    decimal, i = 0, 0
    while(binary != 0):
        dec = binary % 10
        decimal = decimal + dec * pow(2, i)
        binary = binary//10
        i += 1
    return decimal 

def bra_ket(bitstring1,bitstring2):
    """takes two bitstrings 
    returns |bitstring1><bitstring2|"""
    n1 = len(bitstring1)
    N1 = 2**n1
    n2 = len(bitstring2)
    N2 = 2**n2
    output = np.zeros([N1,N2])
    output[binaryToDecimal(bitstring1),binaryToDecimal(bitstring2)] = 1
    return output

def tensor_product(a,b):
    shape_a = np.shape(a)
    shape_b = np.shape(b)
    a = a[:,np.newaxis,:,np.newaxis]
    a = a[:,np.newaxis,:,np.newaxis]*b[np.newaxis,:,np.newaxis,:]
    a.shape = (shape_a[0]*shape_b[0],shape_a[1]*shape_b[1])
    return a

#Main Code, Quantum Simulator with Measurement
def QSimulator(gates, wires, n, size=1):
    """Takes a LIST of gates (e.g. [H,CNOT,DFT]) and a nested list wires (e.g. [[1],[1,3],[1,2,3,4]])
    
    returns list of bitstrings corresponding to measurements of final state
    
    initial state is assumed to be |00...0>"""
    #initialize variables and state
    N = 2**n
    state = np.identity(N)[:,0]
    #use tensor product to make all gates act on all qubits
    modified_gates = []
    for a in range(len(gates)):
        #define current gate and wires we're working with
        current_gate = gates[a]
        current_wires = wires[a]
        #define the N and n of current_gate (not to be confused with N and n of total circuit)
        current_N = len(current_gate)
        current_n = int(np.log2(current_N))
        #initialize modified gate to be zero N by N matrix
        modified_gate = np.zeros([N,N])
        #iterate over all cells in matrix of current_gate
        for i in range(current_N):
            for j in range(current_N):
                #if current_gate[i,j] == 0, no need to proceed
                if current_gate[i,j] == 0:
                    continue
                #i --> column index; j --> row index
                m = 0 #bit index (0<=m<n)
                #construct elementary member of decomposition of modified_gate 
                elementary = [np.identity(2)]*n
                #turn i and j into current_n-bit strings
                rawbitstringi = format(i,'b')
                bitstringi = ('0'*(current_n-len(rawbitstringi)))+rawbitstringi
                rawbitstringj = format(j,'b')
                bitstringj = ('0'*(current_n-len(rawbitstringj)))+rawbitstringj
                #iterate over wires
                for x in current_wires:
                    #implement key pattern in elementary decompositions (KEY PART OF CODE)
                    elementary[x-1] = bra_ket(bitstringi[m],bitstringj[m])
                    m += 1
                #construct modified_gate from each elementary decomposition
                modified_gate += current_gate[i,j]*reduce(tensor_product,elementary)
        #fill list of modified_gates
        modified_gates.append(modified_gate)
    #perform quantum computation with matrix multiplication
    for a in range(len(modified_gates)):
        state = np.dot(modified_gates[a],state)
    #return list of bitstrings corresponding to random measurements of state
    return [measure(state) for i in range(size)]

#Test with 3 qubits; should get 50% '000' and 50% '101'
test = QSimulator([H, CNOT],[[1],[1,3]],3,size=10)
print(test)