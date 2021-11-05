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

#Main Code, Quantum Simulator with Measurement
def QSimulator(gates, wires, n):
    """Takes a LIST of gates (e.g. [H,CNOT,DFT]) and a nested list wires (e.g. [[1],[1,3],[1,2,3,4]])
    
    returns bitstring corresponding to a measurement of final state
    
    initial state is assumed to be |00...0>"""
    #initialize variables and state
    N = 2**n
    state = np.identity(N)[:,0]
    #use tensor product to make all gates act on all qubits
    modified_gates = []
    for a in range(len(gates)):
        current_gate = gates[a]
        current_N = len(current_gate)
        current_n = int(np.log2(current_N))
        modified_gate = np.zeros([N,N])
        for i in range(current_N):
            for j in range(current_N):
                #i --> column index; j --> row index
                current_wires = wires[a]
                elementary = []
                rawbitstringi = format(i,'b')
                bitstringi = '0'*(current_n-len(rawbitstringi))+rawbitstringi
                rawbitstringj = format(j,'b')
                bitstringj = '0'*(current_n-len(rawbitstringj))+rawbitstringj
                for x in range(n):
                    if current_wires != [] and x+1 == current_wires[0]:
                        elementary.append(bra_ket(bitstringi[0],bitstringj[0]))
                        bitstringi = bitstringi[1:]
                        bitstringj = bitstringj[1:]
                        current_wires = current_wires[1:]
                    else:
                        elementary.append(np.identity(2))
                modified_gate += current_gate[i,j]*reduce(np.kron,elementary)
        modified_gates.append(modified_gate)
    for a in range(len(modified_gates)):
        state = np.dot(modified_gates[a],state)
    return measure(state)

#Test with EPR Pair; should get 50% '00' and 50% '11'
test = [QSimulator([H, CNOT],[[1],[1,2]],2) for i in range(10)]
print(test)
