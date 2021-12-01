import math
import numpy as np
from cmath import polar, exp
import matplotlib.pyplot as plt

#Some Useful Functions for Matrix Implementation of QComputing

def tensor_product(a,b):
    """fast tensor product
    
    from https://stackoverflow.com/questions/56067643/speeding-up-kronecker-products-numpy"""
    shape_a = np.shape(a)
    shape_b = np.shape(b)
    a = a[:,np.newaxis,:,np.newaxis]
    a = a[:,np.newaxis,:,np.newaxis]*b[np.newaxis,:,np.newaxis,:]
    a.shape = (shape_a[0]*shape_b[0],shape_a[1]*shape_b[1])
    return a

H = 0.5**0.5 * np.array([[1, 1],\
                         [1, -1]]) 

def HFT(n):
    """Hadamard Transform"""
    if n == 1:
        return H
    return tensor_product(HFT(n-1),H)

def DFT(N):
    """Discrete Fourier Transform"""
    output = np.empty([N,N],dtype=complex)
    omega = exp(2*np.pi*1j/N)
    for i in range(N):
        for j in range(N):
            output[i][j]= omega**(-i*j)
    output *= N**-0.5
    return output

def IDFT(N):
    """Inverse Discrete Fourier Transform"""
    output = np.empty([N,N],dtype=complex)
    omega = exp(2*np.pi*1j/N)
    for i in range(N):
        for j in range(N):
            output[i][j]= omega**(i*j)
    output *= N**-0.5
    return output

def bitwiseXOR(string1, string2):
    """performs bitwise XOR of two strings of the same length"""
    output = ""
    for i in range(len(string1)):
        output += str((int(string1[i])+int(string2[i]))%2)
    return output

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

def Q(f, n_input):
    """generates matrix which implements classical function quantumly
    f should be a function from {0,1}^n_input --> {0,1}^m for some int m"""
    possible_outputs = {}
    for x in range(2**n_input):
        possible_outputs[f(x)] = x
    N_output = len(list(possible_outputs.keys()))
    n_output = math.ceil(math.log2(N_output))
    output_to_bit = {}
    i = 0
    for f_output in list(possible_outputs.keys()):
        output_to_bit[f_output] = i
        i += 1
    output = np.zeros([2**(n_input+n_output),2**(n_input+n_output)])
    for x in range(2**n_input):
        raw_bit_x = format(x,"b")
        bit_x = ('0'*(n_input-len(raw_bit_x)))+raw_bit_x
        classical_output = output_to_bit[f(x)]
        raw_bit_output = format(classical_output,"b")
        bit_output = ('0'*(n_output-len(raw_bit_output)))+raw_bit_output
        for b in range(2**n_output):
            raw_b = format(b,"b")
            bit_b = ('0'*(n_output-len(raw_b)))+raw_b
            modified_output = bitwiseXOR(bit_b,bit_output)
            output[int(bit_x+modified_output,2),int(bit_x+bit_b,2)] = 1
    return output

#"""Asher: Writing a function like this may be cheating a little too much..."""
"""def is_prime(N: int) -> bool:
    if N <= 3:
        return N > 1
    if N % 2 == 0:
        return False
    if N % 3 == 0:
        return False
    i = 5
    while i ** 2 <= N:
        if N % i == 0 or N % (i+2) == 0:
            return False
        i += 6
    return True
"""
def euclideanAlg(a,N): 
    """euclidean algorithm for finding GCD"""
    if (a==0):
        return N
    else:
        return euclideanAlg(N % a, a)

def visualize_f(N,a):
    y =[]

    def f(x):
        return (a**x)%N

    #plot same as to https://qiskit.org/textbook/ch-algorithms/shor.html 
    x = np.arange(N)
    for i in range(N):
        b = f(i)
        y.append(b)
    plt.plot(x, y)
    plt.ylabel("a**x mod N")
    plt.xlabel("x")
    plt.title("Periodic Function in Shor's Alg.")
    plt.show()


#Simon's Period-Finding Alg. (Part 1, Quantum Computation) (returns bitstrings)

def internalQuantPerFind(N, a, size):
    """takes N, a, and integer size
    returns list of bitstrings results from measurement of quantum state
    each bitstring should give a hint as to the period of f(x) = (a**x)%N"""
    n = math.ceil(math.log2(N))
    
    def f(x):
        return (a**x)%N

    compute = Q(f,n)
    n_total = int(math.log2(np.shape(compute)[0]))

    rotate1 = tensor_product(HFT(n),np.identity(2**(n_total-n)))
    rotate2 = tensor_product(IDFT(2**n),np.identity(2**(n_total-n)))
 
    state = np.zeros(2**n_total)
    state[0] = 1
    for gate in [rotate1,compute,rotate2]:
        state = np.dot(gate,state)
    return [measure(state)[0:n] for i in range(size)] #will return a list of bitstrings of length 'size'

#Simon's Period-Finding Alg. (Part 2, Using Bitstrings to Find Period) (returns bitstrings)

def QuantPeriodFinding(N : int, a : int) -> int: #quantum Period finding algorthm Should return r
    """
    QuantPeriodFinding finds the period 'r' of the functon (a^x)Mod(N)
    
    args: 
        N: number to be factored
        a: random guess to help find factors of N
    
    returns: 
        r: is the period of the function (a^x)Mod(N) which will help us find factors of N
    """
    size = 1000
    bitstrings = internalQuantPerFind(N, a, size) 
    base10vals, numbits, fractionalVals, remainderTotals = [], [], [], []

    for i in range(size):
         base10vals.append(int(bitstrings[i], 2)) #converts value to base 10
         numbits.append(len(bitstrings[i])) #stores the number of bits in given string
         fractionalVals.append(base10vals[i]/2**(numbits[i])) #calculates fractional value for each val #/2^num bits

    minremainder, Rforminremainder, indexCount = 10000, 1, 0
    for r in range(2, int(N/2)): 
        remainders = [] 
        for i in range(size):
            temp = fractionalVals[i]*r
            remainders.append(min(temp-int(temp),int(temp)+1-temp))
        remainderTotals.append(np.sum(remainders))
        if(remainderTotals[indexCount] < minremainder):
            minremainder = remainderTotals[indexCount]
            Rforminremainder = r
        indexCount += 1
    #Not sure if this is correct
    if a**(Rforminremainder/4)%N == 1:
        return Rforminremainder/4
    elif a**(Rforminremainder/3)%N == 1:
        return Rforminremainder/3
    elif a**(Rforminremainder/2)%N == 1:
        return Rforminremainder/2
    elif a**(Rforminremainder)%N == 1:
        return Rforminremainder
    else:
        return 1

#Main Code: Shor's Algorithm

#Put the Composite Number Here
def ShorsAlgo(N):
    if (N % 2) == 0:
            return 2, N/2
    while True:
        #1) Pick a random number 1<a<N
        a = np.random.randint(1,N)                                  #possibly find better random num gen

        #2) Compute K = GCD(a,N) using Euclidean Algorithm 
        K = euclideanAlg(a, N)

        #3) Check K
        if K != 1:                                                  #If if K!=1 then it is non trivial(i.e. it is a factor)
            non_trivial_factor = K                                  #WE DID IT, no quantum needed
            return non_trivial_factor

        #4 Use the quantum period-finding subroutine to find r

        r = QuantPeriodFinding(N, a)                                #yet to be made but will return r

        #5 & 6 If r is even and if a^r/2 != -1%N then the factors are as such:
        if ((r % 2) == 0) and (a**(r/2) != (-1)%N):                 #should this be 'and' or 'or'????
            non_trivial_divisor1 = euclideanAlg(a**(r/2) + 1, N)
            non_trivial_divisor2 = euclideanAlg(a**(r/2) - 1, N)
            return non_trivial_divisor1, non_trivial_divisor2

    visualize_f(N,a)


#print(ShorsAlgo(7*13)) #better example


def periodgraph(N : int, a : int) -> int: #quantum Period finding algorthm Should return r
    """
    periodgraph creates a graph the sum of remainders as a function of 'r'
    
    args: 
        N: number to be factored
        a: random guess to help find factors of N
    
    returns: 
        r: Plot sum of remainders as a function of r
    """
    size = 1000
    bitstrings = internalQuantPerFind(N, a, size) 
    base10vals, numbits, fractionalVals, remainderTotals = [], [], [], []

    for i in range(size):
         base10vals.append(int(bitstrings[i], 2)) #converts value to base 10
         numbits.append(len(bitstrings[i])) #stores the number of bits in given string
         fractionalVals.append(base10vals[i]/2**(numbits[i])) #calculates fractional value for each val #/2^num bits

    minremainder, Rforminremainder, indexCount = 10000, 1, 0
    rvals = list(range(2, int(N/2)))
    print(rvals)
    for r in range(2, int(N/2)): 
        remainders = [] 
        for i in range(size):
            temp = fractionalVals[i]*r
            remainders.append(min(temp-int(temp),int(temp)+1-temp)) #Asher: Note --> min(temp-int(temp),int(temp)+1-temp)
        remainderTotals.append(np.sum(remainders))
        if(remainderTotals[indexCount] < minremainder):
            minremainder = remainderTotals[indexCount]
            Rforminremainder = r
        indexCount += 1
    xticks = np.arange(0,N/2+1,2)
    plt.figure(figsize=(12, 6), dpi=80)
    plt.plot(rvals, remainderTotals, '--o')
    plt.xticks(xticks)
    plt.xlabel("r (period value)", fontsize = 12)
    plt.ylabel("total remainder values", fontsize = 12)



periodgraph(91, 15)