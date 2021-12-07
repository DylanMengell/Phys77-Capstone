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
    #generate pdf from quantum state by taking norm squared of each amplitude
    pdf = list(map(lambda complexnum: polar(complexnum)[0]**2,state))
    #run through each element of pdf to determine measurement
    for i in range(len(pdf)):
        if random<pdf[i]:
            rawoutput = format(i,"b")
            return '0'*(n-len(rawoutput))+rawoutput
        random -= pdf[i]

def Q(f, n_input):
    """generates matrix which implements classical function quantumly
    f should be a function from {0,1}^n_input --> {0,1}^m for some int m"""
    #determine min number of bits necessary to represent all outputs
    possible_outputs = {}
    for x in range(2**n_input):
        possible_outputs[f(x)] = x
    N_output = len(list(possible_outputs.keys()))
    n_output = math.ceil(math.log2(N_output))
    output_to_bit = {}
    #associate each output with a corresponding number (equivalent to bitstring)
    i = 0
    for f_output in list(possible_outputs.keys()):
        output_to_bit[f_output] = i
        i += 1
    #generate matrix of correct size with all zeros
    output = np.zeros([2**(n_input+n_output),2**(n_input+n_output)])
    #figure out where to insert 1s in order to make correct matrix
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
            output[int(bit_x+modified_output,2),int(bit_x+bit_b,2)] = 1 #key line: implements qcomputing law for generating matrix of Qf
    return output

def euclideanAlg(a,N): 
    """euclidean algorithm for finding GCD"""
    if a == 0:                                  #Once a = 0, we can return N
        return N
    else:                                       #recursion so that N%a is the new 'a' and the old 'a' becomes the new value for N
        return euclideanAlg(N % a, a)

#Simon's Period-Finding Alg. (Part 1, Quantum Computation) (returns bitstrings)

def internalQuantPerFind(N, a, size):
    """takes N, a, and integer size
    returns list of bitstrings results from measurement of quantum state
    each bitstring should give a hint as to the period of f(x) = (a**x)%N"""
    n = math.ceil(math.log2(N))
    
    def f(x):
        return (a**x)%N

    #generate compute gate
    compute = Q(f,n)
    n_total = int(math.log2(np.shape(compute)[0]))

    #generate two rotations gates
    rotate1 = tensor_product(HFT(n),np.identity(2**(n_total-n)))
    rotate2 = tensor_product(IDFT(2**n),np.identity(2**(n_total-n)))
    
    #execute quantum computation
    state = np.zeros(2**n_total)
    state[0] = 1
    for gate in [rotate1,compute,rotate2]:
        state = np.dot(gate,state)
    
    #measure final state "size" times
    return [measure(state)[0:n] for i in range(size)] #will return a list of bitstrings of length 'size'

#Simon's Period-Finding Alg. (Part 2, Using Bitstrings to Find Period)

def QuantPeriodFinding(N : int, a : int) -> int:
    """
    QuantPeriodFinding finds the period 'r' of the functon (a^x) mod N
    
    args: 
        N: number to be factored
        a: random guess to help find factors of N
    
    returns: 
        r: the period of the function (a^x) mod N, which will help us find factors of N
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
    secondary_checks = [4,3,2]
    for check in secondary_checks:
        if (a**int(Rforminremainder/check))%N == 1:
            return int(Rforminremainder/check)
    if (a**int(Rforminremainder))%N == 1:
        return Rforminremainder
    else:
        return None

#Main Code: Shor's Algorithm

#Put the Composite Number Here
def ShorsAlgo(N):
    #Check is N is even, if so we can skip everything and return 2 & N/2
    if (N % 2) == 0:
            return 2, int(N/2)
    tries = 0  #initialize 'tries' at 0 to keep count of loops
    while tries < int(math.log2(N)):
        #0) Keep track of number of iterations through while loop (if too large, we have confidence N is prime) 
        tries += 1

        #1) Pick a random number 1<a<N
        a = np.random.randint(1,N)
        print("a =", a) 
                                         
        #2) Compute K = GCD(a,N) using Euclidean Algorithm 
        K = euclideanAlg(a, N)

        #3) Check K
        if K != 1:                                                  #If if K!=1 then it is non trivial(i.e. it is a factor)
            non_trivial_factor = K                                  #WE DID IT, no quantum needed
            return non_trivial_factor, int(N/non_trivial_factor)

        #4 Use the quantum period-finding subroutine to find r
        r = QuantPeriodFinding(N, a)  
        print("r =", r)
        if r == None:
            continue #if we fail to find the period of f(x) = (a**x)%N, we just pick a new value of a                        

        #5 If r is even and if (a**(r/2))%N != N-1 then the factors are as such:
        if ((r % 2) == 0) and ((a**int(r/2))%N != N-1):
            """bit_power = format(int(r/2),"b")
            a_power_2 = [a]
            for i in range(len(bit_power)-1):
                a_power_2.append((a_power_2[i]**2)%N)
            result = 1
            for i in range(len(bit_power)):
                if bit_power[::-1][i] == "1":
                    result *= a_power_2[i]
            result %= N #should give a**(r/2)%N
            non_trivial_divisor1 = euclideanAlg(result - 1, N)
            non_trivial_divisor2 = euclideanAlg(result + 1, N)"""
            non_trivial_divisor1 = euclideanAlg(a**int(r/2) - 1, N)
            non_trivial_divisor2 = euclideanAlg(a**int(r/2) + 1, N)
            return non_trivial_divisor1, non_trivial_divisor2
    return "Having run through Shor's algorithm {0} times, we have confidence that {1} is prime.".format(tries, N)

print(ShorsAlgo(89))

#Useful Visualization Functions (used to generate diagrams for presentation)

def visualize_f(N,a):                                                                           #plot to show the periodicity of function & to prove there is a true period
    """Plot same as to https://qiskit.org/textbook/ch-algorithms/shor.html"""
    y = []                                                                                      #create an array for y values
    def f(x):
        return (a**x)%N                                                                         #define a function (a*x)%N to show periodicity
    x = np.arange(N)
    for i in range(N):                                                                          #loop to create Y values
        b = f(i)
        y.append(b)
    plt.plot(x, y)                                                                              #plot
    plt.ylabel(r"${0}^x$ mod ${1}$".format(a,N))
    plt.xlabel(r"$x$")
    plt.title(r"Periodic Function in Shor's Alg: $f(x)={0}^x$ mod ${1}$".format(a,N))
    plt.show()

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
    rvals = list(range(2, int(N/2)-1))
    for r in rvals: 
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
    plt.show()

#periodgraph(91, 19)