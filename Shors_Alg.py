#Shors Alg. Classical part
import numpy as np


def QuantPeriodFinding(): #quantum Period finding algorthm Should return r
    """
    def part4(X,Y): #(Non_trivial_divisor, N)
        p = np.linspace(1,Y)
        for j in p:
            m = np.linspace(1,Y)
            for l in m:
                #remember j and l are m and p
                if ((X**j) == (l * Y) + 1): #ERROR HERE
                    return (j,l)

    p,m = part4(non_trivial_divisor, N)
    #goal find a better a
    def findBettera (p,m,N):
        a = ((m * N) + 1)**p
        return a

    a = findBettera(p,m,N)
    print(a)
    """
    return 0

def euclideanAlg(a,N): #Euclidean algorithm for finding GCD
    if (a==0):
        return N
    else:
        return euclideanAlg(N % a, a)

#Put the Composite Number Here
def ShorsAlgo(N):
    isdone = True
    while not (isdone):
        if (N % 2) == 0:
            return ('The value is even please enter an odd N!')

        #1) Pick a random number 1<a<N
        a = np.random.randint(1,N,1)  #possibly find better random num gen

        #2) Compute K = GCD(a,N) using Euclidean Algorithm 
        K = euclideanAlg(a, N)

        #3) Check K
        if K != 1: #If if K!=1 then it is non trivial(i.e. it is a factor)
            non_trivial_factor = K #WE DID IT, no quantum needed
            return non_trivial_factor

        #4 Use the quantum period-finding subroutine to find r
        #(non_trivial_divisor)**p = m * N +1

        r = QuantPeriodFinding() #yet to be made but will return r

        #5 if r is even and if a^r/2 != -1%N then the factors are as such:
        if ((r % 2) == 0) and (a**(r/2) != (-1)%N): #should this be 'and' or 'or'????
            non_trivial_divisor1 = euclideanAlg(a**(r/2) + 1, N)
            non_trivial_divisor2 = euclideanAlg(a**(r/2) - 1, N)
            return non_trivial_divisor1, non_trivial_divisor2