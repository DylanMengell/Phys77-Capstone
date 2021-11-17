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

        #1 Pick a random number 1<a<N
        non_trivial_divisor = np.random.randint(1,N,1)  #possibly find better random num gen

        #Compute K = GCD(a,N) using Euclidean Algorithm 
        #works!
        K = euclideanAlg(non_trivial_divisor, N)

        #3
        if K != 1:
            non_trivial_factor = K #WE DID IT, no quantum needed
            print('Your GCD is: ', non_trivial_factor)
        elif K == 1:
            print('The values are Coprime. Your GCD is 1, therefore we must continue to step 4.')

        #4 Use the quantum period-finding subroutine to find r
        #(non_trivial_divisor)**p = m * N +1

        #Need to create a way to test out an infinite number of values 
        #for p and m, including fractions

        #leave this up to quauntum computing

        #5 if r is odd, go back to step 1
        if (r % 2) == 0:
            non_trivial_divisor = np.random.randint(0,N,1) 
        #6 if a^(r/2) = -1(mod(N)) go back to step 1

        # 7 other wise 
        # K1 = euclideanAlg((a**(r/2) + 1, N)
        # and 
        # K2 = euclideanAlg((a**(r/2) - 1, N)
        #^^ these are the non_trivial_factors of N
            
        #goal = find non_trivial_factor


        #==============add parts 5 and 6 to end inf loop