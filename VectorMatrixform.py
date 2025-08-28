import numpy as np

class VectorMatrixForm:
    def __init__(self) :
        pass
    def VecCreate (self,func,n,ntot):        #VecCreate stores data of Function(i,j) in a vector
        x = [func(y,x) for x in range(ntot) for y in range(n)]
        a = []
        a.append(x)
        return a

    def MatrixCreate(self,array, n , ntot):
        for i in range (ntot):
            for j in range (n):
                for k in range (ntot):
                    for l in range (n):
                        array[len(array):] = [[l,j ,k,i]]
        x=np.array(array).reshape(n*ntot,n*ntot,4)
        return x

