
import numpy as np
from scipy import linalg

# We know that python uses IEEE double precision by default

eps = 10**(-12)
    
A = np.array([[eps, 2,3],
             [5, 7, 11], 
             [13, 17, 19]], dtype = float)

def gaussian_elim(A): # Gaussian elimination without pivoting
    
    n = len(A)
    
    U = A.copy()   # Don't modify the original matrix
    
    for j in range(n): 
        
        for i in range(j+1, n):
            if U[j, j] == 0:
                raise ZeroDivisionError("Zero pivot! need pivoting.")
            
            m_ij = U[i, j] / U[j, j] # Compute the multiplier
            
            U[i, :] = U[i, :] - m_ij * U[j, :] # Subtract m*row j from row i
            
            # print(U)  # To see every step
            
    return U
             

U = gaussian_elim(A)

print(linalg.det(U), linalg.det(A))