

from scipy import linalg
import matplotlib.pyplot as plt
import numpy as np

N = 101

diag_entries = [2*np.sin(2*np.pi*j/N) for j in range(1, N+1)]

# Create the main diagonal
main_d_matrix = np.diag(diag_entries) 

# Create the upper diagonal
upper_d_matrix = np.diag(np.ones(N-1), 1)

# Create the lower diagonal
lower_d_matrix = np.diag(-1*np.ones(N-1), -1)

# Finish constructing A and add the corner entries
A = main_d_matrix + upper_d_matrix + lower_d_matrix
A[0, N-1] = -1
A[N-1, 0] = 1

# Find eigenvalues
eigs = np.linalg.eigvals(A)

# Extract the real and imag components for plot
x = [z.real for z in eigs]
y = [z.imag for z in eigs]

# Plotting
plt.figure(figsize=(6, 6))
plt.scatter(x, y, color="k", label = "Original A")
plt.axis('equal')
plt.xlabel('Real Axis')
plt.ylabel('Imaginary Axis')

colors = ["C0", "C1", "C2", "C3"]

# Make E:
E = np.random.normal(loc=0, scale=(1/np.sqrt(N)), size=(N, N)) 

epsilons = [1e-12, 1e-9, 1e-6, 1e-3]

for i in range(len(epsilons)):
    A_mod = A+(epsilons[i]*E)
    
    # Same eigenvale process
    eigs = np.linalg.eigvals(A_mod)
    x = [z.real for z in eigs]
    y = [z.imag for z in eigs]
    
    plt.scatter(x, y, 
                s = 10, 
                label = f"Epsilon = {epsilons[i]}", 
                c = colors[i])
    
plt.legend(fontsize = 7)

###########

reals = np.linspace(-2, 2, 100)
imags = np.linspace(-2, 2, 100)

def f(real, imag):
    z = real + 1j*imag
    matrix = z*np.eye(N)-A
    matrix_inverse = np.linalg.inv(matrix)
    return np.linalg.norm(matrix_inverse, ord=2)

real, imag = np.meshgrid(reals, imags)

# Vectorize f so it can be used element by emlement
f_vec = np.vectorize(f)

levels = [eps**(-1) for eps in epsilons]

# Contour plot
CS = plt.contour(real, imag, 
            f_vec(real, imag), 
            levels = levels[::-1],
            colors = colors[::-1])  #levels increasing


######



print("Spectral Radius:", np.max(abs(eigs)))
print("Spectral Abscissa:", np.max(y))

norm = np.linalg.norm(A, ord=2)

print("2-norm:", norm)










