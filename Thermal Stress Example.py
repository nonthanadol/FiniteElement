### FEM Program for Solving System of Linear Spring Problems

import math, numpy, sys

node = 5                                            # Number of Nodes
nele = 4                                            # Number of Elements

n = node                                            # Number of Equations

u = n*[0]                                           # Initial Guess for Solution Vector
F = n*[0]                                           # Initializing of Load Vector

kele  = 0*numpy.identity(2)                         # Element Matrix of Two-Node Line Element
f     = [2*[0] for i in range(nele)]                # Initializing of Element Load Vectors

N   = 100                                           # Maximum Iteration of PCGNE Solver
tol = 1.e-9                                         # Tolerance for PCGNE Solver


###--------------------------------------------------------------------------------------------
### Read Node Data:
###--------------------------------------------------------------------------------------------

R1 = open('thermal_input_node_data.txt')
text = R1.readline()
inode, xc, c1, c2, c3 = [],[],[],[],[]
for line in R1.readlines():
    fields = line.split()
    inode.append(int(fields[0]))
    xc.append(float(fields[1]))
    c1.append(float(fields[2]))
    c2.append(float(fields[3]))
    c3.append(float(fields[4]))
R1.close()

if (node != len(inode)):
    print("Please check parameter 'node' (Number of Nodes) and input_node_data.txt.")
    sys.exit()

print('Node data was read properly.','\n')
#print(inode,ibc,c1,c2)
#sys.exit()


###--------------------------------------------------------------------------------------------
### Read Element Data:
###--------------------------------------------------------------------------------------------

R2 = open("thermal_input_element_data.txt")
text = R2.readline()
iele, n1, n2, A, E, alpha, dT = [],[],[],[],[],[],[]
for line in R2.readlines():
    fields = line.split()
    iele.append(int(fields[0]))
    n1.append(int(fields[1]))
    n2.append(int(fields[2]))
    A.append(float(fields[3]))
    E.append(float(fields[4]))
    alpha.append(float(fields[5]))
    dT.append(float(fields[6]))
R2.close()

if (nele != len(iele)):
    print("Please check parameter 'nele' (Number of Elements) and input_element_data.txt.")
    sys.exit()

print('Element data was read properly.','\n')
#print(iele,n1,n2,k)
#sys.exit()


###--------------------------------------------------------------------------------------------
### Assemble Stiffness Matrix:
###--------------------------------------------------------------------------------------------

def K_element(j,nele,n1,n2,xc,A,E):
    i1 = n1[j] - 1
    i2 = n2[j] - 1
    dl = xc[i2] - xc[i1]
    k  = A[j]*E[j]/dl*1000
    kele[0][0] =  k
    kele[0][1] = -k
    kele[1][0] = -k
    kele[1][1] =  k
    return i1, i2, kele

def assem_K(n,nele,n1,n2,xc,A,E):
    K = 0*numpy.identity(2*n)
    for j in range(nele):
        i1, i2, kele = K_element(j,nele,n1,n2,xc,A,E)
        K[i1][i1] = K[i1][i1] + kele[0][0]
        K[i1][i2] = K[i1][i2] + kele[0][1]
        K[i2][i1] = K[i2][i1] + kele[1][0]
        K[i2][i2] = K[i2][i2] + kele[1][1]
    return K

K = assem_K(n,nele,n1,n2,xc,A,E)
for j in range(n): K[j][n+j] = -1
print('System stiffness matrix was calculated.','\n')
#print(K)
#sys.exit()


###--------------------------------------------------------------------------------------------
### Apply Boundary Conditions:
###--------------------------------------------------------------------------------------------

def apply_bc(n,node,c1,c2,c3,K):
    b = 2*n*[0]
    for j in range(node):
        K[n+j][  j] = c1[j]
        K[n+j][n+j] = c2[j]
        b[n+j]      = c3[j]
    return K, b

K, b = apply_bc(n,node,c1,c2,c3,K)
print('Boundary conditions were specified.','\n')
#print(K)
#sys.exit()


###--------------------------------------------------------------------------------------------
### Thermal Load:
###--------------------------------------------------------------------------------------------

def thermal_load(nele,n1,n2,A,E,alpha,dT):
    for j in range(nele):
        i1 = n1[j] - 1
        i2 = n2[j] - 1
        TL = A[j]*E[j]*alpha[j]*dT[j]/1000
        b[i1] = b[i1] - TL
        b[i2] = b[i2] + TL
    return b

b = thermal_load(nele,n1,n2,A,E,alpha,dT)
#print(b)
#sys.exit()


###--------------------------------------------------------------------------------------------
### Solve The System of Linear Equations:
###--------------------------------------------------------------------------------------------

def pcgne(n,A,b,N,tol):
    x = n*[1]
    for k in range(N+1):
        if (k == 0):                                                # PCGNE Step 1
            r       = b - numpy.dot(A,x)                            # PCGNE Step 2
            norm_r0 = numpy.linalg.norm(r)
            AT      = numpy.transpose(A)
            Ar      = numpy.dot(AT,r)
            M       = numpy.diag(numpy.diagonal(A))
            for i in range(n): 
                if (abs(M[i][i]) < 1.e-6): M[i][i] = 1.
            M1      = numpy.linalg.inv(M)
            z       = numpy.dot(M1,Ar)
            d       = z
        else:
            if (k == 1): zAr = numpy.dot(z,Ar)
            else:        zAr = zAr_new            
            Ad      = numpy.dot(A,d)
            AdAd    = numpy.dot(Ad,Ad)
            alpha   = zAr/AdAd                                      # PCGNE Step 3
            x       = x + numpy.dot(alpha,d)                        # PCGNE Step 4
            r       = r - numpy.dot(alpha,Ad)                       # PCGNE Step 5
            norm_rk = numpy.linalg.norm(r)
            error   = norm_rk/norm_r0                               # PCGNE Step 6
            if (error >= tol and k <  N):                           # PCGNE Step 7-1
                Ar      = numpy.dot(AT,r)
                z       = numpy.dot(M1,Ar)
                zAr_new = numpy.dot(z,Ar)
                beta    = zAr_new/zAr
                d       = z + numpy.dot(beta,d)
            if (error >= tol and k == N):                           # PCGNE Step 7-2
                print('\n Root of system of linear equations cannot be found within',
                      '%2d' % N, 'PCGNE iterations.\n')
            if (error < tol):                                       # PCGNE Step 8
                return x

x = pcgne(2*n,K,b,N,tol)
print('Solution vector was found.','\n')
#print(x)
#sys.exit()


###--------------------------------------------------------------------------------------------
### Extract Node Solution:
###--------------------------------------------------------------------------------------------

for j in range(n):
    u[j] = x[j]
    F[j] = x[n+j]
print('Node solution was extracted.','\n')
#print(F)
#sys.exit()


###--------------------------------------------------------------------------------------------
### Print Node Solution:
###--------------------------------------------------------------------------------------------

W1 = open('thermal_output_node_solution.txt','w')
print('Node Solution:', file = W1)
print('%10s' % 'Node','%20s' % 'Displacement(mm)','%20s' % 'External Force(N)', file = W1)
for i in range(node):
    print('%10d' % (i+1),'%20.6f' % u[i],'%20.6f' % F[i], file = W1)
W1.close()
print('Node solution was written.','\n')
#sys.exit()


###--------------------------------------------------------------------------------------------
### Calculate Element Load Vectors:
###--------------------------------------------------------------------------------------------

def cal_element(nele,n1,n2,u,xc,A,E,alpha,dT):
    for j in range(nele):
        i1, i2, kele = K_element(j,nele,n1,n2,xc,A,E)
        uele = [u[i1],u[i2]]
        f[j] = numpy.dot(kele,uele)
        TL   = A[j]*E[j]*alpha[j]*dT[j]/1000
        f[j][0] = f[j][0] + TL
        f[j][1] = f[j][1] - TL
    return f

f = cal_element(nele,n1,n2,u,xc,A,E,alpha,dT)
print('Element load vectors were calculated.','\n')
#sys.exit()


###--------------------------------------------------------------------------------------------
### Calculation of Normal Stresses:
###--------------------------------------------------------------------------------------------

stress = nele*[0]
for j in range(nele):
    stress[j] = f[j][1]/A[j]

#print(stress)
#sys.exit()


###--------------------------------------------------------------------------------------------
### Print Element Solution:
###--------------------------------------------------------------------------------------------

W2 = open('thermal_output_element_solution.txt','w')
print('Element Solution:', file = W2)
print('%10s' % 'Element','%20s' % 'f1(N)','%20s' % 'f2(N)','%20s' % 'Stress(MPa)', file = W2)
for j in range(nele):
    print('%10d' % (j+1),'%20.6f' % f[j][0],'%20.6f' % f[j][1],'%20.6f' % stress[j], file = W2)
W2.close()
print('Element solution was written.','\n')


