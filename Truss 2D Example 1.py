### FEM Program for Solving System of Linear Spring Problems

import math, numpy, sys

node = 4                                            # Number of Nodes
nele = 5                                            # Number of Elements

n = 2*node                                          # Number of Equations

u = n*[0]                                           # Initial Guess for Solution Vector
F = n*[0]                                           # Initializing of Load Vector

kele  = 0*numpy.identity(4)                         # Element Matrix of Two-Node Line Element
f     = [4*[0] for i in range(nele)]                # Initializing of Element Load Vectors

N   = 100                                           # Maximum Iteration of PCGNE Solver
tol = 1.e-9                                         # Tolerance for PCGNE Solver


###--------------------------------------------------------------------------------------------
### Read Node Data:
###--------------------------------------------------------------------------------------------

R1 = open('truss_input_node_data.txt') #c1=ux ,c2=Fx ,c3=value ,c4=v=uy ,c5=Fy ,c6=value
text = R1.readline()
inode, xc, yc, c1, c2, c3, c4, c5, c6 = [],[],[],[],[],[],[],[],[]
for line in R1.readlines():
    fields = line.split()
    inode.append(int(fields[0]))
    xc.append(float(fields[1]))
    yc.append(float(fields[2]))
    c1.append(float(fields[3]))
    c2.append(float(fields[4]))
    c3.append(float(fields[5]))
    c4.append(float(fields[6]))
    c5.append(float(fields[7]))
    c6.append(float(fields[8]))
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

R2 = open("truss_input_element_data.txt")
text = R2.readline()
iele, n1, n2, A, E = [],[],[],[],[]
for line in R2.readlines():
    fields = line.split()
    iele.append(int(fields[0]))
    n1.append(int(fields[1]))
    n2.append(int(fields[2]))
    A.append(float(fields[3]))
    E.append(float(fields[4]))
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

def K_element(j,nele,n1,n2,xc,yc,A,E):
    i1 = n1[j] - 1
    i2 = n2[j] - 1
    dx = xc[i2] - xc[i1]
    dy = yc[i2] - yc[i1]
    dl = (dx*dx + dy*dy)**0.5
    C  = dx/dl
    S  = dy/dl
    k  = A[j]*E[j]/dl*1000
    kele[0][0] =  k*C*C
    kele[0][1] =  k*C*S
    kele[0][2] = -k*C*C
    kele[0][3] = -k*C*S
    kele[1][0] =  kele[0][1]
    kele[1][1] =  k*S*S
    kele[1][2] = -k*C*S
    kele[1][3] = -k*S*S
    kele[2][0] =  kele[0][2]
    kele[2][1] =  kele[1][2]
    kele[2][2] =  k*C*C
    kele[2][3] =  k*C*S
    kele[3][0] =  kele[0][3]
    kele[3][1] =  kele[1][3]
    kele[3][2] =  kele[2][3]
    kele[3][3] =  k*S*S
    return i1, i2, kele

def assem_K(n,nele,n1,n2,xc,yc,A,E):
    K = 0*numpy.identity(2*n)
    for j in range(nele):
        i1, i2, kele = K_element(j,nele,n1,n2,xc,yc,A,E)
        K[2*i1]  [2*i1]   = K[2*i1]  [2*i1]   + kele[0][0]
        K[2*i1]  [2*i1+1] = K[2*i1]  [2*i1+1] + kele[0][1]
        K[2*i1]  [2*i2]   = K[2*i1]  [2*i2]   + kele[0][2]
        K[2*i1]  [2*i2+1] = K[2*i1]  [2*i2+1] + kele[0][3]
        K[2*i1+1][2*i1]   = K[2*i1+1][2*i1]   + kele[1][0]
        K[2*i1+1][2*i1+1] = K[2*i1+1][2*i1+1] + kele[1][1]
        K[2*i1+1][2*i2]   = K[2*i1+1][2*i2]   + kele[1][2]
        K[2*i1+1][2*i2+1] = K[2*i1+1][2*i2+1] + kele[1][3]
        K[2*i2]  [2*i1]   = K[2*i2]  [2*i1]   + kele[2][0]
        K[2*i2]  [2*i1+1] = K[2*i2]  [2*i1+1] + kele[2][1]
        K[2*i2]  [2*i2]   = K[2*i2]  [2*i2]   + kele[2][2]
        K[2*i2]  [2*i2+1] = K[2*i2]  [2*i2+1] + kele[2][3]
        K[2*i2+1][2*i1]   = K[2*i2+1][2*i1]   + kele[3][0]
        K[2*i2+1][2*i1+1] = K[2*i2+1][2*i1+1] + kele[3][1]
        K[2*i2+1][2*i2]   = K[2*i2+1][2*i2]   + kele[3][2]
        K[2*i2+1][2*i2+1] = K[2*i2+1][2*i2+1] + kele[3][3]
    return K

K = assem_K(n,nele,n1,n2,xc,yc,A,E)
for j in range(n): K[j][n+j] = -1
print('System stiffness matrix was calculated.','\n')
#print(K)
#sys.exit()


###--------------------------------------------------------------------------------------------
### Apply Boundary Conditions:
###--------------------------------------------------------------------------------------------

def apply_bc(n,node,c1,c2,c3,c4,c5,c6,K):
    b = 2*n*[0]
    for j in range(node):
        K[n+2*j][  2*j] = c1[j]
        K[n+2*j][n+2*j] = c2[j]
        b[n+2*j]        = c3[j]
        K[n+2*j+1][  2*j+1] = c4[j]
        K[n+2*j+1][n+2*j+1] = c5[j]
        b[n+2*j+1]          = c6[j]
    return K, b

K, b = apply_bc(n,node,c1,c2,c3,c4,c5,c6,K)
print('Boundary conditions were specified.','\n')
#print(K)
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

W1 = open('truss_output_node_solution.txt','w')
print('Node Solution:', file = W1)
print('%10s' % 'Node','%18s' % 'u(mm)','%18s' % 'v(mm)',
      '%18s' % 'Fx(N)','%18s' % 'Fy(N)', file = W1)
for i in range(node):
    print('%10d' % (i+1),'%18.6f' % u[2*i],'%18.6f' % u[2*i+1],
          '%18.6f' % F[2*i],'%18.6f' % F[2*i+1], file = W1)
W1.close()
print('Node solution was written.','\n')
#sys.exit()


###--------------------------------------------------------------------------------------------
### Calculate Element Load Vectors:
###--------------------------------------------------------------------------------------------

def cal_element(nele,n1,n2,u,xc,yc,A,E):
    for j in range(nele):
        i1, i2, kele = K_element(j,nele,n1,n2,xc,yc,A,E)
        uele = [u[2*i1],u[2*i1+1],u[2*i2],u[2*i2+1]]
        f[j] = numpy.dot(kele,uele)
    return f

f = cal_element(nele,n1,n2,u,xc,yc,A,E)
print('Element load vectors were calculated.','\n')
#sys.exit()


###--------------------------------------------------------------------------------------------
### Calculation of Normal Stresses:
###--------------------------------------------------------------------------------------------

stress = nele*[0]
T      = 0*numpy.identity(4)
for j in range(nele):
    i1 = n1[j] - 1
    i2 = n2[j] - 1
    dx = xc[i2] - xc[i1]
    dy = yc[i2] - yc[i1]
    dl = (dx*dx + dy*dy)**0.5
    C  = dx/dl
    S  = dy/dl
    T[0][0] =  C
    T[0][1] =  S
    T[1][0] = -S
    T[1][1] =  C
    T[2][2] =  C
    T[2][3] =  S
    T[3][2] = -S
    T[3][3] =  C
    uele = [u[2*i1],u[2*i1+1],u[2*i2],u[2*i2+1]]
    transform = numpy.dot(T,uele)
    stress[j] = E[j]*1000/dl*numpy.dot([-1,0,1,0],transform)

#print(stress)
#sys.exit()
    

###--------------------------------------------------------------------------------------------
### Print Element Solution:
###--------------------------------------------------------------------------------------------

W2 = open('truss_output_element_solution.txt','w')
print('Element Solution:', file = W2)
print('%10s' % 'Element','%18s' % 'fx1(N)','%18s' % 'fy1(N)',
      '%18s' % 'fx2(N)','%18s' % 'fy2(N)','%18s' % 'Stress(MPa)', file = W2)
for j in range(nele):
    print('%10d' % (j+1),'%18.6f' % f[j][0],'%18.6f' % f[j][1],
          '%18.6f' % f[j][2],'%18.6f' % f[j][3],'%18.6f' % stress[j], file = W2)
W2.close()
print('Element solution was written.','\n')


