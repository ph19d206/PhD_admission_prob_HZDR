from math import*
import matplotlib.pyplot as plt
import random

#define zero matrix of row by colloum
def zeros(row, col):
    return [[0 for _ in range(col)] for _ in range(row)]




# Matrix multiplication
def dot(Amatrix,Bmatrix):
    A_row=len(Amatrix)
    A_col=len(Amatrix[0])
    B_row=len(Bmatrix)
    B_col=len(Bmatrix[0])
    if A_col==B_row:
        C_Zero= zeros(A_row,B_col)
        for i in range(A_row):
            for j in range(B_col):
                A_dder=0
                for k in range(A_col):
                    A_dder=A_dder+ Amatrix[i][k]*Bmatrix[k][j]
                C_Zero[i][j]=A_dder
    else:
        print("Row colloum not some, recheck the matrix")
    return C_Zero

# N by N Identity matrix
def Identity(D_imension):
    I_m_atrix=zeros(D_imension,D_imension)

    for i in range(D_imension):
        I_m_atrix[i][i]=1
    return I_m_atrix


# scalar multipli of matrix
def scalardot(scalar,matrix):
    
    Copy_matrix= zeros(len(matrix),len(matrix[0]))
    # make a copy of matrix
    for i in range(len(Copy_matrix)):
        for j in range(len(Copy_matrix[0])):
            Copy_matrix[i][j]=matrix[i][j]
    
    for i in range(len(Copy_matrix)):
        for j in range(len(Copy_matrix[0])):
            Copy_matrix[i][j]=scalar*Copy_matrix[i][j]
    return Copy_matrix


# substruct matrix
def submatrix(a,b):
    A_row=len(a)
    A_col=len(a[0])
    B_row=len(b)
    B_col=len(b[0])
    if A_row==B_row and A_col==B_col:
        Substracted_Matrix=zeros(len(a),len(a[0]))
        for i in range(len(a)):
            for j in range(len(a[0])):
                Substracted_Matrix[i][j]=a[i][j]-b[i][j]
        return Substracted_Matrix
    else:
        print("Recheck the substructed matrix.")
        
        
        
# doing complex conjugate
def dagger(a):
    Conjugate_Matrix=zeros(len(a[0]),len(a))
    
    for i in range(len(Conjugate_Matrix)):
            for j in range(len(Conjugate_Matrix[0])):
                Conjugate_Matrix[i][j]=complex.conjugate(complex(a[j][i]))
    return Conjugate_Matrix


# addstruct matrix
def addmatrix(a,b):
    A_row=len(a)
    A_col=len(a[0])
    B_row=len(b)
    B_col=len(b[0])
    if A_row==B_row and A_col==B_col:
        Substracted_Matrix=zeros(len(a),len(a[0]))
        for i in range(len(a)):
            for j in range(len(a[0])):
                Substracted_Matrix[i][j]=a[i][j]+b[i][j]
        return Substracted_Matrix
    else:
        print("Recheck the adding matrices.")
        
        
# Inner product
def inner(a,b):
    Inner_Product_list=dot(dagger(a),b)
    return Inner_Product_list[0][0]
    
        
        
# ITM is the iteration function

def ITM(psi):
    psi=dot(H,psi)
    norm = sqrt(abs(inner(psi,psi)))
    psi=scalardot((1/norm),psi)
    return (psi)



# Energy iteration

def Eit(psi):
    Hpsi= dot(H,psi)
    Hexpectation=inner(psi,Hpsi)
    norm1 = sqrt(abs(inner(psi,psi)))
    energy= Hexpectation/norm1
    energy=energy.real
    return (energy)



# gram smidt of psi0 and psi1 in two dimension
def GramSmidt(psi0,psi1):
    innerproductg =inner(psi0,psi1)
    sqnorm = abs(inner(psi0,psi0))
    psi1=submatrix(psi1,scalardot((innerproductg/sqnorm),psi0))
    psi1=scalardot((1/sqrt(abs(inner(psi1,psi1)))),psi1)
    return psi1
    
    
# parameter value and other initial value
alpha=1.3 # Alpha value
Nmax=100 # order of matrix
H=zeros(Nmax,Nmax)
H1=zeros(Nmax,Nmax)
lam=35
accuracy=10**(-6)


# making hamiltonian matrix
for i in range(Nmax):
    for j in range(Nmax):
        p1=(float((i+1)+(j+1)))**0.52
        p2=log((i+1)+(j+1)+alpha)
        H1[i][j]=sin(p1/p2)
        

H=submatrix(H1,scalardot(lam,Identity(Nmax)))

                                                                                           



# Initial state (random )
psi0=zeros(Nmax,1)
for i in range(len(psi0)):
    psi0[i][0]=random.random() + random.random()*1j
    
psi0=scalardot((1/sqrt(abs(inner(psi0,psi0)))),(psi0))



# Calculating energy at each Iteration and storing in array of E
i=0
N=[]
E=[]


while 1==1:
    N.append(i)
    E.append(Eit(psi0))
    ppsi0=psi0
    psi0=ITM(psi0)
    #print(abs(np.dot(np.matrix.getH(ppsi0-psi0),(ppsi0-psi0))))#testing for state convergance
    if i !=0 and abs(E[i]-E[i-1])<accuracy:
        break
    i = 1+i

    

print("Ground state energy= ",E[-1]+lam)                              





#ploting data of error vs number of iteration for ground state
error=zeros(len(N)-1,1)


NI = zeros(len(N)-1,1)# Iteration number
for i in range(len(N)-1):
    error[i][0] = abs((E[i]-E[i+1])/E[-1])
    NI[i][0]=i+1





# 1st excited state
psi1=zeros(Nmax,1)
for i in range(len(psi1)):
    psi1[i][0]=random.random() + random.random()*1j
    
psi1=scalardot((1/sqrt(abs(inner(psi1,psi1)))),(psi1))


j=0
N1=[]
E1=[]


while 1==1:
    N1.append(j)
    psi1=GramSmidt(psi0,psi1)# doing orthogonal with ground state
    E1.append(Eit(psi1))
    ppsi1=psi1
    psi1=ITM(psi1)# Iteration of 1st energy state
    if j !=0 and abs(E1[j]-E1[j-1])<accuracy:
        break
    j = 1+j

print("1st excited state energy= ",E1[-1]+lam)                          






#ploting data of error vs number of iteration for first excited state
error1=zeros(len(N1)-1,1)

NI1 = zeros(len(N1)-1,1)# Iteration number
for i in range(len(N1)-1):
    error1[i][0] = abs((E1[i]-E1[i+1])/E1[-1])
    NI1[i][0]=i+1



#plots

#plt.plot(NI,error,'.')                                                
#plt.xlabel("number of Iteration")
#plt.ylabel("Error")
#plt.title("error vs iteration for ground state energy")
#plt.savefig("ground state Error vs Iteration")






plt.plot(NI1,error1,'.')
plt.xlabel("number of Iteration")
plt.ylabel("Error")  
plt.title("error vs iteration for First excited state energy")
plt.savefig("excited state Error vs Iteration")
