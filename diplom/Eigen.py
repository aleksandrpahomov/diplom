import matplotlib.pyplot as plt
from scipy.integrate import quad
from math import cos, pi, exp, inf, sin, sqrt,cosh,sinh
from  sympy import *
import numpy as np
from decimal import Decimal, getcontext
getcontext().prec=40

c_par=10000
d=0.025
l=0.5
d_1=0.08
l_1=0.3


M=(d_1**2*l_1)/(d**2*l)
J=(l_1**2/3+d_1**2/16)/l**2
a=l_1/(2*l)
I_par=1

fig, ax = plt.subplots()


def solveFunc(A,B,C,D,beta,x):
    return A * cosh(beta * x) + B * sinh(beta * x) + C * cos(beta * x) + D * sin(beta * x)

def solveFuncFirstPartial(A,B,C,D,beta,x):
    return beta*(A * sinh(beta * x) + B * cosh(beta * x) - C * sin(beta * x) + D * cos(beta * x))

def Norm(A,B,C,D,beta):
    func = lambda x:pow(solveFunc(A,B,C,D,beta,x),2)
    return  sqrt(quad(func,0,1)[0] +M*pow(solveFunc(A,B,C,D,beta,1),2) + M*J*pow(solveFuncFirstPartial(A,B,C,D,beta,1),2)+ M*a*2*(solveFunc(A,B,C,D,beta,1)*solveFuncFirstPartial(A,B,C,D,beta,1)))

def initMatrixForSystem (beta):
    L=beta**4
    matrix= np.array([
        [
            -beta**2*(cosh(beta)+cos(beta))+ L*M*J*beta*(sinh(beta)+sin(beta)) + a*L*M*(cosh(beta)-cos(beta)) - c_par*a*(cosh(beta)-cos(beta)) - a*a*c_par*beta*(sinh(beta) + sin(beta))
            ,
            -beta ** 2 * (sinh(beta) + sin(beta)) + L * M * J * beta * (cosh(beta) - cos(beta)) + a * L * M*(sinh(beta)-sin(beta)) - c_par*a*(sinh(beta)-sin(beta)) - a*a*c_par*beta*(cosh(beta)-cos(beta))
        ],
        [
            -beta ** 3 * (sinh(beta) - sin(beta)) - (L * M -c_par)* (cosh(beta) - cos(beta)) - beta*a * (L * M -c_par) * (sinh(beta) + sin(beta))
            ,
            -beta ** 3 * (cosh(beta) + cos(beta)) - (L * M -c_par) * (sinh(beta) - sin(beta)) - beta * a * (L * M -c_par)* (cosh(beta) - cos(beta))
        ]
    ], float)
    return matrix


def Det(matrix):
    return np.linalg.det(matrix)


def MPD(left:Decimal,right:Decimal):
    eps=0.000001
    matrixLeft = initMatrixForSystem(left)
    matrixRight = initMatrixForSystem(right)
    a = Det(matrixLeft)
    b=Det(matrixRight)
    midDet = 1
    if a*b > 0:
        return -1
    while abs(midDet)>eps:
        matrixLeft = initMatrixForSystem(left)
        matrixRight = initMatrixForSystem(right)
        a = Det(matrixLeft)*exp(-2*left)
        b = Det(matrixRight)*exp(-2*right)
        midpoint = (left + right) / 2.0
        if(abs(left-midpoint)<Decimal('0.000000000000000001')):
           return -1
        matrixMid = initMatrixForSystem(midpoint)
        midDet = Det(matrixMid)*exp(-2*midpoint)
        if abs(midDet) < eps:
           if (abs(a) < eps):
                return a
        if (abs(b) < eps):
            return b
        if  a * midDet <0:
            right = midpoint
        if  b*midDet < 0:
            left = midpoint
    return midpoint



def drawFunc(A,B,C,D,beta):
    norma=Norm(A,B,C,D,beta)
    #print(norma)
    x=[]
    val=[]
    h=0.01
    i=0
    print(A/norma," ",B/norma, " ", C/norma, " ", D/norma)
    while i<1:
        x.append(i)
        val.append(solveFunc(A/norma,B/norma,C/norma,D/norma,beta,i))
        i+=h
    ax.plot(x,val)

currentBeta =0.01
count = 0
betas = []
determinant = []
iter = 0
left =currentBeta
right= currentBeta
i=10
x=[]
y=[]


while count != 5:
    currentBeta+=0.01
    right = currentBeta
    #print('dsfsff =========== ',bet)
    nul = MPD(left, right)
    if(nul!=-1):
        betas.append(nul)
        print(nul)
        A = initMatrixForSystem(nul)
        print(np.linalg.det(A))
        count += 1
    left=right
print(betas)
for b in betas:
    A=initMatrixForSystem(b)
    print(np.linalg.det(A))
    B=np.array([0.0, 0.000001])
    vect = np.linalg.solve(A,B)
    print(vect)
    drawFunc(-vect[0],-vect[1],vect[0],vect[1],b)
b=betas[0]
plt.grid=True
plt.ylim(-2,2)
#plt.title("c= "+str(c_par)+"  "+ "d="+str(d)+"  " + "l= "+ str(l) +"  " +"d1= "+ str(d_1)+ "  "+ "l1= "+str(l_1)+"  ")
plt.legend(['b1='+str('%.3f' % betas[0]),'b2='+str('%.3f' % betas[1]),'b3='+str('%.3f' % betas[2]),'b4='+str( '%.3f' % betas[3]),'b5='+str( '%.3f' % betas[4])])
plt.show()
