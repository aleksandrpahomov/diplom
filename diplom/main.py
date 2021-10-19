import matplotlib.pyplot as plt
from scipy.integrate import quad
from math import cos, pi, exp, inf, sin, sqrt
from decimal import Decimal, getcontext
from math import cos, pi, inf, sin, sqrt,cosh,sinh
from cmath import *
from  sympy import *
import numpy as np

getcontext().prec=40

fig = plt.figure()

"""
beta_values=[0.6915441894531253, 2.38301791429519, 4.987424687072574, 7.981025818548971, 11.076666303282645]
A_koef=[0.733692144174013,0.744583564070441,0.960026437505121,0.988655041866321,
0.995152213388382]
B_koef=[-0.807762180846274,-0.855799387858502,-0.948412398079204,-0.989283131837640 ,
-0.995123061305420]
C_koef=[-0.733692144174013,-0.744583564070441,-0.960026437505121,-0.988655041866321,
-0.995152213388382]
D_koef=[0.807762180846274,0.855799387858502,0.948412398079204,0.989283131837640 ,
0.995123061305420 ]


# c=  1
beta_values=[0.7786181640625006, 2.3830430293083116, 4.9874290977417814, 7.981026398017876, 11.076666426063557]
A_koef=[-0.576640462074326,0.744604201369003,0.960024741053652,0.988654875328396,
0.995152186926061]
B_koef=[0.560203020438132,-0.855798625296382,-0.948410862096128,-0.989282964361043 ,
-0.995123034851658]
C_koef=[0.576640462074326,-0.744604201369003,-0.960024741053652,-0.988654875328396,
-0.995152186926061]
D_koef=[-0.560203020438132,0.855798625296382,0.948410862096128,0.989282964361043 ,
0.995123034851658]


# c=500
beta_values=[2.3733526372909477, 2.9758718597888745, 4.9902230392395825, 7.98135405104595, 11.076734880357794]
A_koef=[-0.726375397284920 ,0.156793099082760,-0.958798286146876,-0.988559109987661,
0.995137367532296 ]
B_koef=[0.844319704752901,-0.139855648162372,0.947287610782835,0.989186667379068 ,
-0.995108220231743 ]
C_koef=[0.726375397284920 ,-0.156793099082760,0.958798286146876,0.988559109987661,
-0.995137367532296 ]
D_koef=[-0.844319704752901,0.139855648162372,-0.947287610782835,-0.989186667379068 ,
0.995108220231743]

"""
# c=  10000
beta_values=[2.377160571180278, 4.9561120951920135, 6.313065616488366, 7.991526915021116,
             11.078192385937829]
A_koef=[0.735959221477949,0.952683429611060,-0.174908574554219,0.984017471587358,
        -0.994790999415068]
B_koef=[-0.851625217867241,-0.940151190706570,0.175617524678446,-0.984627683174865,
        0.994761954550826]
C_koef=[-0.735959221477949,-0.952683429611060,0.174908574554219,-0.984017471587358,
        0.994790999415068]
D_koef=[0.851625217867241,0.940151190706570,-0.175617524678446,0.984627683174865,
        -0.994761954550826]




V_params=[[A_koef[0],B_koef[0],C_koef[0],D_koef[0]],[A_koef[1],B_koef[1],C_koef[1],D_koef[1]],[A_koef[2],B_koef[2],C_koef[2],D_koef[2]],[A_koef[3],B_koef[3],C_koef[3],D_koef[3]],[A_koef[4],B_koef[4],C_koef[4],D_koef[4]]]

a=[[1.0,1.0,1.0,1.0,1.0],
   [1.0,1.0,1.0,1.0,1.0],[1.0,1.0,1.0,1.0,1.0],
   [1.0,1.0,1.0,1.0,1.0],[1.0,1.0,1.0,1.0,1.0]] # a[0][0]= koef a_11 for matrix
omega_1=[beta_values[0]**2,beta_values[1]**2,beta_values[2]**2,beta_values[3]**2,beta_values[4]**2] #array with omega values
mu = 0.25

delta =0.5
p =0.5
alpha = 0.75
d=0.025
l=0.5
d_1=0.08
l_1=0.3
M=float((d_1**2*l_1)/(d**2*l))
J=float((l_1**2/3+d_1**2/16)/(l**2))
S_1=(3.14*d_1*l_1)/(l**2)
S_2=(d_1*l_1)/(d*l)
#M=1
#J=1
#a=1
c_par=10000
a_par=l_1/(2*l)
#I_par=(pi*(d**4))/16
I_par=1
fig, ax = plt.subplots()
def solveFunc(V_koef,beta,x):
    return V_koef[0] * cosh(beta * x) + V_koef[1] * sinh(beta * x) + V_koef[2] * cos(beta * x) + V_koef[3] * sin(beta * x)

def solveFuncFirstPartial(V_koef,beta,x):
    return beta*(V_koef[0]* sinh(beta * x) + V_koef[1] * cosh(beta * x) - V_koef[2] * sin(beta * x) + V_koef[3] * cos(beta * x))

def scalarMult(V1_koef,beta_1,V2_koef,beta_2):
    func = lambda x:solveFunc(V1_koef,beta_1,x)*solveFunc(V2_koef,beta_2,x)
    return float(quad(func,0,1)[0] + M*solveFunc(V1_koef,beta_1,1)*solveFunc(V2_koef,beta_2,1) + M * J * solveFuncFirstPartial(V1_koef,beta_1,1)*solveFuncFirstPartial(V2_koef,beta_2,1)+ M*a_par*(solveFuncFirstPartial(V1_koef,beta_1,1)*solveFunc(V2_koef,beta_2,1)+solveFuncFirstPartial(V2_koef,beta_2,1)*solveFunc(V1_koef,beta_1,1)))


def computeCoefficientsForMatrix(n,m):
    V_n_partial=lambda x: solveFuncFirstPartial(V_params[n],beta_values[n],x)
    V_m_partial= lambda x: solveFuncFirstPartial(V_params[m],beta_values[m],x)
    V_n = lambda  x : solveFunc(V_params[n],beta_values[n],x)
    V_m = lambda  x : solveFunc(V_params[m],beta_values[m],x)
    return scalarMult(V_params[n],beta_values[n],V_params[m],beta_values[m])+S_1*(J*V_n_partial(1)+ a_par*V_n(1))*V_m_partial(1)+S_2*(V_n(1)+a_par*V_n_partial(1))*V_m(1)

print(type(a[0][0]))
print(type(computeCoefficientsForMatrix(0,0)))

for i in range(len(a)):
    for j in range(len(a)):
        a[i][j]=computeCoefficientsForMatrix(i,j)
        print('A[',i+1,', ',j+1,']= ',a[i][j])

func = lambda x: exp(-p * (x ** (1 + delta))) / (x ** alpha)
answer = quad(func, 0.1, 1000)[0]



def P(i,j,sigma,omega_small,mu,E_s,E_c,a0_prev,bigOmega_prev):
    sigma=complex(sigma)
    omega_small = complex(omega_small)
    mu = complex(mu)
    E_s = complex(E_s)
    E_c = complex(E_c)
    a0_prev = complex(a0_prev)
    bigOmega_prev = complex(bigOmega_prev)

    Res = complex(-(sigma + bigOmega_prev) ** 2 + complex(0,1) * a0_prev * a[i-1][j-1] * (sigma + bigOmega_prev) - omega_small ** 2 * (
                1 - mu * E_c - mu * complex(0,1) * E_s))
    return Res;

def solveDop(m,n,sigma,omega_small,mu,E_s,E_c,a0_prev,bigOmega_prev):
    i_new=0
    j_new =0
    zero= complex(1,1)
    M = np.array([[zero,zero,zero,zero],[zero,zero,zero,zero],[zero,zero,zero,zero],[zero,zero,zero,zero]],complex)
    for i in range (5):
        if (i!=m and i_new <4):
            for j in range (5):
                if( j!=n and j_new <4):
                    if(i==j):
                        M[i_new][j_new]= P(i+1,j+1,sigma,omega_small,mu,E_s,E_c,a0_prev,bigOmega_prev)
                    else:
                        M[i_new][j_new] = complex(0,sigma)*a0_prev*complex(a[i][j])
                    #print(i,j)
                    j_new+=1

            i_new+=1
            j_new=0
    #print(M)
    Res=np.linalg.det(M)
    return Res



def Phi_1(sigma,omega_small,mu,E_s,E_c,a0_prev,bigOmega_prev):
    sigma = complex(sigma)
    omega_small = complex(omega_small)
    mu = complex(mu)
    E_s = complex(E_s)
    E_c = complex(E_c)
    a0_prev = complex(a0_prev)
    bigOmega_prev = complex(bigOmega_prev)
    L=complex(0,1)*sigma
    dop_12= solveDop(0,1,sigma,omega_small,mu,E_s,E_c,a0_prev,bigOmega_prev)
    dop_13= solveDop(0,2,sigma,omega_small,mu,E_s,E_c,a0_prev,bigOmega_prev)
    dop_14= solveDop(0,3,sigma,omega_small,mu,E_s,E_c,a0_prev,bigOmega_prev)
    dop_15= solveDop(0,4,sigma,omega_small,mu,E_s,E_c,a0_prev,bigOmega_prev)
    dop_11= solveDop(0,0,sigma,omega_small,mu,E_s,E_c,a0_prev,bigOmega_prev)
    Res=complex((a0_prev*a[0][1]*dop_12 - a0_prev*a[0][2]*dop_13+a0_prev*a[0][3]*dop_14-a0_prev*a[0][4]*dop_15)/(dop_11))
    return Res;


def A0_k(sigma,omega_small,mu,E_s,E_c,Phi_par):
    znam=sqrt(omega_small**2*(1.0-mu*E_c)-Phi_par.real)
    if (znam.is_Mul):
        return 0
    Res = (-(omega_small**2)*mu*E_s+Phi_par.imag)/znam

    return Res
def Omega_big_k(sigma,omega_small,mu,E_s,E_c,Phi):
    znam=sqrt(omega_small**2*(1.0-mu*E_c)-Phi.real)
    if (znam.is_Mul):
        return 0
    Res = -sigma+znam
    return Res


E_c = 0
E_s =0


sigma =-133
L=complex(0,sigma)


Phi_11=complex(0);
Phi_12=complex(0);
Phi_13=complex(0);

a0 =[[0.0],[0.0],[0.0],[0.0],[0.0]]
bigOmega= [[0.0],[0.0],[0.0],[0.0],[0.0]]

sigma_list=[]
while sigma <=0:
    sigma_list.append(sigma)
    if (sigma < -3):
        sigma += 5
    else:
        sigma += 0.1
sigma_list.append(0.0)
for sigma in sigma_list:
    print(sigma)
    ecos = lambda x: (exp(-p * (x ** (1 + delta))) * cos(sigma * x)) / ((x ** alpha) * answer)
    E_c = quad(ecos, 0.01, 500)[0]
    esin = lambda x: (exp(-p * (x ** (1 + delta))) * sin(sigma * x)) / ((x ** alpha) * answer)
    E_s = quad(esin, 0.01, 500)[0]
    Phi_11=Phi_1(sigma,omega_1[0],mu,E_s,E_c,a0[0][len(a0[0])-1],bigOmega[0][len(bigOmega[0])-1])
    Phi_12=Phi_1(sigma,omega_1[1],mu,E_s,E_c,a0[1][len(a0[1])-1],bigOmega[1][len(bigOmega[1])-1])
    Phi_13=Phi_1(sigma,omega_1[2],mu,E_s,E_c,a0[2][len(a0[2])-1],bigOmega[2][len(bigOmega[2])-1])
    Phi_14 =Phi_1(sigma,omega_1[3],mu,E_s,E_c,a0[3][len(a0[3])-1],bigOmega[3][len(bigOmega[3])-1])
    Phi_15 = Phi_1(sigma, omega_1[4], mu, E_s, E_c, a0[4][len(a0[4]) - 1], bigOmega[4][len(bigOmega[4]) - 1])
    a0[0].append((A0_k(sigma,omega_1[0],mu,E_s,E_c,Phi_11)))
    a0[1].append((A0_k(sigma, omega_1[1], mu, E_s, E_c, Phi_12)))
    a0[2].append((A0_k(sigma, omega_1[2], mu, E_s, E_c, Phi_13)))
    a0[3].append((A0_k(sigma, omega_1[3], mu, E_s, E_c, Phi_14)))
    a0[4].append((A0_k(sigma, omega_1[2], mu, E_s, E_c, Phi_15)))
    bigOmega[0].append((Omega_big_k(sigma,omega_1[0],mu,E_s,E_c,Phi_11)))
    bigOmega[1].append((Omega_big_k(sigma, omega_1[1], mu, E_s, E_c, Phi_12)))
    bigOmega[2].append((Omega_big_k(sigma, omega_1[2], mu, E_s, E_c, Phi_13)))
    bigOmega[3].append((Omega_big_k(sigma, omega_1[3], mu, E_s, E_c, Phi_14)))
    bigOmega[4].append((Omega_big_k(sigma, omega_1[4], mu, E_s, E_c, Phi_15)))

    #print(type(A0_k(sigma,omega[0],mu,E_s,E_c,Phi_1)))


for i in range(5):
    a0[i].pop(0)
    zero = bigOmega[i].pop(0)

print(a0[2])
print(bigOmega[2])
plt.plot(a0[0], bigOmega[0],c="red")
plt.plot(a0[1], bigOmega[1],c="green")
plt.plot(a0[2], bigOmega[2],c="orange")
plt.plot(a0[3], bigOmega[3],c="black")
plt.plot(a0[4], bigOmega[4],c="blue")
#plt.title("c= "+str(c_par)+"  "+ "d="+str(d)+"  " + "l= "+ str(l) +"  " +"d1= "+ str(d_1)+ "  "+ "l1= "+str(l_1)+"  ")
plt.ylim(0,125)
plt.xlabel("a0")
plt.ylabel("Omega")
plt.show()

plt.grid(True)



