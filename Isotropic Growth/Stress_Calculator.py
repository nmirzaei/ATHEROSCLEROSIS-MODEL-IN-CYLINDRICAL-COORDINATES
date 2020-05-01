from dolfin import *
import numpy as np
import sympy as sym
from sympy import symbols
from sympy import atan2,Abs
from stress import stress
from GradF import GradF
from MaxTenComp import MaxTenComp
from Reader import Reader
import os
import csv

#The following lines take an already converted Gmsh mesh and import it
mesh = Mesh('Mesh.xml')
######################################################################

# Facet functions
Volume = MeshFunction('size_t' , mesh , 'Mesh_physical_region.xml' )  #saves the interior info of the mesh
bnd_mesh = MeshFunction('size_t', mesh , 'Mesh_facet_region.xml')  #saves the boundary info of the mesh
######################################################################


# Optimization options for the form compiler
parameters['form_compiler']['cpp_optimize'] = True
ffc_options = {'optimize': True}
######################################################################

# define function space
V = VectorFunctionSpace(mesh, 'P', 2)
S = FunctionSpace(mesh,'R',0)
T1 = TensorFunctionSpace(mesh,'DG',1,shape=(3,3))
N = V.dim()
d = mesh.geometry().dim()
######################################################################


# Construct integration measure using these markers
ds = Measure('ds', subdomain_data=bnd_mesh)
dx = Measure('dx', subdomain_data=Volume)
######################################################################


# Define functions
u = Function(V)
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
######################################################################

#defining spatial coordinates for changing to cylindrical coordinates
x = SpatialCoordinate(mesh)
######################################################################



#Defining a function that sets all the irrelevant layer values to zero.
def extract_values(u,cell_function,subdomain_id, V):
    dofmap=V.dofmap()
    mesh = V.mesh()
    for cell in cells(mesh):
        if cell_function[cell.index()]!=subdomain_id:
            dofs = dofmap.cell_dofs(cell.index())
            for dof in dofs:
                u.vector()[dof]=0.0
    return u
######################################################################


######################################################################
# Elasticity parameter for the Hopzafel cylinder (intima)
mu_i, nu_i, beta_i, eta_i, rho_i, phi_i = 27.9, 0.49, 170.88, 263.66, 0.51, 60.3*pi/180
######################################################################

# Elasticity parameter for the Hopzafel cylinder (Media)
mu_m, nu_m, beta_m, eta_m, rho_m, phi_m = 1.27, 0.49, 8.21, 21.60, 0.25, 20.61*pi/180
######################################################################

# Elasticity parameter for the Hopzafel cylinder (adventitia)
mu_a, nu_a, beta_a, eta_a, rho_a, phi_a = 7.56, 0.49, 85.03, 38.57, 0.55, 67*pi/180
######################################################################


#Collagen fibers direction in cylindrical coordinates  (before it was (0,cos,sin) I changed it to (0,sin,cos) which I believe is the right one )
b_i = Expression(('0','sin(phi)','cos(phi)'),degree=0,phi=phi_i)
b_m = Expression(('0','sin(phi)','cos(phi)'),degree=0,phi=phi_m)
b_a = Expression(('0','sin(phi)','cos(phi)'),degree=0,phi=phi_a)
######################################################################


#growth for the media
ginv_m = Identity(3)
G_m = Expression('1', degree=0)
######################################################################

#growth for the adventitia
ginv_a = Identity(3)
G_a = Expression('1', degree=0)
######################################################################


#facet Normals for blood pressure direction
n = Expression(('-1','0','0'),degree=0)
######################################################################


#User input for the desirable stress
response = raw_input('How many t do you want to calculate the stresses?')
M = int(response)


response1 = raw_input('Calculate stress every other -- term:')
Skip = int(response1)

response2 = raw_input('At what step do you want the color map:')
M2 = int(response1)
######################################################################


#looping prameters
counter=0
epsilon, sigma, t = 0.0, 16.0, 0.0
deps, dsig, dt= 0.001, 2.0, 0.001
epsmax, sigmax = 0.001, 16.0
j=0
k=0
ID=0
max11c_i = []
max11t_i = []
max22c_i = []
max22t_i = []
max33c_i = []
max33t_i = []
max11c_m = []
max11t_m = []
max22c_m = []
max22t_m = []
max33c_m = []
max33t_m = []
max11c_a = []
max11t_a = []
max22c_a = []
max22t_a = []
max33c_a = []
max33t_a = []
######################################################################

vtkfile1 = File('local location/Sig_11_total_(%d).pvd' %(M2))
vtkfile2 = File('local location/Sig_22_total_(%d).pvd' %(M2))
vtkfile3 = File('local location/Sig_33_total_(%d).pvd' %(M2))

while counter<=M:


 #Reading the desirable displacement and the growth parameters
 if counter%Skip==0:
    u = Function(V,'local location/u%d.xml' %(counter))
    a1, b1, c1 = Reader(counter)
    if counter!=0:
        #finds the isotropic parameter
        alpha = (-1+((1+a1*epsilon)*(1+b1*epsilon)*(1+c1*epsilon))**(0.3333))/epsilon
        ######################################################################
    else:
        alpha = Constant(0)
 ######################################################################

    #Growth tensor
    EXP = Expression('eps*exp(-a*(x[1]-1.68)*(x[1]-1.68))',degree=0,a=14.2857142857,eps=epsilon)
    PARS2 = as_tensor([[alpha,0,0],[0,alpha,0],[0,0,alpha]])
    g_i = Identity(3)+EXP*PARS2
    G_i = det(g_i)
    ginv_i = inv(g_i)
    ######################################################################

# Kinematics for initma
    II = Identity(3)            # Identity tensor
    F = GradF(u,x)
    J = det(F)
    Fe_i= F*ginv_i
    C_i = F.T*F                   # Right Cauchy-Green tensor
    Ce_i= variable(Fe_i.T*Fe_i)
######################################################################


# Kinematics for media
    Fe_m= F*ginv_m
    C_m = F.T*F                   # Right Cauchy-Green tensor
    Ce_m= variable(Fe_m.T*Fe_m)
######################################################################

# Kinematics for adventitia
    Fe_a= F*ginv_a
    C_a = F.T*F                   # Right Cauchy-Green tensor
    Ce_a= variable(Fe_a.T*Fe_a)
######################################################################


# Invariants of deformation tensors (intima)
    I_i = tr(C_i)
    Ie_i=tr(Ce_i)
    J_i = det(Fe_i)   #Jacobian of F or Fe?
######################################################################


# Invariants of deformation tensors (media)
    I_m = tr(C_m)
    Ie_m=tr(Ce_m)
    J_m = det(Fe_m)   #Jacobian of F or Fe?
######################################################################


 # Invariants of deformation tensors (media)
    I_a = tr(C_a)
    Ie_a=tr(Ce_a)
    J_a = det(Fe_a)   #Jacobian of F or Fe?
######################################################################

 #for hopzafel (intima)
    I4_i = conditional(lt(dot(b_i, Ce_i*b_i),Constant(1)),1,dot(b_i, Ce_i*b_i))
 #for hopzafel (media)
    I4_m = conditional(lt(dot(b_m, Ce_m*b_m),Constant(1)),1,dot(b_m, Ce_m*b_m))
 #for hopzafel (adventitia)
    I4_a = conditional(lt(dot(b_a, Ce_a*b_a),Constant(1)),1,dot(b_a, Ce_a*b_a))
######################################################################

 #Expression needed for stress calculation
    H1_i=beta_i*(rho_i*(I4_i-1)**2+(1-rho_i)*(Ie_i-3)**2)
    H1_m=beta_m*(rho_m*(I4_m-1)**2+(1-rho_m)*(Ie_m-3)**2)
    H1_a=beta_a*(rho_a*(I4_a-1)**2+(1-rho_a)*(Ie_a-3)**2)
    H2_i=rho_i*(I4_i-1)*outer(b_i,b_i)+(1-rho_i)*(Ie_i-3)*II
    H2_m=rho_m*(I4_m-1)*outer(b_m,b_m)+(1-rho_m)*(Ie_m-3)*II
    H2_a=rho_a*(I4_a-1)*outer(b_a,b_a)+(1-rho_a)*(Ie_a-3)*II
 #####################################################################

# #Calculating the stress tensors in each layer
    sig_i=(mu_i/J_i)*Fe_i*Fe_i.T+(4*eta_i/J_i)*exp(H1_i)*Fe_i*H2_i*Fe_i.T+(2*nu_i*mu_i*(J_i-1)/(1-2*nu_i))*II-(mu_i/J_i)*II
    sig_m=(mu_m/J_m)*Fe_m*Fe_m.T+(4*eta_m/J_m)*exp(H1_m)*Fe_m*H2_m*Fe_m.T+(2*nu_m*mu_m*(J_m-1)/(1-2*nu_m))*II-(mu_m/J_m)*II
    sig_a=(mu_a/J_a)*Fe_a*Fe_a.T+(4*eta_a/J_a)*exp(H1_a)*Fe_a*H2_a*Fe_a.T+(2*nu_a*mu_a*(J_a-1)/(1-2*nu_a))*II-(mu_a/J_a)*II
####################################################################

# #projecting the calculated stresses
    Sigma_i = project(sig_i,T1)
    Sigma_m = project(sig_m,T1)
    Sigma_a = project(sig_a,T1)
 #####################################################################

 #Setting all the irrelevent layer values equal to zero
    Sig_i = Function(T1); Sig_i.vector()[:]=Sigma_i.vector()
    Sig_i = extract_values(Sig_i,Volume,7,T1)
    Sig_m = Function(T1); Sig_m.vector()[:]=Sigma_m.vector()
    Sig_m = extract_values(Sig_m,Volume,8,T1)
    Sig_a = Function(T1); Sig_a.vector()[:]=Sigma_a.vector()
    Sig_a = extract_values(Sig_a,Volume,9,T1)
 #######################################################################

    if counter==M2:
        #This function creates a color map of stress on the domain
        stress(u,t,Sig_i,Sig_m,Sig_a,vtkfile1,vtkfile2,vtkfile3)
        ######################################################################

 #Sending to this function for calcualting the maximum tensile and compressive stresses
    i1, i2, i3, i4, i5, i6, m1, m2, m3, m4, m5, m6, a1, a2, a3, a4, a5, a6 = MaxTenComp(u,t,Sig_i,Sig_m,Sig_a)
    max11c_i.append(i1)
    max11t_i.append(i2)
    max22c_i.append(i3)
    max22t_i.append(i4)
    max33c_i.append(i5)
    max33t_i.append(i6)
    max11c_m.append(m1)
    max11t_m.append(m2)
    max22c_m.append(m3)
    max22t_m.append(m4)
    max33c_m.append(m5)
    max33t_m.append(m6)
    max11c_a.append(a1)
    max11t_a.append(a2)
    max22c_a.append(a3)
    max22t_a.append(a4)
    max33c_a.append(a5)
    max33t_a.append(a6)
 #######################################################################


#making sure we start saving results after the pressure is imposed (if there is any pressure). Also making sure we save results every 10 steps
 if sigma<=sigmax:
    sigma+=dsig
 if sigma>sigmax:
    epsilon+=deps
    t+=dt
    counter+=1
    ID+=1
 if ID!=0:
    sigma-=dsig
    ID=0
######################################################################

#printing the current loop variables
 print(epsilon,flush=True)
 print(t,flush=True)
 print(sigma,flush=True)
######################################################################

######################################################################
 k+=1
######################################################################



#Export maximum stress in csv files
import csv

with open('max11c_i.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max11c_i)

with open('max11t_i.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max11t_i)

with open('max22c_i.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max22c_i)

with open('max22t_i.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max22t_i)

with open('max33c_i.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max33c_i)

with open('max33t_i.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max33t_i)

with open('max11c_m.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max11c_m)

with open('max11t_m.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max11t_m)

with open('max22c_m.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max22c_m)

with open('max22t_m.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max22t_m)

with open('max33c_m.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max33c_m)

with open('max33t_m.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max33t_m)

with open('max11c_a.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max11c_a)

with open('max11t_a.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max11t_a)

with open('max22c_a.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max22c_a)

with open('max22t_a.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max22t_a)

with open('max33c_a.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max33c_a)

with open('max33t_a.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(max33t_a)
######################################################################
