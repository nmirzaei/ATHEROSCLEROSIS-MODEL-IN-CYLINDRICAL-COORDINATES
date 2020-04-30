from dolfin import *
import math
import numpy as np
import sympy as sym
from sympy import symbols
from sympy import atan2,Abs
from stress import stress
from Reader import Reader
from GradF import GradF
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
T1 = TensorFunctionSpace(mesh,'DG',0,shape=(3,3))
V1 = VectorFunctionSpace(mesh, 'P', 1, dim=3)
SS = FunctionSpace(mesh,'P',1)
S = FunctionSpace(mesh,'R',0)
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


def extractor(u,cell_function,subdomain_id, V):
    dofmap=V.dofmap()
    mesh = V.mesh()
    for cell in cells(mesh):
        if cell_function[cell.index()]!=subdomain_id:
            dofs = dofmap.cell_dofs(cell.index())
            for dof in dofs:
                u.vector()[dof]=0.0
    return u


def Norm(u):
    a = sqrt(dot(u,u))
    return a
    


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
  
######################################################################
    

#looping prameters
counter=0
epsilon, sigma, t = 0.0, 16.0, 0.0
deps, dsig, dt= 0.001, 2.0, 0.001
epsmax, sigmax = 0.001, 16.0  
j=0
k=0
ID=0
Angles_i = File('Local Drive/Angles_i.pvd')
Angles_m = File('Local Drive/Angles_m.pvd')
Angles_a = File('Local Drive/Angles_a.pvd')
######################################################################


while counter<=M:

 
 #Reading the desirable displacement and growth parameter
 if counter%Skip==0:    
    u = Function(V,'Local Drive/u%d.xml' %(counter))
    a1, b1, c1 = Reader(counter)
    if counter!=0:
        alpha = a1
        beta = b1
        gamma =c1
    else:
        alpha = Constant(0)
        beta = Constant(0)
        gamma = Constant(0)
 ######################################################################

    #Growth tensor
    EXP = Expression('eps*exp(-a*(x[1]-1.68)*(x[1]-1.68))',degree=0,a=14.2857142857,eps=epsilon)
    PARS2 = as_tensor([[alpha,0,0],[0,gamma,0],[0,0,beta]])
    g_i = Identity(3)+EXP*PARS2
    G_i = det(g_i)
    ginv_i = inv(g_i)
    ######################################################################

    #Kinematics for initma
    #U1.interpolate(u)
    II = Identity(3)            # Identity tensor
    F = GradF(u,x)
    J = det(F)
    Fe_i= F*ginv_i
    C_i = F.T*F                   # Right Cauchy-Green tensor
    Ce_i= Fe_i.T*Fe_i
    Ce_i = variable(Ce_i)
    ######################################################################
 

    # Kinematics for media
    Fe_m= F*ginv_m
    C_m = F.T*F                   # Right Cauchy-Green tensor
    Ce_m= Fe_m.T*Fe_m
    Ce_m = variable(Ce_m)
    ######################################################################

    # Kinematics for adventitia
    Fe_a= F*ginv_a
    C_a = F.T*F                   # Right Cauchy-Green tensor
    Ce_a= Fe_a.T*Fe_a
    Ce_a = variable(Ce_a)
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
    
  
    #Calculating the energy    
    Angle_i = dot(F,b_i)
    PHI_i = atan_2(Angle_i[1],Angle_i[2])
    Phi_i = extractor(project(PHI_i,SS),Volume,7,SS)
    Phi_i.rename('Phi_i','Phi_i')
    Angles_i<<(Phi_i,t)
    
    Angle_m = dot(F,b_m)
    PHI_m = atan_2(Angle_m[1],Angle_m[2])
    Phi_m = extractor(project(PHI_m,SS),Volume,8,SS)
    Phi_m.rename('Phi_m','Phi_m')
    Angles_m<<(Phi_m,t)
    
    Angle_a = dot(F,b_a)
    PHI_a = atan_2(Angle_a[1],Angle_a[2])
    Phi_a = extractor(project(PHI_a,SS),Volume,9,SS)
    Phi_a.rename('Phi_a','Phi_a')
    Angles_a<<(Phi_a,t)
    ######################################################################
   

 
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

#updating
 k+=1
######################################################################



