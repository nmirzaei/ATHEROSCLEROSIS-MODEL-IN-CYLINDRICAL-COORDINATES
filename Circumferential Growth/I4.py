#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 10:37:44 2019

@author: root
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 17:00:39 2019

@author: root
"""

from dolfin import *
import numpy as np
import sympy as sym
from sympy import symbols
from sympy import atan2,Abs
from stress import stress
from GradF import GradF
import os

#The following lines take an already converted Gmsh mesh and import it
mesh = Mesh('2D_cylinder.xml')  
######################################################################

# Facet functions
Volume = MeshFunction('size_t' , mesh , '2D_cylinder_physical_region.xml' )  #saves the interior info of the mesh
bnd_mesh = MeshFunction('size_t', mesh , '2D_cylinder_facet_region.xml')  #saves the boundary info of the mesh
######################################################################


# Optimization options for the form compiler
parameters['form_compiler']['cpp_optimize'] = True
ffc_options = {'optimize': True}
######################################################################    

# define function space
V = VectorFunctionSpace(mesh, 'P', 2)
T1 = TensorFunctionSpace(mesh,'DG',0,shape=(3,3))
SS= FunctionSpace(mesh,'P',2)
N = V.dim()
d = mesh.geometry().dim()
######################################################################



# Construct integration measure using these markers
ds = Measure('ds', subdomain_data=bnd_mesh)
dx = Measure('dx', subdomain_data=Volume)
######################################################################


# Define functions
U1 = Function(V)
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
mu_i, nu_i, beta_i, eta_i, rho_i, phi_i = 27.9, 0.49, 170.88, 0, 0.51, 60.3*pi/180

# Elasticity parameter for the Hopzafel cylinder (Media)
mu_m, nu_m, beta_m, eta_m, rho_m, phi_m = 1.27, 0.49, 8.21, 0, 0.25, 20.61*pi/180

# Elasticity parameter for the Hopzafel cylinder (adventitia)
mu_a, nu_a, beta_a, eta_a, rho_a, phi_a = 7.56, 0.49, 85.03, 0, 0.55, 67*pi/180
######################################################################


#Collagen fibers direction in cylindrical coordinates  (before it was (0,cos,sin) I changed it to (0,sin,cos) which I believe is the right one )
b_i = Expression(('0','sin(phi)','cos(phi)'),degree=0,phi=phi_i)
b_m = Expression(('0','sin(phi)','cos(phi)'),degree=0,phi=phi_m)
b_a = Expression(('0','sin(phi)','cos(phi)'),degree=0,phi=phi_a)
######################################################################


#growth function for intima (because x[0]->x->r, x[1]->y->z, x[2]->z->theta we want to not have growth in z hence entery 22 is 1)
ginv_i = Expression((('1','0','0'),('0','1','0'),('0','0','1/(1+epsilon*exp(-a*(x[1]-1.68)*(x[1]-1.68)))')),degree=1,epsilon=0.0,a=14.2857142857)
G_i = Expression('(1+epsilon*exp(-a*(x[1]-1.68)*(x[1]-1.68)))', degree=1, epsilon=0.0,a=14.2857142857)
G1 = Expression('1/(1+epsilon*exp(-a*(x[1]-1.68)*(x[1]-1.68)))', degree=1, epsilon=0.0,a=14.2857142857)
#G_i_inv = project(G1,SS)
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

#radius variable
radius = Expression("x[0]",degree=0)
######################################################################

#looping prameters
counter=0
epsilon, sigma, t = 0.0, 12.0, 0.0
deps, dsig, dt= 0.001, 2.0, 0.001
epsmax, sigmax = 0.001, 0.0  
j=0
k=0
ID=0
vtkfile1 = File('/media/user1/Seagate Backup Plus Drive/Dropbox/Research/Fenics codes Paper 1/Cylindrical_Strip_Paper/Anisotropic_Growth/Theta_Growth_Scaled_I4/Theta_Growth/I4_i.pvd')
vtkfile2 = File('/media/user1/Seagate Backup Plus Drive/Dropbox/Research/Fenics codes Paper 1/Cylindrical_Strip_Paper/Anisotropic_Growth/Theta_Growth_Scaled_I4/Theta_Growth/I4_m.pvd')
vtkfile3 = File('/media/user1/Seagate Backup Plus Drive/Dropbox/Research/Fenics codes Paper 1/Cylindrical_Strip_Paper/Anisotropic_Growth/Theta_Growth_Scaled_I4/Theta_Growth/I4_a.pvd')
######################################################################


while counter<=M:

 
 #Reading the desirable displacement
 if counter%Skip==0:   
    u = Function(V,'/media/user1/Seagate Backup Plus Drive/Dropbox/Research/Fenics codes Paper 1/Cylindrical_Strip_Paper/Anisotropic_Growth/Theta_Growth_Scaled_I4/Theta_Growth/displacement/u%d.xml' %(counter))

 ######################################################################   

# Kinematics for initma
    U1.interpolate(u)
    II = Identity(3)            # Identity tensor
    F = GradF(U1,x)
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
    i44_i = project(I4_i,SS)
    i4_i = extract_values(i44_i,Volume,7,SS)
    i4_i.rename('i4_i','i4_i')
    vtkfile1<<(i4_i,t)
 #for hopzafel (media)
    I4_m = conditional(lt(dot(b_m, Ce_m*b_m),Constant(1)),1,dot(b_m, Ce_m*b_m))
    i44_m = project(I4_m,SS)
    i4_m = extract_values(i44_m,Volume,8,SS)
    i4_m.rename('i4_m','i4_m')
    vtkfile2<<(i4_m,t)
 #for hopzafel (adventitia)
    I4_a = conditional(lt(dot(b_a, Ce_a*b_a),Constant(1)),1,dot(b_a, Ce_a*b_a))
    i44_a = project(I4_a,SS)
    i4_a = extract_values(i44_a,Volume,9,SS)
    i4_a.rename('i4_a','i4_a')
    vtkfile3<<(i4_a,t)
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
 print(epsilon)#,flush=True)
 print(t)#,flush=True)  
 print(sigma)#,flush=True)   
######################################################################

#updating
 ginv_i.epsilon = epsilon
 G_i.epsilon = epsilon
 k+=1




