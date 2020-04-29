from dolfin import *
import numpy 
import sympy as sym
from sympy import symbols
from sympy import atan2,Abs
from animation import animation
from GradF import GradF
from gradv import gradv
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


#Changing the directory to the lustre directory
os.chdir("Location on Lustre storage")
######################################################################


#Auxillary function spaces for projection and plotting purposes
S = FunctionSpace(mesh,'R',0)
V = VectorFunctionSpace(mesh,'P',2)
######################################################################

#Dimension
d = mesh.geometry().dim()
######################################################################


#Save the reference config
file=File('validation/Volume.pvd')
file<<Volume
file=File('validation/bnd.pvd')
file<<bnd_mesh
file=File('validation/mesh.pvd')
file<<mesh
######################################################################

#defining spatial coordinates for changing to cylindrical coordinates 
x = SpatialCoordinate(mesh)
######################################################################
        

# Construct integration measure using these markers
ds = Measure('ds', subdomain_data=bnd_mesh)
dx = Measure('dx', subdomain_data=Volume)
######################################################################


# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)               # Solution function
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


#Collagen fibers direction in cylindrical coordinates  
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


#Files for saving the domain evolution and possibly growth parameter evolution
vtkfile11 = File("Anisotropic/Vol_rec.pvd")
vtkfile22 = File("Anisotropic/bnd_rec.pvd")
vtkfile33 = File("Anisotropic/domain_rec.pvd")
vtkfile55 = File("Anisotropic/deformmed_N.pvd")
#####################################################################

#Material radius
radius = Expression("x[0]",degree=0)
######################################################################


#looping prameters
counter=0
kk = 0.0
epsilon, sigma, t = 0.0, 8.0, 0.0
deps, dsig, dt= 0.001, 8.0, 0.001
epsmax, sigmax = 3.0, 16.0 
ID = 0
######################################################################


while int(epsilon*1000)<=int(epsmax*1000):
 
 #Reading the parameters saved by the Growth-Parameters codes
 if counter!=0:
     a1, b1, c1 = Reader(counter)
 else:
     a1 = Constant(0)
     b1 = Constant(0)
     c1 = Constant(0)
 ######################################################################

 #Growth in the intima with the unknown parameters
 EXP = Expression('eps*exp(-a*(x[1]-1.68)*(x[1]-1.68))',degree=0,a=14.2857142857,eps=epsilon)
 PARS2 = as_tensor([[a1,0,0],[0,c1,0],[0,0,b1]])
 g_i = Identity(3)+EXP*PARS2
 G_i = det(g_i)
 ginv_i = inv(g_i)
 ######################################################################

    
 # Kinematics for initma
 II = Identity(3)            # Identity tensor
 F = GradF(u,x)              # We use the custom Gradient we have defined
 J = det(F)
 Fe_i= F*ginv_i
 Fe_i = variable(Fe_i)
 C_i = F.T*F                 # Right Cauchy-Green tensor
 Ce_i= Fe_i.T*Fe_i
 ######################################################################
 

 # Kinematics for media
 Fe_m= F*ginv_m
 Fe_m = variable(Fe_m)
 C_m = F.T*F                   # Right Cauchy-Green tensor
 Ce_m= Fe_m.T*Fe_m
######################################################################
 

# Kinematics for adventitia
 Fe_a= F*ginv_a
 Fe_a = variable(Fe_a)
 C_a = F.T*F                   # Right Cauchy-Green tensor
 Ce_a= Fe_a.T*Fe_a
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
 
 #printing the current loop variables
 print(epsilon)#,flush=True)
 print(t)#,flush=True)  
 print(sigma)#,flush=True)   
 ######################################################################
 
 
 #Pulling back the normals using Nanson's formula
 NansonOp = J*inv(F).T
 deformed_NN = dot(NansonOp,n)
 NormN = sqrt(dot(deformed_NN,deformed_NN))
 deformed_N = as_vector([deformed_NN[0],deformed_NN[1]])
 ######################################################################


 # Stored strain energy density (neo-Hookean Hopzafel model)
 psi_i = ((mu_i/2)*(Ie_i - 3)+(eta_i/beta_i)*(exp(beta_i*(rho_i*(I4_i-1)**2+(1-rho_i)*(Ie_i-3)**2))-1)+((mu_i*nu_i)/(1-2*nu_i))*(J_i-1)**2- mu_i*ln(J_i))
 psi_m = ((mu_m/2)*(Ie_m - 3)+((mu_m*nu_m)/(1-2*nu_m))*(J_m-1)**2- mu_m*ln(J_m)+(eta_m/beta_m)*(exp(beta_m*(rho_m*(I4_m-1)**2+(1-rho_m)*(Ie_m-3)**2))-1))
 psi_a = ((mu_a/2)*(Ie_a - 3)+((mu_a*nu_a)/(1-2*nu_a))*(J_a-1)**2- mu_a*ln(J_a)+(eta_a/beta_a)*(exp(beta_a*(rho_a*(I4_a-1)**2+(1-rho_a)*(Ie_a-3)**2))-1))
 ######################################################################
 
#Stress variables for the weak form
 TT_i = G_i*diff(psi_i,Fe_i)*(ginv_i.T)
 TT_m = diff(psi_m,Fe_m)*(ginv_m.T)
 TT_a = diff(psi_a,Fe_a)*(ginv_a.T)
######################################################################
 
#Creating the required test functions according to our change of coordinates
 GradV = gradv(v,x)
 testV = as_vector([v[0],0,v[1]])
######################################################################
 
 
#The weak form
 try:
     Pi = (2*pi)*(inner(TT_i,GradV)*radius*dx(7)+inner(TT_m,GradV)*radius*dx(8)+inner(TT_a,GradV)*radius*dx(9)+(sigma)*radius*dot(deformed_N,testV)*ds(1))#-(2*sigma/3)*dot(Non_Loc_G_Correction,testV)*radius*dx(7)
 except:
     Pi = (2*pi)*(inner(TT_i,GradV)*radius*dx(7)+inner(TT_m,GradV)*radius*dx(8)+inner(TT_a,GradV)*radius*dx(9))
######################################################################

#Solve
 solve(Pi == 0, u ,
       solver_parameters={"newton_solver": {"relative_tolerance": 9e-10,   
       "absolute_tolerance": 9e-10,"maximum_iterations": 80}})
######################################################################



 #making sure we start saving results after the pressure is imposed (if there is any pressure). Also making sure we save results every 10 steps    
 if sigma<=sigmax:
     sigma+=dsig
 if sigma>sigmax:
     #Save animation every 10 steps to save space and time
     if counter%10==0:
         animation(u,t,vtkfile11,vtkfile22,vtkfile33)
     ######################################################################

     #Save the displacements for later stress computations
     File("Anisotropic/displacement/u%d.xml" %(counter))<<project(u,V) #Save the displacements for later process
     ######################################################################

     #Update loop info     
     epsilon+=deps
     counter+=1
     t+=dt
     ID+=1
     ######################################################################

 ######################################################################   
 if ID!=0:
     sigma-=dsig
     ID=0
 ######################################################################
     
 ######################################################################
 kk+=1
 ######################################################################


 



