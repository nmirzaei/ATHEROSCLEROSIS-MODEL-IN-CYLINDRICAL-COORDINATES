from dolfin import *
import numpy 
import sympy as sym
from sympy import symbols
from sympy import atan2,Abs
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


#Changing the directory to the lustre directory
os.chdir("/lustre/scratch/nmirzaei")
######################################################################


# define function space     
V1 = VectorElement('P', triangle,1)         #Vector elements for the mixed space
V2 = FiniteElement('R', triangle,0)         #Real line elements for the mixed space
element = MixedElement([V1,V2,V2])
Z = FunctionSpace(mesh,element)             #mixed space
######################################################################

#Auxillary function spaces for projection and plotting purposes
S = FunctionSpace(mesh,'R',0)
V = VectorFunctionSpace(mesh,'P',1)
######################################################################

#Dimension
d = mesh.geometry().dim()
######################################################################


#defining spatial coordinates for changing to cylindrical coordinates 
x = SpatialCoordinate(mesh)
######################################################################
        

# Construct integration measure using these markers
ds = Measure('ds', subdomain_data=bnd_mesh)
dx = Measure('dx', subdomain_data=Volume)
######################################################################


# Define functions
du = TrialFunction(Z)            # Incremental displacement
v  = TestFunction(Z)             # Test function
Sol  = Function(Z)               # Solution function
######################################################################


#Breaking the solution
u , alpha, beta = split(Sol)
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


#Material radius
radius = Expression("x[0]",degree=0)
######################################################################


#looping prameters
counter=0
kk = 0.0
epsilon, sigma, t = 0.0, 16.0, 0.0
deps, dsig, dt= 0.001, 1.0, 0.001
epsmax, sigmax = 3.0, 16.0 
ID = 0
######################################################################


while int(epsilon*1000)<=int(epsmax*1000):
 


 #Growth in the intima with the unknown parameters
 if int(epsilon*1000)!=0:
     #this will enforce the Jacobian constraint
     gamma = (((1+epsilon)/((1+alpha*epsilon)*(1+beta*epsilon))-1)/epsilon)
 else:
     alpha = Constant(0)
     beta = Constant(0)
     gamma = Constant(0)
 ######################################################################

 #Setting up the growth in the intima    
 EXP = Expression('eps*exp(-a*(x[1]-1.68)*(x[1]-1.68))',degree=0,a=14.2857142857,eps=epsilon)     
 PARS2 = as_tensor([[alpha,0,0],[0,gamma,0],[0,0,beta]])
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
 print(epsilon,flush=True)
 print(t,flush=True)  
 print(sigma,flush=True)   
 ######################################################################
 
 
 #Pulling back the normals using Nanson's formula
 NansonOp = J*inv(F).T
 deformed_NN = dot(NansonOp,n)
 NormN = sqrt(dot(deformed_NN,deformed_NN))
 deformed_N = as_vector([deformed_NN[0],deformed_NN[1]])
 ######################################################################


 # Stored strain energy density (neo-Hookean Hopzafel model)
 psi_i = G_i*((mu_i/2)*(Ie_i - 3)+(eta_i/beta_i)*(exp(beta_i*(rho_i*(I4_i-1)**2+(1-rho_i)*(Ie_i-3)**2))-1)+((mu_i*nu_i)/(1-2*nu_i))*(J_i-1)**2- mu_i*ln(J_i))
 psi_m = ((mu_m/2)*(Ie_m - 3)+((mu_m*nu_m)/(1-2*nu_m))*(J_m-1)**2- mu_m*ln(J_m)+(eta_m/beta_m)*(exp(beta_m*(rho_m*(I4_m-1)**2+(1-rho_m)*(Ie_m-3)**2))-1))
 psi_a = ((mu_a/2)*(Ie_a - 3)+((mu_a*nu_a)/(1-2*nu_a))*(J_a-1)**2- mu_a*ln(J_a)+(eta_a/beta_a)*(exp(beta_a*(rho_a*(I4_a-1)**2+(1-rho_a)*(Ie_a-3)**2))-1))
 ######################################################################
 
 
 #The energy form
 try:
     Pi = 2*pi*(psi_i*radius*dx(7)+psi_m*radius*dx(8)+psi_a*radius*dx(9))#+Constant(Pressure_Energy))
 except:
     Pi = 2*pi*(psi_i*radius*dx(7)+psi_m*radius*dx(8)+psi_a*radius*dx(9))
 #####################################################################


 #Gateaux derivatives
 Fac = derivative(Pi, Sol, v)
 Jac = derivative(Fac, Sol, du)
 ######################################################################
 
 #Here you can use either solver. Regular (commented) or the snes.
# #Solve
# solve(Fac == 0, Sol, J=Jac,
#       solver_parameters={"newton_solver": {"relative_tolerance": 9e-10,   
#       "absolute_tolerance": 9e-10,"maximum_iterations": 80}})
########################################################################
# #Gateaux derivatives
# Fac = derivative(Pi, Sol, v)
# Jac = derivative(Fac, Sol, du)
# ######################################################################

# #In here we have no boundary conditions 
 bcs= []
 ######################################################################

 #Defining solver parameters
 snes_solver_parameters = {"nonlinear_solver": "snes",
                       "snes_solver": {"linear_solver": "lu",
                                       "absolute_tolerance":1e-8,
                                       "maximum_iterations": 50,
                                       "report": True,
                                       "error_on_nonconvergence": False}}
 ######################################################################

 #Defining the problem to solve
 problem = NonlinearVariationalProblem(Fac, Sol, bcs, Jac)
 ######################################################################
 
 #solve for Sol
 solver  = NonlinearVariationalSolver(problem)
 solver.parameters.update(snes_solver_parameters)
 info(solver.parameters, True)
 (iter, converged) = solver.solve()
 ######################################################################

 #Break the solution into components
 u, alpha, beta  = split(Sol)
 ######################################################################


 #making sure we start saving results after the pressure is imposed (if there is any pressure). Also making sure we save results every 10 steps    
 if sigma<=sigmax:
     sigma+=dsig
 if sigma>sigmax:
 ######################################################################

     #Saving the growth parameters
     File("Anisotropic/Growth/alpha(%d).xml" %(counter))<<project(alpha,S)
     File("Anisotropic/Growth/beta(%d).xml" %(counter))<<project(beta,S)
     File("Anisotropic/Growth/gamma(%d).xml" %(counter))<<project(gamma,S)
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


 



