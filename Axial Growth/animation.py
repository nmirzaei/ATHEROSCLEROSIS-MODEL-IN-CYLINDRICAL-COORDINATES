from dolfin import *
import numpy as np
import os
#This function creates our animation
def animation(uu,t,vtkfile11,vtkfile22,vtkfile33):

    #changing back to the original directory to read the files
    os.chdir("home directory")
    ######################################################################

    #create the mesh and its volume and boundary
    Domain=Mesh('Mesh.xml')
    Bulk = MeshFunction('size_t' , Domain , 'Mesh_physical_region.xml' )  #saves the interior info of the mesh
    Boundary = MeshFunction('size_t', Domain , 'Mesh_facet_region.xml')
    ##################################################################################################################

    #changing to the scratch directory to save the files
    os.chdir("Lustre directory")
    ######################################################################

    #Create required function spaces
    R = VectorFunctionSpace(Domain, "R", 0)
    V= VectorFunctionSpace(Domain, "P", 1)
    u = Function(V)
    U = as_vector([uu[0],uu[1]])
    uuu = project(U,V)
    u.interpolate(uuu)
    ##################################################################################################################

    #Finding the center of mass for the reference domain
    position = Function(V)
    position.assign(Expression(["x[0]", "x[1]"], element=V.ufl_element()))
    c = TestFunction(R)
    volume = assemble(Constant(1.0)*dx(domain=Domain))
    centroid = assemble(dot(c, position)*dx)
    f = centroid / volume
    f_np = f.get_local()
    ##################################################################################################################

    #Move the mesh
    ALE.move(Domain,u)
    ##################################################################################################################

    #Finding the center of mass for the deformed mesh and move it to the old cener of mass
    new_V = VectorFunctionSpace(Domain, 'P', 1)
    new_R = VectorFunctionSpace(Domain, "R", 0)
    new_position = Function(new_V)
    new_position.assign(Expression(["x[0]", "x[1]"], element=new_V.ufl_element()))
    new_c = TestFunction(new_R)
    new_volume = assemble(Constant(1.0)*dx(domain=Domain))
    new_centroid = assemble(dot(new_c, new_position)*dx)
    new_f = new_centroid / new_volume
    new_f_np = new_f.get_local() # numpy array
    #deviation of the deformed center of mass from the original one
    dev = Constant(f_np - new_f_np)
    #moving the deformed mesh back to (0.5,0.5,0.5) which is the center of mass for the original mesh
    ALE.move(Domain , dev)
    ##################################################################################################################



    #Saving the frames
    vtkfile11<<(Bulk,t)
    vtkfile22<<(Boundary,t)
    vtkfile33<<(Domain,t)
    ##################################################################################################################
