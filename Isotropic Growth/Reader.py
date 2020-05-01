from dolfin import *
import numpy as np
import os

#This function creates our animation
def Reader(counter):
#    os.chdir("/home/1925/Fenics/Isotropic_Growth_1")
    #create the mesh and its volume and boundary
    Domain=Mesh('2D_cylinder.xml')
    Bulk = MeshFunction('size_t' , Domain , '2D_cylinder_physical_region.xml' )  #saves the interior info of the mesh
    Boundary = MeshFunction('size_t', Domain , '2D_cylinder_facet_region.xml') 
    ##################################################################################################################
    
    S = FunctionSpace(Domain,'R',0)
    
    a1 = Function(S,'/media/user1/Seagate Backup Plus Drive/Dropbox/Research/Fenics codes Paper 1/Cylindrical_Strip_Paper/Anisotropic_Growth/Anisotropic Growth Parameters/Growth/alpha(%d).xml' %(counter))
    b1= Function(S,'/media/user1/Seagate Backup Plus Drive/Dropbox/Research/Fenics codes Paper 1/Cylindrical_Strip_Paper/Anisotropic_Growth/Anisotropic Growth Parameters/Growth/beta(%d).xml' %(counter))
    c1 = Function(S,'/media/user1/Seagate Backup Plus Drive/Dropbox/Research/Fenics codes Paper 1/Cylindrical_Strip_Paper/Anisotropic_Growth/Anisotropic Growth Parameters/Growth/gamma(%d).xml' %(counter))
#    os.chdir("/lustre/scratch/nmirzaei")
    
    return a1.vector().get_local()[0], b1.vector().get_local()[0], c1.vector().get_local()[0]
