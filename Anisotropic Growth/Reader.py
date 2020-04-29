from dolfin import *
import numpy as np
import os

#This function reads the gowth parametes saved by our growth parameter finder code 
#The os.chdir and the 
def Reader(counter):

    os.chdir("Location on HPC storage")

    #create the mesh and its volume and boundary
    Domain=Mesh('2D_cylinder.xml')
    Bulk = MeshFunction('size_t' , Domain , '2D_cylinder_physical_region.xml' )  #saves the interior info of the mesh
    Boundary = MeshFunction('size_t', Domain , '2D_cylinder_facet_region.xml') 
    ##################################################################################################################
    
    S = FunctionSpace(Domain,'R',0)
    
    a1 = Function(S,'Location on local Drive/alpha(%d).xml' %(counter))
    b1= Function(S,'Location on local Drive/beta(%d).xml' %(counter))
    c1 = Function(S,'Location on local Drive/gamma(%d).xml' %(counter))
    os.chdir("Location on Lustre storage")
    
    return a1.vector().get_local()[0], b1.vector().get_local()[0], c1.vector().get_local()[0]
