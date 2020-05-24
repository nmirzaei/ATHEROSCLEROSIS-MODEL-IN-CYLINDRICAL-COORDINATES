from dolfin import *
import numpy as np

#this function calculates the maximum tensile and compressive stress
def MaxTenComp(U,t,Sig_i,Sig_m,Sig_a):

    #reading the original configs
    Domain = Mesh('Mesh.xml')
    bulk = MeshFunction('size_t' , Domain , 'Mesh_physical_region.xml' )  #saves the interior info of the Domain
    boundary = MeshFunction('size_t', Domain , 'Mesh_facet_region.xml')
    ###################################################################################

    #allowing extrapolation in case of mesh refinement
    U.set_allow_extrapolation(True)
    ###################################################################################

    #Splitting the tensors into single entries
    [Sig_11_i,Sig_12_i,Sig_13_i,Sig_21_i,Sig_22_i,Sig_23_i,Sig_31_i,Sig_32_i,Sig_33_i] = Sig_i.split(True)
    [Sig_11_m,Sig_12_m,Sig_13_m,Sig_21_m,Sig_22_m,Sig_23_m,Sig_31_m,Sig_32_m,Sig_33_m] = Sig_m.split(True)
    [Sig_11_a,Sig_12_a,Sig_13_a,Sig_21_a,Sig_22_a,Sig_23_a,Sig_31_a,Sig_32_a,Sig_33_a] = Sig_a.split(True)
    ####################################################################################



    #defining function spaces
    VV = VectorFunctionSpace(Domain,'P',1)
    WW = FunctionSpace(Domain,'DG',0)
    TT = TensorFunctionSpace(Domain,'DG',0,shape=(3,3))
    u = Function(VV)
    u.interpolate(U)
    ###################################################################################


    #Plotting according to response
    sigma_11_total = project(Sig_11_i+Sig_11_m+Sig_11_a,WW)
    sigma_22_total = project(Sig_22_i+Sig_22_m+Sig_22_a,WW)
    sigma_33_total = project(Sig_33_i+Sig_33_m+Sig_33_a,WW)
    ###################################################################################


    #Calculates the maximum compressive stress
    def MaxComp(f, subdomains, subd_id):
        '''Minimum of f over subdomains cells marked with subd_id'''
        V = f.function_space()

        dm = V.dofmap()

        subd_dofs = np.unique(np.hstack(
                [dm.cell_dofs(c.index()) for c in SubsetIterator(subdomains, subd_id)]))

        Array= f.vector().get_local()[subd_dofs]
        Array[Array>0]=0
        return np.max(abs(Array))
    ###################################################################################

    #Calculates the maximum tensile stress
    def MaxTen(f, subdomains, subd_id):
        '''Minimum of f over subdomains cells marked with subd_id'''
        V = f.function_space()

        dm = V.dofmap()

        subd_dofs = np.unique(np.hstack(
                [dm.cell_dofs(c.index()) for c in SubsetIterator(subdomains, subd_id)]))

        Array= f.vector().get_local()[subd_dofs]
        Array[Array<0]=0
        return np.max(abs(Array))
    ###################################################################################

    #calculating the center of mass for the reference domain
    R = VectorFunctionSpace(Domain, "R", 0)
    V= VectorFunctionSpace(Domain, "P", 1)
    position = Function(V)
    position.assign(Expression(["x[0]", "x[1]"], element=V.ufl_element()))
    c = TestFunction(R)
    bulk0 = assemble(Constant(1.0)*dx(domain=Domain))
    centroid = assemble(dot(c, position)*dx)
    f = centroid / bulk0
    f_np = f.get_local()
    ###################################################################################

    #Moving the mesh
    ALE.move(Domain,u)
    ###################################################################################

    #finding the new center of mass
    new_V = VectorFunctionSpace(Domain, 'P', 1)
    new_R = VectorFunctionSpace(Domain, "R", 0)
    new_position = Function(new_V)
    new_position.assign(Expression(["x[0]", "x[1]"], element=new_V.ufl_element()))
    new_c = TestFunction(new_R)
    new_bulk = assemble(Constant(1.0)*dx(domain=Domain))
    new_centroid = assemble(dot(new_c, new_position)*dx)
    new_f = new_centroid / new_bulk
    new_f_np = new_f.get_local()
    dev = Constant(f_np - new_f_np)
    ###################################################################################


    #moving the deformed Domain back to the center of mass for the original Domain
    ALE.move(Domain , dev)
    ###################################################################################



    #sending stresses to have their maxima computed
    Max_11c_i = MaxComp(sigma_11_total,bulk,7)
    Max_11t_i = MaxTen(sigma_11_total,bulk,7)
    Max_22c_i = MaxComp(sigma_22_total,bulk,7)
    Max_22t_i = MaxTen(sigma_22_total,bulk,7)
    Max_33c_i = MaxComp(sigma_33_total,bulk,7)
    Max_33t_i = MaxTen(sigma_33_total,bulk,7)
    Max_11c_m = MaxComp(sigma_11_total,bulk,8)
    Max_11t_m = MaxTen(sigma_11_total,bulk,8)
    Max_22c_m = MaxComp(sigma_22_total,bulk,8)
    Max_22t_m = MaxTen(sigma_22_total,bulk,8)
    Max_33c_m = MaxComp(sigma_33_total,bulk,8)
    Max_33t_m = MaxTen(sigma_33_total,bulk,8)
    Max_11c_a = MaxComp(sigma_11_total,bulk,9)
    Max_11t_a = MaxTen(sigma_11_total,bulk,9)
    Max_22c_a = MaxComp(sigma_22_total,bulk,9)
    Max_22t_a = MaxTen(sigma_22_total,bulk,9)
    Max_33c_a = MaxComp(sigma_33_total,bulk,9)
    Max_33t_a = MaxTen(sigma_33_total,bulk,9)
    print "All Done!"
    ###################################################################################

    ###################################################################################
    return Max_11c_i, Max_11t_i, Max_22c_i, Max_22t_i, Max_33c_i, Max_33t_i,Max_11c_m, Max_11t_m, Max_22c_m, Max_22t_m, Max_33c_m, Max_33t_m,Max_11c_a, Max_11t_a, Max_22c_a, Max_22t_a, Max_33c_a, Max_33t_a
    ###################################################################################
