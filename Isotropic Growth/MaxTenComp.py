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
    WW = FunctionSpace(Domain,'DG',1)
    TT = TensorFunctionSpace(Domain,'DG',1,shape=(3,3))
    u = Function(VV)
    u.interpolate(U)
    ###################################################################################


    #Plotting according to response
    sigma_11_total = project(Sig_11_i+Sig_11_m+Sig_11_a,WW)
    sigma_22_total = project(Sig_22_i+Sig_22_m+Sig_22_a,WW)
    sigma_33_total = project(Sig_33_i+Sig_33_m+Sig_33_a,WW)
    ###################################################################################


    #This function finds the maximum compressive stress
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

    #This function finds the maximum tensile stress
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


    #Sending stresses to max finders
    Max_11c_i = MaxComp(Sig_11_i,bulk,7)
    Max_11t_i = MaxTen(Sig_11_i,bulk,7)
    Max_22c_i = MaxComp(Sig_22_i,bulk,7)
    Max_22t_i = MaxTen(Sig_22_i,bulk,7)
    Max_33c_i = MaxComp(Sig_33_i,bulk,7)
    Max_33t_i = MaxTen(Sig_33_i,bulk,7)
    Max_11c_m = MaxComp(Sig_11_m,bulk,8)
    Max_11t_m = MaxTen(Sig_11_m,bulk,8)
    Max_22c_m = MaxComp(Sig_22_m,bulk,8)
    Max_22t_m = MaxTen(Sig_22_m,bulk,8)
    Max_33c_m = MaxComp(Sig_33_m,bulk,8)
    Max_33t_m = MaxTen(Sig_33_m,bulk,8)
    Max_11c_a = MaxComp(Sig_11_a,bulk,9)
    Max_11t_a = MaxTen(Sig_11_a,bulk,9)
    Max_22c_a = MaxComp(Sig_22_a,bulk,9)
    Max_22t_a = MaxTen(Sig_22_a,bulk,9)
    Max_33c_a = MaxComp(Sig_33_a,bulk,9)
    Max_33t_a = MaxTen(Sig_33_a,bulk,9)
    print "All Done!"
    ###################################################################################

    ###################################################################################
    return Max_11c_i, Max_11t_i, Max_22c_i, Max_22t_i, Max_33c_i, Max_33t_i,Max_11c_m, Max_11t_m, Max_22c_m, Max_22t_m, Max_33c_m, Max_33t_m,Max_11c_a, Max_11t_a, Max_22c_a, Max_22t_a, Max_33c_a, Max_33t_a
    ###################################################################################
