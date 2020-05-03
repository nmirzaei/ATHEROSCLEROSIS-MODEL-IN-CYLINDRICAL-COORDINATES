from dolfin import *

#This function calculates the gradient of a vector in cylindrical coordinates
def gradv(v,x):
    return sym(as_tensor([[v[0].dx(0),  v[0].dx(1), 0 ],
                          [v[1].dx(0), v[1].dx(1), 0],
                          [0, 0, v[0]/x[0]]]))

