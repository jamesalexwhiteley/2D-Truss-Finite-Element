
# A basic 2D Truss finite element (stiffness matrix mathod) Python script
# based on a MATLAB script by Angelo Simone, in "An Introduction to the analysis of slender structures", TU Delft , 2011

import numpy as np 

coord = np.array([[20,20], [20,0], [0,20], [0,0]]) # nodal co-ordinates 
# first node is top right, then counterclockwise, members as per order in conn 
conn = np.array([[2,0], [3,0], [0,1], [2,1], [3,1], [2,3]]) # nodal connectivity
dofNode = 2 # DOF's per node
index = np.zeros(4)

E = 10*10^6
A = np.zeros(6)
A[0] = 1
A[1] = np.sqrt(2) / 2 
A[2] = 1
A[3] = np.sqrt(2) / 2 
A[4] = 1
A[5] = 1

AE = A*E # element axial stiffness 

nNodes= np.size(coord , 0)
nElms= np.size(conn , 0)

# array Ku = f
nDofs = dofNode*nNodes  # total numbers of DOF's in system
k = np.zeros((nDofs , nDofs)) # stiffness matrix 
u = np.zeros((nDofs , 1)) # displacement array 
f = np.zeros((nDofs , 1)) # external force matrix 

constrainedDofs = np.array([4,5,6,7]) # constrained DOF's, i.e. boundary conditions
f[1] = 500; # external load at node n 
# assessemble stiffness matrix 
for ind in range(nElms):
    elmConn = conn[ind, :] # element connectivity
    x1 = coord[elmConn[0], 0]
    x2 = coord[elmConn[1], 0]
    y1 = coord[elmConn[0], 1]
    y2 = coord[elmConn[1], 1]
    len= np.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)) # element length 
    
    cos_theta = (x2-x1) / len # direction cosines golbal-local axes
    sin_theta = (y2-y1) / len 

    k_local = AE[ind] / len*np.array([[1,0,-1,0], [0,0,0,0], [-1,0,1,0], [0,0,0,0]]) # standard stiffnes matrix for Truss elements 

    R = np.array([[cos_theta, sin_theta],[-sin_theta, cos_theta]]) # rotation matrix 
    Zero = np.zeros((2,2))
    T = np.bmat([[R,Zero],[Zero,R]]) # translation matrix 

    ke = np.transpose(T)*k_local*T # global stiffness matrix 
    index[0] = dofNode*elmConn[0] # node 1 -> system DOF's 1 and 2 
    index[1] = dofNode*elmConn[0]+1 
    index[2] = dofNode*elmConn[1] # node 2 -> system DOF's 3 and 4
    index[3] = dofNode*elmConn[1]+1 
#   index local k into global k 
    ndof = index.shape[0]
    for i in range(ndof):
        ii = int(index[i])
        for j in range(ndof):
            jj = int(index[j])
            k[ii,jj] = k[ii,jj] + ke[i,j] 

# apply boundary condition and reduce global k

k[constrainedDofs, :] = 0
k[:, constrainedDofs] = 0


# print(k[4,5,6,7])
# k[constrainedDofs, constrainedDofs] = np.identity(constrainedDofs.shape[0])
k[4:8, 4:8] = np.identity(constrainedDofs.shape[0])


a = np.linalg.inv(k)*f


# plot results 
# mag = 400
# hold on
# for e=1:nElements
#     x =coord(conn(e , : ) , 1) 
#     y =coord(conn(e , : ) , 2)
#     u =a(2*conn(e , : ) -1) 
#     v =a(2*conn(e , : ) ) 

#     title( ' Deformed plot ' )
#     axis equal
#     plot (x , y , ' r--o ' )
#     plot (x+mag*u , y+mag*v , ' k-o ')
