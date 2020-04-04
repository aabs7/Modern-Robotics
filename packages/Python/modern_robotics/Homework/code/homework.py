# coding: utf-8
import numpy as np 
import core as c

##Given in question#########################################

W1 = 0.109
W2 = 0.082
L1 = 0.425
L2 = 0.392
H1 = 0.089
H2 = 0.095    #in meters

Tsd = np.array([[0,1,0,-0.5],
                [0,0,-1,0.1],
                [-1,0,0,0.1],
                [0,0,0,1]])

Blist = np.array([[0, 1,  0, W1 + W2, 0,   L1 + L2],
                  [0, 0,  1, H2, -L1 - L2,   0],
                  [0, 0,  1, H2, -L2, 0],
                  [0, 0,  1,  H2, 0, 0],
                  [0,-1,  0, -W2, 0, 0],
                  [0, 0,  1, 0, 0, 0]]).T

M = np.array([[-1, 0,  0, L1 + L2], 
              [ 0, 0,  1, W1 + W2], 
              [ 0, 1, 0,  H1 - H2], 
              [ 0, 0,  0, 1]])

thetalist0 = np.array([3,0.5,-2,2.5,-1.2,-2.3])
eomg = 0.001
ev = 0.0001

#########################################################


###  Function defined below
def IKinBodyIterates(Blist,M,T,thetalist0,eomg,ev):
   
    thetalist = np.array(thetalist0).copy()
    i = 0
    thetamatrix = thetalist
    Tsb = c.FKinBody(M,Blist,thetalist)
    V_b = c.se3ToVec(c.MatrixLog6(np.dot(c.TransInv(Tsb),Tsd)))

    ## Print here ##
    print('Iteration {}:'.format(i))
    print("joint vector : {}".format(thetalist))
    print("SE(3) end−effector config:{}".format(Tsb.reshape(1,np.size(Tsb))))    
    print("error twist V_b : {}".format(V_b))
    print("angular error magnitude ||omega_b||:{}".format(np.linalg.norm(V_b[0:2])))
    print("linear error magnitude ||v_b||:{}".format(np.linalg.norm(V_b[3:5])))
    print("\r\n")

    while True:
        i += 1
        #calculate end effector position after thetalist
        Tsb = c.FKinBody(M,Blist,thetalist)
        # calculate error twist
        ## Tsb* x Tsd = Tbd i.e transformation from b-current position to d-desired position
        # log(Tbd) gives velocity in terms of twist to reach that point 
        ## [Vb] = log(Tsb* x Tsd) and this is se3 matrix
        ##  we need vector 
        ## Vb = se3toVec([Vb])
        V_b = c.se3ToVec(c.MatrixLog6(np.dot(c.TransInv(Tsb),Tsd)))
        e_omg = np.linalg.norm(V_b[0:2])
        e_v = np.linalg.norm(V_b[3:5])
        
        if eomg > e_omg and ev > e_v:
            break

        Jb = c.JacobianBody(Blist,thetalist)
        ## theta(i+1) = theta(i) + Jb* x V_b
        thetalist = thetalist + np.dot(np.linalg.pinv(Jb),V_b)
        
        #save in matrix
        thetamatrix = np.vstack((thetamatrix,thetalist))

        ##Print here
        print('Iteration {}:'.format(i))
        print("joint vector : {}".format(thetalist))
        print("SE(3) end−effector config:{}".format(Tsb.reshape(1,np.size(Tsb))))    
        print("error twist V_b : {}".format(V_b))
        print("angular error magnitude ||omega_b||:{}".format(np.linalg.norm(V_b[0:2])))
        print("linear error magnitude ||v_b||:{}".format(np.linalg.norm(V_b[3:5])))
        print("\r\n")
    
    np.savetxt("iterates.csv", thetamatrix,fmt='%.3f', delimiter=",")
    return (thetalist, thetamatrix)

 
##thetalist, not err, thetamatrix = 
list,matrix = IKinBodyIterates(Blist,M,Tsd,thetalist0,eomg,ev)




