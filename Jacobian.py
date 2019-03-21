#----------------------------------------------------------------------------------------------------------------------#
#FEA Truss Analysis Software
#Created by Brad Young
#Mar 10, 2019
#----------------------------------------------------------------------------------------------------------------------#
import numpy as np
#----------------------------------------------------------------------------------------------------------------------#

def calcJacobian(augMatrix,disp_vector,forces):
    #define maximum acceptable error
    max_error = 0.000001
    #finished sim flag will be set to true when all displacements have less than max error
    finished_sim = False
    #Set the displacement vector to a float
    disp_vector = disp_vector.astype(float)
    #Create an initial displacement for calculating error, updated at end of each iteration
    last_iteration = np.zeros(disp_vector.shape[0])

    # Specify max number of iterations (prevent endless loops if solution doesn't converge)
    for x in range(0,10000):
        #If all errors are < max_error, stop the analysis
        if finished_sim == True:
            break
        # Create array to hold errors less than max_error, reset every iteration
        finished_check = []

        #For each variable (u1,v1,u2,v2,...) perform Jacobi method
        for i in range(0,len(disp_vector)):
            #If the initial value is set to 0, dont modify it (since it is a fixed point)
            if disp_vector[i] == 0:
                pass
            else:
            #Jacobi is split into two parts: See documentation for explanation
                pt1 = np.dot(augMatrix[i][0:i],disp_vector[0:i])
                pt2 = np.dot(augMatrix[i][i+1:augMatrix.shape[0]],disp_vector[i+1:augMatrix.shape[0]])
                disp_vector[i] = float(forces[i] - pt1 - pt2) / augMatrix[i][i]

                #Prevent dividing by zero (OLD CODE)
                # if augMatrix[i][i] == 0:
                #     disp_vector[i] = (forces[i] - pt1 - pt2) / 0.1
                # else:
                #     disp_vector[i] = float(forces[i]-pt1-pt2)/augMatrix[i][i]

        #Calculate error of each entry in displacement matrix
        for z in range(0,len(disp_vector)):
            current_error = abs(last_iteration[z] - disp_vector[z])
            #If the current_error < max_error add entry to finished_check list
            if current_error < max_error:
                finished_check.append(current_error)
            #Set the last_iteration vector to the current displacement vector for use in the next iteration
            last_iteration[z] = disp_vector[z]

        #If the finished_check list has the same number of entries as the columns of the matrix, then the solution has
        #converged and analysis can stop
        if len(finished_check) == disp_vector.shape[0]:
            finished_sim = True

    return(disp_vector)


