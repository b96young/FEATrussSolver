#----------------------------------------------------------------------------------------------------------------------#
#FEA Truss Analysis Software
#Created by Brad Young
#Mar 10, 2019
#----------------------------------------------------------------------------------------------------------------------#
import numpy as np
import math
import matplotlib.pyplot as plt
#----------------------------------------------------------------------------------------------------------------------#
def transformForce(force_list):
    transformed_forces = []

    for node_force in force_list:
        angle = node_force[1]*3.14159/180
        force = []
        force_x = node_force[0]*math.cos(angle)
        force_y = node_force[0]*math.sin(angle)
        force.append(force_x)
        force.append(force_y)
        transformed_forces.append(force)

    return transformed_forces

def diaToArea(diameter):
    area = pow(diameter,2)*3.14159/4
    return area

def augmentMatrix(elements,num_nodes):
    # Function takes the list of element objects and the number of nodes in system

    # Need to double shape of matrix since going from u -> u,v
    augmented_m = np.zeros((num_nodes * 2, num_nodes * 2))

    # Function separates the four quadrants of the individual stiffness matrix, moves them into augmented matrix
    # Detailed explanation in documentation
    # Quadrant labeling system:
    # [1,1,3,3]
    # [1,1,3,3]
    # [2,2,4,4]
    # [2,2,4,4]
    # Python slicing notation: (top:bottom, left:right)

    # Reposition first quadrant
    first_quad = elements.stiff_matrix[0:2, 0:2]
    # need to multiply the node index by 2, since transforming from u -> u,v
    first_quad_top_index = elements.connecting_nodes[0]*2
    # bottom index is 2 below the top index since moving a 2x2 matrix
    first_quad_bottom_index = first_quad_top_index + 2
    # need to multiply the node index by 2, since transforming from u -> u,v
    first_quad_left_index = elements.connecting_nodes[0]*2
    # right index is 2 to the right of the left index since moving a 2x2 matrix
    first_quad_right_index = first_quad_top_index + 2
    # insert the 2x2 matrix into the augmented matrix
    augmented_m[first_quad_top_index:first_quad_bottom_index,first_quad_left_index:first_quad_right_index] = first_quad

    # Reposition second quadrant
    sec_quad = elements.stiff_matrix[2:4, 0:2]
    sec_quad_top_index = elements.connecting_nodes[1]*2
    sec_quad_bottom_index = sec_quad_top_index + 2
    sec_quad_left_index = elements.connecting_nodes[0]*2
    sec_quad_right_index = sec_quad_left_index + 2
    augmented_m[sec_quad_top_index:sec_quad_bottom_index, sec_quad_left_index:sec_quad_right_index] = sec_quad

   # Reposition third quadrant
    third_quad = elements.stiff_matrix[0:2, 2:4]
    third_quad_top_index = elements.connecting_nodes[0] * 2
    third_quad_bottom_index = third_quad_top_index + 2
    third_quad_left_index = elements.connecting_nodes[1] * 2
    third_quad_right_index = third_quad_left_index + 2
    augmented_m[third_quad_top_index:third_quad_bottom_index, third_quad_left_index:third_quad_right_index] = third_quad

    # Reposition fourth quadrant
    four_quad = elements.stiff_matrix[2:4, 2:4]
    four_quad_top_index = elements.connecting_nodes[1] * 2
    four_quad_bottom_index = four_quad_top_index + 2
    four_quad_left_index = elements.connecting_nodes[1] * 2
    four_quad_right_index = four_quad_left_index + 2
    augmented_m[four_quad_top_index:four_quad_bottom_index, four_quad_left_index:four_quad_right_index] = four_quad

    return augmented_m


def combineAugmented(elementList):
    # Performs matrix addition of each element's augmented matrix to form global stiffness matrix

    # Convert list of objects into list of matrices
    mList = []
    for element in elementList:
        mList.append(element.augMatrix)

    # Initialize a zeros matrix of same shape as an augmented matrix of the system
    matrix_shape = mList[0].shape
    global_k = np.zeros(matrix_shape)

    # Add matrices together and return global stiffness matrix
    for each in range(0,len(mList)):
        global_k = np.add(mList[each], global_k)

    return global_k

def flattenList(list):
    # Required to merge 2D lists created in JSON file into a 1D list
    flattened = []
    for node in list:
        for point in node:
            flattened.append(point)
    return flattened

def trimMatrix(initial_disp,matrix,vector=False):
    # Condenses input matrix by removing rows and columns corresponding to 0 entries of the initial displacement vector

    #Flatten the initial displacement vector {flattenList definined in this file}
    initial_disp = flattenList(initial_disp)

    # Since function is used for vectors and matrices, only flatten vectors
    if vector == True:
        matrix = flattenList(matrix)

    # Create a list to hold indexes of rows and columns to delete
    to_delete = []
    # Iterate through the initial displacement vector, append indexes to delete into to_delete list
    for i in range(0,len(initial_disp)):
        if initial_disp[i] == 0:
            to_delete.append(i)

    # Delete columns of matrix / vector corresponding to entries of to_delete list
    matrix = np.delete(matrix,to_delete,0)

    # If the input matrix is not a vector, delete the rows as well
    if vector == False:
        matrix = np.delete(matrix, to_delete, 1)

    return matrix

def recompileMatrix(initial_disp,result_vector):
    #Re-introduces previously deleted 0 entries of initial_disp vector into the final_disp vector

    #Flatten the initial displacement vector {flattenList defined in this file}
    initial_disp = flattenList(initial_disp)

    #Create a list to hold the recompiled result vector
    output = []
    #Create a counter to track which entries from the result_vector have been appended to the output list
    counter = 0
    for i in range(0,len(initial_disp)):
        if initial_disp[i] == 0:
            output.append(0)
        else:
            output.append(result_vector[counter])
            counter = counter +1

    return output


def calcStress(elements,disp_vector):
    #Input a list of elements and the final_disp vector, return stress in each element

    #Create a list to store stresses in each element
    stress = []
    for element in elements:
        disp = []
        #Get order of nodes in element
        start_node = element.connecting_nodes[0]*2
        end_node = element.connecting_nodes[1]*2

        #Slice the corresponding displacements of each node from the displacement vector
        #Extend method appends each entry in the list (as opposed to appending the list)
        disp.extend(disp_vector[start_node:start_node+2])
        disp.extend(disp_vector[end_node:end_node+2])

        #Perform the dot product of the stress transformation matrix and the corresponding displacements,
        #append to stress list after each iteration
        stress.append(np.dot(element.stress_matrix,disp))

    return stress


def reactionForces(global_stiffness,disp_vector):
    #Find the reaction forces at each boundary given the global stiffness matrix and the displacement vector
    #R_force = [K]*[U]

    #Find the indexes of the displacement vector that are non-zero (ie not a boundary)
    #Use to create list of indexes of columns to delete from global stiffness matrix
    to_delete = []
    for i in range(0,len(disp_vector)):
        if disp_vector[i] != 0:
            to_delete.append(i)

    #Delete the columns from the global stiffness matrix corresponding to the to_delete list
    global_stiffness = np.delete(global_stiffness,to_delete,0)

    #Calculate the reaction forces using the dot product
    r_forces = np.dot(global_stiffness,disp_vector)

    #Recompile the r_forces by adding 0s where the columns were deleted
    #Output list holds the final reaction force vector
    output = []
    #Counter tracks which entries from r_forces vector have been appended to output list
    counter = 0
    for i in range(0,len(disp_vector)):
        if i in to_delete:
            output.append(0)
        else:
            output.append(r_forces[counter])
            counter = counter +1

    return output


def printTable(name,vector):
    #Print out final results from analysis

    if name == "disp":
        print("Final Displacements")
        for i in range(0,len(vector)):
            if i % 2 == 0:
                print("u%d" % (i/2), " : ", vector[i])
            else:
                print("v%d" % (i/2), " : ", vector[i])
        print("\n")

    if name == "stress":
        print("Stresses")
        for i in range(0,len(vector)):
            print("Element %d" %i, " : ", vector[i])
        print("\n")

    if name == "reaction":
        print("Reaction Forces")
        for i in range(0,len(vector)):
            if i % 2 == 0:
                print("F%dx" % (i/2), " : ", vector[i])
            else:
                print("F%dy" % math.floor(i/2), " : ", vector[i])
        print("\n")


def plotTruss(node_list,element_list,i_disp,i_force):
    #Plot a schematic view of the analyzed structure using matplotlib

    #Plot the elements used in simulation
    for i in range(0,len(element_list)):
        #Create lists to hold the start and finish point of each element, reset after each iteration for new element
        x = []
        y = []
        #Get coords of starting node and ending node
        node1 = element_list[i].connecting_nodes[0]
        node2 = element_list[i].connecting_nodes[1]

        #Append the coordinates of the start point
        x.append(node_list[node1].node_coords[0])
        y.append(node_list[node1].node_coords[1])
        #Append coordinates of end point
        x.append(node_list[node2].node_coords[0])
        y.append(node_list[node2].node_coords[1])
        #Plot the line corresponding to the element
        plt.plot(x,y)
        #Annotate at the middle of the element the element number
        midpoint = element_list[i].midpoint
        plt.annotate("E%d" % i,(midpoint[0],midpoint[1]))

    #Plot the nodes and forces used in the simulation
    for i in range(0,len(node_list)):
        #Get the x,y coordinate of the node
        x = node_list[i].node_coords[0]
        y = node_list[i].node_coords[1]
        #Annotate at the node location its number
        plt.annotate("N%d" % node_list[i].node_num,(x,y))

        #Get the x,y components of the force at the i'th node
        force = i_force[i]
        #Create lists to hold the start and end points of the force vector
        force_x = []
        force_y = []
        #Append the current node coordinates as the start point of the force vector
        force_x.append(x)
        force_y.append(y)

        #Scale the force using a log function (to maintain some relativity of magnitude)
        #If force is between -1 < F < 1 then dont use log function
        #Append the end x coordinate
        if abs(force[0] + force[1]) > 0 :
            if force[0] > 1:
                force_x.append(math.log(force[0]) + x)
            elif force[1] < 1:
                force_x.append(-1*math.log(-1*force[0]) + x)
            else:
                force_x.append(force[0] + x)
        else:
            force_x.append(x)

        # Append the end y coordinate
        if abs(force[1]) > 0:
            if force[1] > 1:
                force_y.append(math.log(force[1]) + y)
            elif force[1] < -1:
                force_y.append(-1*math.log(-1*force[1]) + y)
            else:
                force_y.append(force[1]+ y)
        else:
            force_y.append(y)

        #Calculate the netforce of the vector, if non-zero then plot
        net_force = math.ceil(pow(pow(force[0],2) + pow(force[1],2),0.5))
        if net_force > 0:
            plt.plot(force_x, force_y, linewidth='5', linestyle='--', marker="P", markersize='15', markevery=[1])
            #Annotate the net force of the vector at the end point of the force
            plt.annotate("F: %d" % (net_force),(force_x[1],force_y[1]))

    #Print the location of boundary conditions, indicate which type
    for t in range(0,len(i_disp)):
        #Get x and y coordinates of the current node
        x = node_list[t].node_coords[0]
        y = node_list[t].node_coords[1]

        #If the initial displacement is set to 0, plot a marker showing which type of boundary condition exists
        if i_disp[t][0] == 0 and i_disp[t][1] == 0:
            #Circle marker indicates a pinned joint
            plt.plot(x,y,"o",markersize=15)
        elif i_disp[t][0] == 0:
            # Horizontal triangle indicates rolling pinned joint restricted to y translation(fixed in x direction)
            plt.plot(x, y, ">",markersize=15)
        elif i_disp[t][1] == 0:
            #Vertical triangle indicates rolling pinned joint restricted to x translation (fixed in y direction)
            plt.plot(x, y, "^",markersize=15)

    #Show the figure
    plt.show()
