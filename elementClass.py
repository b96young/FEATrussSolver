#----------------------------------------------------------------------------------------------------------------------#
#FEA Truss Analysis Software
#Created by Brad Young
#Mar 10, 2019
#----------------------------------------------------------------------------------------------------------------------#
import numpy as np
#----------------------------------------------------------------------------------------------------------------------#

#Create object to hold all element parameters
class Element():
    def stiff_matrix_calc(self):
        #Calculates the individual stiffness matrices in global coordinates
        # l = cos(theta)
        # m = sin(theta)
        l = (self.end_coord[0] - self.start_coord[0])/self.length
        m = (self.end_coord[1] - self.start_coord[1])/self.length

        #Create variable to hold stiff_matrix for the element
        self.stiff_matrix = np.multiply(self.E * self.A / self.length, np.array([ \
            [pow(l, 2), l * m, -pow(l, 2), -l * m], \
            [l * m, pow(m, 2), -l * m, -pow(m, 2)], \
            [-pow(l, 2), -l * m, pow(l, 2), l * m], \
            [-l * m, -pow(m, 2), l * m, pow(m, 2)]]))

        self.strain_matrix = np.array([-l, -m, l, m])
        self.stress_matrix = np.multiply(self.E/self.length, np.array([-l, -m, l, m]))

    def element_midpoint(self):
        #Find the midpoint of the element, used to annotate the plot of the structure
        self.midpoint = []
        x = (self.end_coord[0] - self.start_coord[0]) / 2 + self.start_coord[0]
        y = (self.end_coord[1] - self.start_coord[1]) / 2 + self.start_coord[1]
        self.midpoint.append(x)
        self.midpoint.append(y)

    def __init__(self,E,A,node_list,node_connections):
        #Function called when object is created
        self.E = E
        self.A = A
        self.start_coord = node_list[node_connections[0]].node_coords
        self.end_coord = node_list[node_connections[1]].node_coords
        self.connecting_nodes = node_connections
        self.length = pow((pow((self.end_coord[0] - self.start_coord[0]),2) + pow((self.end_coord[1] - self.start_coord[1]),2)),0.5)
        self.augMatrix = []

        #Calculate stiffness matrix when element object is initialized
        self.stiff_matrix_calc()
        #Calculate element midpoint coordinates when element is initialized
        self.element_midpoint()


