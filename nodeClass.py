#----------------------------------------------------------------------------------------------------------------------#
#FEA Truss Analysis Software
#Created by Brad Young
#Mar 10, 2019
#----------------------------------------------------------------------------------------------------------------------#

class Node():
    #Create object for each node
    def __init__(self,coords,node_num):
        #Holds coordinates of the node
        self.node_coords = coords
        #Holds the number label
        self.node_num = node_num
