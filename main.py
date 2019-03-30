#----------------------------------------------------------------------------------------------------------------------#
#FEA Truss Analysis Software
#Created by Brad Young
#Mar 10, 2019
#------------------------------------------------IMPORTS---------------------------------------------------------------#
#IMPORT other python files from this folder
from projectFunc import *
from elementClass import *
from nodeClass import *
from gaussSeidel import *
from loadJSON import *
#---------------------------READ initial conditions from JSON dictionary using loadJSON.py-----------------------------#
file = 'initialConditions5.json'
IC = loadJSON(file).data["FEA"]

#Initialize variables read in json file (elements, nodes, Youngs modulus, area, initial conditions)
element_params = IC['element_params']
element_connections = element_params["node_connections"]
initial_nodes = IC['nodes']
E = element_params["E"]
D = element_params["D"]
initial_displacement = IC['initial_conds']["displacement"]
initial_force = IC["initial_conds"]["force"]
initial_force = transformForce(initial_force)

#Create list to hold all node objects
nodes = []
for i in range(0,len(initial_nodes)):
    #Node {from nodeClass} object requires x,y coordinates and the current index to create
    nodes.append(Node(initial_nodes[i],i))

#Find total number of nodes, used later to create augmented matrix with correct shape
num_nodes = len(nodes)

#Create a list to hold all element objects
elements = []
for i in range(0,len(element_connections)):
    # For each element, create a Element object as seen in ElementClass.py and append to elements list
    elements.append(Element(E[i],diaToArea(D[i]),nodes,element_connections[i]))
    #Calculate the augmented matrix of the element and set its augMatrix variable {augmentMatrix in projectFunc.py}
    elements[i].augMatrix = augmentMatrix(elements[i], num_nodes)

#Create schematic view of truss to verify correct element, node, and force placements {plotTruss in projectFunc.py}
plotTruss(nodes,elements,initial_displacement,initial_force)

#If truss diagram is accurate, proceed with analysis
proceed = input("Correct truss analysis? Type y to continue")
if proceed != 'y':
    exit()

#Combine all individual element augmented matrices together by matrix addition to form global stiffness matrix
#{combineAugmented in projectFunc.py}
globalStiffnessMatrix = combineAugmented(elements)

#Condense global stiffness matrix to solve for non-zero displacements {using trimMatrix located in projectFunc.py}
trimmedMatrix = trimMatrix(initial_displacement,globalStiffnessMatrix)
#Remove entries in force vector corresponding to 0 displacement entries in displacement vector
forceVector = trimMatrix(initial_displacement,initial_force,vector=True)
#Remove the 0 displacement entries in displacement vector
dispVector = trimMatrix(initial_displacement,initial_displacement,vector=True)

#------------All vectors should now have compatible dimensions for solution using Gauss-Seidel method------------------------#
#Iteratively calculate the matrix solution {uses calcGS function located in gaussSeidel.py}
result_vector = calcGS(trimmedMatrix,dispVector,forceVector)

#Add back 0 displacement entries that were previously deleted {recompileMatrix located in projectFunc.py}
final_displacement = recompileMatrix(initial_displacement,result_vector)

#Call {printTable in projectFunc.py} to print out the final displacements to console
printTable("disp",final_displacement)

#Calculate stresses, then print to console {calcStress in projectFunc.py}
calc_stresses = calcStress(elements,final_displacement)
printTable("stress",calc_stresses)

#Calculate reaction forces, print to console {reactionForces in projectFunc.py}
react_forces = reactionForces(globalStiffnessMatrix,final_displacement)
printTable("reaction",react_forces)

