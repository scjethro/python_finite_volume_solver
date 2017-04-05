import riemannsolver
import matplotlib.pyplot as plt

class Cell():

    def __init__(self):
        self._vol = 0.01
        self._mass = 0
        self._mom = 0
        self._energy = 0
        self._density = 0
        self._velocity = 0
        self._pressure = 0
        self._name = 0

        self._rightCell = None
        self._surfaceArea = 1



    # @brief calculate the primitive variables all at once to save time in the loop

    def calcPrims(self):
        self._density = self._mass/self._vol
        self._velocity = self._mom/self._mass
        self._pressure = (gamma - 1)*(self._energy/self._vol - 0.5*self._density*self._velocity**2)

gamma = 5/3
timestep = 0.001

# hydroprogram

cellArray = []
n = 1000

# @brief calculate the flux through all the surfaces
#
# made into a function taking everything needed
# use the variable test to define 1 or 0 and thus whether you should use + or - on the flux
# after updating variables it calculates the primitive variables again, a step I forgot previously
#

def fluxCalc(cell, test, massFlux, momFlux, energyFlux):
    area = cell._surfaceArea

    if (test == 1):
        cell._mass -= massFlux * area * timestep
        cell._mom -= momFlux * area * timestep
        cell._energy -= energyFlux * area * timestep

    elif(test == 0):
        cell._mass += massFlux * area * timestep
        cell._mom += momFlux * area * timestep
        cell._energy += energyFlux * area * timestep

    cell.calcPrims()

# look over the number of cells and populate the list

for i in range(0,n):
    cell = Cell()
    cell._name = i
    print(i)

    if(i < n/2):
        cell._mass = 0.01
        cell._energy = 0.01/(gamma - 1)

    else:
        cell._mass = 0.00125
        cell._energy = 0.001/(gamma - 1)

    # using the defined method to calculate the initial primitive variables
    cell.calcPrims()
    cellArray.append(cell)

# loop over all the cells again and set all the left neighbours to the cells on the right
# I tried to do this within the loop above but couldn't get it to work
for i in range(0, len(cellArray)):
    cellArray[i-1]._rightCell = cellArray[i]

# create the Riemann solver from bwvdnbro/python_finite_volume_solver
solver = riemannsolver.RiemannSolver(gamma)

# define a maximum number of timesteps
maxT = 200;

# loop over the timesteps to repeatedly calculate the Riemann problem and work out the results
for step in range(0, maxT):
    for cell in cellArray:
        # per cell we find the right hand cell
        rightCell = cell._rightCell

        # using the given Riemann solver we solve the Riemann equation
        densSol, velSol, pressSol, _ = solver.solve(cell._density,cell._velocity,cell._pressure,rightCell._density, rightCell._velocity, rightCell._pressure)

        # calculate the fluxes as given in the equations we derived in class, again thanks to bwvdnbro
        massFlux = densSol*velSol
        momFlux = densSol*velSol**2+pressSol
        energyFlux = (pressSol*gamma/(gamma-1)+0.5*densSol*velSol**2)*velSol

        # use the defined function above to calculate the fluxes
        # everything is taken care of in the function
        fluxCalc(cell, 1,massFlux,momFlux,energyFlux)
        fluxCalc(rightCell, 0,massFlux,momFlux,energyFlux)

# loops to just extract the density and position in an easy manner for plotting
x = [i for i in range(0,n)]
density = [cell._density for cell in cellArray]

# plot the results as a scatter plot to see the output
# plt.scatter(density,x)
# plt.show()