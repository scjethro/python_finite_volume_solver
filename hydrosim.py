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

def fluxCalc(cell, test, massFlux, momFlux, energyFlux):
    area = cell._surfaceArea

    if (test == 1):
        cell._mass -= massFlux * area * timestep
        cell._mom -= momFlux * area * timestep
        cell._energy -= energyFlux * area * timestep
        cell.calcPrims()

    else:
        cell._mass += massFlux * area * timestep
        cell._mom += momFlux * area * timestep
        cell._energy += energyFlux * area * timestep
        cell.calcPrims()

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
    cell.calcPrims()
    cellArray.append(cell)

for i in range(0, len(cellArray)):
    cellArray[i-1]._rightCell = cellArray[i]

solver = riemannsolver.RiemannSolver(gamma)

for step in range(0, 200):
    for cell in cellArray:
        rightCell = cell._rightCell

        densSol, velSol, pressSol, _ = solver.solve(cell._density,cell._velocity,cell._pressure,rightCell._density, rightCell._velocity, rightCell._pressure)

        massFlux = densSol*velSol
        momFlux = densSol*velSol**2+pressSol
        energyFlux = (pressSol*gamma/(gamma-1)+0.5*densSol*velSol**2)*velSol

        fluxCalc(cell, 1,massFlux,momFlux,energyFlux)
        fluxCalc(rightCell, 0,massFlux,momFlux,energyFlux)

x = [i for i in range(0,n)]
density = [cell._density for cell in cellArray]

# print("x len", len(x))
# print("density len", len(density))

# plt.scatter(density,x)
# plt.show()