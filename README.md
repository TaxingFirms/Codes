# Codes
The main file is running.jl

Objects
- Taxes
- Parameters
- Equilibrium
- Firm Problem: value function, policy functions

ValueFunctionIteration
INPUT = equilibrium, taxes, parameters, !firmproblem
OUTPUT = firmproblem

Distribution
INPUT = firmproblem, parameters, !equilibrium
OUTPUT = equilibrium.(stationary distribution)

FreeEntry
INPUT = firmproblem, taxes, parameteres, !equilibrium 
OUTPUT = equilibrium.wage

MktClearing
INPUT = firmproblem, !equilibrium
OUTPUT = equilibrium.(mass of entrants)


