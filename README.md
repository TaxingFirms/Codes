# Codes
The main file is running.jl

Objects
- Taxes
- Parameters
- Equilibrium
- Firm Problem: value function, policy functions

Firms.jl
INPUT = equilibrium, taxes, parameters, !firmproblem
OUTPUT = firmproblem

Distribution.jl
INPUT = firmproblem, parameters, !equilibrium
OUTPUT = equilibrium.(stationary distribution)

FreeEntry.jl
INPUT = firmproblem, taxes, parameteres, !equilibrium 
OUTPUT = equilibrium.wage

Aggregation.jl
INPUT = firmproblem, !equilibrium
OUTPUT = equilibrium.(mass of entrants,momemnts, aggregates)


