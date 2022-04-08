function __init__()\
    Unitful.register(HydrationOfCement)
end

#using PeriodicTable, UnitfulMoles, Unitful
using CSV, DataFrames

file_location = joinpath("src", "cemdata18_std_thermo_properties.txt")
d = DataFrame(CSV.File(file_location))

function Avrami(u, p, t)
    k = p[1]
    n = p[2]
    if 0 <= u < 1
        return k/n * log(1-u) * (-log(1-u))^(1-n)
    end
end

function Jander(u, p, t)
    if 0 < u
        return (p * (1 - u)^(2/3)) / (1- (1 - u)^(1/3))
    end
end

function Power_law(α, k , n)
    if 0 <= α <= 1
        return k * (1 - α)^n
    end
end

function Parrot_Killoh(α, k, n)
#using DifferentialEquations
#p=...
#u0=...
#prob = ODEProblem(Avrami,u₀,tspan,p)
#sol = solve(prob)
end

struct Molecule
    name::String
    chemicalFormula::AbstractString
end

struct Anhydre{Molecule}
    name::String
    chemicalFormula::AbstractString
end

function read_molecule(fich)
    #TODO
end

"Defines a cool function. Returns some stuff"
function coolFunc(x::Int)
    return x+2
end
 
"""
Defines an even cooler function. ``LaTeX``.
 
 
### Returns
 * Markdown works in here
"""
function coolFunc2(y::Int)
    return y+3
end
