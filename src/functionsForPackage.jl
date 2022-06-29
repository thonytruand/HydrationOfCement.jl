function __init__()\
    Unitful.register(HydrationOfCement)
end

#using PeriodicTable, UnitfulMoles, Unitful
using CSV, DataFrames, DifferentialEquations, Plots, PhysicalConstants.CODATA2018, LsqFit
using Catalyst, Latexify
using PyCall
chemicals = pyimport("chemicals")

file_location = joinpath("src", "cemdata18_std_thermo_properties.txt")
d = DataFrame(CSV.File(file_location))



chemicals.periodic_table.Na
chemicals.serialize_formula("Pd(NH3)4+3")
chemicals.atoms_to_Hill({"H": 5, "C": 2, "Br": 1})
chemicals.atoms_to_Hill(Dict("H" => 5, "C" => 2, "Br" => 1))
chemicals.molecular_weight(Dict("H" => 5, "C" => 2, "Br" => 1))
chemicals.Hfs("101-81-5")
chemicals.reaction.S0s("101-81-5")
chemicals.periodic_table.Na.CAS
chemicals.periodic_table


Avrami_u(u, p) = 0 <= u < 1 ? p[1]/p[2] * (1-u) * (-log(1-u))^(1-p[2]) : 0

function Avrami(u, p, t)
    return Avrami_u(u, p)
end

Jander_u(u, p) = 0 < u < 1 ? (p[1] * (1 - u)^(2/3)) / (1- (1 - u)^(1/3)) : 0

function Jander(u, p, t)
    return Jander_u(u, p)
end

Power_law_u(u, p) = 0 <= u <= 1 ? p[1] * (1 - u)^p[2] : 0

function Power_law(u, p, t)
    return Power_law_u(u, p)
end

function g(u, p, t)
    return min(Avrami(u, [p[1], p[2]], t), Jander(u, p[3], t), Power_law(u, [p[4], p[5]], t))
end

function h(u, p, t)
    return g(u[1], p, t)* Arrhenius(u[2])
end

function Arrhenius(T)
    Eₐ = 40 #u"kJ/mol" #Valeur moyenne des énergies d'activation pour les composants du ciment (Lothenbach et al., 2008)
    T₀ = 293 #u"K"
    return exp(-Eₐ/R*(1/T - 1/T₀))
end


function Parrot_Killoh(α₀, anhydre)
    K1 = Dict("C3S"=>1.5, "C2S"=>0.5, "C3A"=>1, "C4AF"=>0.37)
    N1 = Dict("C3S"=>0.7, "C2S"=>1., "C3A"=>0.85, "C4AF"=>0.7)
    K2 = Dict("C3S"=>0.05, "C2S"=>0.006, "C3A"=>0.04, "C4AF"=>0.015)
    K3 = Dict("C3S"=>1.1, "C2S"=>0.2, "C3A"=>1., "C4AF"=>0.4)
    N3 = Dict("C3S"=>3.3, "C2S"=>5.0, "C3A"=>3.2, "C4AF"=>3.7)

    tspan = (0.0,100.0)
    p_Avrami = [K1[anhydre], N1[anhydre]]
    prob_Avrami = ODEProblem(Avrami,α₀,tspan,p_Avrami)
    p_Jander = K2[anhydre]
    prob_Jander = ODEProblem(Jander,α₀,tspan,p_Jander)
    p_Power_law = [K3[anhydre], N3[anhydre]]
    prob_Power_law = ODEProblem(Power_law,α₀,tspan,p_Power_law)
    prob = ODEProblem(g,α₀,tspan,[K1[anhydre], N1[anhydre], K2[anhydre], K3[anhydre], (N3[anhydre])])
    return solve(prob)
end

function display_Parrot_Killoh()
    plot()
    for elem in ["C3S", "C2S", "C3A", "C4AF"]
        sol_tmp = Parrot_Killoh(1e-1, elem)
        plot!(sol_tmp,linewidth=5,title="Hydration degree",xaxis="Time (t)",yaxis="α(t)",label=elem)
    end
    return current()
end

gr()
display_Parrot_Killoh()

######################
# TEST Catalyst
######################
function Phreeqc_kinetics(m_over_m0, k, A, p, q, n)
    Ω = 0.1 #Saturation degree (for mineral species)
    return (1-Ω^p)^q * A / V * m_over_m0^n * k
end

#Dans les formules d'hydratation des anhydres, on regarde l'évolution de degré d'hydratation
#Ce degré est défini comme suit : α = 1 - m/m0 où m est la quantité de l'anhydre considéré
#et m0 sa quantité initiale
function Avrami_f(α, k, n)
    return k/n * (1-α) * (-log(1-α))^(1-n)
end

function Jander_f(α, k)
    return (k * (1 - α)^(2/3)) / (1- (1 - α)^(1/3))

Power_law_f(u, p) = p[1] * (1 - u)^p[2]

#Catalyst
rn = @reaction_network begin
    f(1-C₃A,p₁), C₃A --> CaO + H₂O + SiO₂
    p₂, CaO + H₂O + SiO₂ --> C₃A
    #hill(X,v,K,n), ∅ --> X
    #v*(X^(2/3)) / (1 - X^(1/3)), 
    #g(X,p,t), X --> ∅
    #k₁/k₂ * (1-X), Y --> ∅
end p₁ p₂ #v K n #k₁ k₂

osys  = convert(ODESystem, rn)

latexify(osys, starred=true)

pmap  = (:p₁ => 1.5e-5, :p₂ => 1.e-11)
u₀map = [:C₃A => 1., :CaO => 0., :H₂O => 1., :SiO₂ => 0.]

# time interval to solve on
tspan = (0., 10000.)

# create the ODEProblem we want to solve
oprob = ODEProblem(rn, u₀map, tspan, pmap)

sol = solve(oprob, Tsit5(), saveat=10.)
plot(sol)

###Ex2
rs = @reaction_network begin
    c1, S + E --> SE
    c2, SE --> S + E
    c3, SE --> P + E
  end c1 c2 c3
  p = (:c1 => 0.00166, :c2 => 0.0001, :c3 => 0.1) 
  tspan = (0., 100.)
  u0 = [:S => 301., :E => 100., :SE => 0., :P => 0.]  
  
  # solve ODEs
  oprob = ODEProblem(rs, u0, tspan, p)
  osol  = solve(oprob, Tsit5())

plot(osol)

######################
# FIN TEST
######################




begin
    sol = solve(  ODEProblem( f, y₀, (0, 10.0), (a, b, forcing(c)) ) );
    f(y,(a, b, forcing(c)),t) = a - b*y + forcing(c)(t)
    forcing(c)  =  t->c*t
    y₀ = 2.0
end

function trace_fig()
    K1 = Dict("C3S"=>1.5, "C2S"=>0.5, "C3A"=>1, "C4AF"=>0.37)
    N1 = Dict("C3S"=>0.7, "C2S"=>1., "C3A"=>0.85, "C4AF"=>0.7)
    K2 = Dict("C3S"=>0.05, "C2S"=>0.006, "C3A"=>0.04, "C4AF"=>0.015)
    hydration_degree = [0.01:0.01:0.99]
    p = [K1["C3A"], N1["C3A"], K2["C3A"]]
    plot(hydration_degree, Avrami_u(hydration_degree, p))
end
trace_fig()

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
