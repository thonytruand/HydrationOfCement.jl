using CSV

txt = CSV.File("src/StandardThermodynamicProperties.txt"; header=7, ignoreemptyrows=true, missingstring=["-"], comment="#", delim=" ")
#Free Enthalpy of formation, Enthalpy of formation, Entropy, Heat capacity (as a funciton of temperature), Molar volume
DGf, DHf, DS, Cp, Vm = Dict(), Dict(), Dict(), Dict(), Dict()
for elem in txt
    DGf[elem.phases] = elem.FormationFreeEnthalpy
    DHf[elem.phases] = elem.FormationEnthalpy
    DS[elem.phases] = elem.Entropy
    Cp[elem.phases] = T -> elem.a0 + elem.a1 * T + elem.a2 * T^(-2) + elem.a3 * T^(-0.5)
    Vm[elem.phases] = elem.MolarVolume
end
#---------------------------------------------
function read_StandardThermodynamicProperties(StandardThermodynamicProperties_file)
    
end

