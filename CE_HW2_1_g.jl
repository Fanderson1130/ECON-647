using DataFrames
using DataArrays
using Optim
using QuantEcon: gridmake!
using QuantEcon: fix
using QuantEcon: ckron
using Interpolations

cd("$(homedir())/Desktop/School/Econ 647 - Applied Computational/Problem Sets/PS2/Data")
pwd()
df = readtable("Data_2.csv")

#Time Values
t_list = [0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 30.0]

#Setup R matrix
r_matrix = convert(DataArray, df[:2:10])

#Set R values to first row
r_list = convert(Array, r_matrix[1,1:9])/100

#set r0
r0 = r_list[1]

#Gamma
function gamma(kappa_r , sigma_r)
    return sqrt(kappa_r^2+2*sigma_r^2)
end

factors = (1 + t_list *. r_list)
