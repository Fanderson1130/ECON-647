#Problem Set 2 Problem One Bond
#Data
using DataFrames
using DataArrays
using Optim
using QuantEcon: gridmake!
using QuantEcon: fix
using QuantEcon: ckron

cd("$(homedir())/Desktop/School/Econ 647 - Applied Computational/Problem Sets/PS2/Data")
pwd()
df = readtable("Data_2.csv")

#Time Values
t_list = [0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 30.0]

#Setup R matrix
r_matrix = convert(DataArray, df[:2:10])

#Set R values to first row
r_list = convert(Array, r_matrix[1,1:9])

#set r0
r0 = r_list[1]

#Gamma
function gamma(kappa_r , sigma_r)
    return sqrt(kappa_r^2+2*sigma_r^2)
end

function MSE_function(opt)
    kappa_r, theta_r, sigma_r = opt
    g = gamma(kappa_r , sigma_r)
    fx = []
    for i in 1:length(r_list)
        tau = t_list[i]
        forward_rates = ((kappa_r*theta_r*(e^(g*tau)-1))/
        ((g+kappa_r)*(exp(g*tau)-1)+2*g))+
        r0*((4*g^2*e^(g*tau))/((g+kappa_r)*(exp(g*tau)-1)+2*g)^2)
        R = r_list[i]
        sq_dif = (R - forward_rates)^2
        push!(fx, sq_dif)
    end
    MSE = sum(fx) / length(r_list)
    return MSE
end

function CIR_calibration()
    res = optimize(MSE_function, [0.0,0.0,0.0,0.0])
    return Optim.minimizer(res)
end

CIR_calibration()