#Packages
using DataFrames
using DataArrays

#Load Data & Convert to Array
cd("$(homedir())/Desktop/School/Econ 647 - Applied Computational/Problem Sets/PS2/Data")
pwd()
df = readtable("Data_2.csv")
Rates = convert(Array, df[2:10])
tau = [0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 30.0]

A(tau, kappa_r, alpha_r, sigma_r) = (2*sqrt(kappa_r^2+2*sigma_r^2)*exp((sqrt(kappa_r^2+2*sigma_r^2)+kappa_r)*tau/2) /
((sqrt(kappa_r^2+2*sigma_r^2)+kappa_r)*(exp(sqrt(kappa_r^2+2*sigma_r^2)*tau)-1)+ 2*sqrt(kappa_r^2+2*sigma_r^2)))^(2*kappa_r*alpha_r/sigma_r^2)

B(tau, kappa_r, alpha_r, sigma_r) = (2*(exp(sqrt(kappa_r^2+2*sigma_r^2)*tau)-1) /
((sqrt(kappa_r^2+2*sigma_r^2)+kappa_r)*(exp(sqrt(kappa_r^2+2*sigma_r^2)*tau)-1)+ 2*sqrt(kappa_r^2+2*sigma_r^2)))

Z(tau, kappa_r, alpha_r, sigma_r, r) = A(tau, kappa_r, alpha_r, sigma_r) * exp(-B(tau, kappa_r, alpha_r, sigma_r)*r)

ROR(tau, kappa_r, alpha_r, sigma_r, r) = -100*log(Z(tau, r, kappa_r, alpha_r, sigma_r))/tau


#Objective Function (squared residuals)
function Sum_sq(Rates, tau, r, kappa_r, theta_r, sigma_r)
    
    lastDate, lastTau = size(Rates)
    
    resid = copy(Rates)
    
    for i in 1:lastDate
        for k in 1:lastTau
            resid[i,k] = (Rates[i,k] -  ROR(tau[k], r[i], kappa_r, theta_r, sigma_r))^2
        end
    end
    
    return sum(resid)
    
end


#optimization for part A
using Optim

resid_a = Array(Float64, 5, 4)

for j in 1:5
    opt_a(y) = Sum_sq(Rates[j,:]', tau, y[1], y[2], y[3], y[4])
    resid_a[j, :] = Optim.minimizer(optimize(opt_a, ones(4), BFGS()))
end

r = resid_a[:,1]
kappa_r = resid_a[:,2]
alpha_r = resid_a[:,3]
sigma_r = resid_a[:,4]

println("r = ", r)
println("kappa_r = ", kappa_r)
println("alpha_r = ", alpha_r)
println("sigma_r = ", sigma_r)

println("kappa_r*alpha_r = ", kappa_r.*alpha_r)




#Part B

F_b(y) = Sum_sq(Rates, tau, [y[1], y[2], y[3], y[4], y[5]], y[6], y[7], y[8]) 
resid_b = Optim.minimizer(optimize(F_b, ones(8), BFGS()))

r = resid_b[1:5]
kappa_r = resid_b[6]
alpha_r = resid_b[7] 
sigma_r = resid_b[8]

println("r = ", r)
println("kappa_r = ", kappa_r)
println("alpha_r = ", alpha_r)
println("sigma_r = ", sigma_r)

println("kappa_r * alpha_r = ", kappa_r.*alpha_r)