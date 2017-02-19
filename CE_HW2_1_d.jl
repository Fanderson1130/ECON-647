using Optim

#Gamma
function gamma(kappa_r , sigma_r)
    return sqrt(kappa_r^2+2*sigma_r^2)
end


#Big A
function b1(opt)
    kappa_r, theta_r, sigma_r, r = opt
    tau = tau_list
    g = gamma(kappa_r , sigma_r)
    return ((2*g*e^((g+kappa_r)*(tau/2)))/((g+kappa_r)*
    (e^(g*tau)-1)+2*g))^((2*kappa_r*theta_r)/(sigma_r^2))
end

#big B
function b2(opt)
    kappa_r, theta_r, sigma_r, r = opt
    tau = tau_list
    g = gamma(kappa_r , sigma_r)
    return ((2*(e^(g*tau)-1))/((g+kappa_r)*(exp(g*tau)-1)+
    2*g))
end

#Now For Z
function Z(opt)
    kappa_r, theta_r, sigma_r, r = opt
    tau = tau_list
    g = gamma(kappa_r , sigma_r)
    b_1 = b1(opt)
    b_2 = b2(opt)
    return b_1 * e^((-1)*b_2*r)
end

function R_est(opt)
    kappa_r, theta_r, sigma_r, r = opt
    tau = tau_list
    z_1 = Z(opt)
    return (-100)*log(z_1)/tau
end

function MSE_function(opt)
    kappa_r, theta_r, sigma_r, r = opt
    tau = tau_list
    ROR = R_est(opt)
    R = R_list
    MSE = sum((R - ROR) ^ 2) / length(R)
    return MSE
end

function calibration(tau, R)
    tau, R = tau_list, R_list
    res = optimize(MSE_function, zeros(4))
    return Optim.minimizer(res)
end