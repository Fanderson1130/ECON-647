#Gamma
function gamma(kappa_r , sigma_r)
    return sqrt(kappa_r^2+2*sigma_r^2)
end

#Big A
function A(params)
    params = kappa_r, theta_r, sigma_r, r, tau
    g = gamma(kappa_r , sigma_r)
    return ((2*g*e^((g+kappa_r)*(tau/2)))/((g+kappa_r)*
    (e^(g*tau)-1)+2*g))^((2*kappa_r*theta_r)/(sigma_r^2))
end

#Big B 
function B(params)
    params = kappa_r, theta_r, sigma_r, r, tau
    g = gamma(kappa_r , sigma_r)
    return ((2*(e^(g*tau)-1))/((g+kappa_r)*(exp(g*tau)-1)+
    2*g))
end

#Now for Z
function Z(params)
    params = kappa_r, theta_r, sigma_r, r, tau
    b_1 = A(params)
    b_2 = B(params)
    return b_1 * e^((-1)*b_2*r)
end


#doin it the python way
function CIR_fwd_r(opt)
    kappa_r, theta_r, sigma_r = opt
    tau = tau_list
    g = gamma(kappa_r , sigma_r)
    return ((kappa_r*theta_r*(e^(g*tau)-1))/
    ((g+kappa_r)*(exp(g*tau)-1)+2*g))+
    r0*((4*g^2*e^(g*tau))/((g+kappa_r)*(exp(g*tau)-1)+2*g)^2)
end

function MSE_function(opt)
    kappa_r, theta_r, sigma_r = opt
    forward_rates = CIR_fwd_r(opt)
    MSE = sum((R - forward_rates) ^ 2) / length(R)
    # print opt, MSE
    return MSE
end


tau_list = 1.0
R = 4.44

#doin it the python way
function CIR_fwd_r(opt)
    kappa_r, theta_r, sigma_r = opt
    tau = tau_list
    g = gamma(kappa_r , sigma_r)
    return ((kappa_r*theta_r*(e^(g*tau)-1))/
    ((g+kappa_r)*(exp(g*tau)-1)+2*g))+
    r0*((4*g^2*e^(g*tau))/((g+kappa_r)*(exp(g*tau)-1)+2*g)^2)
end

function MSE_function(opt)
    kappa_r, theta_r, sigma_r = opt
    forward_rates = CIR_fwd_r(opt)
    MSE = sum((R - forward_rates) ^ 2) / length(R)
    # print opt, MSE
    return MSE
end

function CIR_calibration(tau, R)
    tau = tau_list
    R = R
    res = optimize(MSE_function, [0.0,0.0,0.0,0.0])
    return Optim.minimizer(res)
end


#Modified Python Approach






