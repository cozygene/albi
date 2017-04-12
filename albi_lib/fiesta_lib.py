import numpy.linalg
from numpy import *
from scipy.stats import norm
import scipy
import sys
import functools

import albi_lib

def binary_robbins_monro_general(function, C, alpha, start, tau, iterations, use_convergence_criterion=False):
    current_endpoint = start

    tau_1 = tau
    mtag = 1
    beta = 1/(mtag*norm.pdf(norm.ppf(alpha)))

    v = (beta**2)*(tau_1**2)
    
    points = []
    const = norm.pdf(norm.ppf(alpha))
    norm_alpha = norm.ppf(alpha)
    
    k = (sqrt(2*pi)*2)/(norm_alpha*exp((-norm_alpha**2)/2))    
    
    jumps_sign = []
    
    for i in range(1, iterations):
        if (use_convergence_criterion and i > 500):
            positive_jumps = sum(jumps_sign[-100:])
            if positive_jumps < 90 and positive_jumps > 10:
                return points[-1]
            
        y = function(current_endpoint)

        t = norm_alpha/sqrt(1+v)
        b = norm.cdf(t)
        c = (v/sqrt(1+v))*scipy.stats.norm.pdf(t)

        jumps_sign.append(y)
        
        v -= c**2/(b*(1-b))
        current_endpoint -= c*(y-b)/(beta*b*(1-b))

        current_endpoint = max(0, current_endpoint)
        current_endpoint = min(current_endpoint, 1)
        
        mtag = abs(k*(current_endpoint - C))
        if mtag == 0:
            mtag = 1        
        
        beta = 1/(mtag*const)

        points.append(current_endpoint)
        
    return points[-1]

def eigenvector_derivative(h2, H2, eigen_values, eigenvectors_as_X):
    calc = albi_lib.OnlyEigenvectorsDerivativeSignCalculator(h2_values=[h2], 
                                                             H2_values=[H2], 
                                                             kinship_eigenvalues=eigen_values, 
                                                             eigenvectors_as_X=eigenvectors_as_X, 
                                                             REML=True)

    return calc.get_derivative_signs(scipy.randn(1, 1, len(eigen_values)))[0, 0, 0]

def binary_robbins_monro_ci(estimated_heritability, eigen_values, eigenvectors_as_X, alpha, start, tau, iterations, use_convergence_criterion=False):
    def invert_quantile(x):
        return eigenvector_derivative(x, estimated_heritability, eigen_values, eigenvectors_as_X) > 0

    return binary_robbins_monro_general(invert_quantile, estimated_heritability, alpha, start, tau, iterations, use_convergence_criterion=False)    

def binary_robbins_monro_boundary_lower(eigen_values, eigenvectors_as_X, alpha, start, tau, iterations, use_convergence_criterion=False):
    def invert_cdf_0(x):
        return eigenvector_derivative(0, x, eigen_values, eigenvectors_as_X) < 0

    return binary_robbins_monro_general(invert_cdf_0, 0, alpha, start, tau, iterations, use_convergence_criterion=False)


def binary_robbins_monro_boundary_upper(eigen_values, eigenvectors_as_X, alpha, start, tau, iterations, use_convergence_criterion=False):
    def invert_cdf_1(x):
        return eigenvector_derivative(1, x, eigen_values, eigenvectors_as_X) < 0

    return binary_robbins_monro_general(invert_cdf_1, 1, alpha, start, tau, iterations, use_convergence_criterion=False)

 
def get_upper_endpoint(eigen_values, eigenvectors_as_X, alpha, tau, h_hat, s, t, s_tag, t_tag, s_tagaim, t_tagaim, iterations, use_convergence_criterion=False):
    start = min(h_hat*1.5, 1)

    # Case 1: s < t
    if s < t:
        if h_hat == 1 or h_hat > t_tag:
            return 1
        if h_hat == 0:
            return s
        
        upper_endpoint_half_alpha = binary_robbins_monro_ci(h_hat, eigen_values, eigenvectors_as_X, 1-alpha/2, start=start, tau=tau, iterations=iterations,
                                                  use_convergence_criterion=use_convergence_criterion)

        if t > upper_endpoint_half_alpha:
            return upper_endpoint_half_alpha

        upper_endpoint_alpha = binary_robbins_monro_ci(h_hat, eigen_values, eigenvectors_as_X, 1-alpha, start=start, tau=tau, iterations=iterations,
                                                  use_convergence_criterion=use_convergence_criterion)
            
        if t < upper_endpoint_alpha:
            return upper_endpoint_alpha
                
        return t

    # Case 2: s > t but s_tagaim < t_tagaim
    elif s_tagaim < t_tagaim:
        delta = (s + t) / 2

        if h_hat == 1 or h_hat > t_tag:
            return 1

        upper_endpoint_alpha = binary_robbins_monro_ci(h_hat, eigen_values, eigenvectors_as_X, 1-alpha, start=start, tau=tau, iterations=iterations,
                                                  use_convergence_criterion=use_convergence_criterion)

        if upper_endpoint_alpha > delta:
            return upper_endpoint_alpha

        return delta

    # Case 3: s < t and s_tagaim < t_tagaim
    else: 
        if h_hat == 1 or h_hat > t_tag:
            return 1

        upper_endpoint_alpha = binary_robbins_monro_ci(h_hat, eigen_values, eigenvectors_as_X, 1-alpha, start=start, tau=tau, iterations=iterations,
                                                  use_convergence_criterion=use_convergence_criterion)

        return upper_endpoint_alpha


def get_lower_endpoint(eigen_values, eigenvectors_as_X, alpha, tau, h_hat, s, t, s_tag, t_tag, s_tagaim, t_tagaim, iterations, use_convergence_criterion=False):
    start = max(h_hat*0.7, 0)
    
    # Case 1: s < t
    if s < t:
        if h_hat == 0 or h_hat < s_tag: 
            return 0
        if h_hat == 1:
            return t
        
        
        lower_endpoint_half_alpha = binary_robbins_monro_ci(h_hat, eigen_values, eigenvectors_as_X, alpha/2, start=start, tau=tau, iterations=iterations,
                                                           use_convergence_criterion=use_convergence_criterion)
            
        if s < lower_endpoint_half_alpha:
            return lower_endpoint_half_alpha
        
        lower_endpoint_alpha = binary_robbins_monro_ci(h_hat, eigen_values, eigenvectors_as_X, alpha, start=start, tau=tau, iterations=iterations,
                                                      use_convergence_criterion=use_convergence_criterion)

        if s > lower_endpoint_alpha:
            return lower_endpoint_alpha
        
        return s
    
    # Case 2: s > t but s_tagaim < t_tagaim
    elif s_tagaim < t_tagaim:
        if h_hat == 0 or h_hat < s_tag: 
            return 0
        
        delta = (s + t) / 2

        if h_hat == 1:
            return delta
        
        start = max(h_hat*0.7, 0)
        
        lower_endpoint_alpha = binary_robbins_monro_ci(h_hat, eigen_values, eigenvectors_as_X, alpha, start=start, tau=tau, iterations=iterations,
                                                      use_convergence_criterion=use_convergence_criterion)

        if lower_endpoint_alpha <= delta:
            return lower_endpoint_alpha

        return delta

    # Case 3: s < t and s_tagaim < t_tagaim
    else: 
        if h_hat == 0 or h_hat < s_tag: 
            return 0

        lower_endpoint_alpha = binary_robbins_monro_ci(h_hat, eigen_values, eigenvectors_as_X, alpha, start=start, tau=tau, iterations=iterations,
                                                    use_convergence_criterion=use_convergence_criterion)            
        return lower_endpoint_alpha
        


def calculate_cis_eigenvectors(kinship_eigenvalues, eigenvectors_as_X, h_hat_values, iterations, alpha, tau=0.4, use_convergence_criterion=False):
    s = binary_robbins_monro_ci(0, kinship_eigenvalues, eigenvectors_as_X, 1-alpha/2, start=0.3, tau=tau, iterations=iterations, use_convergence_criterion=use_convergence_criterion)
    t = binary_robbins_monro_ci(1, kinship_eigenvalues, eigenvectors_as_X, alpha/2, start=0.7, tau=tau, iterations=iterations, use_convergence_criterion=use_convergence_criterion)

    s_tag = binary_robbins_monro_boundary_lower(kinship_eigenvalues, eigenvectors_as_X, 1-alpha, start=0.3, tau=tau, iterations=iterations, use_convergence_criterion=use_convergence_criterion)
    t_tag = binary_robbins_monro_boundary_upper(kinship_eigenvalues, eigenvectors_as_X, alpha, start=0.7, tau=tau, iterations=iterations, use_convergence_criterion=use_convergence_criterion)

    if (s > t):
        s_tagaim = binary_robbins_monro_ci(0, kinship_eigenvalues, eigenvectors_as_X, 1-alpha, start=0.3, tau=tau, iterations=iterations, use_convergence_criterion=use_convergence_criterion)
        t_tagaim = binary_robbins_monro_ci(1, kinship_eigenvalues, eigenvectors_as_X, alpha, start=0.7, tau=tau, iterations=iterations, use_convergence_criterion=use_convergence_criterion)
    else:
        s_tagaim = None
        t_tagaim = None

    print s, t, s_tag, t_tag, s_tagaim, t_tagaim

    upper_endpoints = [get_upper_endpoint(kinship_eigenvalues, eigenvectors_as_X, alpha, tau, h_hat, s, t, s_tag, t_tag, s_tagaim, t_tagaim, iterations, use_convergence_criterion) for h_hat in h_hat_values]
    lower_endpoints = [get_lower_endpoint(kinship_eigenvalues, eigenvectors_as_X, alpha, tau, h_hat, s, t, s_tag, t_tag, s_tagaim, t_tagaim, iterations, use_convergence_criterion) for h_hat in h_hat_values]

    return array([lower_endpoints, upper_endpoints]).T
