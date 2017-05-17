import numpy.linalg
from numpy import *
from scipy.stats import norm
import scipy
import sys
import functools

import progress_bar
import albi_lib

def binary_robbins_monro_general(function, C, alpha, start):
    global rm_tau, rm_iterations, rm_use_convergence_criterion

    current_endpoint = start

    mtag = 1
    beta = 1/(mtag*norm.pdf(norm.ppf(alpha)))

    v = (beta**2)*(rm_tau**2)
    
    points = []
    const = norm.pdf(norm.ppf(alpha))
    norm_alpha = norm.ppf(alpha)
    
    k = (sqrt(2*pi)*2)/(norm_alpha*exp((-norm_alpha**2)/2))    
    
    jumps_sign = []
    
    for i in range(1, rm_iterations):
        if (rm_use_convergence_criterion and i > 500):
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

def eigenvector_derivative(h2, H2, kinship_eigenvalues, eigenvectors_as_X):
    calc = albi_lib.OnlyEigenvectorsDerivativeSignCalculator(h2_values=[h2], 
                                                             H2_values=[H2], 
                                                             kinship_eigenvalues=kinship_eigenvalues, 
                                                             eigenvectors_as_X=eigenvectors_as_X, 
                                                             REML=True)

    return calc.get_derivative_signs(scipy.randn(1, 1, len(kinship_eigenvalues)))[0, 0, 0]

def general_derivative(h2, H2, kinship_eigenvalues, kinship_eigenvectors, covariates):
    calc = albi_lib.GeneralDerivativeSignCalculator(h2_values=[h2], 
                                                    H2_values=[H2], 
                                                    kinship_eigenvalues=kinship_eigenvalues,
                                                    kinship_eigenvectors=kinship_eigenvectors, 
                                                    covariates=covariates,
                                                    REML=True)

    return calc.get_derivative_signs(scipy.randn(1, len(kinship_eigenvalues)))[0, 0]


def binary_robbins_monro_ci(derivative_function, estimated_heritability, alpha, start):
    def invert_quantile(x):
        return derivative_function(x, estimated_heritability) > 0

    return binary_robbins_monro_general(invert_quantile, estimated_heritability, alpha, start)    

def binary_robbins_monro_boundary_lower(derivative_function, alpha, start):
    def invert_cdf_0(x):
        return derivative_function(0, x) < 0

    return binary_robbins_monro_general(invert_cdf_0, 0, alpha, start)


def binary_robbins_monro_boundary_upper(derivative_function, alpha, start):
    def invert_cdf_1(x):
        return derivative_function(1, x) < 0

    return binary_robbins_monro_general(invert_cdf_1, 1, alpha, start)

 
def get_upper_endpoint(derivative_function, alpha, h_hat, s, t, s_tag, t_tag, s_tagaim, t_tagaim):
    start = min(h_hat*1.5, 1)

    # Case 1: s < t
    if s < t:
        if h_hat == 1 or h_hat > t_tag:
            return 1
        if h_hat == 0:
            return s
        
        upper_endpoint_half_alpha = binary_robbins_monro_ci(derivative_function, h_hat, 1-alpha/2, start=start)

        if t > upper_endpoint_half_alpha:
            return upper_endpoint_half_alpha

        upper_endpoint_alpha = binary_robbins_monro_ci(derivative_function, h_hat, 1-alpha, start=start)
            
        if t < upper_endpoint_alpha:
            return upper_endpoint_alpha
                
        return t

    # Case 2: s > t but s_tagaim < t_tagaim
    elif s_tagaim < t_tagaim:
        delta = (s + t) / 2

        if h_hat == 1 or h_hat > t_tag:
            return 1

        upper_endpoint_alpha = binary_robbins_monro_ci(derivative_function, h_hat, 1-alpha, start=start)

        if upper_endpoint_alpha > delta:
            return upper_endpoint_alpha

        return delta

    # Case 3: s < t and s_tagaim < t_tagaim
    else: 
        if h_hat == 1 or h_hat > t_tag:
            return 1

        upper_endpoint_alpha = binary_robbins_monro_ci(derivative_function, h_hat, 1-alpha, start=start)

        return upper_endpoint_alpha


def get_lower_endpoint(derivative_function, alpha, h_hat, s, t, s_tag, t_tag, s_tagaim, t_tagaim):
    start = max(h_hat*0.7, 0)
    
    # Case 1: s < t
    if s < t:
        if h_hat == 0 or h_hat < s_tag: 
            return 0
        if h_hat == 1:
            return t
        
        
        lower_endpoint_half_alpha = binary_robbins_monro_ci(derivative_function, h_hat, alpha/2, start=start)
            
        if s < lower_endpoint_half_alpha:
            return lower_endpoint_half_alpha
        
        lower_endpoint_alpha = binary_robbins_monro_ci(derivative_function, h_hat, alpha, start=start)

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
        
        lower_endpoint_alpha = binary_robbins_monro_ci(derivative_function, h_hat, alpha, start=start)

        if lower_endpoint_alpha <= delta:
            return lower_endpoint_alpha

        return delta

    # Case 3: s < t and s_tagaim < t_tagaim
    else: 
        if h_hat == 0 or h_hat < s_tag: 
            return 0

        lower_endpoint_alpha = binary_robbins_monro_ci(derivative_function, h_hat, alpha, start=start)
        return lower_endpoint_alpha
        


def calculate_cis_internal(h_hat_values, derivative_function, alpha, use_progress_bar):
    s = binary_robbins_monro_ci(derivative_function, 0, 1-alpha/2, start=0.3)
    t = binary_robbins_monro_ci(derivative_function, 1, alpha/2, start=0.7)

    s_tag = binary_robbins_monro_boundary_lower(derivative_function, 1-alpha, start=0.3)
    t_tag = binary_robbins_monro_boundary_upper(derivative_function, alpha, start=0.7)

    if (s > t):
        s_tagaim = binary_robbins_monro_ci(derivative_function, 0, 1-alpha, start=0.3)
        t_tagaim = binary_robbins_monro_ci(derivative_function, 1, alpha, start=0.7)
    else:
        s_tagaim = None
        t_tagaim = None

    lower_endpoints = []
    upper_endpoints = []

    if use_progress_bar:
        it = progress_bar.ProgressBarIter(len(h_hat_values))
    else:
        it = range(len(h_hat_values))

    n = 0
    for i in it:
        h_hat = h_hat_values[n]        
        n += 1
        upper_endpoints.append(get_upper_endpoint(derivative_function, alpha, h_hat, s, t, s_tag, t_tag, s_tagaim, t_tagaim))
        lower_endpoints.append(get_lower_endpoint(derivative_function, alpha, h_hat, s, t, s_tag, t_tag, s_tagaim, t_tagaim))

    return array([lower_endpoints, upper_endpoints]).T

def calculate_cis_eigenvectors(h_hat_values, kinship_eigenvalues, eigenvectors_as_X, iterations, alpha, tau=0.4, use_convergence_criterion=False, use_progress_bar=True):
    global rm_tau, rm_iterations, rm_use_convergence_criterion
    rm_tau, rm_iterations, rm_use_convergence_criterion = tau, iterations, use_convergence_criterion

    derivative_function = functools.partial(eigenvector_derivative, kinship_eigenvalues=kinship_eigenvalues, eigenvectors_as_X=eigenvectors_as_X)
    return calculate_cis_internal(h_hat_values, derivative_function, alpha, use_progress_bar)

def calculate_cis_general(h_hat_values, kinship_eigenvalues, kinship_eigenvectors, covariates, iterations, alpha, tau=0.4, use_convergence_criterion=False, use_progress_bar=True):
    global rm_tau, rm_iterations, rm_use_convergence_criterion
    rm_tau, rm_iterations, rm_use_convergence_criterion = tau, iterations, use_convergence_criterion

    derivative_function = functools.partial(general_derivative, kinship_eigenvalues=kinship_eigenvalues, kinship_eigenvectors=kinship_eigenvectors, covariates=covariates)
    return calculate_cis_internal(h_hat_values, derivative_function, alpha, use_progress_bar)

