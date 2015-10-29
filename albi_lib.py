from numpy import *

import numpy.linalg, numpy.random
import socket, os.path, cPickle, csv, bisect, time, tempfile, subprocess

# Use progress bar is available
itrange = locals().get('trange', range)

def weights_zero_derivative(h2_values, H2_values, kinship_eigenvalues, 
                            eigenvectors_as_X=[-1], REML=True):
    """
    Calculates the weights in the expression of the derivative of the likelihood, in the case that
    the fixed effects X are eigenvectors of the kinship matrix. Each weight is parametrized by 
    three parameters: (i) h^2 - the true value of heritability; (ii) H^2 - the estimated value; 
    (iii) - an index.

    Arguments:
        h2_values - a vector of size N, of all possible values of h^2
        H2_values - a vector of size M, of all possible values of H^2
        kinship_eigenvalues - A vector of size K of the eigenvalues of the kinship matrix, in decreasing order.
        eigenvectors_as_X - A list of indices, of which eigenvectors of the kinship matrix are fixed effects.
        REML - True is REML, False if ML.

    Returns:
        A matrix of size N x M x K, of the weights for the cartesian product of h2_values, H2_values and an index.
    """
    if eigenvectors_as_X is None:
        eigenvectors_as_X = []
    if isinstance(h2_values, int) or isinstance(h2_values, float):
        h2_values = [h2_values]
    if isinstance(H2_values, int) or isinstance(H2_values, float):
        H2_values = [H2_values]
    n_samples = len(kinship_eigenvalues)

    H2_values = reshape(H2_values, (1, len(H2_values), 1))
    h2_values = reshape(h2_values, (len(h2_values), 1, 1))
    kinship_eigenvalues = reshape(kinship_eigenvalues, (1, 1, n_samples))

    projection = ones((1, 1, n_samples))
    projection[0, 0, eigenvectors_as_X] = 0
    
    ds = (kinship_eigenvalues - 1) / (H2_values * (kinship_eigenvalues - 1) + 1)
    denom = n_samples

    if REML:
        ds = projection * (kinship_eigenvalues - 1) / (H2_values * (kinship_eigenvalues-1) + 1)
        denom = n_samples - len(eigenvectors_as_X)

    return projection * (h2_values * (kinship_eigenvalues - 1) + 1) / \
                        (H2_values * (kinship_eigenvalues - 1) + 1) \
                                      * (ds - sum(ds, 2)[:, :, newaxis] / denom)


def calculate_probability_intervals(true_h2, interval_boundaries, kinship_eigenvalues, 
                                    eigenvectors_as_X=[-1], REML=True, monte_carlo_size=1, n_chunks=100, seed=0):
    """
    Calculates the probability of 

    Arguments:
        h2_values - a vector of size N, of all possible values of h^2
        H2_values - a vector of size M, of all possible values of H^2
        kinship_eigenvalues - A vector of size K of the eigenvalues of the kinship matrix, in decreasing order.
        eigenvectors_as_X - A list of indices, of which eigenvectors of the kinship matrix are fixed effects.
        REML - True is REML, False if ML.

    Returns:
        A matrix of size N x M x K, of the weights for the cartesian product of h2_values, H2_values and an index.
    """
    n_samples = len(kinship_eigenvalues)
    n_intervals = len(interval_boundaries)-1

    # Size: n_true_h2 X n_grid X n_individuals
    weights = weights_zero_derivative(true_h2, interval_boundaries, kinship_eigenvalues, eigenvectors_as_X=eigenvectors_as_X, REML=REML)
    
    rng = numpy.random.RandomState(seed)
    prob = zeros((len(true_h2), n_intervals+2))

    for i in itrange(n_chunks):
        # Avoid replicating weights across h2's, so that each sample is independent
        us = rng.normal(size=(monte_carlo_size, len(true_h2), 1, n_samples))
        dotproducts = sum((us**2) * weights[newaxis, :, :, :], 3)     
        
        #prob[:, 1:-1] += mean((dotproducts[:, :, :-1] >= 0) & (dotproducts[:, :, 1:] <= 0), 0)
        #prob[:, 0] += mean(dotproducts[:, :, 0] <= 0, 0)        
        #prob[:, -1] += mean(dotproducts[:, :, -1] >= 0, 0) 
        
        # Size: monte_carlo_size X len(true_h2)
        hit_zero = (dotproducts[:, :, 0] <= 0)
        hit_one = (dotproducts[:, :, -1] >= 0)
        hit_boundary = hit_zero | hit_one

        # Size: len(true_h2)
        prob[:, 0] += mean(hit_zero, 0)        
        prob[:, -1] += mean(hit_one, 0)        
        prob[:, 1:-1] += mean(((dotproducts[:, :, :-1] >= 0) & (dotproducts[:, :, 1:] <= 0)) & ~hit_boundary[:, :, newaxis], 0)
        
    prob /= n_chunks       # Average across chunks
    prob /= sum(prob, 1)[:, newaxis]   # Normalize so sum is 1, in case of multiple local maxima
    return prob



def get_approx_inv_cdf(true_h2s, cdf, alpha, include_one=True):
    """
    Find a value c that satisfies Pr(X <= c) >= alpha.

    Delicate because of the discontinuity points at 0,1. More docs here!!!
    """
    epsilon = 1e-10

    # max{x | P(X <= x) == p_0}
    L = bisect.bisect_right(cdf, cdf[0])-1

    # min{x | P(X <= x) == 1-p_1}
    R = bisect.bisect_left(cdf, cdf[-2])

    # Now, the open interval (true_h2s[L], true_h2s[R]) maps 1-1 to the open interval (p0, 1-p1)
    if (cdf[0] < alpha < cdf[-2]):
        subrange = true_h2s[L:R+1]
        subcdf = cdf[L:R+1]
        idx = bisect.bisect_left(subcdf, alpha)
        if subcdf[idx] == alpha: # If we hit exactly right, simplify things and just return that
            return subrange[idx]
        else:
            # Otherwise, the required value is explicity inside the interval (subcdf[idx-1], subcdf[idx]),
            # so to get the right value, interpolate linearly
            proportion = (alpha - subcdf[idx-1]) / (subcdf[idx] - subcdf[idx-1])
            value = subrange[idx-1] + proportion * (subrange[idx] - subrange[idx-1])
            return value

    # If it's not in that open interval, this either alpha <= p0 or alpha >= 1-p1
    elif (alpha <= cdf[0]):
        if not numpy.isclose(alpha, cdf[0]):   # Should have been (alpha < cdf[0]). Impossible
            #return -inf
            return 0.0
        else:  # alpha == p0, and we want to take the largest of these values
            return true_h2s[L]
    elif (alpha >= cdf[-2]):
        if alpha >= cdf[-1] or numpy.isclose(alpha, cdf[-1]):   
            return 1.0
        else:
            if include_one:
                return 1.0
            else:
                return 1.0 - epsilon
    else:
        assert "Should not happen"

def get_approx_cdf(true_h2s, cdf, real_value):
    """
    Find P(X <= real_value)
    """
    assert 0 <= real_value <= 1
    if real_value == 0.0:
        return cdf[0]
    elif real_value == 1:
        return 1.0
    else:
        idx = bisect.bisect_left(true_h2s, real_value)
        if true_h2s[idx] == real_value:
            return cdf[idx]
        else:        
            proportion = (real_value - true_h2s[idx-1]) / (true_h2s[idx] - true_h2s[idx-1])
            return cdf[idx-1] + (cdf[idx] - cdf[idx-1]) * proportion 

def get_probabilty_in_closed_interval(true_h2s, cdf, interval):
    """
    Pr(X in [interval[0], interval[1]])
    """
    p = get_approx_cdf(true_h2s, cdf, interval[1]) - get_approx_cdf(true_h2s, cdf, interval[0])
    if interval[0] == 0.0:
        p += cdf[0]
    return p

def build_heritability_cis(all_distributions, true_h2s, ar_h2s, values, gamma, prob0, prob1, seed=0, use_randomized_cis=False):  

    def build_ars(true_h2s, distribution, real_value, gamma, verbose=False):
        numerical_error = 1e-5
        epsilon = 1e-10   # Small value that is just above 0 but still lower than all actual values probably
        alpha = 1-gamma
        cdf = cumsum(distribution)

        # Pr(\hat{h2} = 0)
        p0 = distribution[0]

        # Pr(\hat{h2} = 1)
        p1 = distribution[-1]

        # Pr(0 < \hat{h2} <= h2)
        pl = get_approx_cdf(true_h2s, cdf, real_value) - p0

        # Pr(h2 < \hat{h2} < 1)
        pr = 1 - p1 - p0 - pl
        
        if verbose: print p0, pl, pr, p1

        low1 = 0.0
        high1 = get_approx_inv_cdf(true_h2s, cdf, 1-alpha, True)
        #coverage1 = get_probabilty_in_closed_interval(true_h2s, cdf, [low1, high1])

        low2 = get_approx_inv_cdf(true_h2s, cdf, alpha/2, True)
        high2 = get_approx_inv_cdf(true_h2s, cdf, 1 -  alpha/2, True)
        #coverage2 = get_probabilty_in_closed_interval(true_h2s, cdf, [low2, high2])        

        low3 = get_approx_inv_cdf(true_h2s, cdf, alpha, True)
        high3 = 1.0
        #coverage3 = get_probabilty_in_closed_interval(true_h2s, cdf, [low3, high3])        

        #return [low1, high1, coverage1], [low2, high2, coverage2], [low3, high3, coverage3]
        return [low1, high1], [low2, high2], [low3, high3]
    
    accept_regions = zeros((len(true_h2s), 3, 3))
    accept_regions[:,:,:2] = array([build_ars(ar_h2s, all_distributions[i,:], h2, gamma, False) for i,h2 in enumerate(true_h2s)])

    # Regularize for noise:
    # 1. Make sure h^2 is included 
    accept_regions[:,0,1] = maximum(accept_regions[:,0,1], true_h2s)
    accept_regions[:,1,0] = minimum(accept_regions[:,1,0], true_h2s)
    accept_regions[:,1,1] = maximum(accept_regions[:,1,1], true_h2s)
    accept_regions[:,2,0] = minimum(accept_regions[:,2,0], true_h2s)

    # 2. Make sure type II boundaries are over type I and III where relevant
    accept_regions[:,1,1] = maximum(accept_regions[:,1,1], accept_regions[:,0,1])
    accept_regions[:,1,0] = minimum(accept_regions[:,1,0], accept_regions[:,2,0])

    # 3. Make sure they are monotone
    accept_regions[:,0,1] = maximum.accumulate(accept_regions[:,0,1])
    accept_regions[:,1,0] = maximum.accumulate(accept_regions[:,1,0])
    accept_regions[:,1,1] = maximum.accumulate(accept_regions[:,1,1])
    accept_regions[:,2,0] = maximum.accumulate(accept_regions[:,2,0])

    # Calculate coverages
    for i,h2 in enumerate(true_h2s):
        for j in range(3):
            cdf = cumsum(all_distributions[i,:])  
            accept_regions[i,j,2] = get_probabilty_in_closed_interval(ar_h2s, cdf, accept_regions[i,j,:2])

    #return accept_regions

    # Find the first place type I covers 1
    s = where(numpy.isclose(accept_regions[:,0,2], 1))[0][0]

    # Find the last place type III covers 1
    t = where(numpy.isclose(accept_regions[:,2,2], 1))[0][-1]
    d = int((s+t)/2)

    # Find the regions where type II covers exactly 
    w = where(numpy.isclose(accept_regions[:,1,2], gamma))[0]

    # If type II covers a range of h^2, then by definition type I covers before and type III after
    if len(w):
        regions = vstack([accept_regions[:w[0],     0,  :], 
                          accept_regions[w[0]:w[-1], 1, :],
                          accept_regions[w[-1]:    , 2, :]])
    else:
        # Otherwise, either type I and III cover everything (and we're good)
        # or we resort to randomized CIs. In both cases we do this:
        regions = vstack([accept_regions[:d, 0, :],                               
                          accept_regions[d:, 2, :]])   
    
    starts = maximum.accumulate(regions[:,0])
    ends = minimum.accumulate(regions[:,1][::-1])[::-1]

    rng = numpy.random.RandomState(seed)

    # This is relevant only is w is empty, s < t, and use_randomized_cis=True
    u = 1 - (1-gamma) / prob0
    u[:d] = 1
    u[t:] = 0
    u = minimum.accumulate(maximum(0, u))  # denoise and make monotone

    l = 1 - (1-gamma) / prob1
    l[:s] = 0
    l[d:] = 1
    l = maximum.accumulate(maximum(0, l))  # denoise and make monotone

    #print s,t, d, w

    cis = zeros((len(values), 2))
    for i,estimated_h in enumerate(values):
        if use_randomized_cis and s < t:
            rand_value = rng.rand()
            if estimated_h == 0.0:                
                cis[i, 0] = 0.0                
                cis[i, 1] = true_h2s[where(rand_value < u)[0][-1]]
                continue

            elif estimated_h == 1.0:                
                cis[i, 0] = true_h2s[where(rand_value < l)[0][0]]
                cis[i, 1] = 1.0
                continue   

        #print estimated_h, ends, bisect.bisect_left(ends, estimated_h)
        idx = bisect.bisect_left(ends, estimated_h)
        if idx == 0:
            cis[i,0] = 0.0
        else:
            proportion = (estimated_h - ends[idx-1]) / (ends[idx] - ends[idx-1])
            cis[i,0] = true_h2s[idx-1] + (true_h2s[idx] - true_h2s[idx-1]) * proportion 
        
        idx = bisect.bisect_right(starts, estimated_h)
        if idx == len(true_h2s):
            cis[i,1] = 1.0
        else:
            proportion = (estimated_h - starts[idx-1]) / (starts[idx] - starts[idx-1])
            cis[i,1] = true_h2s[idx-1] + (true_h2s[idx] - true_h2s[idx-1]) * proportion 
    return cis#, u, l
       


def example():
    pass