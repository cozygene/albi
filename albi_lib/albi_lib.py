from numpy import * # for: newaxis, isclose, seterr, random, array, maximum, reshape, ones
import numpy.linalg
import bisect


# Ignore divide by 0 warnings
seterr(divide='ignore', invalid='ignore')

class DerivativeSignCalculator(object):
    """
    An abstract class to calculate the sign of the derivative of the likelihood
    under various settings.
    """
    def get_derivative_signs(us):
        raise NotImplementedError("Should be overrided")

class OnlyEigenvectorsDerivativeSignCalculator(DerivativeSignCalculator):
    def __init__(self, h2_values, H2_values, kinship_eigenvalues, eigenvectors_as_X=[-1], REML=True):
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

        Calculates the weights -
            A matrix of size N x M x K, of the weights for the cartesian product of h2_values, H2_values and an index.
        """
        # Set the shapes accordingly, so vectorization will work
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

        # Calculate weights
        projection = ones((1, 1, n_samples))
        if len(eigenvectors_as_X):
            projection[0, 0, eigenvectors_as_X] = 0


        if REML:
            ds = projection * (kinship_eigenvalues - 1) / (H2_values * (kinship_eigenvalues-1) + 1)                
            denom = n_samples - len(eigenvectors_as_X)
        else:
            ds = (kinship_eigenvalues - 1) / (H2_values * (kinship_eigenvalues - 1) + 1)
            denom = n_samples


        self.weights = projection * (h2_values * (kinship_eigenvalues - 1) + 1) / \
                                    (H2_values * (kinship_eigenvalues - 1) + 1) \
                                     * (ds - nansum(ds, 2)[:, :, newaxis] / denom)

        self.weights = nan_to_num(self.weights)


    def get_derivative_signs(self, us):
        """
        Get the signs of the derivative for the supplied vectors.

        Arguments:
            us - a matrix of size T X N X K, of T X N vectors (T vectors per N of all possible values of h^2)

        Returns:
            A matrix of size T X N X M of the derivative signs of each vector at the M given H2 points
        """
        assert shape(us)[1] == shape(self.weights)[0]
        assert shape(us)[2] == shape(self.weights)[2]

        us = us[:,:,newaxis,:]
        dotproducts = sum((us**2) * self.weights[newaxis, :, :, :], 3)   

        return sign(dotproducts)


# Retain for backward compatibility
def weights_zero_derivative(h2_values, H2_values, kinship_eigenvalues, 
                            eigenvectors_as_X=[-1], REML=True):
    return OnlyEigenvectorsDerivativeSignCalculator(h2_values, H2_values, kinship_eigenvalues, eigenvectors_as_X=[-1], REML=True).weights
    

class GeneralDerivativeSignCalculator(DerivativeSignCalculator):
    def __init__(self, h2_values, H2_values, kinship_eigenvalues, kinship_eigenvectors, covariates, REML=True):
        """
        Calculate quantities to be used later when calculating derivative signs.

        Arguments:
            h2_values - a vector of size N, of all possible values of h^2
            H2_values - a vector of size M, of a grid of possible values of H^2
            kinship_eigenvalues - A vector of size K of the eigenvalues of the kinship matrix, in decreasing order.
            kinship_eigenvectors - A matrix of size K x K whose columns are the eigenvectors of the kinship matrix, corresponding to the given eigenvalues.
            covariates - a matrix of size K x P, of P covariates to be used.
            REML - True is REML, False if ML.
        """
        n_samples = len(kinship_eigenvalues)
        n_intervals = len(H2_values)-1
        n_covariates = shape(covariates)[1]

        self.h2_values = array(h2_values)
        self.H2_values = array(H2_values)
        
        # Make sure the eigenvalues are in decreasing order and nonzero
        kinship_eigenvalues = maximum(kinship_eigenvalues, 1e-10)
        reverse_sorted_indices = argsort(kinship_eigenvalues)[::-1]
        self.kinship_eigenvalues = array(kinship_eigenvalues)[reverse_sorted_indices]
        self.kinship_eigenvectors = array(kinship_eigenvectors)[:, reverse_sorted_indices]

        rotated_X = dot(self.kinship_eigenvectors.T, covariates)

        # Size: M X K X P
        self.X_XtdXi_matrices = zeros((len(self.H2_values), n_samples, n_covariates))
        for iH, H2 in enumerate(self.H2_values): 
            self.X_XtdXi_matrices[iH, :, :] = dot(rotated_X, linalg.inv(dot(rotated_X.T * (1.0/(H2*(self.kinship_eigenvalues-1) + 1)), rotated_X)))

        # Size: M X P X K
        self.Xt_Dv_matrices = zeros((len(self.H2_values), n_covariates, n_samples))
        for iH, H2 in enumerate(self.H2_values): 
            self.Xt_Dv_matrices[iH, :, :] = rotated_X.T * (1.0/(H2*(self.kinship_eigenvalues-1) + 1))

        first_logdets = zeros(len(self.H2_values))
        second_logdets = zeros(len(self.H2_values))
        self.B_multipliers = zeros(len(self.H2_values))
        for iH, H2 in enumerate(self.H2_values): 
            first_logdet = nansum((self.kinship_eigenvalues-1) / (H2*(self.kinship_eigenvalues-1) + 1))
            second_logdet = nansum(linalg.inv(dot(rotated_X.T * (1.0/(H2*(self.kinship_eigenvalues-1) + 1)), rotated_X))  *
                                     dot(rotated_X.T * (-(self.kinship_eigenvalues-1)/((H2*(self.kinship_eigenvalues-1) + 1)**2)), rotated_X))
            if REML:
                self.B_multipliers[iH] = (first_logdet + second_logdet) * (1.0/(n_samples - n_covariates))
            else:
                self.B_multipliers[iH] = (first_logdet) * (1.0/n_samples)        

    def get_derivative_signs(self, us):
        """
        Get the signs of the derivative for the supplied vectors.

        Arguments:
            us - a matrix of size N X K, of N vectors, one per value of h2

        Returns:
            A matrix of size N X M of the derivative signs of each vector at the M given H2 points
        """
        n_samples = len(self.kinship_eigenvalues)

        assert shape(us)[0] == len(self.h2_values)
        assert shape(us)[1] == n_samples

        vs = ((self.h2_values[:, newaxis] * (self.kinship_eigenvalues[newaxis, :] - 1) + 1)**0.5) * us

        # Size N X M X P
        w1s = tensordot(vs, self.Xt_Dv_matrices, axes=([1, 2]))

        # Size N X M X K
        w2s = zeros([len(self.h2_values), len(self.H2_values), n_samples])
        for ih, h2 in enumerate(self.h2_values):
            w2s[ih, :, :] = sum(self.X_XtdXi_matrices * w1s[ih, :, newaxis, :], 2)

        # Size N X M X K
        wss = vs[:, newaxis, :] - w2s
        Ass = sum(((self.kinship_eigenvalues[newaxis, newaxis, :] - 1)) / ((self.H2_values[newaxis, :, newaxis] * (self.kinship_eigenvalues[newaxis, newaxis, :] - 1) + 1)**2) * (wss**2), 2)
        Bss = sum(1.0/(self.H2_values[newaxis, :, newaxis] * (self.kinship_eigenvalues[newaxis, newaxis, :] - 1) + 1) * (wss**2), 2)

        dotproducts = Ass - self.B_multipliers[newaxis, :] * Bss

        return sign(dotproducts)

def estimate_distributions_eigenvectors(h2_values, H2_values, kinship_eigenvalues, 
                                        n_random_samples=100, eigenvectors_as_X=[-1], REML=True, seed=0):
    """
    Across a grid of possible estimated values H^2, approximately calculate the probability of either evaluating a boundary 
    estimate (for the boundaries of the grid) or the the estimate falling between each grid points. The probability is 
    defined for a true value of heritability h^2. The probability is estimated with a parametric bootstrap.

    Limited to cases where the covariates are eigenvectors of the kinship matrix; but in these cases, calculation is faster.

    Arguments:
        h2_values - a vector of size N, of all possible values of h^2
        H2_values - a vector of size M, of a grid of possible values of H^2
        kinship_eigenvalues - A vector of size K of the eigenvalues of the kinship matrix, in decreasing order.
        n_random_samples - The number of random samples to use for the parameteric bootstrap. Can be an int or an iterable with the __len__ func implemented
        eigenvectors_as_X - A list of indices, of which eigenvectors of the kinship matrix are fixed effects.
        REML - True is REML, False if ML.
        seed - A seed for the random generator used for the random samples.

    Returns:
        A matrix of size N x (M + 1) of the probabilities, where:
        - The cell at index (i, 0) is the probability of the estimate being smaller than the smallest grid point;
        - The cell at index (i, M) is the probability of the estimate being larger than the largest grid point;
        - The cell at index (i, j), for 0<j<M, is the probability of estimate being between the (j-1)-th and j-th grid points;

          all of the above are for the i-th h2 value.
    """
    monte_carlo_size = 1
    n_samples = len(kinship_eigenvalues)
    n_intervals = len(H2_values)-1

    # Make sure the eigenvalues are in decreasing order and nonzero
    kinship_eigenvalues = array(sorted(maximum(kinship_eigenvalues, 1e-10)))[::-1]

    # Size: N X M X K
    sign_calculator = OnlyEigenvectorsDerivativeSignCalculator(h2_values, H2_values, kinship_eigenvalues, eigenvectors_as_X=eigenvectors_as_X, REML=REML)
    
    rng = random.RandomState(seed)
    prob = zeros((len(h2_values), n_intervals+2))

    if not hasattr(n_random_samples, '__iter__'):  # assumes that n_random_samples is int
        n_random_samples = range(n_random_samples)

    for i in n_random_samples:
        # Avoid replicating weights across h2's, so that each sample is independent
        us = rng.normal(size=(monte_carlo_size, len(h2_values), n_samples))
        dotproducts = sign_calculator.get_derivative_signs(us)
        
        # Size: monte_carlo_size X N
        hit_zero = (dotproducts[:, :, 0] <= 0)
        hit_one = (dotproducts[:, :, -1] >= 0)
        hit_boundary = hit_zero | hit_one

        prob[:, 0] += mean(hit_zero, 0)        
        prob[:, -1] += mean(hit_one, 0)        
        prob[:, 1:-1] += mean(((dotproducts[:, :, :-1] >= 0) & (dotproducts[:, :, 1:] <= 0)) & ~hit_boundary[:, :, newaxis], 0)
        
    prob /= len(n_random_samples)       # Average across chunks
    prob /= sum(prob, 1)[:, newaxis]   
    return prob

def estimate_distributions_general(h2_values, H2_values, kinship_eigenvalues, kinship_eigenvectors, covariates,
                                   n_random_samples=100, REML=True, seed=0):
    """
    Across a grid of possible estimated values H^2, approximately calculate the probability of either evaluating a boundary 
    estimate (for the boundaries of the grid) or the the estimate falling between each grid points. The probability is 
    defined for a true value of heritability h^2. The probability is estimated with a parametric bootstrap.

    Arguments:
        h2_values - a vector of size N, of all possible values of h^2
        H2_values - a vector of size M, of a grid of possible values of H^2
        kinship_eigenvalues - A vector of size K of the eigenvalues of the kinship matrix, in decreasing order.
        kinship_eigenvectors - A matrix of size K x K whose columns are the eigenvectors of the kinship matrix, corresponding to the given eigenvalues.
        covariates - a matrix of size K x P, of P covariates to be used.
        n_random_samples - The number of random samples to use for the parameteric bootstrap. Can be an int or an iterable with the __len__ func implemented
        REML - True is REML, False if ML.
        seed - A seed for the random generator used for the random samples.

    Returns:
        A matrix of size N x (M + 1) of the probabilities, where:
        - The cell at index (i, 0) is the probability of the estimate being smaller than the smallest grid point;
        - The cell at index (i, M) is the probability of the estimate being larger than the largest grid point;
        - The cell at index (i, j), for 0<j<M, is the probability of estimate being between the (j-1)-th and j-th grid points;

          all of the above are for the i-th h2 value.
    """
    n_samples = len(kinship_eigenvalues)
    n_intervals = len(H2_values)-1
    n_covariates = shape(covariates)[1]

    rng = random.RandomState(seed)

    h2_values = array(h2_values)
    H2_values = array(H2_values)
   
    prob = zeros((len(h2_values), n_intervals+2))

    if not hasattr(n_random_samples, '__iter__'):  # assumes that n_random_samples is int
        n_random_samples = range(n_random_samples)

    sign_calculator = GeneralDerivativeSignCalculator(h2_values, H2_values, kinship_eigenvalues, kinship_eigenvectors, covariates, REML=REML)

    # Main loop
    for i in n_random_samples:        
        us = rng.normal(size=(len(h2_values), n_samples))         
        dotproducts = sign_calculator.get_derivative_signs(us)

        hit_zero = (dotproducts[:, 0] <= 0)
        hit_one = (dotproducts[:, -1] >= 0)
        hit_boundary = hit_zero | hit_one

        prob[:, 0] += hit_zero        
        prob[:, -1] += hit_one
        prob[:, 1:-1] += ((dotproducts[:, :-1] >= 0) & (dotproducts[:, 1:] <= 0)) & ~hit_boundary[:, newaxis]
        
    prob /= len(n_random_samples)       # Average across chunks
    prob /= sum(prob, 1)[:, newaxis]   

    return prob

def quantiles(h2_values, cdf, beta):
    """
    Calculate the beta-th quantile of the given cumulative distribution function (cdf) of \hat{h}^2. Namely,
    find the value q satisfying:

    q = inf{h^2 | beta <= Pr(\hat{h}^2 <= h^2)}

    Since there are discontinuity points at 0,1, this translates to the following function:

    - If beta <= Pr(\hat{h}^2 = 0), q = 0
    - If beta >= 1 - Pr(\hat{h}^2 = 1), q = 1
    - Other, it is the inverse of the CDF at the range between these points

    Arguments:
        h2_values - a vector of size N, of the values of h^2 where the given cdf is evaulated. Assumed to 
                    start at 0 and end at 1.
        cdf - a vector of size N+1, The CDF of \hat{h}^2, evaluated at the points of h2_values. 
              To account for the discontinuity points, it is required that:
                 - cdf[0] = Pr(\hat{h}^2 = 0)
                 - cdf[i] = Pr(\hat{h}^2 < h2_values[i-1]), for i = 1..N-1
                 - cdf[N] = Pr(\hat{h}^2 <= 1) == 1

              It follows that cdf[-2] == 1 - Pr(\hat{h}^2 = 1).

        beta - The required quantile.

    Returns:
        q - The required h^2 value corresponding to the beta-th quantile
    """
    assert h2_values[0] == 0 and h2_values[-1] == 1

    # max{h^2 | P(\hat{h}^2 <= h^2) == p_0}
    L = bisect.bisect_right(cdf, cdf[0])-1

    # min{h^2 | P(\hat{h}^2 <= h^2) == 1-p_1}
    R = bisect.bisect_left(cdf, cdf[-2])

    # Now, the open interval (h2_values[L], h2_values[R]) maps 1-1 to the open interval (p0, 1-p1)
    if (cdf[0] < beta < cdf[-2]):
        subrange = h2_values[L:R+1]
        subcdf = cdf[L:R+1]
        idx = bisect.bisect_left(subcdf, beta)
        if subcdf[idx] == beta: # If we hit exactly right, simplify things and just return that
            return subrange[idx]
        else:
            # Otherwise, the required value is explicity inside the interval (subcdf[idx-1], subcdf[idx]),
            # so to get the right value, interpolate linearly
            proportion = (beta - subcdf[idx-1]) / (subcdf[idx] - subcdf[idx-1])
            value = subrange[idx-1] + proportion * (subrange[idx] - subrange[idx-1])
            return value

    # If it's not in that open interval, this either beta <= p0 or beta >= 1-p1
    elif (beta <= cdf[0]) or isclose(beta, cdf[0]):
        return 0.0        
    elif (beta >= cdf[-2]):
        return 1.0

    else:
        assert "Should not happen"

def get_approx_cdf(h2_values, cdf, required_h2_value):
    """
    Find Pr(\hat{h}^2 <= required_h2_value).

    Arguments:
        h2_values - a vector of size N, of the values of h^2 where the given cdf is evaulated. Assumed to 
                    start at 0 and end at 1.
        cdf - a vector of size N+1, The CDF of \hat{h}^2, evaluated at the points of h2_values. 
              To account for the discontinuity points, it is required that:
                 - cdf[0] = Pr(\hat{h}^2 = 0)
                 - cdf[i] = Pr(\hat{h}^2 < h2_values[i-1]), for i = 1..N-1
                 - cdf[N] = Pr(\hat{h}^2 <= 1) == 1

              It follows that cdf[-2] == 1 - Pr(\hat{h}^2 = 1).

        required_h2_value - The value at which to evaulate the cdf.

    Returns:
        CDF at point required_h2_value.
    """
    assert 0 <= required_h2_value <= 1
    if required_h2_value == 0.0:
        return cdf[0]
    elif required_h2_value == 1:
        return 1.0
    else:
        idx = bisect.bisect_left(h2_values, required_h2_value)
        # If we hit exactly right, simplify things and just return that
        if h2_values[idx] == required_h2_value:
            return cdf[idx]
        # Otherwise, the required value is explicity inside the interval (cdf[idx-1], cdf[idx]),
        # so to get the right value, interpolate linearly
        else:        
            proportion = (required_h2_value - h2_values[idx-1]) / (h2_values[idx] - h2_values[idx-1])
            return cdf[idx-1] + (cdf[idx] - cdf[idx-1]) * proportion 

def get_probabilty_in_closed_interval(h2_values, cdf, interval):
    """
    Find Pr(interval[0] <= \hat{h}^2 <= interval[1]).

    Arguments:
        h2_values - a vector of size N, of the values of h^2 where the given cdf is evaulated. Assumed to 
                    start at 0 and end at 1.
        cdf - a vector of size N+1, The CDF of \hat{h}^2, evaluated at the points of h2_values. 
              To account for the discontinuity points, it is required that:
                 - cdf[0] = Pr(\hat{h}^2 = 0)
                 - cdf[i] = Pr(\hat{h}^2 < h2_values[i-1]), for i = 1..N-1
                 - cdf[N] = Pr(\hat{h}^2 <= 1) == 1

              It follows that cdf[-2] == 1 - Pr(\hat{h}^2 = 1).

        interval - The interval for which the probability should be evaluated.

    Returns:
        The probability of the interval.
    """
    p = get_approx_cdf(h2_values, cdf, interval[1]) - get_approx_cdf(h2_values, cdf, interval[0])
    if interval[0] == 0.0:
        p += cdf[0]
    return p

def build_heritability_cis(h2_values, H2_values, all_distributions, estimated_values, 
                           confidence=0.95, use_randomized_cis=False, seed=0): 
    """
    Build confidence intervals for a set of estimated values, given estimator distributions.

    Arguments:
        h2_values - a vector of size N, of values of true h^2 values at which the estimator
                    distribution is given.
        H2_values - a vector of size M, of a grid of values of estimated values at which 
                    the estimator distribution is given, for each true value h^2.
        all_distributions - a matrix of size N x (M + 1), where the i-th row is the estimator
                            distribution of the i-th h^2 value (the output of estimate_distributions)        
        estimated_values - A vector of size P, of estimated heritability values, for which we wish to 
                           calculate confidence intervals.
        confidence - The required confidence level of the CIs we wish to construct (e.g., 95%).
        seed - A seed for the random generator used for the randomized CIs, if needed.
        use_randomized_cis - Should we use randomized CIs, if needed.

    Returns:
        A matrix of size (P X 2), of the confidence intervals for each required estimate value.
    """                          
    alpha = 1-confidence

    # Calculate acceptance regions
    accept_regions = zeros((len(h2_values), 3, 3))
    for i,h2 in enumerate(h2_values):
        cdf = cumsum(all_distributions[i,:])        
        
        accept_regions[i,0,0] = 0.0
        accept_regions[i,0,1] = quantiles(H2_values, cdf, 1-alpha)

        accept_regions[i,1,0] = quantiles(H2_values, cdf, alpha/2)
        accept_regions[i,1,1] = quantiles(H2_values, cdf, 1-alpha/2)

        accept_regions[i,2,0] = quantiles(H2_values, cdf, alpha)
        accept_regions[i,2,1] = 1.0
    
    # Regularize for noise:
    # 1. Make sure h^2 is included 
    accept_regions[:,0,1] = maximum(accept_regions[:,0,1], h2_values)
    accept_regions[:,1,0] = minimum(accept_regions[:,1,0], h2_values)
    accept_regions[:,1,1] = maximum(accept_regions[:,1,1], h2_values)
    accept_regions[:,2,0] = minimum(accept_regions[:,2,0], h2_values)

    # 2. Make sure type II boundaries are over type I and III where relevant
    accept_regions[:,1,1] = maximum(accept_regions[:,1,1], accept_regions[:,0,1])
    accept_regions[:,1,0] = minimum(accept_regions[:,1,0], accept_regions[:,2,0])

    # 3. Make sure they are monotone
    accept_regions[:,0,1] = maximum.accumulate(accept_regions[:,0,1])
    accept_regions[:,1,0] = maximum.accumulate(accept_regions[:,1,0])
    accept_regions[:,1,1] = maximum.accumulate(accept_regions[:,1,1])
    accept_regions[:,2,0] = maximum.accumulate(accept_regions[:,2,0])

    # Calculate coverages
    for i,h2 in enumerate(h2_values):
        for j in range(3):
            cdf = cumsum(all_distributions[i,:])  
            accept_regions[i,j,2] = get_probabilty_in_closed_interval(H2_values, cdf, accept_regions[i,j,:2])

    # Find the first place type I covers 1
    s = where(isclose(accept_regions[:,0,2], 1))[0][0]

    # Find the last place type III covers 1
    t = where(isclose(accept_regions[:,2,2], 1))[0][-1]
    d = int((s+t)/2)

    # Find the regions where type II covers exactly 
    w = where(isclose(accept_regions[:,1,2], confidence))[0]

    # If type II covers a range of h^2, then by definition type I covers before and type III after
    if len(w):
        regions = vstack([accept_regions[:w[0],      0, :], 
                          accept_regions[w[0]:w[-1], 1, :],
                          accept_regions[w[-1]:    , 2, :]])
    else:
        # Otherwise, either type I and III cover everything (and we're good)
        # or we resort to randomized CIs. In both cases we do this:
        regions = vstack([accept_regions[:d, 0, :],                               
                          accept_regions[d:, 2, :]])   
    
    starts = maximum.accumulate(regions[:,0])
    ends = minimum.accumulate(regions[:,1][::-1])[::-1]

    rng = random.RandomState(seed)

    # This is relevant only is w is empty, s < t, and use_randomized_cis=True
    u = 1 - (1-confidence) / all_distributions[:, 0]
    u[:d] = 1
    u[t:] = 0
    u = minimum.accumulate(maximum(0, u))  # denoise and make monotone

    l = 1 - (1-confidence) / all_distributions[:, -1]
    l[:s] = 0
    l[d:] = 1
    l = maximum.accumulate(maximum(0, l))  # denoise and make monotone

    # Calculate confidence intervals
    cis = zeros((len(estimated_values), 2))
    for i, estimated_h in enumerate(estimated_values):
        if use_randomized_cis and s < t:
            rand_value = rng.rand()
            if estimated_h == 0.0:                
                cis[i, 0] = 0.0                
                cis[i, 1] = h2_values[where(rand_value < u)[0][-1]]
                continue

            elif estimated_h == 1.0:                
                cis[i, 0] = h2_values[where(rand_value < l)[0][0]]
                cis[i, 1] = 1.0
                continue   

        idx = bisect.bisect_left(ends, estimated_h)
        if idx == 0:
            cis[i,0] = 0.0
        else:
            proportion = (estimated_h - ends[idx-1]) / (ends[idx] - ends[idx-1])
            cis[i,0] = h2_values[idx-1] + (h2_values[idx] - h2_values[idx-1]) * proportion 
        
        idx = bisect.bisect_right(starts, estimated_h)
        if idx == len(h2_values):
            cis[i,1] = 1.0
        else:
            proportion = (estimated_h - starts[idx-1]) / (starts[idx] - starts[idx-1])
            cis[i,1] = h2_values[idx-1] + (h2_values[idx] - h2_values[idx-1]) * proportion 

    return cis

