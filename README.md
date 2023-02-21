# maximum-likelihood-estimation

<p>

<img src="https://drive.google.com/uc?id=1QDiXNZgvR6hAlW0uuG8plWsbgWwEl2xG" width="30%" />

<img src="https://drive.google.com/uc?id=1pimWlSOK-dbf2LvJ2AiCdvDlapkyJKYg" width="30%" />
  
<img src="https://drive.google.com/uc?id=1dkPAkVpjOuUHCwT6J_5pKsjT5Wird0Fg" width="30%" />

</p>

The basic problem is to determine the best failure distribution implied by n failure times comprised in a sample. In all cases the sample is assumed to be a simple random (or probability) sample where the failure or repaire times are independent observations from a common population.

There are two general approaches to fitting reliability distributions to failure data. The first, and usually preferred, method is to fit a theoretical distribution, such as exponential, Weibull, normal or lognormal distribution. The second is to derive, directly from the collected data, an empirical reliability function or hazard rate function.

Once one or more distributions have been identified, the next step is to estimate the parameters of the distribution.

Maximum likelihood estimation (MLE) is a method of estimating the parameters of an assumed probability distribution. This is achieved by maximizing a likelihood function so that the observed data is most probable. The point in the parameter space that maximizes the likelihood function is called the maximum likelihood estimate.

The MLE for the two-parameter Weibull distribution must be computed numerically. The Newton-Raphson method for solving a nonlinear equation is used.
