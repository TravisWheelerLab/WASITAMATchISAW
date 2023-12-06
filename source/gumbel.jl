# thank you to https://djalil.chafai.net/blog/2018/02/18/gumbel-fit-with-julia/

# Pkg.add("Roots") # https://github.com/JuliaMath/Roots.jl
using Roots
""" 
Computes approximately the Maximum Likelihood Estimator (MLE) of 
the position and scale parameters of a Gumbel distribution with
cumulative distribution function x-> exp(-exp(-(x-mu)/sigma)).
Has mean mu+eulergamma*sigma and variance sigma^2*pi^2/6.
The MLE for mu is a function of data and the MLE for sigma.
The MLE for sigma is the root of a nonlinear function of data,
computed numerically with function fzero() from package Roots.

Reference: page 24 of book ISBN 1-86094-224-5.
Samuel Kotz and Saralees Nadarajah
Extreme value distributions.Theory and applications. 
Imperial College Press, London, 2000. viii+187 pp. 
https://ams.org/mathscinet-getitem?mr=1892574 

``no distribution should be stated without an explanation of how the
parameters are estimated even at the risk that the methods used will 
not stand up to the present rigorous requirements of mathematically 
minded statisticians''. E. J. Gumbel (1958).
"""
function gumbel_fit(data)
  f(x) = x - mean(data) + sum(data.*exp.(-data/x)) / sum(exp.(-data/x)) 
  sigma = fzero(f,sqrt(6)*std(data)/pi) # moment estimator as initial value
  mu = - sigma * log(sum(exp.(-data/sigma))/length(data))
  return mu , sigma
end #function

import Distributions: fit
function fit(::Type{Gumbel}, data)
    position, scale = gumbel_fit(data)
    Gumbel(position, scale)
end