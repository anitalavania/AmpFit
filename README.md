# AmpFit

This is a package that can be used to perform any three-body amplitude fit. I developed it as part of one of my PhD projects. 

## Optimization details
- TMinuit package from ROOT is used to implement Migrad optimization algorithm.
- Migrad calculates second derivative Hessian matrix and used inverse of it as dynamic learning rate. This Hessian matrix can be understood as curvature of the slope.

x_new = x_old - H⁻¹ * ∇f(x_old)


## Particle physics details
- Breit-Wigner, K-matrix, LASS and Flatte formulations are available for resonance descriptions.





(Go ahead, open the black box! :D)
