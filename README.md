# AMS_shape_orient_matlab
This repository contains Matlab code that calculates snow aggregate distribution functions, bulk fall speed quantities, and self-collection quantities used in the Journal of the Atmospheric Science (JAS) article entitled: "How Snow Aggregate Shape and Orientation Variability Affects Fall Speed and Self-Collection Rates." The relevant mathematical details that describe the H-function distribution function forms and bulk fall speed quantities and self-collection rate equations can be found in the Supplementary Materials section of that article. Additional Matlab code is provided for calculating Fox's H-function through numerical integration of the H-function Mellin-Barnes integral representation (see Supplementary Materials section for more details). This function (Fox_H.m) comes from the article: "Performance Analysis of M-QAM Multihop Relaying over mmWave Weibull Fading Channels", Soulimani et al. 2016" but is modified slightly to incorporate array inputs. Bulk fall speed quantities (i.e., number-, mass-, and reflectivity-weighted fall speeds) can be computed exactly using the integration properties of the H-function; that is, integral quantities are explicitly described in terms of gamma function ratios. Self-collection rates are estimated by sampling from the individual probability density functions for size and shape. H-function distribution functions are computed for variables of the form: xi = alpha*a^zeta_a * phi_ba^zeta_ba * phi_ca^zeta_ca, where xi is a power-law transformed variable of 'a' (aggregate ellipsoid size), phi_ba and phi_ca (ellipsoid aspect ratios: phi_ba = b/a and phi_ca = c/a, where a>=b>=c). These new distributions are computed from both the H-function definition (solid lines) and from sampling each individual distribution of a~Gamma(nu,a_n) and (phi_ba,phi_ca) ~ BivBeta(a_ba,b_ba,b_cb).


Additional Matlab files can be found on the Mathworks file exchange website:

File 1: pinky.m
Description: Used to discretize and sample from bivariate beta distribution function. 
Purpose: This allows for a quick and easy way to simulate estimates of self-collection rates and H-function distributions.
https://www.mathworks.com/matlabcentral/fileexchange/35797-generate-random-numbers-from-a-2d-discrete-distribution

File 2: gammaz.m
Description: the gamma function for complex inputs. Notice that this function on the file exchange was renamed to "gamma.m". However, the Fox_H.m file renames this file as "gammaz.m" in order to preserve the stanard matlab function: gamma.m. 
Purpose: Allows for the Fox_H.m to numerically integrate in complex s plane.
https://www.mathworks.com/matlabcentral/fileexchange/3572-gamma
