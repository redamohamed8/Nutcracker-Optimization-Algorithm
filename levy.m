%_________________________________________________________________________%
%  Nutcracker Optimization Algorithm (NOA) source codes demo 1.0               %
%                                                                         %
%  Developed in MATLAB R2019A                                      %
%                                                                         %
%  Author and programmer: Reda Mohamed (E-mail: redamoh@zu.edu.eg) & Mohamed Abdel-Basset (E-mail: mohamedbasset@ieee.org)                              %
%                                                                         %
%   Main paper: Abdel-Basset, M., Mohamed, R.                                    %
%               Nutcracker optimizer,                         %
%               Knowledge-Based Systems, in press,              %
%               DOI: https://doi.org/10.1016/j.knosys.2022.110248   %
%                                                                         %
%_________________________________________________________________________%

function [z] = levy(n,m,beta)

    num = gamma(1+beta)*sin(pi*beta/2); % used for Numerator 
    
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2); % used for Denominator

    sigma_u = (num/den)^(1/beta);% Standard deviation

    u = random('Normal',0,sigma_u,n,m); 
    
    v = random('Normal',0,1,n,m);

    z =u./(abs(v).^(1/beta));

  
  end

