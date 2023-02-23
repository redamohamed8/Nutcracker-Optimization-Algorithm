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

clear all
clc
SearchAgents_no=25; % Number of search agents
Max_iteration=50000; % Maximum number of Function evaluations
RUN_NO=30; %% Number of independent runs

for i=1:23 %% Test functions

  for j=1:RUN_NO
   [lb,ub,dim,fobj]=Get_Functions_details(i);
   [Best_score,Best_pos,Convergence_curve]=NOA(SearchAgents_no,Max_iteration,ub,lb,dim,fobj);
   fitness(1,j)=Best_score;
  end
  fprintf(['Function_ID\t',num2str(i),'\tAverage Fitness:',num2str(mean(fitness(1,:)),20),'\n']);

  %% Drawing Convergence Curve %%
  figure(i)
  h=semilogy(Convergence_curve,'-<','MarkerSize',8,'LineWidth',2,'Color','r');
  h.MarkerIndices = 1000:1000:Max_iteration;
  xlabel('Iteration');
  ylabel('Best Fitness obtained so-far');
  axis tight
  grid off
  box on
  legend({'NOA'});
end