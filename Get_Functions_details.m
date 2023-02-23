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
% lb is the lower bound: lb=[lb_1,lb_2,...,lb_d]
% up is the uppper bound: ub=[ub_1,ub_2,...,ub_d]
% dim is the number of variables (dimension of the problem)

function [lb,ub,dim,fobj] = Get_Functions_details(F)


switch F
    case 1
        fobj = @F1;
        dim=100;
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        
    case 2
        fobj = @F2;
                dim=100;

        lb=-10*ones(1,dim);
        ub=10*ones(1,dim);
        
    case 3
        fobj = @F3;
         dim=100;
        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
       
        
    case 4
        fobj = @F4;
                dim=100;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
        
    case 5
        fobj = @F5;
        dim=100;
        lb=-30*ones(1,dim);;
        ub=30*ones(1,dim);;
        
        
    case 6
        fobj = @F6;
        dim=100;
        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
        
    case 7
        fobj = @F7;
                dim=100;

        lb=-1.28*ones(1,dim);;
        ub=1.28*ones(1,dim);;
        
    case 8
        fobj = @F8;
                dim=100;
        lb=-500*ones(1,dim);;
        ub=500*ones(1,dim);;
        
    case 9
        fobj = @F9;
        dim=100;

        lb=-5.12*ones(1,dim);;
        ub=5.12*ones(1,dim);;
        
    case 10
        fobj = @F10;
          dim=100;
        lb=-32*ones(1,dim);;
        ub=32*ones(1,dim);;
        
    case 11
        fobj = @F11;
         dim=100;
        lb=-600*ones(1,dim);;
        ub=600*ones(1,dim);;
        
    case 12
        fobj = @F12;
        dim=100;
        lb=-50*ones(1,dim);;
        ub=50*ones(1,dim);;
        
    case 13
        fobj = @F13;
        dim=100;
        lb=-50*ones(1,dim);;
        ub=50*ones(1,dim);;
        
    case 14
        fobj = @F14;
        dim=2;
        lb=-65.536*ones(1,dim);;
        ub=65.536*ones(1,dim);;
        
    case 15
        fobj = @F15;
        dim=4;
        lb=-5*ones(1,dim);;
        ub=5*ones(1,dim);;
        
    case 16
        fobj = @F16;
        dim=2;
        lb=-5*ones(1,dim);;
        ub=5*ones(1,dim);;
        
    case 17
        fobj = @F17;
       dim=2;
        lb=[-5,0];
        ub=[10,15];
        
    case 18
        fobj = @F18;
dim=2;
        lb=-2*ones(1,dim);;
        ub=2*ones(1,dim);;
        
    case 19
        fobj = @F19;
        dim=3;
        lb=0*ones(1,dim);;
        ub=1*ones(1,dim);;
        
        
    case 20
        fobj = @F20;
        dim=6;
        lb=0*ones(1,dim);;
        ub=1*ones(1,dim);;
             
        
    case 21
        fobj = @F21;
        dim=4;
        lb=0*ones(1,dim);;
        ub=10*ones(1,dim);;
            
        
    case 22
        fobj = @F22;
        dim=4;
        lb=0*ones(1,dim);;
        ub=10*ones(1,dim);;
            
        
    case 23
        fobj = @F23;
         dim=4;   
        lb=0*ones(1,dim);
        ub=10*ones(1,dim);
    end
end

% F1

function o = F1(x)
o=sum(x.^2);
end

% F2

function o = F2(x)
o=sum(abs(x))+prod(abs(x));
end

% F3

function o = F3(x)
dim=size(x,2);
o=0;
for i=1:dim
    o=o+sum(x(1:i))^2;
end
end

% F4

function o = F4(x)
o=max(abs(x));
end

% F5

function o = F5(x)
dim=size(x,2);
o=sum(100*(x(2:dim)-(x(1:dim-1).^2)).^2+(x(1:dim-1)-1).^2);
end

% F6

function o = F6(x)
o=sum(abs((x+.5)).^2);
end

% F7

function o = F7(x)
dim=size(x,2);
o=sum([1:dim].*(x.^4))+rand;
end

% F8

function o = F8(x)
o=sum(-x.*sin(sqrt(abs(x))));
end

% F9

function o = F9(x)
dim=size(x,2);
o=sum(x.^2-10*cos(2*pi.*x))+10*dim;
end

% F10

function o = F10(x)
dim=size(x,2);
o=-20*exp(-.2*sqrt(sum(x.^2)/dim))-exp(sum(cos(2*pi.*x))/dim)+20+exp(1);
end

% F11

function o = F11(x)
dim=size(x,2);
o=sum(x.^2)/4000-prod(cos(x./sqrt([1:dim])))+1;
end

% F12

function o = F12(x)
dim=size(x,2);
o=(pi/dim)*(10*((sin(pi*(1+(x(1)+1)/4)))^2)+sum((((x(1:dim-1)+1)./4).^2).*...
(1+10.*((sin(pi.*(1+(x(2:dim)+1)./4)))).^2))+((x(dim)+1)/4)^2)+sum(Ufun(x,10,100,4));
end

% F13

function o = F13(x)
dim=size(x,2);
o=.1*((sin(3*pi*x(1)))^2+sum((x(1:dim-1)-1).^2.*(1+(sin(3.*pi.*x(2:dim))).^2))+...
((x(dim)-1)^2)*(1+(sin(2*pi*x(dim)))^2))+sum(Ufun(x,5,100,4));
end

% F14

function o = F14(x)
aS=[-32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32;,...
-32 -32 -32 -32 -32 -16 -16 -16 -16 -16 0 0 0 0 0 16 16 16 16 16 32 32 32 32 32];

for j=1:25
    bS(j)=sum((x'-aS(:,j)).^6);
end
o=(1/500+sum(1./([1:25]+bS))).^(-1);
end

% F15

function o = F15(x)
aK=[.1957 .1947 .1735 .16 .0844 .0627 .0456 .0342 .0323 .0235 .0246];
bK=[.25 .5 1 2 4 6 8 10 12 14 16];bK=1./bK;
o=sum((aK-((x(1).*(bK.^2+x(2).*bK))./(bK.^2+x(3).*bK+x(4)))).^2);
end

% F16

function o = F16(x)
o=4*(x(1)^2)-2.1*(x(1)^4)+(x(1)^6)/3+x(1)*x(2)-4*(x(2)^2)+4*(x(2)^4);
end

% F17

function o = F17(x)
o=(x(2)-(x(1)^2)*5.1/(4*(pi^2))+5/pi*x(1)-6)^2+10*(1-1/(8*pi))*cos(x(1))+10;
end

% F18

function o = F18(x)
o=(1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*(x(1)^2)-14*x(2)+6*x(1)*x(2)+3*x(2)^2))*...
    (30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*(x(1)^2)+48*x(2)-36*x(1)*x(2)+27*(x(2)^2)));
end

% F19

function o = F19(x)
 aH=[3 10 30;.1 10 35;3 10 30;.1 10 35];cH=[1 1.2 3 3.2];
 pH=[.3689 .117 .2673;.4699 .4387 .747;.1091 .8732 .5547;.03815 .5743 .8828];
 o=0;
for i=1:4
    o=o-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
end
end

% F20

function o = F20(x)
aH=[10 3 17 3.5 1.7 8;.05 10 17 .1 8 14;3 3.5 1.7 10 17 8;17 8 .05 10 .1 14];
cH=[1 1.2 3 3.2];
pH=[.1312 .1696 .5569 .0124 .8283 .5886;.2329 .4135 .8307 .3736 .1004 .9991;...
.2348 .1415 .3522 .2883 .3047 .6650;.4047 .8828 .8732 .5743 .1091 .0381];
o=0;
for i=1:4
    o=o-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
end
end

% F21

function o = F21(x)
aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];

o=0;
for i=1:5
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

% F22

function o = F22(x)
 aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
 cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];
 o=0;
for i=1:7
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

% F23

function o = F23(x)
 aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
 cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];

 o=0;
for i=1:10
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

function o=Ufun(x,a,k,m)
     o=k.*((x-a).^m).*(x>a)+k.*((-x-a).^m).*(x<(-a));
end

function y = F24(x)
% 
% Beale function.
% Matlab Code by A. Hedar (Sep. 29, 2005).
% The number of variables n = 2.
% 
y = (1.5-x(1)*(1-x(2)))^2+(2.25-x(1)*(1-x(2)^2))^2+(2.625-x(1)*(1-x(2)^3))^2;
end

function [y] = F25(xx)


x1 = xx(1);
x2 = xx(2);

term1 = (x1 + 2*x2 - 7)^2;
term2 = (2*x1 + x2 - 5)^2;

y = term1 + term2;

end

function [y] = F26(xx)

% % matya
x1 = xx(1);
x2 = xx(2);

term1 = 0.26 * (x1^2 + x2^2);
term2 = -0.48*x1*x2;

y = term1 + term2;

end
