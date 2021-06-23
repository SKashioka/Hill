%% Solve Hill Equation:
% Reference 宇宙機の相対運動_松永_東工大
% 論文: W. H. Clohessy and R. S. Wiltshire :A Terminal Guidance System for Satellite Rendezvous, Aerospace Science, Vol.29, 1960.
% 河野功 : 宇宙機の相対運動ダイナミクスと FF 誘導制 御への応用第 52 回宇宙科学技術連合講演会, 3C16，2008

%%
clc;
clear;


        
%%
x0= [0; 1000; 1000; 0; 1; -1];
tspan = [.0 12*60];
[t,x] = ode45(@(t,x) Hill_eq(t,x), tspan, x0);
plot(t,x(1), '-o');


%%

function Phi = Phi_eq(t,n)
    c = @(x_rad) cos(x_rad);
    s = @(x_rad) sin(x_rad);
    Phi = [4-3*c(n*t)      0   0       s(n*t)          2*(1-c(n*t))    0;
           6*(s(n*t)-n*t)  1   0       -2*(1-c(n*t))   4*s(n*t)-3*n*t  0;
           0               0   c(n*t)  0               0               s(n*t);
           3*s(n*t)        0   0       c(n*t)          2*s(n*t)        0;
           -6*(1 - c(n*t)) 0   0       -2*(s(n*t))     -3+4*c(n*t)     0;
           0               0   -s(n*t) 0               0               c(n*t)];
end
        
function dydt= Hill_eq(t,x)
    n = 0.02;
    dydt = Phi_eq(t,n) * x;
end

