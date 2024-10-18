%% AER 722 Project 1 | Sharvani Yadav, Alexia Economou, Daniel Mielnik
% Clean stuff
clear all;
clc;

syms y

% Constants
function eta = feta(y, s)
    eta = y/s;
end

function alpha_eta = falpha_eta(eta)
    alpha_eta = 5-3*eta;
end

function cl_alpha = fcl_alpha(eta)
    cl_alpha = 2*pi*sqrt(1-eta^2);
end

function GJ_eta = fGJ_eta(k, eta)
    GJ_eta = 8500*(1-k*eta);
end

function error = ferror(q1, q2)
    error = ((q1 - q2)/q1)*100;
end

function qd = div_p(n, y, s, GJ_eta, c, ec, cl_alpha, alpha)
    for i = 1:n
        fi = i*(y^i);
        fi_prime = diff(fi);
        for j = 1:n
            fj = j*(y^j);
            fj_prime = diff(fj);
            
            E_ij(y) = GJ_eta*fi_prime*fj_prime;
            E(i,j) = int(E_ij,0,s);
            
            K_ij(y) = -((c^2)*ec*cl_alpha*(fi*fj));
            K(i,j) = int(K_ij, 0, s);
        end
        F(i,j) = q*int((c^2)*ec*cl_alpha*alpha*fi, 0, s);
        theta = inv(E(i,j)+q*K(i,j))*F(i,:)
        

    end
    
    vpa(E);
    vpa(K);
    
    qd = -eig(double(E),double(K));
    qd = min(qd);
end

cl_beta = 0;
cm_beta = 0;

% Part A:
% Constants
k = 0.25;
c1 = 0.35;
c2 = 0.4;
s = 2;

eta = feta(y, s)
GJ_eta = fGJ_eta(k, eta)
cl_alpha = fcl_alpha(eta)
alpha_eta = falpha_eta(eta)


n = 3;

slope = ((c1/2)-c2)/(s-0)
c(y) = slope*y + (c1 + c2)

ec = c1-0.25*c

for n = 1:10
    qd1(n) = div_p(n, y, s, GJ_eta, c, ec, cl_alpha, alpha_eta);
    qd2 = div_p(n+1, y, s, GJ_eta, c, ec, cl_alpha, alpha_eta);
    
    error = ferror(qd1(n), qd2);
    if error <= 0.1
        break
    end
end

n
error
qd1
%bkjkjkj

% Divergence Dynamic Pressure VS Mode Graph
figure(1)
plot(1:n, qd1)
title('Divergence Dynamic Pressure VS Mode')
xlabel('Mode')
ylabel('Divergence Dynamic Pressure')

y_wt = 2;
K = sqrt((0.5*qd1(n)*ec*(c^2)*cl_alpha)/GJ_eta)
K = subs(K, y, 2)
