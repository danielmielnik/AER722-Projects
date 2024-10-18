%% AER 722 Project 1 | Sharvani Yadav, Alexia Economou, Daniel Mielnik
% Clean stuff
clear all;
clc;

syms y

%% Constants
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

function [qd, E, K, theta] = div_p(n, y, s, GJ_eta, c, ec, cl_alpha, alpha, q, eta)
    for i = 1:n
        fi(i,:) = i*(y^i);
        fi_prime = diff(fi(i,:));
        for j = 1:n
            fj = j*(y^j);
            fj_prime = diff(fj);
            
            E_ij(y) = GJ_eta*fi_prime*fj_prime;
            E(i,j) = int(E_ij,0,s);
            
            K_ij(y) = -((c^2)*ec*cl_alpha*(fi(i,:)*fj));
            K(i,j) = int(K_ij, 0, s);
        end
        F(i,:) = q*int((c^2)*ec*cl_alpha*alpha*fi(i,:), 0, s);
        %theta(i,:) = ((E(i,j)+q*K(i,j))^-1)*q*F(i,:);
        theta(i,:) = sum((((E(1:i,1:i)+q*K(1:i,1:i))^-1)*F(1:i,:)).*fi(1:i,:)); 
        
        c
        GJ_eta
        ec
        cl_alpha

    end
   
    vpa(E)
    vpa(K)
    subs(K, eta, 1)
    K_numeric = double(subs(K, eta, 1))
    
    %qd = -eig(double(E),double(K));
    qd = -eig(double(E),double(K));
    qd = min(qd);

    

    E = E(n, n);
    K = K(n, n);
    theta;

end

cl_beta = 0;
cm_beta = 0;
rho = 1.225;

%% Part A:
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

slope = ((c1/2)-c2)/(s-0)+c2
c(y) = slope*y + (c1 + c2)

ec = c1-0.25*c

for n = 1:10
    [qd1(n), E, K, fi] = div_p(n, y, s, GJ_eta, c, ec, cl_alpha, alpha_eta, 0, 0);
    qd2 = div_p(n+1, y, s, GJ_eta, c, ec, cl_alpha, alpha_eta, 0, 0);
    
    error = ferror(qd1(n), qd2);
    if error <= 0.1
        break
    end
end

n
error
qd1

% Divergence Dynamic Pressure VS Mode Graph
figure(1)
plot(1:n, qd1)
title('Divergence Dynamic Pressure VS Mode')
xlabel('Mode')
ylabel('Divergence Dynamic Pressure')


[qd1(n), E, K, theta] = div_p(n, y, s, GJ_eta, c, ec, cl_alpha, alpha_eta, (qd1(n)/2), 0);

max_deflection = vpa(max(subs(theta, y, 2)), 4)

%% Part B
L = 700;
V = 70;
%s/c > 3;
%maxdeflection < 1 deg
M_max = 300;
V_div = 150;

syms y k s eta

s = 1.7;
k = 0.8;

c1 = 0.3;
c2 = 0.375;

c_mean = ((c1+c2)+(c1+c1/2))/2

slope = ((c1/2)-c2)/(s-0)+c2
c(y) = slope*y + (c1 + c2)

ec = c1-0.25*c

sc_ratio = s/c_mean

int(GJ_eta*c)

GJ_eta = 8500*(1-k*eta);

c = subs(c, y, s*eta)

w_req = vpa(int(GJ_eta*c, eta, 0, 1))

cl_want = 700/(0.5*rho*(V^2)*s)

eta = max(vpa(solve(2*pi*sqrt(1-eta^2) == cl_want,eta)))

GJ_eta = 8500*(1-k*eta);

cl_alpha = 2*pi*sqrt(1-eta^2)

alpha_eta = 5-3*eta;

c = (111 * eta) / 200 + 3/4;

%c = double(subs(c, eta, 1))

[qd1, E, K, theta] = div_p(n, y, s, GJ_eta, c, ec, cl_alpha, alpha_eta, (qd1(n)/2), eta);

max_deflection = vpa(max(subs(theta, y, 2)), 4)

V_div = sqrt((2*qd1)/rho)

M_bend = 0.5*L*(s/2)

%M_bend = 

%hahahahah hehehehe heheh tehehehehe bahahahahah

%% Part C
syms y

[qd1, E, K, theta] = div_p(n, y, s, GJ_eta, c, ec, cl_alpha, alpha_eta, (qd1), eta);


