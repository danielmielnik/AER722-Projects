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
    alpha_eta = (5-3*eta)*(pi/180);
end

function cl_alpha = fcl_alpha(eta)
    cl_alpha = 2*pi*sqrt(1-eta.^2);
end

function GJ_eta = fGJ_eta(k, eta)
    GJ_eta = 8500*(1-k*eta);
end

function error = ferror(q1, q2)
    error = ((q1 - q2)/q1)*100;
end

function [qd, E, K, F, theta] = div_p(n, y, s, GJ_eta, c, ec, cl_alpha, alpha, q, eta)
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
        theta(i,:) = sum((((E(1:i,1:i)+q*K(1:i,1:i))^-1)*F(1:i,:)).*fi(1:i,:));

    end
   
    qd = -eig(double(E),double(K));
    qd = min(qd);

end

rho = 1.225;

%% Part A:
% Constants
k = 0.25;
c1 = 0.35;
c2 = 0.4;
s = 2;

eta = feta(y, s);
GJ_eta = fGJ_eta(k, eta);
cl_alpha = fcl_alpha(eta);
alpha_eta = falpha_eta(eta);


n = 3;

slope = ((c1/2)-c2)/(s-0)+c2;
c(y) = slope*y + (c1 + c2);

ec = c2 - 0.25*c;

for n = 1:10
    [qd1(n), E, K, F, fi] = div_p(n, y, s, GJ_eta, c, ec, cl_alpha, alpha_eta, 0, 0);
    qd2 = div_p(n+1, y, s, GJ_eta, c, ec, cl_alpha, alpha_eta, 0, 0);
    error = ferror(qd1(n), qd2);
    if error <= 0.1
        break
    end
end

vpa(E, 3)
vpa(K, 3)
vpa(qd1(n), 3)

% Divergence Dynamic Pressure VS Mode Graph
figure(1)
plot(1:n, qd1)
title('Divergence Dynamic Pressure VS Mode')
xlabel('Mode')
ylabel('Divergence Dynamic Pressure (Pa)')


[qd1(n), E, K, F, theta] = div_p(n, y, s, GJ_eta, c, ec, cl_alpha, alpha_eta, (qd1(n)/2), 0);


vpa(F, 3)
vpa(theta, 3)
theta_total = vpa(theta(1,:)+theta(2,:)+theta(3,:)+theta(4,:))
max_deflection = double(subs(theta_total, y, 2))

%% Part B
syms y k s eta
L = 700;
V = 70;
M_max = 300;
V_div = 150;

s = 1.7;
k = 0.8;

c1 = 0.3;
c2 = 0.375;

c_mean = ((c1+c2)+(c1+c1/2))/2;

slope = ((c1/2)-c2)/(s-0)+c2;
c(y) = slope*y + (c1 + c2);
ec = c1-0.25*c;
sc_ratio = s/c_mean;

S = int(c, y, 0, s);

vpa(int(int(GJ_eta*c, eta, 0, 1), y, 0, s))
GJ_eta = 8500*(1-k*eta);

c = subs(c, y, s*eta);

w_req = vpa(int(GJ_eta*c, eta, 0, 1));

cl_want = 700/(0.5*rho*(V^2)*S);

eta = max(vpa(solve(2*pi*sqrt(1-eta^2) == cl_want,eta)));
GJ_eta = 8500*(1-k*eta);

cl_alpha = 2*pi*sqrt(1-eta^2);
alpha_eta = 5-3*eta;

c = (111 * eta) / 200 + 3/4;

[qd1, E, K, F, theta] = div_p(n, y, s, GJ_eta, c, ec, cl_alpha, alpha_eta, (qd1(n)/2), eta);

theta_total = vpa(theta(1,:)+theta(2,:)+theta(3,:)+theta(4,:))
max_deflection = vpa(subs(theta_total, y, s), 4)

V_div = sqrt((2*qd1)/rho)

M_bend = 0.5*L*(s/2)

%% Part C
y_values = 0:0.1:1.7;
theta_total = theta(1,:)+theta(2,:)+theta(3,:)+theta(4,:);
theta_values = double(subs(theta_total, y, y_values));

figure(2)
fplot(theta_total, [0 s])
title('Wing Twist Distribution')
xlabel('Wing Span (m)')
ylabel('Wing Twist (deg)')

%% Part D
syms y_values

q = 0.5*rho*(V^2);

eta = feta(y, s);
alpha_values = (5 - 2*eta)*(pi/180);
cl_alpha = fcl_alpha(eta);
S = int(c, y, 0, s);

L = q*S*cl_alpha*(max(theta_values)-alpha_values);
L_rigid = q*S*cl_alpha*max(theta_values);

vpa(L)
vpa(L_rigid)

figure(3)
hold on
fplot(L,[0 s])
fplot(L_rigid,[0 s])
hold off
title('Lift Distribution')
legend('Aeroelastic Wing', 'Rigid Wing')
xlabel('Wing Span (m)')
ylabel('Lift (N)')

%% Part E
L_total = int(L, y, 0, s)
L_rigid_total = int(L_rigid, y, 0, s)
PercentChangeInLift = (abs(L_total-L_rigid_total)/L_rigid_total)*100

M_total = L_total*(1/3)*s
M_rigid_total = L_rigid_total*(1/3)*s
PercentChangeInMoment = (abs(M_total-M_rigid_total)/M_rigid_total)*100

M = L*y
M_rigid = L_rigid*y

figure(4)
hold on
fplot(M,[0 s])
fplot(M_rigid,[0 s])
hold off
title('Bending Moment Distribution')
legend('Aeroelastic Wing', 'Rigid Wing')
xlabel('Wing Span (m)')
ylabel('Moment (Nm)')

%% Part F
ca = 0.3*c;
y1 = 0.5*s;
y2 = 0.8*s;
V = 70;

phi = acos(2*ca-1)
cl_beta = 2*(pi-phi+sin(phi))
cm_beta = (-sin(phi)*(1-cos(phi)))/2

f = y*(4*s-2*y);
f_prime = diff(f);

A = int(c*cl_beta*y, y, y1, y2);
B = int((c^2)*(ec*cl_beta+cm_beta)*f, y, y1, y2);
C = int(c*cl_alpha*f*y, y, 0, s);
D = int((c^2)*ec*cl_alpha*f*f, y, 0, s);
E = int(GJ_eta*f_prime*f_prime, y, 0, s);
F = int(c*cl_alpha*(y^2), y, 0, s);
G = int((c^2)*ec*cl_alpha*f*y, y, 0, s);

q = (A*E)/(A*D-B*C)

rolling_power = V*(((A*(E/q)-D)+B*C)/(C*G+F*((E/q)-D)))

%% Part G
ca_sizes = [2*ca, 0.5*ca]; % Larger, Smaller

phi = acos(2*ca_sizes-1)
cl_beta = 2*(pi-phi+sin(phi))

A = int(c*cl_beta*y, y, y1, y2);
B = int((c^2)*(ec*cl_beta+cm_beta)*f, y, y1, y2);

rolling_power_sizes = V*(((A*(E/q)-D)+B*C)/(C*G+F*((E/q)-D)))