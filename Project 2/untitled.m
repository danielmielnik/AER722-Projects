%% AER 722 Project 2 | Sharvani Yadav, Alexia Economou, Daniel Mielnik

%% Constants
S = 1; % m
c = 0.5; % m
b = 0.5*c; % m
k_1 = 5; % kN/m
k_theta1 = 500; % Nm/rad
k_2 = 1; % kN/m
c_theta1 = 0; % Nms/rad
c_2 = 0; % Ns/m
m = 5; % kg
I_CG = 0.05; % kgm^2
x_g = 0.15; % m
m_1 = 2; % kg
x_m = 0.15; % m
U_max = 100; % m/s
rho = 1.225; % kg/m^3

%% Question 2 | Aeroelastic Equations and Critical Speeds
syms U y lambda

i = 1;
U_est = 1;

while i == 0
    U_est = U_est + 100;
    
    fu = 2*pi*0.5*U^2*S*c;
    
    M = [m+m_1, m*(x_g-b)+m_1*x_m; m*(x_g-b)+m_1*x_m, m*(x_g-b)^2+m_1*x_m^2+I_CG];
    
    B_s = [c_2, c_2*(b/2); c_2*(b/2), c_2*(b/2)^2];
    B_a = [1, b/2; 0, (b/2)^2];
    B_bar = B_s + U*fu*B_a;
    
    E = [k_1+k_2, k_2*(b/2)-k_1*b; k_2*(b/2)-k_1*b, k_2*(b/2)^2-k_1*b^2];
    
    K = [0, 1; 0, 0];
    K_bar = U^2*fu*K;
    
    J = [0, 0; 0, 0];
    
    A = [M, B_bar; J, M];
    C = [J, E+K_bar; -M, J];
    
    lambda = eig(C, -A);

    U_est

    if max(lambda > 0)
        %max(lambda) > 0
        U_flutter = U_est
        i = 1;
        break
    else
        U_est = U_est + 1;
    end

    if U_est > U_max
        break
    end
end

U_est = U_est + 100;

fu = 2*pi*0.5*S*c;

M = [m+m_1, m*(x_g-b)+m_1*x_m; m*(x_g-b)+m_1*x_m, m*(x_g-b)^2+m_1*x_m^2+I_CG];

B_s = [c_2, c_2*(b/2); c_2*(b/2), c_2*(b/2)^2];
B_a = [1, b/2; 0, (b/2)^2];
B_bar = B_s + U*fu*B_a;

E = [k_1+k_2, k_2*(b/2)-k_1*b; k_2*(b/2)-k_1*b, k_2*(b/2)^2-k_1*b^2];

K = [0, 1; 0, 0];
K_bar = U^2*fu*K;

J = [0, 0; 0, 0];

A = [M, B_bar; J, M];
C = [J, E+K_bar; -M, J];

A_v = double(subs(A, U, U_max));
C_v = double(subs(C, U, U_max));

lambdas = eig(C_v, -A_v)

if(lambdas>0) 
    f= 4
end

%Charmatr=[(M(1,1)*lambda.^2+(B_s(1,1)+U*B_bar(1,1))*lambda+E(1,1)), (M(1,2)*lambda.^2+(B_s(1,2)+U*B_bar(1,2))*lambda+E(1,2)+U.^2*K_bar(1,2)); (M(2,1)*lambda.^2+(B_s(2,1)+U*B_bar(2,1))*lambda+E(2,1)), (M(2,2)*lambda.^2+(B_s(2,2)+U*B_bar(2,2))*lambda+E(2,2)+U.^2*K_bar(2,2))]
%CharmatrVar = subs(Charmatr,[M(:,:),E(:,:), B_s(:,:), B_a(:,:), K(:,:), B_bar(:,:), K_bar(:,:)], [M(:,:), E(:,:), B_s(:,:), B_a(:,:), K(:,:), B_bar(:,:), K_bar(:,:)])

CharMat = [M(1,1)*lambda.^2 + U*B_bar(1,1)*lambda + E(1,1), M(1,2)*lambda.^2 + U*B_bar(1,2)*lambda + U.^2*K_bar(1,2); M(2,1)*lambda.^2 + U*B_bar(2,1)*lambda,  M(2,2)*lambda.^2 + U*B_bar(2,2)*lambda + U.^2*K_bar(2,2) + E(2,2)]

CharEqn = det(CharMat)
coll=collect(CharEqn,lambda)

%syms U y
eqn = M*y^2 + B_bar*y + K_bar*U^2 + E;

%cf = coeffs(collect(det(eqn), y), y);
cf = vpa(fliplr(coeffs((CharEqn),lambda)),4);
p0 = cf(1,1)
p1 = cf(1,2)
p2 = cf(1,3)
p3 = cf(1,4)
p4 = cf(1,5)

T3 = vpa(solve(p1*p2*p3 - p1^2*p4 - p0*p3^2 == 0, U),3)
P0 = vpa(solve(p0 == 0, U), 3);
Flutter = vpa(min(T3(T3>0)),5)
Divergence = vpa(min(P0(P0>0)),5)