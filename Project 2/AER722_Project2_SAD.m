%% AER 722 Project 2 | Sharvani Yadav, Alexia Economou, Daniel Mielnik
clear all;
clc;

%% Constants
S = 1; % m
c = 0.5; % m
b = 0.5*c; % m
k1 = 5000; % kN/m
k_theta1 = 500; % Nm/rad
k2 = 1000; % kN/m
c_theta1 = 0; % Nms/rad
c2 = 0; % Ns/m
m = 5; % kg
I_CG = 0.05; % kgm^2
x_g = 0.15; % m
m1 = 2; % kg
x_m = 0.15; % m
rho = 1.225; % kg/m^3
U_max = 100; % m/s

syms U lambda
q = 0.5*rho*U^2;

%% Question 2
M = [m+m1, m*(x_g-b)-m1*(b-x_m); m*(x_g-b)-m1*(b-x_m), m*(x_g-b)^2+m1*(b-x_m)^2+I_CG];

B_s = [c2, c2*(b/2); c2*(b/2), c2*(b/2)^2+c_theta1];
%B_a = ((2*pi*q*S*c)/U)*[1, b/2; -b/2, 0];
B_a = [1, b/2; -b/2, 0];
B_bar_a = pi*rho*c*S*B_a;
B_bar = B_s + U*B_bar_a;

E = [(k1+k2), k2*(b/2)-k1*b; k2*(b/2)-k1*b, k2*(b/2)^2+k1*b^2+k_theta1];

%K = (2*pi*q*S*c)*[0, 1; 0, -b/2];
K = [0, 1; 0, -b/2];
Kb = pi*c*S*rho*K;
K_bar = pi*rho*S*Kb;

J = [0, 0; 0, 0];

A = [M, B_bar; J, M];
C = [J, E+U^2*K_bar; -M, J];

Asub = double(subs(A, U, U_max));
Csub = double(subs(C, U, U_max));

eigenVal = eig(Csub, -Asub);

CharMatrix = [(M(1,1)*lambda^2+(B_s(1,1)+U*B_bar_a(1,1))*lambda+E(1,1)), (M(1,2)*lambda^2+(B_s(1,2)+U*B_bar_a(1,2))*lambda+E(1,2)+U^2*Kb(1,2)); (M(2,1)*lambda^2+(B_s(2,1)+U*B_bar_a(2,1))*lambda+E(2,1)), (M(2,2)*lambda^2+(B_s(2,2)+U*B_bar_a(2,2))*lambda+E(2,2)+U^2*Kb(2,2))];
CharEqn = det(CharMatrix);
eigCollect = collect(CharEqn, lambda);

Cf = vpa(fliplr(coeffs((CharEqn),lambda)),4);
p0 = Cf(1);
p1 = Cf(2);
p2 = Cf(3);
p3 = Cf(4);
p4 = Cf(5);

T3 = p1*p2*p3 - p1^2*p4 - p0*p3^2;
T3 = vpa(solve(T3==0,U),3);

U_Critical = vpa(min(T3(T3>0)), 5)


%% Question 3
U_Critical = 70;

realParts = [];
imagParts = [];

for U_new = 0:1:(1.2*U_Critical)
    Asub = double(subs(A, U, U_new));
    Csub = double(subs(C, U, U_new));
    
    eigenVal = eig(Csub, -Asub);

    realParts = [realParts; real(eigenVal(1)), real(eigenVal(2)), real(eigenVal(3)), real(eigenVal(4))];
    imagParts = [imagParts; imag(eigenVal(1)), imag(eigenVal(2)), imag(eigenVal(3)), imag(eigenVal(4))];

end

% Plotting the results
figure;
hold on;
plot(0:1:(1.2 * U_Critical), realParts);
plot(0:1:(1.2 * U_Critical), imagParts);
xlabel('Critical Speed (m/s)');
ylabel('Eigenvalue Parts');
title('Variation of Eigenvalues with Airspeed');
legend('Real Parts 1', 'Real Parts 2', 'Real Parts 3', 'Real Parts 4', 'Imaginary Parts 1', 'Imaginary Parts 2', 'Imaginary Parts 3', 'Imaginary Parts 4');
grid on;
hold off;


%% Question 4
k1_og = 5000;
k_theta1_og = 500;
k2_og = 1000;

% Spring Constant k1
U_Crit_k1 = [];

for i = 0:1:5

    k1 = k1_og*i;
    k_theta1 = k_theta1_og;
    k2 = k2_og;

    E = [(k1+k2), k2*(b/2)-k1*b; k2*(b/2)-k1*b, k2*(b/2)^2+k1*b^2+k_theta1];

    CharMatrix = [(M(1,1)*lambda^2+(B_s(1,1)+U*B_bar_a(1,1))*lambda+E(1,1)), (M(1,2)*lambda^2+(B_s(1,2)+U*B_bar_a(1,2))*lambda+E(1,2)+U^2*Kb(1,2)); (M(2,1)*lambda^2+(B_s(2,1)+U*B_bar_a(2,1))*lambda+E(2,1)), (M(2,2)*lambda^2+(B_s(2,2)+U*B_bar_a(2,2))*lambda+E(2,2)+U^2*Kb(2,2))];
    CharEqn = det(CharMatrix);

    Cf = vpa(fliplr(coeffs((CharEqn),lambda)),4);
    p0 = Cf(1);
    p1 = Cf(2);
    p2 = Cf(3);
    p3 = Cf(4);
    p4 = Cf(5);
    
    T3 = p1*p2*p3 - p1^2*p4 - p0*p3^2;
    T3 = vpa(solve(T3==0,U),3);
    
    U_Crit_k1 = [U_Crit_k1, vpa(min(T3(T3>0)), 5);];

end

spring_const = 0:k1_og:5*k1_og;

figure;
hold on;
plot(spring_const, U_Crit_k1, '-*');
xlabel('Spring Constant k_1 (N/m)');
ylabel('Critical Speed (m/s)');
title('Critical Speed Vs Spring Constant k_1');
grid on;

for i = 1:length(U_Crit_k1)
    text(0 + (i-1)*k1_og, U_Crit_k1(i), num2str(double(U_Crit_k1(i))), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

hold off;

% Spring Constant k2
U_Crit_k2 = [];

for i = 0:1:5

    k1 = k1_og;
    k_theta1 = k_theta1_og;
    k2 = k2_og*i;

    E = [(k1+k2), k2*(b/2)-k1*b; k2*(b/2)-k1*b, k2*(b/2)^2+k1*b^2+k_theta1];

    CharMatrix = [(M(1,1)*lambda^2+(B_s(1,1)+U*B_bar_a(1,1))*lambda+E(1,1)), (M(1,2)*lambda^2+(B_s(1,2)+U*B_bar_a(1,2))*lambda+E(1,2)+U^2*Kb(1,2)); (M(2,1)*lambda^2+(B_s(2,1)+U*B_bar_a(2,1))*lambda+E(2,1)), (M(2,2)*lambda^2+(B_s(2,2)+U*B_bar_a(2,2))*lambda+E(2,2)+U^2*Kb(2,2))];
    CharEqn = det(CharMatrix);

    Cf = vpa(fliplr(coeffs((CharEqn),lambda)),4);
    p0 = Cf(1);
    p1 = Cf(2);
    p2 = Cf(3);
    p3 = Cf(4);
    p4 = Cf(5);
    
    T3 = p1*p2*p3 - p1^2*p4 - p0*p3^2;
    T3 = vpa(solve(T3==0,U),3);
    
    U_Crit_k2 = [U_Crit_k2, vpa(min(T3(T3>0)), 5);];

end

spring_const = 0:k2_og:5*k2_og;

figure;
hold on;
plot(spring_const, U_Crit_k2, '-*');
xlabel('Spring Constant k_2 (N/m)');
ylabel('Critical Speed (m/s)');
title('Critical Speed Vs Spring Constant k_2');
ax = gca;
chart = ax.Children(1);
grid on;

for i = 1:length(U_Crit_k2)
    text(0 + (i-1)*k2_og, U_Crit_k2(i), num2str(double(U_Crit_k2(i))), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

hold off;

% Spring Constant k theata
U_Crit_ktheta = [];

for i = 0:1:5

    k1 = k1_og;
    k_theta1 = k_theta1_og*i;
    k2 = k2_og;

    E = [(k1+k2), k2*(b/2)-k1*b; k2*(b/2)-k1*b, k2*(b/2)^2+k1*b^2+k_theta1];

    CharMatrix = [(M(1,1)*lambda^2+(B_s(1,1)+U*B_bar_a(1,1))*lambda+E(1,1)), (M(1,2)*lambda^2+(B_s(1,2)+U*B_bar_a(1,2))*lambda+E(1,2)+U^2*Kb(1,2)); (M(2,1)*lambda^2+(B_s(2,1)+U*B_bar_a(2,1))*lambda+E(2,1)), (M(2,2)*lambda^2+(B_s(2,2)+U*B_bar_a(2,2))*lambda+E(2,2)+U^2*Kb(2,2))];
    CharEqn = det(CharMatrix);

    Cf = vpa(fliplr(coeffs((CharEqn),lambda)),4);
    p0 = Cf(1);
    p1 = Cf(2);
    p2 = Cf(3);
    p3 = Cf(4);
    p4 = Cf(5);
    
    T3 = p1*p2*p3 - p1^2*p4 - p0*p3^2;
    T3 = vpa(solve(T3==0,U),3);
    
    U_Crit_ktheta = [U_Crit_ktheta, vpa(min(T3(T3>0)), 5);];

end

spring_const = 0:k_theta1_og:5*k_theta1_og;

figure;
hold on;
plot(spring_const, U_Crit_ktheta, '-*');
xlabel('Spring Constant k theta (N/m)');
ylabel('Critical Speed (m/s)');
title('Critical Speed Vs Spring Constant k theta');
grid on;

for i = 1:length(U_Crit_ktheta)
    text(0 + (i-1)*k_theta1_og, U_Crit_ktheta(i), num2str(double(U_Crit_ktheta(i))), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

hold off;

%% Question 5
k1_vals = 1000:500:9000;
k2_vals = 1000:500:9000;
k_theta1_vals = 0:50:700;

U_max = 0;

% Loop through each combination of x, y, and z
for k1 = k1_vals
    for k2 = k2_vals
        if k1 + k2 <= 10000
            for k_theta1 = k_theta1_vals
                E = [(k1+k2), k2*(b/2)-k1*b; k2*(b/2)-k1*b, k2*(b/2)^2+k1*b^2+k_theta1];

                CharMatrix = [(M(1,1)*lambda^2+(B_s(1,1)+U*B_bar_a(1,1))*lambda+E(1,1)), (M(1,2)*lambda^2+(B_s(1,2)+U*B_bar_a(1,2))*lambda+E(1,2)+U^2*Kb(1,2)); (M(2,1)*lambda^2+(B_s(2,1)+U*B_bar_a(2,1))*lambda+E(2,1)), (M(2,2)*lambda^2+(B_s(2,2)+U*B_bar_a(2,2))*lambda+E(2,2)+U^2*Kb(2,2))];
                CharEqn = det(CharMatrix);
            
                Cf = vpa(fliplr(coeffs((CharEqn),lambda)),4);
                p0 = Cf(1);
                p1 = Cf(2);
                p2 = Cf(3);
                p3 = Cf(4);
                p4 = Cf(5);
                
                T3 = p1*p2*p3 - p1^2*p4 - p0*p3^2;
                T3 = vpa(solve(T3==0,U),3);
                
                U_Crit = vpa(min(T3(T3>0)), 5);

                if any(U_Crit > U_max) && all(U_Crit < 200)
                    U_max = U_Crit;
                    k1_max = k1;
                    k2_max = k2;
                    k_theta1_max = k_theta1;
                end
            end
        end
    end
end

% Valid combinations
U_max
k1_max
k2_max
k_theta1_max


%% Question 6
k1 = k1_og;
k_theta1 = k_theta1_og;
k2 = k2_og;

U_Crit_xm = [];

for i = 0:0.05:c

    x_m = i;

    M = [m+m1, m*(x_g-b)-m1*(b-x_m); m*(x_g-b)-m1*(b-x_m), m*(x_g-b)^2+m1*(b-x_m)^2+I_CG];

    E = [(k1+k2), k2*(b/2)-k1*b; k2*(b/2)-k1*b, k2*(b/2)^2+k1*b^2+k_theta1];
    
    CharMatrix = [(M(1,1)*lambda^2+(B_s(1,1)+U*B_bar_a(1,1))*lambda+E(1,1)), (M(1,2)*lambda^2+(B_s(1,2)+U*B_bar_a(1,2))*lambda+E(1,2)+U^2*Kb(1,2)); (M(2,1)*lambda^2+(B_s(2,1)+U*B_bar_a(2,1))*lambda+E(2,1)), (M(2,2)*lambda^2+(B_s(2,2)+U*B_bar_a(2,2))*lambda+E(2,2)+U^2*Kb(2,2))];
    CharEqn = det(CharMatrix);

    Cf = vpa(fliplr(coeffs((CharEqn),lambda)),4);
    p0 = Cf(1);
    p1 = Cf(2);
    p2 = Cf(3);
    p3 = Cf(4);
    p4 = Cf(5);
    
    T3 = p1*p2*p3 - p1^2*p4 - p0*p3^2;
    T3 = vpa(solve(T3==0,U),3);

    U_Crit_xm = [U_Crit_xm, vpa(min(T3(T3>0)), 5);];

end

figure;
hold on;
plot(0.05:0.05:c, U_Crit_xm, '-*');
xlabel('Point Mass Location (m)');
ylabel('Critical Speed (m/s)');
title('Critical Speed Vs Point Mass Location');
grid on;

for i = 1:length(U_Crit_xm)
    text(0 + (i-1)*0.05, U_Crit_xm(i), num2str(double(U_Crit_xm(i))), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

hold off;


%% Question 7
x_m = 0.15;

% Mechanical Damping c2
U_Crit_c2 = [];

for i = 0:1:20
    c2 = i;

    M = [m+m1, m*(x_g-b)-m1*(b-x_m); m*(x_g-b)-m1*(b-x_m), m*(x_g-b)^2+m1*(b-x_m)^2+I_CG];
    
    E = [(k1+k2), k2*(b/2)-k1*b; k2*(b/2)-k1*b, k2*(b/2)^2+k1*b^2+k_theta1];
    
    B_s = [c2, c2*(b/2); c2*(b/2), c2*(b/2)^2+c_theta1];
    B_bar = B_s + U*B_bar_a;
    
    CharMatrix = [(M(1,1)*lambda^2+(B_s(1,1)+U*B_bar_a(1,1))*lambda+E(1,1)), (M(1,2)*lambda^2+(B_s(1,2)+U*B_bar_a(1,2))*lambda+E(1,2)+U^2*Kb(1,2)); (M(2,1)*lambda^2+(B_s(2,1)+U*B_bar_a(2,1))*lambda+E(2,1)), (M(2,2)*lambda^2+(B_s(2,2)+U*B_bar_a(2,2))*lambda+E(2,2)+U^2*Kb(2,2))];
    CharEqn = det(CharMatrix);

    Cf = vpa(fliplr(coeffs((CharEqn),lambda)),4);
    p0 = Cf(1);
    p1 = Cf(2);
    p2 = Cf(3);
    p3 = Cf(4);
    p4 = Cf(5);
    
    T3 = p1*p2*p3 - p1^2*p4 - p0*p3^2;
    T3 = vpa(solve(T3==0,U),3);

    U_Crit_c2 = [U_Crit_c2, vpa(min(T3(T3>0)), 5);];

end

figure;
hold on;
plot(0:1:20, U_Crit_c2, '-*');
xlabel('Mechanical Damping c_2 (Ns/m)');
ylabel('Critical Speed (m/s)');
title('Critical Speed Vs Mechanical Damping c_2');
grid on;

for i = 1:2:length(U_Crit_c2)
    text(0 + (i-1)*1, U_Crit_c2(i), num2str(double(U_Crit_c2(i))), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

hold off;

% Mechanical Damping c_theta1
c2 = 0;

U_Crit_ctheta1 = [];

for i = 0:0.035:0.7
    c_theta1 = i;

    M = [m+m1, m*(x_g-b)-m1*(b-x_m); m*(x_g-b)-m1*(b-x_m), m*(x_g-b)^2+m1*(b-x_m)^2+I_CG];
    
    E = [(k1+k2), k2*(b/2)-k1*b; k2*(b/2)-k1*b, k2*(b/2)^2+k1*b^2+k_theta1];
    
    B_s = [c2, c2*(b/2); c2*(b/2), c2*(b/2)^2+c_theta1];
    B_bar = B_s + U*B_bar_a;
    
    CharMatrix = [(M(1,1)*lambda^2+(B_s(1,1)+U*B_bar_a(1,1))*lambda+E(1,1)), (M(1,2)*lambda^2+(B_s(1,2)+U*B_bar_a(1,2))*lambda+E(1,2)+U^2*Kb(1,2)); (M(2,1)*lambda^2+(B_s(2,1)+U*B_bar_a(2,1))*lambda+E(2,1)), (M(2,2)*lambda^2+(B_s(2,2)+U*B_bar_a(2,2))*lambda+E(2,2)+U^2*Kb(2,2))];
    CharEqn = det(CharMatrix);

    Cf = vpa(fliplr(coeffs((CharEqn),lambda)),4);
    p0 = Cf(1);
    p1 = Cf(2);
    p2 = Cf(3);
    p3 = Cf(4);
    p4 = Cf(5);
    
    T3 = p1*p2*p3 - p1^2*p4 - p0*p3^2;
    T3 = vpa(solve(T3==0,U),3);

    U_Crit_ctheta1 = [U_Crit_ctheta1, vpa(min(T3(T3>0)), 5);];

end

figure;
hold on;
plot(0:0.035:0.7, U_Crit_ctheta1, '-*');
xlabel('Mechanical Damping c theta (Nms/rad)');
ylabel('Critical Speed (m/s)');
title('Critical Speed Vs Mechanical Damping c theta');
grid on;

for i = 1:2:length(U_Crit_ctheta1)
    text(0 + (i-1)*0.035, U_Crit_ctheta1(i), num2str(double(U_Crit_ctheta1(i))), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

hold off;