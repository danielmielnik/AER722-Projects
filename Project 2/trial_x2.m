clc
  %bound to change
syms k1 k2 k_theta U;
syms lambda;

x_m=0.15;
x_g=0.15;
%
m=5;
m1=2;
s=1;
c=0.5;
I_CG=0.05;
rho=1.225;
c_theta1=0;
c2=0;
b=c/2;
%l1=b-x_1
%l2=b-x_2
%lg=(x_g-b)
%lc=(x_c-b)
%Ib=I_CG+m*(lg)^2;
syms M E B_s B_a K B_bar_a K_bar_ [2 2]
Mv=[m-m1 m*(x_g-b)-m1*x_m; m*(x_g-b)-m1*x_m m*(x_g-b)^2-m1*x_m^2+I_CG]
Ev=[(k1+k2) k2*(b/2)-k1*b; k2*(b/2)-k1*b k2*(b/2)^2-k1*b^2+k_theta]
Bsv=[c2 c2*(b/2); c2*(b/2) c2*(b/2)^2]
Bav=[1 b/2; b/2 (b/2)^2]
Kv=[0 1; 0 b/2]
Bab=pi*c*s*rho*Bav
Kb=pi*c*s*rho*Kv
%Matrix Solution method
Z=zeros([2,2])

B=Bsv +U*Bab
A=[Mv, B; Z, Mv]
Avar=double(subs(A,[k1, k2,k_theta,U],[5000,1000,500,100]))
EK=Ev+U^2*Kb
C=[Z,EK;-Mv, Z]
Cvar=double(subs(C,[k1, k2, k_theta,U],[5000,1000,500,100]))
EigenVal=eig(Cvar,-Avar)
if(EigenVal>0) 
    f= 4
end

Charmatr=[(M(1,1)*lambda^2+(B_s(1,1)+U*B_bar_a(1,1))*lambda+E(1,1)), (M(1,2)*lambda^2+(B_s(1,2)+U*B_bar_a(1,2))*lambda+E(1,2)+U^2*K_bar_(1,2)); (M(2,1)*lambda^2+(B_s(2,1)+U*B_bar_a(2,1))*lambda+E(2,1)), (M(2,2)*lambda^2+(B_s(2,2)+U*B_bar_a(2,2))*lambda+E(2,2)+U^2*K_bar_(2,2))]
CharmatrVar=subs(Charmatr,[M(:,:),E(:,:), B_s(:,:), B_a(:,:), K(:,:), B_bar_a(:,:), K_bar_(:,:)], [Mv(:,:), Ev(:,:), Bsv(:,:), Bav(:,:), Kv(:,:), Bab(:,:), Kb(:,:)])

CharEqn = det(CharmatrVar)
coll=collect(CharEqn,lambda)

CharEqnP = subs(coll,[k1, k2,k_theta],[5000,1000,500])

Coef=vpa(fliplr(coeffs((CharEqnP),lambda)),4);
P4 = vpa(Coef(1,1),3);
P3 = vpa(Coef(1,2),3);
P2= vpa(Coef(1,3),3);
P1=vpa(Coef(1,4),3);
P0=vpa(Coef(1,5),3)
T=vpa(P1*P2*P3-P1^2*P4-P0*P3^2,3);
T = vpa(solve(T==0,U),3)
P0 = vpa(solve(P0 == 0, U), 3)

Flutter = vpa(min(T(T>0)),5)
Divergence = vpa(min(P0(P0>0)),5)

% clc;
% clear;
% 
% % Parameters (adjustable)
% syms k1 k2 k_theta U lambda;
% x_m = 0.15; x_g = 0.15; 
% m = 5; m1 = 2; s = 1; 
% c = 0.5; I_CG = 0.05; 
% rho = 1.225; b = c / 2;
% 
% % Mass matrix
% Mv = [m + m1, m * (x_g - b) + m1 * x_m; 
%       m * (x_g - b) + m1 * x_m, m * (x_g - b)^2 + m1 * x_m^2 + I_CG];
% 
% % Elastic and damping matrices
% Ev = [(k1 + k2), k2 * (b / 2) - k1 * b; 
%       k2 * (b / 2) - k1 * b, k2 * (b / 2)^2 - k1 * b^2 + k_theta];
% 
% Bsv = [0, 0; 0, 0]; % Structural damping is zero in the base case
% Bav = [1, b / 2; b / 2, (b / 2)^2];
% Kv = [0, 1; 0, b / 2];
% Bab = pi * c * s * rho * Bav;
% Kb = pi * c * s * rho * Kv;
% 
% % System matrices
% Z = zeros(2, 2);
% B = Bsv + U * Bab;
% A = [Mv, B; Z, Mv];
% EK = Ev + U^2 * Kb;
% C = [Z, EK; -Mv, Z];
% 
% % Substitution for numerical evaluation
% k1_val = 5000; k2_val = 1000; k_theta_val = 500;
% Umax = 100;
% Cvar = double(subs(C, [k1, k2, k_theta], [k1_val, k2_val, k_theta_val]));
% Avar = double(subs(A, [k1, k2, k_theta], [k1_val, k2_val, k_theta_val]));
% 
% % Eigenvalue analysis
% U_vals = linspace(1, Umax, 500); % Airspeed range
% real_parts = zeros(length(U_vals), 4);
% imag_parts = zeros(length(U_vals), 4);
% 
% for i = 1:length(U_vals)
%     Cvar_i = double(subs(Cvar, U, U_vals(i)));
%     Avar_i = double(subs(Avar, U, U_vals(i)));
%     EigenVal = eig(Cvar_i, -Avar_i);
%     real_parts(i, :) = real(EigenVal);
%     imag_parts(i, :) = imag(EigenVal);
% end
% 
% % Plot eigenvalue trends
% figure;
% subplot(2, 1, 1);
% plot(U_vals, real_parts, 'LineWidth', 1.5);
% title('Real Parts of Eigenvalues vs Airspeed');
% xlabel('Airspeed (m/s)');
% ylabel('Real(Eigenvalues)');
% grid on;
% 
% subplot(2, 1, 2);
% plot(U_vals, imag_parts, 'LineWidth', 1.5);
% title('Imaginary Parts of Eigenvalues vs Airspeed');
% xlabel('Airspeed (m/s)');
% ylabel('Imaginary(Eigenvalues)');
% grid on;
% 
% % Characteristic equation determinant
% Charmatr = [Mv(1, 1) * lambda^2 + (Bsv(1, 1) + U * Bab(1, 1)) * lambda + Ev(1, 1), ...
%             Mv(1, 2) * lambda^2 + (Bsv(1, 2) + U * Bab(1, 2)) * lambda + Ev(1, 2) + U^2 * Kb(1, 2);
%             Mv(2, 1) * lambda^2 + (Bsv(2, 1) + U * Bab(2, 1)) * lambda + Ev(2, 1), ...
%             Mv(2, 2) * lambda^2 + (Bsv(2, 2) + U * Bab(2, 2)) * lambda + Ev(2, 2) + U^2 * Kb(2, 2)];
% CharEqn = det(Charmatr);
% 
% % Numerical evaluation of determinant
% det_vals = zeros(size(U_vals));
% for i = 1:length(U_vals)
%     det_vals(i) = double(subs(CharEqn, [k1, k2, k_theta, U], [k1_val, k2_val, k_theta_val, U_vals(i)]));
% end
% 
% % Plot determinant
% figure;
% plot(U_vals, det_vals, 'LineWidth', 1.5);
% title('Determinant of Characteristic Matrix vs Airspeed');
% xlabel('Airspeed (m/s)');
% ylabel('Determinant');
% grid on;
% 
% % Finding critical speeds
% % Flutter speed occurs when the real part of an eigenvalue transitions from negative to positive.
% flutter_speed = min(U_vals(max(real_parts, [], 2) > 0));
% 
% % Divergence speed occurs when det(CharEqn) = 0 for the lowest speed.
% divergence_speed = min(U_vals(det_vals < 1e-6 & det_vals > -1e-6)); % Threshold for near-zero determinant
% 
% % Display results
% fprintf('Flutter Speed: %.3f m/s\n', flutter_speed);
% fprintf('Divergence Speed: %.3f m/s\n', divergence_speed);
