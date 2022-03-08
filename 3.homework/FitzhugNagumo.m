clear all

%parameters, i = {P,Q,R,S,T}
theta_i = [-1/3*pi,-1/12*pi,0,1/12*pi,1/2*pi];
a_i     = [1.2,-5.0,30,-7.5,0.75];
b_i     = [0.25,0.1,0.1,0.1,0.4];


w  = 2*pi/72;
f2 = 0.23;
A  = 0.15;

timelimit = [0 10];
initial_cond = [-1 0 0];


ECG_mod(theta_i,a_i,b_i,w,f2,A,timelimit,initial_cond);


function [t,z] =ECG_mod(theta_i,a_i,b_i,w,f2,A,timelimit,initial_cond)

a = (1- sqrt(S(1)^2+S(2)^2));
theta = atan2(S(2),S(1));
d_theta = rem(theta - theta_i,2*pi);
z0 = A *sin(2*pi*f2*timelimit);
dECG = @(t,S) [a*S(1)-w*S(2);a*S(2)+w*S(1);
    -(a_i*d_theta*exp(-(d_theta^2/2*b_i^2))-S(3)+z0)]; 

[t,S_sol] = ode45(dECG,timelimit,initial_cond);

figure
plot3(S_sol(:,1),S_sol(:,2),S_sol(:,3))
xlabel('x')
ylabel('y')
zlabel('z')
title('trajectory of the output')
end
function 
