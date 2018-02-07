clc
clear all;

%Declaring symbolic variables
syms l1 l2 I1 I2 theta1(t) theta2(t) m1 m2 a1 a2;
syms t1s t2s t3s td1s td2s td3s;
syms g t;

%Position of M1
x1 = l1 * cos(theta1(t));
y1 = l1 * sin(theta1(t));

X1 = [x1;y1];

%Position of M2
x2 = (a1 * cos(theta1(t))) + (l2 * cos(theta1(t) + theta2(t)));
y2 = (a1*sin(theta1(t))) + (l2*sin(theta1(t)+theta2(t)));

X2 = [x2;y2];

%Velocity
V1 = diff(X1, t);
V2 = diff(X2, t);

%Velocity squares (dot product)
v1 = simplify(V1.' * V1);
v2 = simplify(V2.' * V2);

%Calculating Kinetic energies
ke1 = (0.5 * m1 * v1) + (0.5 * I1 * (diff(theta1(t),t))^2);
ke2 = (0.5 * m2 * v2) + (0.5 * I2 * (diff(theta1(t), t) + diff(theta2(t), t))^2);

KE = simplify(ke1 + ke2);

%Calculating Potential energies
pe1 = m1 * g * y1;
pe2 = m2 * g * y2;

PE = simplify(pe1 + pe2);

%Performing Euler-Lagrange equation

L = KE - PE;
L = simplify(L);

%Substituting the temporary variables first and then taking derivatives wrt
%time
params = {theta1(t), theta2(t), diff(theta1(t),t), diff(theta2(t),t)};
params_sub = {t1s, t2s, td1s, td2s};

t1 = diff(subs(diff(subs( L, params, params_sub), td1s),params_sub, params))...
    - subs( diff(subs( L, params, params_sub), t1s),params_sub, params);
t1 = simplify(t1);

t2 = diff(subs(diff(subs( L, params, params_sub), td2s),params_sub, params))...
    - subs( diff(subs( L, params, params_sub), t2s),params_sub, params);
t2 = simplify(t2);

disp("Torques");

T = [t1;t2];

disp(T);

%Calculating numeric value using arbitary values
disp("Numeric value")
subs(T, {l1, l2, a1, a2, theta1, theta2, m1, m2, I1, I2}, {5, 7, 10, 14, -pi/2, 0, 5, 6, 8,9})