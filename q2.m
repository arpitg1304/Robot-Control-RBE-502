clc

%Declaring symbolic variables
syms l1 l2 l3;
syms theta1(t) theta2(t) theta3(t);
syms m1 m2 m3;
syms t g;
syms t1s t2s t3s td1s td2s td3s;

%Position of M1
x1 = 0;
y1 = 0;
z1 = l1;
X1 = [x1;y1;z1];

%Position of Point M2
x2 = l2 * cos(theta2(t)) * cos(theta1(t));
y2 = l2 * cos(theta2(t)) * sin(theta1(t));
z2 = l1 + (l2 * sin(theta2(t)));
X2 = [x2;y2;z2];

%Position of point M3
x3 = ((l2 * cos(theta2(t))) + l3 * cos(theta2(t) + theta3(t))) * cos(theta1(t));
y3 = ((l2 * cos(theta2(t))) + l3 * cos(theta2(t) + theta3(t))) * sin(theta1(t));
z3 = l1 + (l2 * sin(theta2(t))) + (l3 * sin(theta2(t) + theta3(t)));
X3 = [x3;y3;z3];

%Velocity
V1 = diff(X1, t);
V2 = diff(X2, t);
V3 = diff(X3, t);

%Calculating Potential Energies
pe1 = m1 * g * z1;
pe2 = m2 * g * z2;
pe3 = m3 * g * z3;
PE = simplify(pe1 + pe2 + pe3);

%Velocity squares (dot product)
v1 = simplify(V1.' * V1);
v2 = simplify(V2.' * V2);
v3 = simplify(V3.' * V3);

%Calculating Kinetic Energies
ke1 = 0.5 * m1 * v1;
ke2 = 0.5 * m2 * v2;
ke3 = 0.5 * m3 * v3;
KE = simplify(ke1 + ke2 + ke3);


%Performing Euler-Lagrange equation

L = (KE - PE);
L = simplify(L);

%Substituting the temporary variables first and then taking derivatives wrt
%time
params = {theta1(t), theta2(t), theta3(t), diff(theta1(t),t), diff(theta2(t),t), diff(theta3(t),t)};
params_sub = {t1s, t2s, t3s, td1s, td2s, td3s};

%Calculating Torques
t1 = diff(subs(diff(subs( L, params, params_sub), td1s),params_sub, params))...
    - subs( diff(subs( L, params, params_sub), t1s),params_sub, params);
t1 = simplify(t1);

t2 = diff(subs(diff(subs( L, params, params_sub), td2s),params_sub, params))...
    - subs( diff(subs( L, params, params_sub), t2s),params_sub, params);
t2 = simplify(t2);

t3 = diff(subs(diff(subs( L, params, params_sub), td3s),params_sub, params))...
    - subs( diff(subs( L, params, params_sub), t3s),params_sub, params);
t3 = simplify(t3);

disp("Torques");

T = [t1;t2;t3];

disp(T);