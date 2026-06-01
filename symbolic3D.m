%% derive_3D_brachiator_symbolic.m
% Symbolic derivation of the 3-D brachiating robot model
% Requires: Symbolic Math Toolbox

clear; clc;

%% ------------------------------------------------------------------------
% 1) Symbolic variables
%% ------------------------------------------------------------------------
syms th1 phi th2 real
syms dth1 dphi dth2 real
syms ddth1 ddphi ddth2 real

syms m1 m2 real positive
syms l1 l2 lc1 lc2 real positive
syms I1 I2 real positive
syms g real positive

q   = [th1;  phi;  th2];
dq  = [dth1; dphi; dth2];
ddq = [ddth1; ddphi; ddth2];

%% ------------------------------------------------------------------------
% 2) Shorthand trig terms
%% ------------------------------------------------------------------------
s1  = sin(th1);
c1  = cos(th1);

s2  = sin(th2);
c2  = cos(th2);

s12 = sin(th1 + th2);
c12 = cos(th1 + th2);

%% ------------------------------------------------------------------------
% 3) Centre-of-mass positions in the world frame
%    These follow the corrected kinematics:
%    - link 1 COM depends on (th1, phi)
%    - link 2 COM is obtained from end of link 1 + offset of link 2 COM
%% ------------------------------------------------------------------------
r1 = lc1 * [ ...
    s1*cos(phi);
    s1*sin(phi);
   -c1 ];

r2 = [ ...
    (l1*s1 + lc2*c2)*cos(phi);
    (l1*s1 + lc2*c2)*sin(phi);
   -l1*c1 + lc2*s2 ];

r1 = simplify(r1);
r2 = simplify(r2);

disp('COM position of link 1, r1(q) = ');
pretty(r1)

disp('COM position of link 2, r2(q) = ');
pretty(r2)

%% ------------------------------------------------------------------------
% 4) Translational Jacobians
%% ------------------------------------------------------------------------
Jv1 = simplify(jacobian(r1, q));
Jv2 = simplify(jacobian(r2, q));

disp('Translational Jacobian of link 1, Jv1 = ');
pretty(Jv1)

disp('Translational Jacobian of link 2, Jv2 = ');
pretty(Jv2)

%% ------------------------------------------------------------------------
% 5) Rotational contribution (as used in your thesis model)
%    These matrices match the structure in your text.
%    If later you use full CAD-based inertia tensors, this block can be
%    replaced by a Jw-based derivation.
%% ------------------------------------------------------------------------
Mrot1 = [ ...
    1, 0, 0;
    0, sin(th1)^2, 0;
    0, 0, 0 ];

Mrot2 = [ ...
    1, 0, 1;
    0, sin(th1 + th2)^2, 0;
    1, 0, 1 ];

%% ------------------------------------------------------------------------
% 6) Inertia matrix M(q)
%% ------------------------------------------------------------------------
M = simplify( ...
    m1*(Jv1.'*Jv1) + ...
    m2*(Jv2.'*Jv2) + ...
    I1*Mrot1 + ...
    I2*Mrot2 );

disp('Inertia matrix M(q) = ');
pretty(M)

%% ------------------------------------------------------------------------
% 7) Potential energy and gravity vector
%% ------------------------------------------------------------------------
V = simplify( ...
    m1*g*r1(3) + ...
    m2*g*r2(3) );

G = simplify(jacobian(V, q).');

disp('Potential energy V(q) = ');
pretty(V)

disp('Gravity vector G(q) = ');
pretty(G)

%% ------------------------------------------------------------------------
% 8) Coriolis / centrifugal vector via Christoffel symbols
%    Cvec(q,dq) such that:
%    M(q) ddq + Cvec(q,dq) + G(q) = tau_gen
%% ------------------------------------------------------------------------
Cvec = christoffelVector(M, q, dq);

disp('Coriolis / centrifugal vector Cvec(q,dq) = ');
pretty(Cvec)

%% ------------------------------------------------------------------------
% 9) Full manipulator equation
%% ------------------------------------------------------------------------
tau_gen = sym('tau_gen', [3 1], 'real');

EOM = simplify(M*ddq + Cvec + G - tau_gen);

disp('Manipulator equation residual M*ddq + Cvec + G - tau_gen = ');
pretty(EOM)

%% ------------------------------------------------------------------------
% 10) Actuation map
%     Only yaw (phi) and elbow (theta2) are actuated.
%% ------------------------------------------------------------------------
B = [0 0;
     1 0;
     0 1];

u = sym('u', [2 1], 'real');
tau_gen_from_actuators = B*u;

disp('Actuation matrix B = ');
disp(B)

disp('Generalized torque from actuators B*u = ');
pretty(tau_gen_from_actuators)

%% ------------------------------------------------------------------------
% 11) Optional: export symbolic expressions to MATLAB functions
%     Uncomment if you want numerical function handles automatically.
%% ------------------------------------------------------------------------
% matlabFunction(M, G, Cvec, Jv1, Jv2, ...
%     'File', 'brachiator3D_symbolic_dynamics', ...
%     'Vars', {q, dq, m1, m2, l1, l2, lc1, lc2, I1, I2, g});

%% ------------------------------------------------------------------------
% 12) Optional: save results to MAT-file
%% ------------------------------------------------------------------------
save('brachiator3D_symbolic_model.mat', ...
    'q','dq','ddq','r1','r2','Jv1','Jv2','M','V','G','Cvec','B');

fprintf('\nSymbolic derivation complete. Results saved to brachiator3D_symbolic_model.mat\n');

%% ========================================================================
% Local function: Christoffel-based Coriolis vector
%% ========================================================================
function Cvec = christoffelVector(M, q, dq)
    % Returns the Coriolis/centrifugal vector Cvec such that
    % M(q) qdd + Cvec(q,qd) + G(q) = tau

    n = length(q);
    Cvec = sym(zeros(n,1));

    for k = 1:n
        expr = sym(0);
        for i = 1:n
            for j = 1:n
                Gamma_kij = 1/2 * ( ...
                    diff(M(k,j), q(i)) + ...
                    diff(M(k,i), q(j)) - ...
                    diff(M(i,j), q(k)) );
                expr = expr + Gamma_kij * dq(i) * dq(j);
            end
        end
        Cvec(k) = simplify(expr);
    end
end
