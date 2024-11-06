function [tau1x, tau1y, tau2] = control_law(theta1x, theta1y, theta2, dtheta1x, dtheta1y, dtheta2, ...
                                            x1_ref, y1_ref, x2_ref, dtheta1x_ref, dtheta1y_ref, dtheta2_ref)

    % Parameters for the robot
    m1 = 4.8; % Mass of link 1
    m2 = 4.8; % Mass of link 2
    l1 = 1;   % Length of link 1
    l2 = 1;   % Length of link 2
    I1 = (1/3) * m1 * l1^2; % Inertia of link 1
    I2 = (1/3) * m2 * l2^2; % Inertia of link 2
    lc1 = 0.5 * l1; % Center of mass of link 1
    lc2 = 0.5 * l2; % Center of mass of link 2
    g = 9.81; % Gravity

    % PD Control gains (separate for each axis)
    Kp1x = 180; Kd1x = 20; % Gains for theta1x
    Kp1y = 2.337e-7 ; Kd1y = 9.8403e-5; % Gains for theta1y (different from theta1x)
    Kp2 = 224.831; Kd2 = 0.0094105;   % Gains for theta2
    
    % Compute the errors between the reference and current positions
    e1x = x1_ref - theta1x; % Error in theta1x
    e1y = y1_ref - theta1y; % Error in theta1y
    e2 = x2_ref - theta2;   % Error in theta2
    
    % Derivative errors (rate of error)
    de1x = dtheta1x_ref - dtheta1x; % Error in theta1x velocity
    de1y = dtheta1y_ref - dtheta1y; % Error in theta1y velocity
    de2 = dtheta2_ref - dtheta2;    % Error in theta2 velocity
    
    % Control inputs (PD controllers with separate gains)
    u1x = Kp1x * e1x + Kd1x * de1x; % PD control for theta1x
    u1y = Kp1y * e1y + Kd1y * de1y; % PD control for theta1y
    u2 = Kp2 * e2 + Kd2 * de2;      % PD control for theta2
    
    %% Mass Matrix M(theta) for the fully actuated system

    % Inertia terms for theta1x, theta1y, theta2
    m11x = m2 * l1^2 * cos(theta1x) + m2 * lc1^2 * cos(theta1x)^2 + I1;
    m12x = 0;
    m13 = 0;
    m14 = l1 * lc2 * m2 * cos(theta2) * cos(theta1x);  % Coupling between theta1x and theta2

    m21y = m1 * lc1^2 * cos(theta1y)^2 + I1;  % Inertia term for theta1y
    m22y = 0;

    m44 = m2 * lc2^2 * cos(theta2)^2 + I2;  % Inertia term for theta2
    
    % Assemble the mass matrix (fully actuated)
    M = [m11x, m12x, 0, m14;    % Theta1x dynamics
         0, m21y, 0, 0;         % Theta1y dynamics (no coupling)
         0, 0, I1, 0;           % Theta1z (if needed, currently no dynamics)
         m14, 0, 0, m44];       % Theta2 dynamics with coupling

    %% Coriolis Matrix C(theta, dtheta) for the fully actuated system

    % Coriolis terms for theta1x, theta1y, theta2
    c11x = -dtheta1x * (m2 * cos(theta1x) * sin(theta1x) * l1^2 + m1 * cos(theta1x) * sin(theta1x) * lc1^2);
    c14 = -dtheta2 * l1 * lc2 * m2 * cos(theta1x) * sin(theta2);  % Coupling between theta1x and theta2
    
    c21y = -dtheta1y * lc1^2 * m1 * cos(theta1y) * sin(theta1y);  % Coriolis for theta1y

    c41 = -dtheta1x * l1 * lc2 * m2 * cos(theta2) * sin(theta1x);  % Coupling between theta2 and theta1x
    c44 = -dtheta2 * lc2^2 * m2 * cos(theta2) * sin(theta2);       % Coriolis for theta2
    
    % Assemble the Coriolis matrix (fully actuated)
    C = [c11x, 0, 0, c14;   % Coriolis for theta1x
         0, c21y, 0, 0;     % Coriolis for theta1y
         0, 0, 0, 0;        % Theta1z (if needed, currently no dynamics)
         c41, 0, 0, c44];   % Coriolis for theta2
    
    %% Gravity Vector G(theta) for the fully actuated system

    G1x = -g * m2 * (lc2 * sin(theta2 + theta1x) + l1 * sin(theta1x)) - g * lc1 * m1 * sin(theta1x);  % Gravity for theta1x
    G1y = -g * lc1 * m1 * sin(theta1y);   % Gravity for theta1y
    G2 = -g * lc2 * m2 * sin(theta2 + theta1x);  % Gravity for theta2
    
    % Assemble the gravity vector (fully actuated)
    G = [G1x; G1y; 0; G2];  % Note: Theta1z is not considered here (can be added if needed)
    
    %% Compute generalized gravity force/torque (𝜏_g)
    tau_g = G;  % Since the gravitational forces have already been calculated as G

    %% Matrix B (for uncertainties handling)
    % Assuming B is a diagonal matrix that defines actuation on the joints.
    %B = diag([0, 1, 0, 1]);  % This could vary depending on how actuators affect the system.
    % Introduce adaptive gains based on error magnitude
    uncertain_gain_theta1x = 1 + 0.1 * abs(e1x + de1x);
    uncertain_gain_theta1y = 1 + 0.1 * abs(e1y + de1y);  % Gain adapts based on error for theta1y
    uncertain_gain_theta2 = 1 + 0.1 * abs(e2 + de2);     % Gain adapts based on error for theta2

    % Define B with variable entries for uncertain actuated joints
    B = diag([0, uncertain_gain_theta1y, 0, uncertain_gain_theta2]);

    %% Control Law with Matrix B and 𝜏_g
    % We need to account for the uncertainty in the system model using the pseudo-inverse of B.
    B_inv = pinv(B);  % Compute the pseudo-inverse of B.

    %% Compute the torques (PFL control law with 𝜏_g)
    % Tau = B_inv * (M * control_input + C * dtheta + G - tau_g)
    dtheta = [dtheta1x; dtheta1y; 0; dtheta2]; % Velocity vector for all joints
    
    % Control inputs for theta1x, theta1y, and theta2
    control_input = [u1x; u1y; 0; u2];
    
    % Compute the torques, factoring in matrix B and the generalized gravitational force
    tau = B_inv * (M * control_input + C * dtheta + G - tau_g);
    
    % Extract torques tau1x, tau1y, and tau2
    tau1x = tau(1); % Torque for theta1x
    tau1y = tau(2); % Torque for theta1y
    tau2 = tau(4);  % Torque for theta2
end
