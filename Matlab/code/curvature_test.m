clear; clc;

T = readmatrix('csvfiles\SiliconeRubber\SiliconeRubber_L0.108_N50_r0.009_Tt1248.50.csv');
x = T(:,1);
z = T(:,3);  % Ensure z is the correct lateral coordinate

% Ensure column vectors
x = x(:)*100;
z = z(:)*100;

% Create the reference path
refPath = referencePathFrenet([x z]);

% Create a Frenet trajectory generator
connector = trajectoryGeneratorFrenet(refPath);

% Define initial and terminal states
initState = [0 0 0 0 0 0];  
termState = [0.01*100 0 0 0 0 0]; 

% Generate the trajectory
[~, trajGlobal] = connect(connector, initState, termState, 10);

% Visualization
show(refPath);
hold on;
axis equal;
plot(trajGlobal.Trajectory(:,1), trajGlobal.Trajectory(:,2), 'b');
legend(["Waypoints", "Reference Path", "Trajectory"]);

figure(2)
kappa = curvature(refPath, z);
plot(kappa)