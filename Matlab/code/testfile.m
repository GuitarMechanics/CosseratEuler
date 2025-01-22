hat=@(y)[0,-y(3),y(2);y(3),0,-y(1);-y(2),y(1),0];
global p R j n m v u q w vs us vt ut qt wt vst ust vh uh vsh ush qh wh nLL mLL x y z X Y Z  %Make vars available in whole program
%Parameters
force = 10;
L = 0.6;                       %Length (before strain)
N = 50;                        %Spatial resolution
E = 207e6;                     %Young's modulus
r = 0.01;                     %Cross-section radius
rt1 = [0.01;0;0];
rt2 = [0;0.01;0];
rho = 8000;                    %Density
% g = [-9.81;0;0];               %Gravity
g = [0;0;0];               %Gravity(ignored)
Bse = zeros(3);                %Material damping coefficients - shear and extension
Bbt = 1e-6*eye(3);             %Material damping coefficients - bending and torsion
C = 0.03*eye(3);               %Square-law-drag damping coefficients
dt = 0.015;                    %Time step
alpha = -0.2;                  %BDF-alpha parameter
STEPS = 300;                   %Number of timesteps to completion
vstar = @(s)[0;0;1];           %Value of v when static and absent loading
ustar = @(s)[0;0;0];           %Precurvature
vsstar = @(s)[0;0;1]
usstar = @(s)[0;0;0]
%Boundary Conditions
for i = 1 : STEPS
    p{i,1} = [0;0;0];          %Clamped base
    R{i,1} = eye(3);
    q{i,1} = [0;0;0];
    w{i,1} = [0;0;0];
end
nL = 0.0*g;                    %Start with a weight hung at the tip
mL = [0;0;0];

%Dependent Parameter Calculations
A = pi*r^2;                                 %Cross-sectional area
J = diag([pi*r^4/4  pi*r^4/4  pi*r^4/2]);   %Inertia
G = E/( 2*(1+0.3) );                        %Shear modulus
Kse = diag([G*A, G*A, E*A]);                %Stiffness matrix - shear and extension
Kbt = diag([E*J(1,1), E*J(2,2), G*J(3,3)]); %Stiffness matrix - bending and torsion
ds = L/(N-1);                               %Grid distance (before strain)
c0 = (1.5 + alpha) / ( dt*(1+alpha) );      %BDF-alpha coefficients
c1 = -2/dt;
c2 = (0.5 + alpha) / ( dt*(1+alpha) );
d1 = alpha / (1+alpha);

%Main Simulation
i = 1;
% fsolve(@staticIVP, zeros(6,1)); %Solve static BVP w/ shooting method
% applyStaticBDFalpha();
v = Kse\R'*n + vstar(ds*(j-1));
vh{i+1,j} = (c1+c2)*v{i,j};
disp(vh)