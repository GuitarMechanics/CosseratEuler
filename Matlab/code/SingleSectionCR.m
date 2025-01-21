function SingleSectionCR(Force)
    %clear all;
    clc;
    hat=@(y)[0,-y(3),y(2);y(3),0,-y(1);-y(2),y(1),0];
    global p R j n m v u q w vs us vt ut qt wt vst ust vh uh vsh ush qh wh nLL mLL x y z X Y Z  %Make vars available in whole program
    %Parameters
    L = 0.6;                       %Length (before strain)
    N = 100;                        %Spatial resolution
    E = 207e6;                     %Young's modulus
    r = 0.01;                     %Cross-section radius
    rt1 = [0.01;0;0];
    rt2 = [0;0.01;0];
    rho = 8000;                    %Density
    g = [-9.81;0;0];               %Gravity
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
    fsolve(@staticIVP, zeros(6,1)); %Solve static BVP w/ shooting method
    applyStaticBDFalpha();
%       visualize1();
    
    for i = 2 : STEPS

       if i < 5
            Tt1 = 0;
            Tt2 = 0;
        else
            Tt1 = Force;
            Tt2 = 0;
        end 

        fsolve(@dynamicIVP, [n{i-1,1}; m{i-1,1}]); %Solve semi-discretized PDE w/ shooting
        applyDynamicBDFalpha();
%          visualize();
    end

    
    
    %Function Definitions
    function applyStaticBDFalpha()
        for j = 1 : N-1
            vh{i+1,j} = (c1+c2)*v{i,j};
            uh{i+1,j} = (c1+c2)*u{i,j};
            vsh{i+1,j} = (c1+c2)*vs{i,j};
            ush{i+1,j} = (c1+c2)*us{i,j};
            qh{i+1,j} = [0;0;0];
            wh{i+1,j} = [0;0;0];
            q{i,j} = [0;0;0];
            w{i,j} = [0;0;0];
        end
    end

    function applyDynamicBDFalpha()
        for j = 1 : N-1        
            vh{i+1,j} = c1*v{i,j} + c2*v{i-1,j} + d1*vt{i,j};
            uh{i+1,j} = c1*u{i,j} + c2*u{i-1,j} + d1*ut{i,j};
            vsh{i+1,j} = c1*vs{i,j} + c2*vs{i-1,j} + d1*vst{i,j};
            ush{i+1,j} = c1*us{i,j} + c2*us{i-1,j} + d1*ust{i,j};
            qh{i+1,j} = c1*q{i,j} + c2*q{i-1,j} + d1*qt{i,j};
            wh{i+1,j} = c1*w{i,j} + c2*w{i-1,j} + d1*wt{i,j};
        end
    end

    function E = staticIVP(G)
        n{i,1} = G(1:3);
        m{i,1} = G(4:6);

        %Euler's method
        for j = 1 : N-1
            [ps, Rs, ns, ms, us{i,j}, vs{i,j} ,v{i,j}, u{i,j}] = staticODE(p{i,j},R{i,j},n{i,j},m{i,j});
            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} + ds*ns;
            m{i,j+1} = m{i,j} + ds*ms;
        end
        E = [ n{i,N} - nL;  m{i,N} - mL ];
    end

    function E = dynamicIVP(G)
        n{i,1} = G(1:3);
        m{i,1} = G(4:6);

        %Euler's method
        for j = 1 : N-1
            [ps, Rs, ns, ms, qs, ws, vs{i,j}, us{i,j},...
                 v{i,j}, u{i,j}, vt{i,j}, ut{i,j},...
                 qt{i,j}, wt{i,j},vst{i,j}, ust{i,j}] = dynamicODE(p{i,j},R{i,j},n{i,j},m{i,j},q{i,j},w{i,j});
            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} + ds*ns;
            m{i,j+1} = m{i,j} + ds*ms;
            q{i,j+1} = q{i,j} + ds*qs;
            w{i,j+1} = w{i,j} + ds*ws;
            
        end
        E = [n{i,N} - nLL ;  m{i,N} - mLL];
        
    end

    function [ps, Rs, ns, ms, vs, us, v, u] = staticODE(p,R,n,m)
        v = Kse\R'*n + vstar(ds*(j-1));
        u = Kbt\R'*m + ustar(ds*(j-1));
        
        ptsb = hat(u)*rt1+v;
        Tt = 0;
        At = -Tt/norm(ptsb)^3*hat(ptsb)*hat(ptsb);
        B = hat(rt1)*At;
        Gt = -At*hat(rt1);
        H =-B*hat(rt1);
        a = At*(hat(u)*ptsb);
        b = hat(rt1)*a;
        d = Kse*vsstar(ds*(j-1))-hat(u)*Kse*(v-vstar(ds*(j-1)))-a-R'*rho*A*g;
        c = Kbt*usstar(ds*(j-1))-hat(u)*Kbt*(u-ustar(ds*(j-1)))-hat(v)*Kse*(v-vstar(ds*(j-1)))-b;

        
        Mat = [Kse+At,Gt;B,Kbt+H];
        vs = 1/det(Mat)*((Kbt+H)*d-Gt*c);
        us = 1/det(Mat)*(-B*d+(Kse+At)*c);
        
         
        
        ps = R*v;
        Rs = R*hat(u);
        ns = -rho*A*g;
        ms = -hat(ps)*n;
    end

    function [ps,Rs,ns,ms,qs,ws,vs,us,v,u,vt,ut,qt,wt,vst,ust] = dynamicODE(p,R,n,m,q,w)
        v = (Kse + c0*Bse)\(R'*n + Kse*vstar(ds*(j-1)) - Bse*vh{i,j});
        u = (Kbt + c0*Bbt)\(R'*m + Kbt*ustar(ds*(j-1)) - Bbt*uh{i,j});
        vt = c0*v + vh{i,j};
        ut = c0*u + uh{i,j};
        qt = c0*q + qh{i,j};
        wt = c0*w + wh{i,j};
                
        ptsb1 = hat(u)*rt1+v;
        At1 = -Tt1/norm(ptsb1)^3*hat(ptsb1)*hat(ptsb1);
        Gt1 = -At1*hat(rt1);
        H1 =hat(rt1)*Gt1;
        a1 = At1*(hat(u)*ptsb1);
        b1 = hat(rt1)*a1;
        
        
        ptsb2 = hat(u)*rt2+v;
        At2 = -Tt2/norm(ptsb2)^3*hat(ptsb2)*hat(ptsb2);
        Gt2 = -At2*hat(rt2);
        H2 =hat(rt2)*Gt2;
        a2 = At2*(hat(u)*ptsb2);
        b2 = hat(rt2)*a2;
        
        a = a1 + a2;
        b = b1 + b2;
        At = At1 + At2;
        Gt = Gt1 + Gt2;
        H = H1 + H2;
               
        
        LamdaN = -a + rho*A*(hat(w)*q + qt)+ C*q.*abs(q)- R'*rho*A*g;
        LamdaM = -b + rho*(hat(w)*J*w + J*wt) -hat (v)*(Kse*(v-vstar(ds*(j-1)))+Bse*vt);
        GammaV = hat(u)*(Kse*(v-vstar(ds*(j-1)))+Bse*vt)-Kse*vsstar(ds*(j-1))+Bse*vsh{i,j};
        GammaU = hat(u)*(Kbt*(u-ustar(ds*(j-1)))+Bbt*ut)-Kbt*usstar(ds*(j-1))+Bbt*ush{i,j};

        
        
        Mat = [Kse+c0*Bse+At,Gt;Gt',Kbt+c0*Bbt+H];
        
        us = 1/det(Mat)*(-Gt'*(-GammaV+LamdaN)+(Kse+c0*Bse+At)*(-GammaU+LamdaM));
        vs = 1/det(Mat)*((Kbt+c0*Bbt+H)*(-GammaV+LamdaN)-Gt*(-GammaU+LamdaM));
        vst = c0*vs + vsh{i,j};
        ust = c0*us + ush{i,j};
        
        
        
        pts1 = R*hat(u)*rt1+R*v;
        pts2 = R*hat(u)*rt2+R*v;

        nLL = -Tt1*pts1/norm(pts1)-Tt2*pts2/norm(pts2);
        mLL = -Tt1*hat(R*rt1)*pts1/norm(pts1)-Tt2*hat(R*rt2)*pts2/norm(pts2);
        
        ps = R*v;
        Rs = R*hat(u);
        ns = rho*A*R*(hat(w)*q + qt) -R*(a+At*vs+Gt*us);
        ms = -R*(b+Gt'*vs+H*us)+rho*R*(hat(w)*J*w + J*wt) - hat(ps)*n;
        qs = vt - hat(u)*q + hat(w)*v;
        ws = ut - hat(u)*w;
    end

   function visualize1()
        for j = 1 : N,  x(j) = p{i,j}(1);  y(j) = p{i,j}(2); z(j) = p{i,j}(3);   end
        figure (2)
        fig = plot3(z,y,x); axis([-0.05*L 1.1*L  -0.1*L 0.1*L -0.05*L 0.1*L]);
        filename = sprintf('SingleSectionCR_L%.2f_N%d_r%.4f_Tt1%.2f_Tt2%.2f.csv', L, N, r, Tt1, Tt2);
        csvwrite(filename, [z' y' x']);
        xlabel('z (m)');  ylabel('y (m)'); zlabel('x (m)')
        hold on; grid on;  drawnow;  pause(0.05);
        fprintf(filename)
    end


    function visualize()
        if rem(i,10) == 0; 
        for j = 1 : N,  x(j) = p{i,j}(1);  y(j) = p{i,j}(2); z(j) = p{i,j}(3);   end
        figure (1)
        fig = plot3(z,y,x); axis([-0.05*L 1.1*L  -0.1*L 0.1*L -0.05*L 0.1*L]);
        xlabel('z (m)');  ylabel('y (m)'); zlabel('x (m)')
        hold on; grid on;  drawnow;  pause(0.05);
        end
    end

    for i = 1 : STEPS, U(i)=i*dt; X(i)=p{i,N}(1); Y(i)=p{i,N}(2); Z(i)=p{i,N}(3); end
%     figure (1)    
%     subplot(2,1,1)
%     plot(U,X);
%     xlabel('t (s)');  ylabel('x (m)'); title('Tip Displacement - X Component');
%     subplot(2,1,2)
%     plot(U,Y);
%     xlabel('t (s)');  ylabel('y (m)'); title('Tip Displacement - Y Component');
    visualize1()
%     saveas(gcf, '../results/SingleSectionCR.png')

end

