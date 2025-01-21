function ConcentricTubeCR
clear all;
    hat=@(y)[0,-y(3),y(2);y(3),0,-y(1);-y(2),y(1),0];
    global p R j n m v u q w vs us vt ut qt wt vst ust vh uh vsh ush qh wh nLL mLL n2 m2 nstatic2 mstatic2 pt x z %Make vars available in whole program
    %Parameters
    L = 0.2;                                %Length (before strain)
    N1 = 20;                                %Spatial resolution
    E1 = 207e9;                             %Young's modulus-first  segment
    E2 = 207e9;                             %Young's modulus-Second  segment 
    r1 = 0.001;                             %Cross-section radius-first  segment
    r2 = 0.0008;                            %Cross-section radius-Second  segment
    rt11 = [0.01;0;0];                      % tendon offset for the first section (near the base)
    rt12 = [0;0.01;0];                      % tendon offset for the first section (near the base)
    rt21 = [0.005;0;0];                     % tendon offset for the second section (near the tip)
    rt22 = [0;0.005;0];                     % tendon offset for the second section (near the tip)
    rho1 = 8000;                            %Density--First Segment
    rho2 = 8000;                            %Density--Second Segment
    rho3 = rho1 + (rho2 - rho1)*(r2/r1)^2;
    g = [-9.81;0;0];                        %Gravity
    Bse1 = zeros(3);                        %Material damping coefficients - shear and extension-First segment
    Bse2 = zeros(3);                        %Material damping coefficients - shear and extension-Second segment 
    Bse3 = Bse1;
    Bbt1 = 1e-6*eye(3);                     %Material damping coefficients - bending and torsion--First segment
    Bbt2 = 1e-6*eye(3);                     %Material damping coefficients - bending and torsion--Second segment
    Bbt3 = 1e-6*eye(3);
    C1 = 0.03*eye(3);                       %Square-law-drag damping coefficients-First segment
    C2 = 0.03*eye(3);                       %Square-law-drag damping coefficients-Second segment
    C3 = 0.03*eye(3);
    dt = 0.01;                              %Time step
    alpha = -0.2;                           %BDF-alpha parameter
    STEPS1 = 75;                            %Number of timesteps for first segment deflection
    STEPS2 = 20;                            %Number of timesteps for second segment deflection
    STEPS3 = 250;
    vstar = @(s)[0;0;1];                    %Value of v when static and absent loading
    ustar = @(s)[0;0;0];                    %Precurvature
    vsstar = @(s)[0;0;0]
    usstar = @(s)[0;0;0]
    %Boundary Conditions
    for i = 1 : STEPS1 + STEPS2 + STEPS3
        p{i,1} = [0;0;0];                   %Clamped base
        R{i,1} = eye(3);
        q{i,1} = [0;0;0];
        w{i,1} = [0;0;0];
    end
  
    nL = 0.0*g;                             %Start with a weight hung at the tip
    mL = [0;0;0];

    %Dependent Parameter Calculations First Segment
    A1 = pi*(r1^2-r2^2);                                    %Cross-sectional area
    A2 = pi*r2^2;                                           %Cross-sectional area
    A3 = A1 + A2;
    J1 = diag([pi*(r1^4/4-r2^4/4)  pi*(r1^4/4-r2^4/4)  pi*(r1^4/2-r2^4/4)]);    %Inertia
    J2 = diag([pi*r2^4/4  pi*r2^4/4  pi*r2^4/2]);                               %Inertia
    J3 = J1 + J2;
    G1 = E1/( 2*(1+0.3) );                                                      %Shear modulus
    G2 = E2/( 2*(1+0.3) );                                                      %Shear modulus
    Kse1 = diag([G1*A1, G1*A1, E1*A1]);                                         %Stiffness matrix - shear and extension
    Kse2 = diag([G2*A2, G2*A2, E2*A2]);                                         %Stiffness matrix - shear and extension
    Kse3 = Kse1;
    Kbt1 = diag([E1*J1(1,1), E1*J1(2,2), G1*J1(3,3)]);                          %Stiffness matrix - bending and torsion
    Kbt2 = diag([E2*J2(1,1), E2*J2(2,2), G2*J2(3,3)]);                          %Stiffness matrix - bending and torsion
    Kbt3 = Kbt1;
    ds = L/(N1-1);                                                              %Grid distance (before strain)
    c0 = (1.5 + alpha) / ( dt*(1+alpha) );                                      %BDF-alpha coefficients
    c1 = -2/dt;
    c2 = (0.5 + alpha) / ( dt*(1+alpha) );
    d1 = alpha / (1+alpha);
    
    %Main Simulation
    % Static Solution
    i = 1;
    Kse=Kse3; Kbt=Kbt3; A=A3; C=C3; rt=rt11; rho=rho3;
    fsolve(@staticIVP, zeros(6,1));             %Solve static BVP w/ shooting method
    applyStaticBDFalpha();                      % Caclulate Time derivatives based on BDF-Alpha
    NN=N1;
     U(i)=i*dt;X(i)=p{i,N1}(1); Y(i)=p{i,N1}(2); Z(i)=p{i,N1}(3);
%     visualize();                               % Plot Funcion

   

    %Dynamic Solution----------------------------------------
    for i = 2 : STEPS1
         Tt11 = 5;
         Tt12 = 2;
        
        fsolve(@dynamicIVP, [n{i-1,1}; m{i-1,1}]); %Solve semi-discretized PDE w/ shooting for first interval
        applyDynamicBDFalpha();                    % Caclulate Time derivatives based on BDF-Alpha
        NN=N1;
        U(i)=i*dt; X(i)=p{i,N1}(1); Y(i)=p{i,N1}(2); Z(i)=p{i,N1}(3);
%         visualize();                               % Plot Funcion
    end
    
    
    % Second Segment elongation
    for i = STEPS1+1 : STEPS1+STEPS2

            N2= N1+ i- STEPS1;
            ds = L/(N2-1);                               %Grid distance (before strain)
            fsolve(@staticIVPP, [n{i-1,1}; m{i-1,1}]); %Solve static BVP w/ shooting method
            applyStaticBDFalphaa();
            L = L + 0.01;
            U(i)=i*dt; X(i)=p{i,N2}(1); Y(i)=p{i,N2}(2); Z(i)=p{i,N2}(3);
%             visualize2();
    end
      
    for i = STEPS1+STEPS2+1 : STEPS1+STEPS2+STEPS3
         Tt11 = 5;
         Tt12 = 2;
         Tt21 = 5;
         Tt22 = 2;
               
            fsolve(@dynamicIVPP, [n{i-1,1}; m{i-1,1}]); %Solve semi-discretized PDE w/ shooting for second interval
            applyDynamicBDFalphaa()
            U(i)=i*dt; X(i)=p{i,N2}(1); Y(i)=p{i,N2}(2); Z(i)=p{i,N2}(3);
%             visualize2();
    end
    
   %------------------------------------- 
    
    %Function Definitions
    function applyStaticBDFalpha()
        for j = 1 : N1-1
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
    function applyStaticBDFalphaa()
        for j = 1 : N1-1
            vh{i+1,j} = (c1+c2)*v{i,j};
            uh{i+1,j} = (c1+c2)*u{i,j};
            vsh{i+1,j} = (c1+c2)*vs{i,j};
            ush{i+1,j} = (c1+c2)*us{i,j};
            qh{i+1,j} = [0;0;0];
            wh{i+1,j} = [0;0;0];
            q{i,j} = [0;0;0];
            w{i,j} = [0;0;0];
        end
        for j = N1 : N2-1
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
        for j = 1 : N1-1       
            vh{i+1,j} = c1*v{i,j} + c2*v{i-1,j} + d1*vt{i,j};
            uh{i+1,j} = c1*u{i,j} + c2*u{i-1,j} + d1*ut{i,j};
            vsh{i+1,j} = c1*vs{i,j} + c2*vs{i-1,j} + d1*vst{i,j};
            ush{i+1,j} = c1*us{i,j} + c2*us{i-1,j} + d1*ust{i,j};
            qh{i+1,j} = c1*q{i,j} + c2*q{i-1,j} + d1*qt{i,j};
            wh{i+1,j} = c1*w{i,j} + c2*w{i-1,j} + d1*wt{i,j};
        end
    end

    function applyDynamicBDFalphaa()
        for j = 1 : N2-1       
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
        for j = 1 : N1-1
            [ps, Rs, ns, ms, us{i,j}, vs{i,j}, v{i,j}, u{i,j}] = staticODE(p{i,j},R{i,j},n{i,j},m{i,j});
            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} + ds*ns;
            m{i,j+1} = m{i,j} + ds*ms;
        end
        E = [ n{i,N1} - nL;  m{i,N1} - mL ];
    end

    function E = staticIVPP(G)
        n{i,1} = G(1:3);
        m{i,1} = G(4:6);

       for j = 1 : N1-2
            [ps, Rs, ns, ms, us{i,j}, vs{i,j}, v{i,j}, u{i,j}] = staticODE1(p{i,j},R{i,j},n{i,j},m{i,j});
            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} + ds*ns;
            m{i,j+1} = m{i,j} + ds*ms;
       end   
       for j = N1-1 : N1-1
            [ps, Rs, ns, ms, us{i,j}, vs{i,j}, v{i,j}, u{i,j}] = staticODEN2(p{i,j},R{i,j},n{i,j},m{i,j});
            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} - nstatic2 + ds*ns;
            m{i,j+1} = m{i,j} - mstatic2 + ds*ms;
        end
        %Euler's method
        for j = N1 : N2-1
            [ps, Rs, ns, ms, us{i,j}, vs{i,j}, v{i,j}, u{i,j}] = staticODE2(p{i,j},R{i,j},n{i,j},m{i,j});
            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} + ds*ns;
            m{i,j+1} = m{i,j} + ds*ms;
        end
        E = [ n{i,N2} - nL;  m{i,N2} - mL ];
    end


    function E = dynamicIVP(G)
        n{i,1} = G(1:3);
        m{i,1} = G(4:6);

        %Euler's method
        for j = 1 : N1 - 1
            [ps, Rs, ns, ms, qs, ws, vs{i,j}, us{i,j},...
                 v{i,j}, u{i,j}, vt{i,j}, ut{i,j},...
                 qt{i,j}, wt{i,j},vst{i,j}, ust{i,j}] = ...
                dynamicODE(p{i,j},R{i,j},n{i,j},m{i,j},q{i,j},w{i,j});
            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} + ds*ns;
            m{i,j+1} = m{i,j} + ds*ms;
            q{i,j+1} = q{i,j} + ds*qs;
            w{i,j+1} = w{i,j} + ds*ws;
            
        end
        
        E = [n{i,N1} - nLL ;  m{i,N1} - mLL];
        
    end
 
    function E = dynamicIVPP(G)
        n{i,1} = G(1:3);
        m{i,1} = G(4:6);

        %Euler's method
        for j = 1 : N1-2
            [ps, Rs, ns, ms, qs, ws, vs{i,j}, us{i,j},...
                 v{i,j}, u{i,j}, vt{i,j}, ut{i,j},...
                 qt{i,j}, wt{i,j},vst{i,j}, ust{i,j}] = dynamicODE1(p{i,j},R{i,j},n{i,j},m{i,j},q{i,j},w{i,j});
            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} + ds*ns;
            m{i,j+1} = m{i,j} + ds*ms;
            q{i,j+1} = q{i,j} + ds*qs;
            w{i,j+1} = w{i,j} + ds*ws;
            
        end
       for j = N1-1 :N1-1 
        [ps, Rs, ns, ms, qs, ws, vs{i,j}, us{i,j},...
                 v{i,j}, u{i,j}, vt{i,j}, ut{i,j},...
                 qt{i,j}, wt{i,j},vst{i,j}, ust{i,j}] = ...
                dynamicODEN2(p{i,j},R{i,j},n{i,j},m{i,j},q{i,j},w{i,j});
            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} - n2 + ds*ns;
            m{i,j+1} = m{i,j} - m2 + ds*ms;
            q{i,j+1} = q{i,j} + ds*qs;
            w{i,j+1} = w{i,j} + ds*ws;
            
        end
        for j = N1 : N2 - 1
        [ps, Rs, ns, ms, qs, ws, vs{i,j}, us{i,j},...
                 v{i,j}, u{i,j}, vt{i,j}, ut{i,j},...
                 qt{i,j}, wt{i,j},vst{i,j}, ust{i,j}] = ...
                dynamicODE2(p{i,j},R{i,j},n{i,j},m{i,j},q{i,j},w{i,j});
            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} + ds*ns;
            m{i,j+1} = m{i,j} + ds*ms;
            q{i,j+1} = q{i,j} + ds*qs;
            w{i,j+1} = w{i,j} + ds*ws;
            
        end
        
        
        
        E = [n{i,N2} - nLL ;  m{i,N2} - mLL];
         
    end



    function [ps, Rs, ns, ms,  us, vs, v, u] = staticODE(p,R,n,m)
        v = Kse3\R'*n + vstar(ds*(j-1));
        u = Kbt3\R'*m + ustar(ds*(j-1));
        
        ptsb = hat(u)*rt+v;
        Tt = 0;
        At = -Tt/norm(ptsb)^3*hat(ptsb)*hat(ptsb);
        B = -hat(rt)*At;
        Gt = -At*hat(rt);
        H =B*hat(rt);
        a = At*(hat(u)*ptsb);
        b = hat(rt)*a;
        d = Kse*vsstar(ds*(j-1))-hat(u)*Kse*(v-vstar(ds*(j-1)))-a-R'*rho*A*g;
        c = Kbt*usstar(ds*(j-1))-hat(u)*Kbt*(u-ustar(ds*(j-1)))-hat(v)*Kse*(v-vstar(ds*(j-1)))-b;

        
        Mat = [Kse+At,Gt;B,Kbt+H];
        vs = 1/det(Mat)*((Kbt+H)*d-Gt*c);
        us = 1/det(Mat)*(-B*d+(Kse+At)*c);
        
        pts = R*hat(u)*rt+R*v;
        
        
        ps = R*v;
        Rs = R*hat(u);
        ns = -rho3*A3*g-R*(a+At*vs+Gt*us);
        ms = -R*(b+B*vs+H*us)-hat(ps)*n;
    end

    function [ps, Rs, ns, ms,  us, vs, v, u] = staticODE1(p,R,n,m)
        v = Kse3\R'*n + vstar(ds*(j-1));
        u = Kbt3\R'*m + ustar(ds*(j-1));
        
        Tt11 = 5;
        Tt12 = 2;
        ptsb1 = hat(u)*rt11+v;
        At1 = -Tt11/norm(ptsb1)^3*hat(ptsb1)*hat(ptsb1);
        B1 = -hat(rt11)*At1;
        Gt1 = -At1*hat(rt11);
        H1 =B1*hat(rt11);
        a1 = At1*(hat(u)*ptsb1);
        b1 = hat(rt11)*a1;
        
        ptsb2 = hat(u)*rt12+v;
        At2 = -Tt12/norm(ptsb2)^3*hat(ptsb2)*hat(ptsb2);
        B2 = hat(rt12)*At2;
        Gt2 = -At2*hat(rt12);
        H2 =-B2*hat(rt12);
        a2 = At2*(hat(u)*ptsb2);
        b2 = hat(rt12)*a2;
        
        a = a1 + a2;
        b = b1 + b2;
        H = H1 + H2;
        Gt = Gt1 + Gt2;
        B = B1 + B2;
        At = At1 + At2;
        
                
        d = Kse3*vsstar(ds*(j-1))-hat(u)*Kse3*(v-vstar(ds*(j-1)))-a-R'*rho3*A3*g;
        c = Kbt3*usstar(ds*(j-1))-hat(u)*Kbt3*(u-ustar(ds*(j-1)))-hat(v)*Kse3*(v-vstar(ds*(j-1)))-b;

        
        Mat = [Kse3+At,Gt;B,Kbt3+H];
        vs = 1/det(Mat)*((Kbt3+H)*d-Gt*c);
        us = 1/det(Mat)*(-B*d+(Kse3+At)*c);
         
        ps = R*v;
        Rs = R*hat(u);
        ns = -rho3*A3*g-R*(a+At*vs+Gt*us);
        ms = -R*(b+B*vs+H*us)-hat(ps)*n;
    end

   function [ps, Rs, ns, ms,  us, vs, v, u] = staticODEN2(p,R,n,m)
        v = Kse3\R'*n + vstar(ds*(j-1));
        u = Kbt3\R'*m + ustar(ds*(j-1));
        
        
        Tt11 = 5;
        Tt12 = 2;
        ptsb1 = hat(u)*rt11+v;
        At1 = -Tt11/norm(ptsb1)^3*hat(ptsb1)*hat(ptsb1);
        B1 = hat(rt11)*At1;
        Gt1 = -At1*hat(rt11);
        H1 =-B1*hat(rt11);
        a1 = At1*(hat(u)*ptsb1);
        b1 = hat(rt11)*a1;
        
        ptsb2 = hat(u)*rt12+v;
        At2 = -Tt12/norm(ptsb2)^3*hat(ptsb2)*hat(ptsb2);
        B2 = hat(rt12)*At2;
        Gt2 = -At2*hat(rt12);
        H2 =-B2*hat(rt12);
        a2 = At2*(hat(u)*ptsb2);
        b2 = hat(rt12)*a2;
        
        a = a1 + a2;
        b = b1 + b2;
        H = H1 + H2;
        Gt = Gt1 + Gt2;
        B = B1 + B2;
        At = At1 + At2;
        
                
        d = Kse3*vsstar(ds*(j-1))-hat(u)*Kse3*(v-vstar(ds*(j-1)))-a-R'*rho3*A3*g;
        c = Kbt3*usstar(ds*(j-1))-hat(u)*Kbt3*(u-ustar(ds*(j-1)))-hat(v)*Kse3*(v-vstar(ds*(j-1)))-b;

        
        Mat = [Kse3+At,Gt;B,Kbt3+H];
        vs = 1/det(Mat)*((Kbt3+H)*d-Gt*c);
        us = 1/det(Mat)*(-B*d+(Kse3+At)*c);
        
        pts1 = R*hat(u)*rt11+R*v;
        pts2 = R*hat(u)*rt12+R*v;
        
        nstatic2 = -Tt11*pts1/norm(pts1)  -Tt12*pts2/norm(pts2);
        mstatic2 = -Tt11*hat(R*rt11)*pts1/norm(pts1) -Tt12*hat(R*rt12)*pts2/norm(pts2);
        
        ps = R*v;
        Rs = R*hat(u);
        ns = -rho3*A3*g-R*(a+At*vs+Gt*us);
        ms = -R*(b+B*vs+H*us)-hat(ps)*n;
   end    
   

   function [ps, Rs, ns, ms,  us, vs, v, u] = staticODE2(p,R,n,m)
        v = Kse2\R'*n + vstar(ds*(j-1));
        u = Kbt2\R'*m + ustar(ds*(j-1));
        
        ptsb = hat(u)*rt11+v;
        Tt = 0;
        At = -Tt/norm(ptsb)^3*hat(ptsb)*hat(ptsb);
        B = hat(rt11)*At;
        Gt = -At*hat(rt11);
        H =-B*hat(rt11);
        a = At*(hat(u)*ptsb);
        b = hat(rt11)*a;
        d = Kse2*vsstar(ds*(j-1))-hat(u)*Kse2*(v-vstar(ds*(j-1)))-a-R'*rho2*A2*g;
        c = Kbt2*usstar(ds*(j-1))-hat(u)*Kbt2*(u-ustar(ds*(j-1)))-hat(v)*Kse2*(v-vstar(ds*(j-1)))-b;

        
        Mat = [Kse2+At,Gt;B,Kbt2+H];
        vs = 1/det(Mat)*((Kbt2+H)*d-Gt*c);
        us = 1/det(Mat)*(-B*d+(Kse2+At)*c);
         
        ps = R*v;
        Rs = R*hat(u);
        ns = -rho3*A3*g-R*(a+At*vs+Gt*us);
        ms = -R*(b+B*vs+H*us)-hat(ps)*n;
    end




    function [ps,Rs,ns,ms,qs,ws,vs, us,v,u,vt,ut,qt,wt,vst,ust] = dynamicODE(p,R,n,m,q,w)
        v = (Kse3 + c0*Bse3)\(R'*n + Kse3*vstar(ds*(j-1)) - Bse3*vh{i,j});
        u = (Kbt3 + c0*Bbt1)\(R'*m + Kbt3*ustar(ds*(j-1)) - Bbt1*uh{i,j});
        vt = c0*v + vh{i,j};
        ut = c0*u + uh{i,j};
        qt = c0*q + qh{i,j};
        wt = c0*w + wh{i,j};
        
        ptsb1 = hat(u)*rt11+v;
        At1 = -Tt11/norm(ptsb1)^3*hat(ptsb1)*hat(ptsb1);
        Gt1 = -At1*hat(rt11);
        H1 =hat(rt11)*Gt1;
        a1 = At1*(hat(u)*ptsb1);
        b1 = hat(rt11)*a1;
        
        
        ptsb2 = hat(u)*rt12+v;
        At2 = -Tt12/norm(ptsb2)^3*hat(ptsb2)*hat(ptsb2);
        Gt2 = -At2*hat(rt12);
        H2 =hat(rt12)*Gt2;
        a2 = At2*(hat(u)*ptsb2);
        b2 = hat(rt12)*a2;
        
        a = a1 + a2;
        b = b1 + b2;
        At = At1 + At2;
        Gt = Gt1 + Gt2;
        H = H1 + H2;
        
        
        
        
        LamdaN = -a + rho3*A3*(hat(w)*q + qt)+ C*q.*abs(q)- R'*rho3*A3*g;
        LamdaM = -b + rho3*(hat(w)*J3*w + J3*wt) -hat (v)*(Kse3*(v-vstar(ds*(j-1)))+Bse3*vt);
        GammaV = hat(u)*(Kse3*(v-vstar(ds*(j-1)))+Bse3*vt)-Kse3*vsstar(ds*(j-1))+Bse3*vsh{i,j};
        GammaU = hat(u)*(Kbt3*(u-ustar(ds*(j-1)))+Bbt3*ut)-Kbt3*usstar(ds*(j-1))+Bbt3*ush{i,j};
        

        
        Mat = [Kse3+c0*Bse3+At,Gt;Gt',Kbt3+c0*Bbt3+H];
        
        us = 1/det(Mat)*(-Gt'*(-GammaV+LamdaN)+(Kse3+c0*Bse1+At)*(-GammaU+LamdaM));
        vs = 1/det(Mat)*((Kbt3+c0*Bbt1+H)*(-GammaV+LamdaN)-Gt*(-GammaU+LamdaM));
        vst = c0*vs + vsh{i,j};
        ust = c0*us + ush{i,j};
                

        pts1 = R*hat(u)*rt11+R*v;
        pts2 = R*hat(u)*rt12+R*v;
       
        nLL = -Tt11*pts1/norm(pts1)  -Tt12*pts2/norm(pts2);
        mLL = -Tt11*hat(R*rt11)*pts1/norm(pts1) -Tt12*hat(R*rt12)*pts2/norm(pts2);
        
        ps = R*v;
        Rs = R*hat(u);
        ns = rho3*A3*R*(hat(w)*q + qt) -R*(a+At*vs+Gt*us)+R*C*q.*abs(q) - rho3*A3*g;
        ms = -R*(b+Gt'*vs+H*us)+rho3*R*(hat(w)*J3*w + J3*wt) - hat(ps)*n;
        qs = vt - hat(u)*q + hat(w)*v;
        ws = ut - hat(u)*w;
    end
    
    function [ps,Rs,ns,ms,qs,ws,vs, us,v,u,vt,ut,qt,wt,vst,ust] = dynamicODE1(~,R,n,m,q,w)
        v = (Kse3 + c0*Bse3)\(R'*n + Kse3*vstar(ds*(j-1)) - Bse3*vh{i,j});
        u = (Kbt3 + c0*Bbt1)\(R'*m + Kbt3*ustar(ds*(j-1)) - Bbt1*uh{i,j});
        vt = c0*v + vh{i,j};
        ut = c0*u + uh{i,j};
        qt = c0*q + qh{i,j};
        wt = c0*w + wh{i,j};
                
        ptsb1 = hat(u)*rt11+v;
        At1 = -Tt11/norm(ptsb1)^3*hat(ptsb1)*hat(ptsb1);
        Gt1 = -At1*hat(rt11);
        H1 =hat(rt11)*Gt1;
        a1 = At1*(hat(u)*ptsb1);
        b1 = hat(rt11)*a1;
        
        
        ptsb2 = hat(u)*rt12+v;
        At2 = -Tt12/norm(ptsb2)^3*hat(ptsb2)*hat(ptsb2);
        Gt2 = -At2*hat(rt12);
        H2 =hat(rt12)*Gt2;
        a2 = At2*(hat(u)*ptsb2);
        b2 = hat(rt12)*a2;
        
        ptsb3 = hat(u)*rt21+v;
        At3 = -Tt21/norm(ptsb3)^3*hat(ptsb3)*hat(ptsb3);
        Gt3 = -At3*hat(rt21);
        H3 =hat(rt21)*Gt3;
        a3 = At3*(hat(u)*ptsb3);
        b3 = hat(rt21)*a3;
        
        
        ptsb4 = hat(u)*rt22+v;
        At4 = -Tt22/norm(ptsb4)^3*hat(ptsb4)*hat(ptsb4);
        Gt4 = -At4*hat(rt22);
        H4 =hat(rt22)*Gt4;
        a4 = At4*(hat(u)*ptsb4);
        b4 = hat(rt22)*a4; 
        
        a = a1 + a2 + a3 + a4;
        b = b1 + b2 + b3 + b4;
        At = At1 + At2 + At3 + At4;
        Gt = Gt1 + Gt2 + Gt3 + Gt4;
        H = H1 + H2 + H3 + H4;
        
              
        LamdaN = -a + rho3*A3*(hat(w)*q + qt)+ C*q.*abs(q)- R'*rho3*A3*g;
        LamdaM = -b + rho3*(hat(w)*J3*w + J3*wt) -hat (v)*(Kse3*(v-vstar(ds*(j-1)))+Bse1*vt);
        GammaV = hat(u)*(Kse3*(v-vstar(ds*(j-1)))+Bse3*vt)-Kse3*vsstar(ds*(j-1))+Bse3*vsh{i,j};
        GammaU = hat(u)*(Kbt3*(u-ustar(ds*(j-1)))+Bbt3*ut)-Kbt3*usstar(ds*(j-1))+Bbt3*ush{i,j};
        

        
        Mat = [Kse3+c0*Bse3+At,Gt;Gt',Kbt3+c0*Bbt1+H];
        
        us = 1/det(Mat)*(-Gt'*(-GammaV+LamdaN)+(Kse3+c0*Bse1+At)*(-GammaU+LamdaM));
        vs = 1/det(Mat)*((Kbt3+c0*Bbt1+H)*(-GammaV+LamdaN)-Gt*(-GammaU+LamdaM));
        vst = c0*vs + vsh{i,j};
        ust = c0*us + ush{i,j};

        
        
        ps = R*v;
        Rs = R*hat(u);
        ns = rho3*A3*R*(hat(w)*q + qt) -R*(a+At*vs+Gt*us)+R*C*q.*abs(q) - rho3*A3*g;
        ms = -R*(b+Gt'*vs+H*us)+rho3*R*(hat(w)*J3*w + J3*wt) - hat(ps)*n;
        qs = vt - hat(u)*q + hat(w)*v;
        ws = ut - hat(u)*w;
    end

   function [ps,Rs,ns,ms,qs,ws,vs, us,v,u,vt,ut,qt,wt,vst,ust] = dynamicODEN2(p,R,n,m,q,w)
        v = (Kse3 + c0*Bse3)\(R'*n + Kse3*vstar(ds*(j-1)) - Bse3*vh{i,j});
        u = (Kbt3 + c0*Bbt1)\(R'*m + Kbt3*ustar(ds*(j-1)) - Bbt1*uh{i,j});
        vt = c0*v + vh{i,j};
        ut = c0*u + uh{i,j};
        qt = c0*q + qh{i,j};
        wt = c0*w + wh{i,j};
                
        ptsb1 = hat(u)*rt11+v;
        At1 = -Tt11/norm(ptsb1)^3*hat(ptsb1)*hat(ptsb1);
        Gt1 = -At1*hat(rt11);
        H1 =hat(rt11)*Gt1;
        a1 = At1*(hat(u)*ptsb1);
        b1 = hat(rt11)*a1;
        
        
        ptsb2 = hat(u)*rt12+v;
        At2 = -Tt12/norm(ptsb2)^3*hat(ptsb2)*hat(ptsb2);
        Gt2 = -At2*hat(rt12);
        H2 =hat(rt12)*Gt2;
        a2 = At2*(hat(u)*ptsb2);
        b2 = hat(rt12)*a2;
        
        ptsb3 = hat(u)*rt21+v;
        At3 = -Tt21/norm(ptsb3)^3*hat(ptsb3)*hat(ptsb3);
        Gt3 = -At3*hat(rt21);
        H3 =hat(rt21)*Gt3;
        a3 = At3*(hat(u)*ptsb3);
        b3 = hat(rt21)*a3;
        
        
        ptsb4 = hat(u)*rt22+v;
        At4 = -Tt22/norm(ptsb4)^3*hat(ptsb4)*hat(ptsb4);
        Gt4 = -At4*hat(rt22);
        H4 =hat(rt22)*Gt4;
        a4 = At4*(hat(u)*ptsb4);
        b4 = hat(rt22)*a4; 
        
        a = a1 + a2 + a3 + a4;
        b = b1 + b2 + b3 + b4;
        At = At1 + At2 + At3 + At4;
        Gt = Gt1 + Gt2 + Gt3 + Gt4;
        H = H1 + H2 + H3 + H4;
        
              
        LamdaN = -a + rho3*A3*(hat(w)*q + qt)+ C*q.*abs(q)- R'*rho3*A3*g;
        LamdaM = -b + rho3*(hat(w)*J3*w + J3*wt) -hat (v)*(Kse3*(v-vstar(ds*(j-1)))+Bse1*vt);
        GammaV = hat(u)*(Kse3*(v-vstar(ds*(j-1)))+Bse3*vt)-Kse3*vsstar(ds*(j-1))+Bse3*vsh{i,j};
        GammaU = hat(u)*(Kbt3*(u-ustar(ds*(j-1)))+Bbt3*ut)-Kbt3*usstar(ds*(j-1))+Bbt3*ush{i,j};
        
        
        Mat = [Kse3+c0*Bse1+At,Gt;Gt',Kbt3+c0*Bbt1+H];
        
        us = 1/det(Mat)*(-Gt'*(-GammaV+LamdaN)+(Kse3+c0*Bse1+At)*(-GammaU+LamdaM));
        vs = 1/det(Mat)*((Kbt3+c0*Bbt1+H)*(-GammaV+LamdaN)-Gt*(-GammaU+LamdaM));
        vst = c0*vs + vsh{i,j};
        ust = c0*us + ush{i,j};

        pts1 = R*hat(u)*rt11+R*v;
        pts2 = R*hat(u)*rt12+R*v;
        
        n2 = -Tt11*pts1/norm(pts1)  -Tt12*pts2/norm(pts2);
        m2 = -Tt11*hat(R*rt11)*pts1/norm(pts1) -Tt12*hat(R*rt12)*pts2/norm(pts2);
        
        ps = R*v;
        Rs = R*hat(u);
        ns = rho3*A3*R*(hat(w)*q + qt) -R*(a+At*vs+Gt*us)+R*C*q.*abs(q) - rho3*A3*g;
        ms = -R*(b+Gt'*vs+H*us)+rho3*R*(hat(w)*J3*w + J3*wt) - hat(ps)*n;
        qs = vt - hat(u)*q + hat(w)*v;
        ws = ut - hat(u)*w;
   end

   function [ps,Rs,ns,ms,qs,ws,vs, us,v,u,vt,ut,qt,wt,vst,ust] = dynamicODE2(p,R,n,m,q,w)
        v = (Kse2 + c0*Bse2)\(R'*n + Kse2*vstar(ds*(j-1)) - Bse2*vh{i,j});
        u = (Kbt2 + c0*Bbt2)\(R'*m + Kbt2*ustar(ds*(j-1)) - Bbt2*uh{i,j});
        vt = c0*v + vh{i,j};
        ut = c0*u + uh{i,j};
        qt = c0*q + qh{i,j};
        wt = c0*w + wh{i,j};
                
        
        ptsb3 = hat(u)*rt21+v;
        At3 = -Tt21/norm(ptsb3)^3*hat(ptsb3)*hat(ptsb3);
        Gt3 = -At3*hat(rt21);
        H3 =hat(rt21)*Gt3;
        a3 = At3*(hat(u)*ptsb3);
        b3 = hat(rt21)*a3;
        
        
        ptsb4 = hat(u)*rt22+v;
        At4 = -Tt22/norm(ptsb4)^3*hat(ptsb4)*hat(ptsb4);
        Gt4 = -At4*hat(rt22);
        H4 =hat(rt22)*Gt4;
        a4 = At4*(hat(u)*ptsb4);
        b4 = hat(rt22)*a4; 
        
        a = a3 + a4;
        b = b3 + b4;
        At = At3 + At4;
        Gt = Gt3 + Gt4;
        H = H3 + H4;
        
              
        LamdaN = -a + rho2*A2*(hat(w)*q + qt)+ C2*q.*abs(q)- R'*rho2*A2*g;
        LamdaM = -b + rho2*(hat(w)*J2*w + J2*wt) -hat (v)*(Kse2*(v-vstar(ds*(j-1)))+Bse2*vt);
        GammaV = hat(u)*(Kse2*(v-vstar(ds*(j-1)))+Bse2*vt)-Kse2*vsstar(ds*(j-1))+Bse2*vsh{i,j};
        GammaU = hat(u)*(Kbt2*(u-ustar(ds*(j-1)))+Bbt2*ut)-Kbt2*usstar(ds*(j-1))+Bbt2*ush{i,j};
        

        
        Mat = [Kse2+c0*Bse2+At,Gt;Gt',Kbt2+c0*Bbt2+H];
        
        us = 1/det(Mat)*(-Gt'*(-GammaV+LamdaN)+(Kse2+c0*Bse2+At)*(-GammaU+LamdaM));
        vs = 1/det(Mat)*((Kbt2+c0*Bbt2+H)*(-GammaV+LamdaN)-Gt*(-GammaU+LamdaM));
        vst = c0*vs + vsh{i,j};
        ust = c0*us + ush{i,j};

        pts3 = R*hat(u)*rt21+R*v;
        pts4 = R*hat(u)*rt22+R*v;
        
        nLL = -Tt21*pts3/norm(pts3)  -Tt22*pts4/norm(pts4);
        mLL = -Tt21*hat(R*rt21)*pts3/norm(pts3) -Tt22*hat(R*rt22)*pts4/norm(pts4);
        
        ps = R*v;
        Rs = R*hat(u);
        ns = rho2*A2*R*(hat(w)*q + qt) -R*(a+At*vs+Gt*us)+R*C2*q.*abs(q) - rho2*A2*g;
        ms = -R*(b+Gt'*vs+H*us)+rho2*R*(hat(w)*J2*w + J2*wt) - hat(ps)*n;
        qs = vt - hat(u)*q + hat(w)*v;
        ws = ut - hat(u)*w;
    end

      
    function visualize()
        for j = 1 : N1,  x(j) = p{i,j}(1);  y(j) = p{i,j}(2); z(j) = p{i,j}(3);   end
        figure (1)
        fig = plot3(z,y,x); axis([-0.05*L 1.1*L  -0.1*L 0.1*L -0.05*L 0.2*L]);
        xlabel('z (m)');  ylabel('y (m)'); zlabel('x (m)')
        hold on; grid on;  drawnow;  pause(0.05);
      end
    function visualize2()
        for j = 1 : N2,  x(j) = p{i,j}(1);  y(j) = p{i,j}(2); z(j) = p{i,j}(3);   end
        figure (1)
        fig = plot3(z,y,x); axis([-0.05*L 1.1*L  -0.1*L 0.1*L -0.05*L 0.2*L]);
        xlabel('z (m)');  ylabel('y (m)'); zlabel('x (m)')
        hold on; grid on;  drawnow;  pause(0.05);
      end
    
    subplot(2,1,1)
    plot(U,X);
    xlabel('t (s)');  ylabel('x (m)'); title('Tip Displacement - X Component');
    subplot(2,1,2)
    plot(U,Y);
    xlabel('t (s)');  ylabel('y (m)'); title('Tip Displacement - Y Component');
    saveas(gcf, '../results/ConcentricTubeCR.png')
 end

