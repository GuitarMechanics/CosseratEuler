function CoManipulativeTypeCCR
    clear all
    hat=@(y)[0,-y(3),y(2);y(3),0,-y(1);-y(2),y(1),0];
    hat2=@(y)[y(3,2);y(1,3);y(2,1)];
    global p R j n m v u q w vs us vt ut qt wt vst ust vh uh vsh ush qh wh nLL mLL pt x y z x1 z1%Make vars available in whole program
    global p2 R2 n2 m2 v2 u2 q2 w2 vs2  us2 vst2 ust2 vt2 ut2 qt2 wt2 vh2 uh2 vsh2 ush2 qh2 wh2 nLL2 mLL2 pt2 xx yy zz Nc Mc q3 w3 xxx zzz
        %Parameters
    L = 0.4;                 %Length (before strain)
    L1 = 0.2062;
    L2 = 0.4;
    N = 40;                  %Spatial resolution
    N1 = 2;
    N2 = 40;
    E = 207e9;               %Young's modulus
    E1 = 207e11;             %Young's modulus-Object
    r = 0.001;               %Cross-section radius
    rt11 = [0;0;0.01];
    rt12 = [0;0.01;0];
    rt21 = [0;0;0.01];
    rt22 = [0;0.01;0];
    rho = 8000;             %Density
    rho1 = 500;             %Density-Object
    g = [-9.81;0;0];        %Gravity
    Bse = zeros(3);         %Material damping coefficients - shear and extension
    Bbt = 1e-6*eye(3);      %Material damping coefficients - bending and torsion
    C = 0.03*eye(3);        %Square-law-drag damping coefficients
    dt = 0.02;              %Time step
    alpha = -0.2;           %BDF-alpha parameter
    STEPS = 100;            %Number of timesteps to completion
    vstar = @(s)[1;0;0];      % Value of v when static and absent loading
    ustar = @(s)[0;-2.62;0];  % Precurvature
    vsstar = @(s)[0;0;0];
    usstar = @(s)[0;0;0];
    vstar1 = @(s)[0;0;1];     % Value of v when static and absent loading
    ustar1 = @(s)[0;0;0];     % Precurvature
    vsstar1 = @(s)[0;0;0];
    usstar1 = @(s)[0;0;0];
    
    %Boundary Conditions
    for i = 1 : STEPS
        p{i,1} = [0;0;-0.3]; %Clamped base
        R{i,1} = eye(3);
        q{i,1} = [0;0;0];
        w{i,1} = [0;0;0];
    end
    vstar2 = @(s)[1;0;0]; %Value of v when static and absent loading
    ustar2 = @(s)[0;2.62;0]; %Precurvature
    vsstar2 = @(s)[0;0;0];
    usstar2 = @(s)[0;0;0];
    %Boundary Conditions
    for i = 1 : STEPS
        p2{i,1} = [0;0;0.3]; %Clamped base
        R2{i,1} = eye(3);
        q2{i,1} = [0;0;0];
        w2{i,1} = [0;0;0];
    end
    
    
    
    
    nL = 0.0*g;  %Start with a weight hung at the tip
    mL = [0;0;0];

    %Dependent Parameter Calculations
    A = pi*r^2;                                 %Cross-sectional area
    J = diag([pi*r^4/4  pi*r^4/4  pi*r^4/2]);   %Inertia
    G = E/( 2*(1+0.3) );                        %Shear modulus
    G1 = E1/( 2*(1+0.3) );  
%   G = 36600; 
    Kse = diag([G*A, G*A, E*A]);                %Stiffness matrix - shear and extension
    Kse1 = diag([G1*A, G1*A, E1*A]);
    Kbt = diag([E*J(1,1), E*J(2,2), G*J(3,3)]); %Stiffness matrix - bending and torsion
    Kbt1 = diag([E1*J(1,1), E1*J(2,2), G1*J(3,3)]);
    ds = L/(N-1);                               %Grid distance (before strain)
    ds1 = L1/2/(N1-1);
    ds2 = L2/(N2-1); 
    c0 = (1.5 + alpha) / ( dt*(1+alpha) );      %BDF-alpha coefficients
    c1 = -2/dt;
    c2 = (0.5 + alpha) / ( dt*(1+alpha) );
    d1 = alpha / (1+alpha);
    
    %Main Simulation
    theta = pi/2.872;
    Rrot = [cos(theta) 0 -sin(theta);0 1 0; sin(theta) 0 cos(theta)];
    i = 1;
    fsolve(@staticIVP, zeros(6,1)); %Solve static BVP w/ shooting method
    applyStaticBDFalpha();
    fsolve(@staticIVP2, zeros(6,1));
    applyStaticBDFalpha2();
%     visualize1();
  
    Nc{1} = [0;0;0];
    Mc{1} = [0;0;0];
    for i = 2 : STEPS
        Tt=0;
        Tt11=0;
        Tt12=0;
        Tt21=4;
        Tt22=1;
        
        fsolve(@Errorfc,[Nc{i-1} ; Mc{i-1}]);
%         visualize();
    end
%         
 %             
    function E2 = Errorfc(G)
        Nc{i} = G(1:3);
        Mc{i} = G(4:6);
       
        fsolve(@dynamicIVP, [n{i-1,1}; m{i-1,1}]); %Solve semi-discretized PDE w/ shooting
        applyDynamicBDFalpha();
        fsolve(@dynamicIVP2, [n2{i-1,1}; m2{i-1,1}]); %Solve semi-discretized PDE w/ shooting
        applyDynamicBDFalpha2();
        E2 = [(p{i,N+N1} - p2{i,N2}); hat2((R2{i,N2}*R{i,N+N1}' - R2{1,N2}*R{1,N+N1}')-(R2{i,N2}*R{i,N+N1}' - R2{1,N2}*R{1,N+N1}')')];
     end
    
    %Function Definitions
    function applyStaticBDFalpha()
        for j = 1 : N+N1-1
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
    function applyStaticBDFalpha2()
        for j = 1 : N2-1
            vh2{i+1,j} = (c1+c2)*v2{i,j};
            uh2{i+1,j} = (c1+c2)*u2{i,j};
            vsh2{i+1,j} = (c1+c2)*vs2{i,j};
            ush2{i+1,j} = (c1+c2)*us2{i,j};
            qh2{i+1,j} = [0;0;0];
            wh2{i+1,j} = [0;0;0];
            q2{i,j} = [0;0;0];
            w2{i,j} = [0;0;0];
        end
    end

    function applyDynamicBDFalpha()
        for j = 1 : N+N1-1        
            vh{i+1,j} = c1*v{i,j} + c2*v{i-1,j} + d1*vt{i,j};
            uh{i+1,j} = c1*u{i,j} + c2*u{i-1,j} + d1*ut{i,j};
            vsh{i+1,j} = c1*vs{i,j} + c2*vs{i-1,j} + d1*vst{i,j};
            ush{i+1,j} = c1*us{i,j} + c2*us{i-1,j} + d1*ust{i,j};
            qh{i+1,j} = c1*q{i,j} + c2*q{i-1,j} + d1*qt{i,j};
            wh{i+1,j} = c1*w{i,j} + c2*w{i-1,j} + d1*wt{i,j};
        end
    end
    function applyDynamicBDFalpha2()
        for j = 1 : N2-1        
            vh2{i+1,j} = c1*v2{i,j} + c2*v2{i-1,j} + d1*vt2{i,j};
            uh2{i+1,j} = c1*u2{i,j} + c2*u2{i-1,j} + d1*ut2{i,j};
            vsh2{i+1,j} = c1*vs2{i,j} + c2*vs2{i-1,j} + d1*vst2{i,j};
            ush2{i+1,j} = c1*us2{i,j} + c2*us2{i-1,j} + d1*ust2{i,j};
            qh2{i+1,j} = c1*q2{i,j} + c2*q2{i-1,j} + d1*qt2{i,j};
            wh2{i+1,j} = c1*w2{i,j} + c2*w2{i-1,j} + d1*wt2{i,j};
        end
    end

    function E = staticIVP(G)
        n{i,1} = G(1:3);
        m{i,1} = G(4:6);

        %Euler's method
        for j = 1 : N-2
            [ps, Rs, ns, ms, us{i,j}, vs{i,j} ,v{i,j}, u{i,j}] = staticODE(p{i,j},R{i,j},n{i,j},m{i,j});
            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} + ds*ns;
            m{i,j+1} = m{i,j} + ds*ms;
                        
        end
        for j = N-1 : N-1
            [ps, Rs, ns, ms, us{i,j}, vs{i,j} ,v{i,j}, u{i,j}] = staticODE(p{i,j},R{i,j},n{i,j},m{i,j});
            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} + ds*ns ;
            m{i,j+1} = m{i,j} + ds*ms;
                        
        end
        for j = N : N+N1-1
            [ps, Rs, ns, ms, us{i,j}, vs{i,j} ,v{i,j}, u{i,j}] = staticODE1(p{i,j},R{i,j},n{i,j},m{i,j});
            p{i,j+1} = p{i,j} + ds1*ps;
            R{i,j+1} = (R{i,j} + ds1*Rs) ;
            n{i,j+1} = n{i,j} + ds1*ns;
            m{i,j+1} = m{i,j} + ds1*ms;
                        
        end
        
        E = [ n{i,N+N1} - nL;  m{i,N+N1} - mL ];
    end
    function E = staticIVP2(G)
        n2{i,1} = G(1:3);
        m2{i,1} = G(4:6);

        %Euler's method
        for j = 1 : N2-1
            [ps2, Rs2, ns2, ms2, us2{i,j}, vs2{i,j}, v2{i,j}, u2{i,j}] = staticODE2(p2{i,j},R2{i,j},n2{i,j},m2{i,j});
            p2{i,j+1} = p2{i,j} + ds2*ps2;
            R2{i,j+1} = R2{i,j} + ds2*Rs2;
            n2{i,j+1} = n2{i,j} + ds2*ns2;
            m2{i,j+1} = m2{i,j} + ds2*ms2;
        end
        E = [ n2{i,N2} - nL;  m2{i,N2} - mL ];
    end

    function E = dynamicIVP(G)
        n{i,1} = G(1:3);
        m{i,1} = G(4:6);

        %Euler's method
        for j = 1 : N-2
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
        for j = N-1 : N-1
            [ps, Rs, ns, ms, qs, ws, vs{i,j}, us{i,j},...
                 v{i,j}, u{i,j}, vt{i,j}, ut{i,j},...
                 qt{i,j}, wt{i,j},vst{i,j}, ust{i,j}] = dynamicODEN2(p{i,j},R{i,j},n{i,j},m{i,j},q{i,j},w{i,j});
            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} - nLL + ds*ns;
            m{i,j+1} = m{i,j} - mLL + ds*ms;
            q{i,j+1} = q{i,j} + ds*qs;
            w{i,j+1} = w{i,j} + ds*ws;
            
            
        end
        for j = N : N+N1-1
            [ps, Rs, ns, ms, qs, ws, vs{i,j}, us{i,j},...
                 v{i,j}, u{i,j}, vt{i,j}, ut{i,j},...
                 qt{i,j}, wt{i,j},vst{i,j}, ust{i,j}] = dynamicODE2(p{i,j},R{i,j},n{i,j},m{i,j},q{i,j},w{i,j});
            p{i,j+1} = p{i,j} + ds1*ps;
            R{i,j+1} = R{i,j} + ds1*Rs;
            n{i,j+1} = n{i,j} + ds1*ns;
            m{i,j+1} = m{i,j} + ds1*ms;
            q{i,j+1} = q{i,j} + ds1*qs;
            w{i,j+1} = w{i,j} + ds1*ws;
            
            
        end
        
        
        E = [n{i,N+N1} - Nc{i} ;  m{i,N+N1} - Mc{i}];
        
    end
    function E = dynamicIVP2(G)
        n2{i,1} = G(1:3);
        m2{i,1} = G(4:6);

        %Euler's method
        for j = 1 : N2-2
            [ps2, Rs2, ns2, ms2, qs2, ws2, vs2{i,j}, us2{i,j},...
                 v2{i,j}, u2{i,j}, vt2{i,j}, ut2{i,j},...
                 qt2{i,j}, wt2{i,j},vst2{i,j}, ust2{i,j}] = dynamicODE3(p2{i,j},R2{i,j},n2{i,j},m2{i,j},q2{i,j},w2{i,j});
            p2{i,j+1} = p2{i,j} + ds2*ps2;
            R2{i,j+1} = R2{i,j} + ds2*Rs2;
            n2{i,j+1} = n2{i,j} + ds2*ns2;
            m2{i,j+1} = m2{i,j} + ds2*ms2;
            q2{i,j+1} = q2{i,j} + ds2*qs2;
            w2{i,j+1} = w2{i,j} + ds2*ws2;
            
            
        end
        for j = N2-1 : N2-1
            [ps2, Rs2, ns2, ms2, qs2, ws2, vs2{i,j}, us2{i,j},...
                 v2{i,j}, u2{i,j}, vt2{i,j}, ut2{i,j},...
                 qt2{i,j}, wt2{i,j},vst2{i,j}, ust2{i,j}] = dynamicODE3(p2{i,j},R2{i,j},n2{i,j},m2{i,j},q2{i,j},w2{i,j});
            p2{i,j+1} = p2{i,j} + ds2*ps2;
            R2{i,j+1} = R2{i,j} + ds2*Rs2;
            n2{i,j+1} = n2{i,j} + ds2*ns2;
            m2{i,j+1} = m2{i,j} + ds2*ms2;
            q2{i,j+1} = q2{i,j} + ds2*qs2;
            w2{i,j+1} = w2{i,j} + ds2*ws2;
            
            
        end
      E = [n2{i,N2} - nLL2 + Nc{i} ;  m2{i,N2} - mLL2 + Mc{i} ];
        
    end



    function [ps, Rs, ns, ms, vs, us, v, u] = staticODE(p,R,n,m)
        
        v = Kse\R'*n + vstar(ds*(j-1));
        u = Kbt\R'*m + ustar(ds*(j-1));
        
        Tt11 = 0;
        Tt12 = 0;
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
        At = At1 + At2;
        B = B1 + B2;
        Gt = Gt1 + Gt2;
        H = H1 +H2;
        
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

    function [ps, Rs, ns, ms, vs, us, v, u] = staticODE1(p,R,n,m)
        R1 = R * [cos(theta) 0 -sin(theta);0 1 0; sin(theta) 0 cos(theta)]';
        v = Kse1\R1'*n + vstar1(ds1*(j-1));
        u = Kbt1\R1'*m + ustar1(ds1*(j-1));
        
        Tt11 = 0;
        Tt12 = 0;
        
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
        At = At1 + At2;
        B = B1 + B2;
        Gt = Gt1 + Gt2;
        H = H1 +H2;
        
        
        d = Kse1*vsstar1(ds1*(j-1))-hat(u)*Kse1*(v-vstar1(ds1*(j-1)))-a-R1'*rho1*A*g;
        c = Kbt1*usstar1(ds1*(j-1))-hat(u)*Kbt1*(u-ustar1(ds1*(j-1)))-hat(v)*Kse1*(v-vstar1(ds1*(j-1)))-b;
        
        
        Mat = [Kse1+At,Gt;B,Kbt1+H];
        vs = 1/det(Mat)*((Kbt1+H)*d-Gt*c);
        us = 1/det(Mat)*(-B*d+(Kse1+At)*c);
        
         
        
        ps = R1*v;
        Rs = R1*hat(u);
        ns = -rho1*A*g;
        ms = -hat(ps)*n;
    end
    function [ps2, Rs2, ns2, ms2,  vs2, us2, v2, u2] = staticODE2(p2,R2,n2,m2)
        v2 = Kse\R2'*n2 + vstar2(ds2*(j-1));
        u2 = Kbt\R2'*m2 + ustar2(ds2*(j-1));

        Tt21 = 0;
        Tt22 = 0;
        
        ptsb3 = hat(u2)*rt21+v2;
        At3 = -Tt21/norm(ptsb3)^3*hat(ptsb3)*hat(ptsb3);
        B3 = hat(rt21)*At3;
        Gt3 = -At3*hat(rt21);
        H3 =-B3*hat(rt21);
        a3 = At3*(hat(u2)*ptsb3);
        b3 = hat(rt21)*a3;
        
        ptsb4 = hat(u2)*rt22+v2;
        At4 = -Tt22/norm(ptsb4)^3*hat(ptsb4)*hat(ptsb4);
        B4 = hat(rt22)*At4;
        Gt4 = -At4*hat(rt22);
        H4 =-B4*hat(rt22);
        a4 = At4*(hat(u2)*ptsb4);
        b4 = hat(rt22)*a4;
        
        a = a3 + a4;
        b = b3 + b4;
        At = At3 + At4;
        B = B3 + B4;
        Gt = Gt3 + Gt4;
        H = H3 +H4;
        
        d = Kse*vsstar2(ds2*(j-1))-hat(u2)*Kse*(v2-vstar2(ds2*(j-1)))-a-R2'*rho*A*g;
        c = Kbt*usstar2(ds2*(j-1))-hat(u2)*Kbt*(u2-ustar2(ds2*(j-1)))-hat(v2)*Kse*(v2-vstar2(ds2*(j-1)))-b;

        
        Mat = [Kse+At,Gt;B,Kbt+H];
        vs2 = 1/det(Mat)*((Kbt+H)*d-Gt*c);
        us2 = 1/det(Mat)*(-B*d+(Kse+At)*c);
        
        ps2 = R2*v2;
        Rs2 = R2*hat(u2);
        ns2 = -rho*A*g;
        ms2 = -hat(ps2)*n2;
    end




    function [ps,Rs,ns,ms,qs,ws,vs, us,v,u,vt,ut,qt,wt,vst,ust] = dynamicODE1(p,R,n,m,q,w)
        v = (Kse + c0*Bse)\(R'*n + Kse*vstar(ds*(j-1)) - Bse*vh{i,j});
        u = (Kbt + c0*Bbt)\(R'*m + Kbt*ustar(ds*(j-1)) - Bbt*uh{i,j});
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
        
        
        LamdaN = -a + rho*A*(hat(w)*q + qt)+ C*q.*abs(q)- R'*rho*A*g;
        LamdaM = -b + rho*(hat(w)*J*w + J*wt) -hat (v)*(Kse*(v-vstar(ds*(j-1)))+Bse*vt);
        GammaV = hat(u)*(R'*n)-Kse*vsstar(ds*(j-1))+Bse*vsh{i,j};
        GammaU = hat(u)*(R'*m)-Kbt*usstar(ds*(j-1))+Bbt*ush{i,j};

        
        Mat = [Kse+c0*Bse+At,Gt;Gt',Kbt+c0*Bbt+H];
        
        us = 1/det(Mat)*(-Gt'*(-GammaV+LamdaN)+(Kse+c0*Bse+At)*(-GammaU+LamdaM));
        vs = 1/det(Mat)*((Kbt+c0*Bbt+H)*(-GammaV+LamdaN)-Gt*(-GammaU+LamdaM));
        vst = c0*vs + vsh{i,j};
        ust = c0*us + ush{i,j};
        
        ps = R*v;
        Rs = R*hat(u);
        ns = rho*A*R*(hat(w)*q + qt) -R*(a+At*vs+Gt*us)+R*C*q.*abs(q) - rho*A*g;
        ms = -R*(b+Gt'*vs+H*us)+rho*R*(hat(w)*J*w + J*wt) - hat(ps)*n;
        qs = vt - hat(u)*q + hat(w)*v;
        ws = ut - hat(u)*w;
    end
    function [ps,Rs,ns,ms,qs,ws,vs, us,v,u,vt,ut,qt,wt,vst,ust] = dynamicODEN2(p,R,n,m,q,w)
        v = (Kse + c0*Bse)\(R'*n + Kse*vstar(ds*(j-1)) - Bse*vh{i,j});
        u = (Kbt + c0*Bbt)\(R'*m + Kbt*ustar(ds*(j-1)) - Bbt*uh{i,j});
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
        
        
        LamdaN = -a + rho*A*(hat(w)*q + qt)+ C*q.*abs(q)- R'*rho*A*g;
        LamdaM = -b + rho*(hat(w)*J*w + J*wt) -hat (v)*(Kse*(v-vstar(ds*(j-1)))+Bse*vt);
        GammaV = hat(u)*(R'*n)-Kse*vsstar(ds*(j-1))+Bse*vsh{i,j};
        GammaU = hat(u)*(R'*m)-Kbt*usstar(ds*(j-1))+Bbt*ush{i,j};

        
        Mat = [Kse+c0*Bse+At,Gt;Gt',Kbt+c0*Bbt+H];
        
        us = 1/det(Mat)*(-Gt'*(-GammaV+LamdaN)+(Kse+c0*Bse+At)*(-GammaU+LamdaM));
        vs = 1/det(Mat)*((Kbt+c0*Bbt+H)*(-GammaV+LamdaN)-Gt*(-GammaU+LamdaM));
        vst = c0*vs + vsh{i,j};
        ust = c0*us + ush{i,j};

        pts1 = R*hat(u)*rt11+R*v;
        pts2 = R*hat(u)*rt12+R*v;
        nLL = -Tt11*pts1/norm(pts1)-Tt12*pts2/norm(pts2);
        mLL = -Tt11*hat(R*rt11)*pts1/norm(pts1)-Tt12*hat(R*rt12)*pts2/norm(pts2);
        
        ps = R*v;
        Rs = R*hat(u);
        ns = rho*A*R*(hat(w)*q + qt) -R*(a+At*vs+Gt*us)+R*C*q.*abs(q) - rho*A*g;
        ms = -R*(b+Gt'*vs+H*us)+rho*R*(hat(w)*J*w + J*wt) - hat(ps)*n;
        qs = vt - hat(u)*q + hat(w)*v;
        ws = ut - hat(u)*w;
    end
    function [ps,Rs,ns,ms,qs,ws,vs, us,v,u,vt,ut,qt,wt,vst,ust] = dynamicODE2(p,R,n,m,q,w)
        R1 = R * Rrot';
        v = (Kse1 + c0*Bse)\(R1'*n + Kse1*vstar1(ds1*(j-1))- Bse*vh{i,j} );
        u = (Kbt1 + c0*Bbt)\(R1'*m + Kbt1*ustar1(ds1*(j-1)) - Bbt*uh{i,j});
        vt = c0*v + vh{i,j};
        ut = c0*u + uh{i,j};
        qt = c0*q + qh{i,j};
        wt = c0*w + wh{i,j};
                
        ptsb = hat(u)*rt11+v;
         Ttt = 0;
        At = -Ttt/norm(ptsb)^3*hat(ptsb)*hat(ptsb);
        Gt = -At*hat(rt11);
        H =hat(rt11)*Gt;
        a = At*(hat(u)*ptsb);
        b = hat(rt11)*a;
        LamdaN = -a + rho1*A*(hat(w)*q + qt)+ C*q.*abs(q)- R1'*rho1*A*g;
        LamdaM = -b + rho1*(hat(w)*J*w + J*wt) -hat (v)*(Kse*(v-vstar1(ds1*(j-1)))+Bse*vt);
        GammaV = hat(u)*(R1'*n)-Kse1*vsstar1(ds1*(j-1))+Bse*vsh{i,j};
        GammaU = hat(u)*(R1'*m)-Kbt1*usstar1(ds1*(j-1))+Bbt*ush{i,j};

        
        Mat = [Kse+c0*Bse+At,Gt;Gt',Kbt+c0*Bbt+H];
        
        us = 1/det(Mat)*(-Gt'*(-GammaV+LamdaN)+(Kse1+c0*Bse+At)*(-GammaU+LamdaM));
        vs = 1/det(Mat)*((Kbt1+c0*Bbt+H)*(-GammaV+LamdaN)-Gt*(-GammaU+LamdaM));
        vst = c0*vs + vsh{i,j};
        ust = c0*us + ush{i,j};

      
        ps = R1*v;
        Rs = R1*hat(u);
        ns = rho1*A*R1*(hat(w)*q + qt) -R1*(a+At*vs+Gt*us)+R1*C*q.*abs(q) - rho1*A*g;
        ms = -R1*(b+Gt'*vs+H*us)+rho1*R1*(hat(w)*J*w + J*wt) - hat(ps)*n;
        qs = vt - hat(u)*q + hat(w)*v;
        ws = ut - hat(u)*w;
    end



    function [ps2,Rs2,ns2,ms2,qs2,ws2,vs2, us2,v2,u2,vt2,ut2,qt2,wt2,vst2,ust2] = dynamicODE3(p2,R2,n2,m2,q2,w2)
        v2 = (Kse + c0*Bse)\(R2'*n2 + Kse*vstar2(ds2*(j-1)) - Bse*vh2{i,j});
        u2 = (Kbt + c0*Bbt)\(R2'*m2 + Kbt*ustar2(ds2*(j-1)) - Bbt*uh2{i,j});
        vt2 = c0*v2 + vh2{i,j};
        ut2 = c0*u2 + uh2{i,j};
        qt2 = c0*q2 + qh2{i,j};
        wt2 = c0*w2 + wh2{i,j};
                
        ptsb3 = hat(u2)*rt21+v2;
        At3 = -Tt21/norm(ptsb3)^3*hat(ptsb3)*hat(ptsb3);
        Gt3 = -At3*hat(rt21);
        H3 =hat(rt21)*Gt3;
        a3 = At3*(hat(u2)*ptsb3);
        b3 = hat(rt21)*a3;
        
        
        ptsb4 = hat(u2)*rt22+v2;
        At4 = -Tt22/norm(ptsb4)^3*hat(ptsb4)*hat(ptsb4);
        Gt4 = -At4*hat(rt22);
        H4 =hat(rt22)*Gt4;
        a4 = At4*(hat(u2)*ptsb4);
        b4 = hat(rt22)*a4; 
        
        a = a3 + a4;
        b = b3 + b4;
        At = At3 + At4;
        Gt = Gt3 + Gt4;
        H = H3 + H4;
        
              
        LamdaN = -a + rho*A*(hat(w2)*q2 + qt2)+ C*q2.*abs(q2)- R2'*rho*A*g;
        LamdaM = -b + rho*(hat(w2)*J*w2 + J*wt2) -hat (v2)*(Kse*(v2-vstar2(ds2*(j-1)))+Bse*vt2);
        GammaV = hat(u2)*(R2'*n2)-Kse*vsstar2(ds2*(j-1))+Bse*vsh2{i,j};
        GammaU = hat(u2)*(R2'*m2)-Kbt*usstar2(ds2*(j-1))+Bbt*ush2{i,j};

        
        Mat = [Kse+c0*Bse+At,Gt;Gt',Kbt+c0*Bbt+H];
        
        us2 = 1/det(Mat)*(-Gt'*(-GammaV+LamdaN)+(Kse+c0*Bse+At)*(-GammaU+LamdaM));
        vs2 = 1/det(Mat)*((Kbt+c0*Bbt+H)*(-GammaV+LamdaN)-Gt*(-GammaU+LamdaM));
        vst2 = c0*vs2 + vsh2{i,j};
        ust2 = c0*us2 + ush2{i,j};

        pts3 = R2*hat(u2)*rt21+R2*v2;
        pts4 = R2*hat(u2)*rt22+R2*v2;
           
        
        nLL2 = -Tt21*pts3/norm(pts3)-Tt22*pts4/norm(pts4);
        mLL2 = -Tt21*hat(R2*rt21)*pts3/norm(pts3)-Tt22*hat(R2*rt22)*pts4/norm(pts4);
        
        ps2 = R2*v2;
        Rs2 = R2*hat(u2);
        ns2 = rho*A*R2*(hat(w2)*q2 + qt2) -R2*(a+At*vs2+Gt*us2)+R2*C*q2.*abs(q2) - rho*A*g;
        ms2 = -R2*(b+Gt'*vs2+H*us2)+rho*R2*(hat(w2)*J*w2 + J*wt2) - hat(ps2)*n2;
        qs2 = vt2 - hat(u2)*q2 + hat(w2)*v2;
        ws2 = ut2 - hat(u2)*w2;
    end

    function visualize1()
        for j = 1 : N+N1,  x(j) = p{i,j}(1); y(j) = p{i,j}(2);   z(j) = p{i,j}(3);   end
        for j = 1 : N2,  xx(j) = p2{i,j}(1); yy(j) = p2{i,j}(2);  zz(j) = p2{i,j}(3);   end
        plot3(z,y,x,zz,yy,xx);  
        str = sprintf('Iteration %f', i);
        title(str);  xlabel('z (m)');  ylabel('y (m)'); zlabel('x (m)');
        hold on; grid on;  drawnow;  pause(0.05); 
    end

    function visualize()
        for j = 1 : N+N1,  x(j) = p{i,j}(1);  y(j) = p{i,j}(2); z(j) = p{i,j}(3);   end
        for j = 1 : N2,  xx(j) = p2{i,j}(1);  yy(j) = p2{i,j}(2); zz(j) = p2{i,j}(3);   end
         plot3(z,y,x,zz,yy,xx);  
        str = sprintf('Iteration %f', i);
        title(str);  xlabel('z (m)');  ylabel('y (m)'); zlabel('x (m)');
        hold on; grid on;  drawnow; pause(0.05); 
    end

    for i = 1 : STEPS, U(i) = i*dt; X(i)= (p{i,N}(3) + p2{i,N2}(3))/2-(p{1,N}(3) + p2{1,N2}(3))/2; y(i)= (p{i,N}(2) + p2{i,N2}(2))/2-(p{1,N}(2) + p2{1,N2}(2))/2; Z(i)= (p{i,N}(1) + p2{i,N2}(1))/2-(p{1,N}(1) + p2{1,N2}(1))/2;end
    subplot(2,1,1)
    plot(U,X);
    xlabel('t (s)');  ylabel('x (m)'); title('Tip Displacement - X Component');
    subplot(2,1,2)
    plot(U,y);
    xlabel('t (s)');  ylabel('y (m)'); title('Tip Displacement - Y Component');
    saveas(gcf, '../results/CoManipulativeTypeCCR.png')
end

