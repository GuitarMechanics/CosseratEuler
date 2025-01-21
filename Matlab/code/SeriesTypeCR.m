function SeriesTypeCR
    clear all
    hat=@(y)[0,-y(3),y(2);y(3),0,-y(1);-y(2),y(1),0];
    global p R j n m v u q w vs us vt ut qt wt vst ust vh uh vsh ush qh wh nLL mLL n2 m2 pt x z %Make vars available in whole program
    %Parameters
    L = 0.4;                    %Length (before strain)
    N = 40;                     %Spatial resolution
    E = 207e9;                  %Young's modulus
    E2 = 157e9;
    r = 0.001;
    r2 = 0.0009;                %Cross-section radius
    rt11 = [0.01;0;0];          % tendon offset for the first section (near the base)
    rt12 = [0;0.01;0];          % tendon offset for the first section (near the base)
    rt21 = [0.005;0;0];         % tendon offset for the second section (near the tip)
    rt22 = [0;0.005;0];         % tendon offset for the second section (near the tip)
    rho = 8000;
    rho2 = 7500;                %Density
    g = [-9.81;0;0];            %Gravity
    Bse = zeros(3);
    Bse2 = zeros(3);            %Material damping coefficients - shear and extension
    Bbt = 1e-6*eye(3);          %Material damping coefficients - bending and torsion
    Bbt2 = 1e-6*eye(3);
    C = 0.03*eye(3);
    C2 = 0.03*eye(3);           %Square-law-drag damping coefficients
    dt = 0.015;                 %Time step
    alpha = -0.2;               %BDF-alpha parameter
    STEPS = 500;                %Number of timesteps to completion
    vstar = @(s)[0;0;1];        %Value of v when static and absent loading
    ustar = @(s)[0;0;0];        %Precurvature
    vsstar = @(s)[0;0;1]
    usstar = @(s)[0;0;0]
    %Boundary Conditions
    for i = 1 : STEPS
        p{i,1} = [0;0;0];       %Clamped base
        R{i,1} = eye(3);
        q{i,1} = [0;0;0];
        w{i,1} = [0;0;0];
    end
    nL = 0.0*g;                 %Start with a weight hung at the tip
    mL = [0;0;0];

    %Dependent Parameter Calculations
    A = pi*r^2;                                         %Cross-sectional area
    A2 = pi*r2^2;                                       %Cross-sectional area
    J = diag([pi*r^4/4  pi*r^4/4  pi*r^4/2]);           %Inertia
    J2 = diag([pi*r2^4/4  pi*r2^4/4  pi*r2^4/2]);       %Inertia
    G = E/( 2*(1+0.3) );                                %Shear modulus
    G2 = E2/( 2*(1+0.3) );                              %Shear modulus
    Kse = diag([G*A, G*A, E*A]);                        %Stiffness matrix - shear and extension
    Kse2 = diag([G2*A2, G2*A2, E2*A2]);                 %Stiffness matrix - shear and extension
    Kbt = diag([E*J(1,1), E*J(2,2), G*J(3,3)]);         %Stiffness matrix - bending and torsion
    Kbt2 = diag([E2*J2(1,1), E2*J2(2,2), G2*J2(3,3)]);  %Stiffness matrix - bending and torsion
    ds = L/(N-1);                                       %Grid distance (before strain)
    c0 = (1.5 + alpha) / ( dt*(1+alpha) );              %BDF-alpha coefficients
    c1 = -2/dt;
    c2 = (0.5 + alpha) / ( dt*(1+alpha) );
    d1 = alpha / (1+alpha);
    
    %Main Simulation
    i = 1;
    fsolve(@staticIVP, zeros(6,1));                     %Solve static BVP w/ shooting method
    applyStaticBDFalpha();
    visualize();
    
    for i = 2 : STEPS
        if i < 5
            Tt11 = 0;   
            Tt12 = 0;
            Tt21 = 0;   
            Tt22 = 0;
        else
            if i < STEPS/2
                Tt11 = 5;   
                Tt12 = 2;
                Tt21 = 0;  
                Tt22 = 0;
            else
                Tt11 = 5;   
                Tt12 = 2;
                Tt21 = 5;  
                Tt22 = 2;
            end
                
        end
            
        fsolve(@dynamicIVP, [n{i-1,1}; m{i-1,1}]); %Solve semi-discretized PDE w/ shooting
        applyDynamicBDFalpha();
        visualize();
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
        for j = 1 : N/2 - 2
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
       for j = N/2-1 :N/2-1 
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
        for j = N/2 : N - 1
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
        
        
        
        E = [n{i,N} - nLL ;  m{i,N} - mLL];
        
    end

    function [ps, Rs, ns, ms, vs, us, v, u] = staticODE(p,R,n,m)
        v = Kse\R'*n + vstar(ds*(j-1));
        u = Kbt\R'*m + ustar(ds*(j-1));
        
        ptsb = hat(u)*rt11+v;
        Tt = 0;
        At = -Tt/norm(ptsb)^3*hat(ptsb)*hat(ptsb);
        B = hat(rt11)*At;
        Gt = -At*hat(rt11);
        H =-B*hat(rt11);
        a = At*(hat(u)*ptsb);
        b = hat(rt11)*a;
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

    function [ps,Rs,ns,ms,qs,ws,vs, us,v,u,vt,ut,qt,wt,vst,ust] = dynamicODE(p,R,n,m,q,w)
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
        
              
        LamdaN = -a + rho*A*(hat(w)*q + qt)+ C*q.*abs(q)- R'*rho*A*g;
        LamdaM = -b + rho*(hat(w)*J*w + J*wt) -hat (v)*(Kse*(v-vstar(ds*(j-1)))+Bse*vt);
        GammaV = hat(u)*(Kse*(v-vstar(ds*(j-1)))+Bse*vt)-Kse*vsstar(ds*(j-1))+Bse*vsh{i,j};
        GammaU = hat(u)*(Kbt*(u-ustar(ds*(j-1)))+Bbt*ut)-Kbt*usstar(ds*(j-1))+Bbt*ush{i,j};


        
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
        
               
        LamdaN = -a + rho*A*(hat(w)*q + qt)+ C*q.*abs(q)- R'*rho*A*g;
        LamdaM = -b + rho*(hat(w)*J*w + J*wt) -hat (v)*(Kse*(v-vstar(ds*(j-1)))+Bse*vt);
        GammaV = hat(u)*(Kse*(v-vstar(ds*(j-1)))+Bse*vt)-Kse*vsstar(ds*(j-1))+Bse*vsh{i,j};
        GammaU = hat(u)*(Kbt*(u-ustar(ds*(j-1)))+Bbt*ut)-Kbt*usstar(ds*(j-1))+Bbt*ush{i,j};


        
        Mat = [Kse+c0*Bse+At,Gt;Gt',Kbt+c0*Bbt+H];
        
        us = 1/det(Mat)*(-Gt'*(-GammaV+LamdaN)+(Kse+c0*Bse+At)*(-GammaU+LamdaM));
        vs = 1/det(Mat)*((Kbt+c0*Bbt+H)*(-GammaV+LamdaN)-Gt*(-GammaU+LamdaM));
        vst = c0*vs + vsh{i,j};
        ust = c0*us + ush{i,j};

        pts1 = R*hat(u)*rt11+R*v;
        pts2 = R*hat(u)*rt12+R*v;
             
        
        n2 = -Tt11*pts1/norm(pts1)  -Tt12*pts2/norm(pts2);
        m2 = -Tt11*hat(R*rt11)*pts1/norm(pts1) -Tt12*hat(R*rt12)*pts2/norm(pts2);
        
        ps = R*v;
        Rs = R*hat(u);
        ns = rho*A*R*(hat(w)*q + qt) -R*(a+At*vs+Gt*us)+R*C*q.*abs(q) - rho*A*g;
        ms = -R*(b+Gt'*vs+H*us)+rho*R*(hat(w)*J*w + J*wt) - hat(ps)*n;
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
        LamdaM = -b + rho2*(hat(w)*J*w + J*wt) -hat (v)*(Kse2*(v-vstar(ds*(j-1)))+Bse2*vt);
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
        mLL = -Tt21*hat(R*rt21)*pts3/norm(pts3)  -Tt22*hat(R*rt22)*pts4/norm(pts4);
        
        ps = R*v;
        Rs = R*hat(u);
        ns = rho2*A2*R*(hat(w)*q + qt) -R*(a+At*vs+Gt*us)+R*C2*q.*abs(q) - rho2*A2*g;
        ms = -R*(b+Gt'*vs+H*us)+rho2*R*(hat(w)*J2*w + J2*wt) - hat(ps)*n;
        qs = vt - hat(u)*q + hat(w)*v;
        ws = ut - hat(u)*w;
   end

      function visualize()
        if rem(i,10) == 0; 
        for j = 1 : N,  x(j) = p{i,j}(1);  y(j) = p{i,j}(2); z(j) = p{i,j}(3);   end
        figure (1)
        fig = plot3(z,y,x); axis([-0.05*L 1.1*L  -0.1*L 0.1*L -0.05*L 0.2*L]);
        xlabel('z (m)');  ylabel('y (m)'); zlabel('x (m)')
        hold on; grid on;  drawnow;  pause(0.05);
        end
      end

   
   for i = 1 : STEPS, U(i)=i*dt; X(i)=p{i,N}(1); Y(i)=p{i,N}(2); Z(i)=p{i,N}(3); end
    subplot(2,1,1)
    plot(U,X);
    xlabel('t (s)');  ylabel('x (m)'); title('Tip Displacement - X Component');
    subplot(2,1,2)
    plot(U,Y);
    xlabel('t (s)');  ylabel('y (m)'); title('Tip Displacement - Y Component');
    saveas(gcf, '../results/SeriesTypeCR.png')
end

