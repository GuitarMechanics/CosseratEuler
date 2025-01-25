%This is the main running file of the capsule. To run what you whant, just select the targeted configuration and uncommand the others. 

% Single section continuum robot:                                       SingleSectionCR
% Multiple section Series type Continuum robot:                         SeriesTypeCR 
% Tendon bent Concentric tube:                                          ConcentricTubeCR
% parallel/series Supportive type cooperative continuum robot           SupportiveTypeCR
% Comanipulative type cooperaive continuum robot:                       CoManipulativeTypeCR 

% THE CODE HAS BEEN SET ON SINGLE SECTION CR BY DEFAULT

% generlize >> rl/EI (temporary)

clear all; clc;
E = 207e9; %Pa<<springsteel
% E = 70e9; %Pa<<NiTi
A = pi*(0.01)^2; %m^2
momentI = (pi*(0.1)^4)/64
unit_EI = E*momentI
unit_EA = E*A
force_over_EA = ([1, 2, 3, 4, 5, 6,7,8,9,10]+10) * 1e-3;
% force_over_EA = [7e-05,8e-05];
force_REI_ratio = linspace(5,20,10)*1e-4;
lengthvalues = linspace(0.05,0.3,11);
% lengthvalues = [0.8 1.0 1.2,1.4,1.6,1.8,2.0];
% for lengthindex = 1:length(lengthvalues)
for lengthindex = 1:length(lengthvalues)
   clc;
    forceratio = (lengthvalues(lengthindex))*0.01/unit_EI;
    lengthval = lengthvalues(lengthindex);

    resolution = 50;
    timestepcount = 300;
    for forceindex = 1:length(force_REI_ratio)
        force = force_REI_ratio(forceindex) / forceratio;
        try
            SingleSectionCR(force,lengthval,resolution, timestepcount)
        catch
            resolution = resolution + 50;
            try
                SingleSectionCR(force,lengthval,resolution, timestepcount)
            catch
                resolution = resolution + 50;
                try
                    SingleSectionCR(force,lengthval,resolution, timestepcount)
                catch
                    break %% 해상도 상향 2트해도 안되면 다음 길이로 넘어감
                end
            end
        end
    end
end


% for i = 0.05:0.15:0.8
%     clear all; clc;
%     SingleSectionCR(i)
% end

%SeriesTypeCR
%ConcentricTubeCR
%SupportiveTypeCCR
%CoManipulativeTypeCCR

        
            
            
            
            
            

  
            
        