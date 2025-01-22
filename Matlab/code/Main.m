%This is the main running file of the capsule. To run what you whant, just select the targeted configuration and uncommand the others. 

% Single section continuum robot:                                       SingleSectionCR
% Multiple section Series type Continuum robot:                         SeriesTypeCR 
% Tendon bent Concentric tube:                                          ConcentricTubeCR
% parallel/series Supportive type cooperative continuum robot           SupportiveTypeCR
% Comanipulative type cooperaive continuum robot:                       CoManipulativeTypeCR 

% THE CODE HAS BEEN SET ON SINGLE SECTION CR BY DEFAULT

% generlize >> t/EI (temporary)

clear all; clc;
E = 207e6; %Pa
A = pi*(0.1)^2; %m^2
unit_EA = E*A
force_over_EA = [0.5e-05, 1e-05, 1.5e-05, 2e-05, 2.5e-05, 3e-05, 3.5e-05, 4e-05];
resolution = 50;
length = [0.2, 0,4, 0.6, 0.8, 1];
for i = 1:8
    force = force_over_EA(i)*unit_EA;
    try
        SingleSectionCR(force,length,resolution, 300)
    catch
        resolution = resolution + 50;
        SingleSectionCR(force,length,resolution, 300)
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

        
            
            
            
            
            

  
            
        