%This is the main running file of the capsule. To run what you whant, just select the targeted configuration and uncommand the others. 

% Single section continuum robot:                                       SingleSectionCR
% Multiple section Series type Continuum robot:                         SeriesTypeCR 
% Tendon bent Concentric tube:                                          ConcentricTubeCR
% parallel/series Supportive type cooperative continuum robot           SupportiveTypeCR
% Comanipulative type cooperaive continuum robot:                       CoManipulativeTypeCR 

% THE CODE HAS BEEN SET ON SINGLE SECTION CR BY DEFAULT

clear all; clc;
forces = [70,80,90,100,110,120,130,140];
resolution = 1000;
length = 1;
for i = 1:8
    SingleSectionCR(forces(i),length,resolution)
end


% for i = 0.05:0.15:0.8
%     clear all; clc;
%     SingleSectionCR(i)
% end

%SeriesTypeCR
%ConcentricTubeCR
%SupportiveTypeCCR
%CoManipulativeTypeCCR

        
            
            
            
            
            

  
            
        