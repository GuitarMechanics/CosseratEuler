%This is the main running file of the capsule. To run what you whant, just select the targeted configuration and uncommand the others. 

% Single section continuum robot:                                       SingleSectionCR
% Multiple section Series type Continuum robot:                         SeriesTypeCR 
% Tendon bent Concentric tube:                                          ConcentricTubeCR
% parallel/series Supportive type cooperative continuum robot           SupportiveTypeCR
% Comanipulative type cooperaive continuum robot:                       CoManipulativeTypeCR 

% THE CODE HAS BEEN SET ON SINGLE SECTION CR BY DEFAULT

for i = 0:0.1:0.5
    clear all; clc;
    SingleSectionCR(i);
end

%SeriesTypeCR
%ConcentricTubeCR
%SupportiveTypeCCR
%CoManipulativeTypeCCR

        
            
            
            
            
            

  
            
        