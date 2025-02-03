%This is the main running file of the capsule. To run what you whant, just select the targeted configuration and uncommand the others. 

% Single section continuum robot:                                       SingleSectionCR
% Multiple section Series type Continuum robot:                         SeriesTypeCR 
% Tendon bent Concentric tube:                                          ConcentricTubeCR
% parallel/series Supportive type cooperative continuum robot           SupportiveTypeCR
% Comanipulative type cooperaive continuum robot:                       CoManipulativeTypeCR 

% THE CODE HAS BEEN SET ON SINGLE SECTION CR BY DEFAULT

% generlize >> rl/EI (temporary)

%최적화 시도해볼것:
%스레드 생성을 메인으로 빼기(계속 만들었다 지웠다 하다보니 발적화되는거같다)
%clc 적절히 배치해서 콘솔 메모리 확보
%컨피겨레이션 피겨 hold 지워서 렉안걸리게하기기
%fsolve 루프 조기탈출설정

% 기타 사항:
% 치수 변경(얇게)
% l/r 비율이 100 이상이 되는 경우를 어케 처리할 지 생각해볼 것(Python/test_codes.py 참고)
% fsolve 뻑나면 exitcode 띄우고 코드정지
% 실리콘러버는 포아송비가 다름-->SingleSectionCR에서 G의 값을 조정
% 실리콘러버 댐핑때문인지 수렴이 안됨 --> 타임스텝 dt 1/10으로 적용, 댐핑 줄이기
clear; clc;
% E = 200e9; %Pa<<sus304
% E = 80e9; %Pa<<NiTi, Austenite
E = 50e6; % sliicone rubber
% force_over_EA = [7e-05,8e-05];
force_REI_ratio = linspace(5,15,6); % Niti, SUS
% force_REI_ratio = linspace(5,15,6)*10; % Silicone 
lengthvalues = linspace(0.01,0.5,11);
r_values = linspace(0.001,0.01,10);
delete(gcp); % Close any existing parallel pool
parpool("threads",6); % Start parallel pool
try
    clf(1)
    clf(2)
catch
    fprintf('No figure to clear')
end
for radindex = 1:length(r_values)
    A = pi * (r_values(radindex)^2);
    momentI = (pi*r_values(radindex)^4)/64;
    unit_EI = E*momentI;

    for lengthindex = 1:length(lengthvalues)
        forceratio = (lengthvalues(lengthindex))*r_values(radindex)/unit_EI;
        lengthval = lengthvalues(lengthindex);
        l_r_ratio = lengthval / r_values(radindex)
        resolution = 50;
        timestepcount = 1500;
        
        if l_r_ratio > 50 || l_r_ratio < 10
            fprintf(sprintf('Loop skipped by the l_r_ratio %0.2f',l_r_ratio))
            continue
        end
            for forceindex = 1:length(force_REI_ratio)
                force = force_REI_ratio(forceindex) / forceratio;
                clc;
                SingleSectionCR(force,lengthval,resolution, timestepcount,r_values(radindex),E)
                clear global
                clf(1)
            end
        clf(2)
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

        
            
            
            
            
            

  
            
        