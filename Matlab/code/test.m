clc;
force_REI_ratio = linspace(5,15,4);
lengthvalues = linspace(0.05,0.5,4);
r_values = linspace(0.001,0.01,4);
for radindex = 1:length(r_values)
    for lengthindex = 1:length(lengthvalues)
        l_r_ratio = lengthval / r_values(radindex);
        if l_r_ratio > 40 || l_r_ratio < 10
            fprintf(sprintf('Loop skipped by the l_r_ratio %0.3f, l = %0.2f r = %0.2f\n',l_r_ratio,lengthvalues(lengthindex),r_values(radindex)))
            continue
        end
        
        for forceindex = 1:length(force_REI_ratio)
            fprintf(sprintf('l_r_ratio %0.3f, l = %0.2f r = %0.3f\n',l_r_ratio,lengthvalues(lengthindex),r_values(radindex)))
        end


    end
end
