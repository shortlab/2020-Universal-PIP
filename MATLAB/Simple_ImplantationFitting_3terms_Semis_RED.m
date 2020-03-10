ImplantationData_All = importdata('Annihilation_Coords_1000000_1mm_diameter_10uCi_22Na_Hist_PEN_SEMIS_2.txt');
Material_Properties = readtable('Material_Data_PositronImplantation_Semis.csv');


% Note: tried 10, 50, 100, 200 iterations and only saw difference from 10
% to 50, thus choosing 100 is certainly safe

%%% USE DISTANCE 

num_materials = 1+max(ImplantationData_All(:,1));
expo_fits_3terms_distance_semi=zeros(num_materials,5);
residuals_3terms_distance_semi=zeros(num_materials,4000);
% giving the following criteria: sse rsquare dfe adjrsquare rmse
fit_goodness_3terms_distance_semi = zeros(num_materials,5);

ft = fittype('a*exp(b*x) + c*exp(d*x) + (1-(a+c))*exp(f*x)'); 
      
for mat_num = [1:num_materials]
    mat_num
    ImplantationData = ImplantationData_All((ImplantationData_All(:,1) == (mat_num-1)), :);
    z_norm = ImplantationData(1:4000,4)/max(ImplantationData(1:4000,4));
      
    distance = [0:3999];
    distance = distance(:);
        
    fit_goodness = inf;
        
    for rep_one = 1:200
        options = fitoptions(ft);
        options.Lower = [0,-inf,0,-inf,-inf]; 
        options.Upper = [1.0,0,1.0,0,0];      
        if (rep_one > 1) | (mat_num > 1)
            options.StartPoint = [previous_fit_z.a previous_fit_z.b previous_fit_z.c previous_fit_z.d previous_fit_z.f];
        end       
        [fit_z, gof] = fit(distance,z_norm,ft,options);
        if gof.rmse < fit_goodness
            previous_fit_z = fit_z;
            fit_goodness = gof.rmse;
            fit_goodness_all = gof;
           
        end
    end
   
   bf = previous_fit_z;
   best_rmse = fit_goodness;
   total = bf.a*exp(distance*bf.b) + bf.c*exp(distance*bf.d) + (1-(bf.a+bf.c))*exp(distance*bf.f); 
   fit_goodness_3terms_distance_semi(mat_num, :) = [fit_goodness_all.sse fit_goodness_all.rsquare fit_goodness_all.dfe fit_goodness_all.adjrsquare fit_goodness_all.rmse];
    
   expo_fits_3terms_distance_semi(mat_num,:) = [bf.a bf.b bf.c bf.d bf.f];
   residuals_3terms_distance_semi(mat_num,:) = z_norm - total;
    
   
end



%%% USE DISTANCE/DENSITY

ft = fittype('a*exp(b*x) + c*exp(d*x) + (1-(a+c))*exp(f*x)'); 
num_materials = 1+max(ImplantationData_All(:,1));
expo_fits_3terms_distancedensity_semi=zeros(num_materials,5);
residuals_3terms_distancedensity_semi=zeros(num_materials,4000);
% giving the following criteria: sse rsquare dfe adjrsquare rmse
fit_goodness_3terms_distancedensity_semi = zeros(num_materials,5);

      
for mat_num = [1:num_materials]
    mat_num
    ImplantationData = ImplantationData_All((ImplantationData_All(:,1) == (mat_num-1)), :);
    z_norm = ImplantationData(1:4000,4)/max(ImplantationData(1:4000,4));
      
    distance = [0:3999];
    density = table2array((Material_Properties(mat_num, "density_gcm3")));
    distance = distance(:)/density;
        
    fit_goodness = inf;
        
    for rep_one = 1:200
        options = fitoptions(ft);
        options.Lower = [0,-inf,0,-inf,-inf]; 
        options.Upper = [1.0,0,1.0,0,0];      
        if (rep_one > 1) | (mat_num > 1)
            options.StartPoint = [previous_fit_z.a previous_fit_z.b previous_fit_z.c previous_fit_z.d previous_fit_z.f];
        end       
        [fit_z, gof] = fit(distance,z_norm,ft,options);
        if gof.rmse < fit_goodness
            previous_fit_z = fit_z;
            fit_goodness = gof.rmse;
            fit_goodness_all = gof;
           
        end
    end
   
   bf = previous_fit_z;
   best_rmse = fit_goodness;
   total = bf.a*exp(distance*bf.b) + bf.c*exp(distance*bf.d) + (1-(bf.a+bf.c))*exp(distance*bf.f); 
   fit_goodness_3terms_distancedensity_semi(mat_num,:) = [fit_goodness_all.sse fit_goodness_all.rsquare fit_goodness_all.dfe fit_goodness_all.adjrsquare fit_goodness_all.rmse];
    
   expo_fits_3terms_distancedensity_semi(mat_num,:) = [bf.a bf.b bf.c bf.d bf.f];
   residuals_3terms_distancedensity_semi(mat_num,:) = z_norm - total;
    
   
end

writematrix(expo_fits_3terms_distance_semi, "normalizedFitCoef_3TermDistance_Semi_RED.csv")
writematrix(expo_fits_3terms_distancedensity_semi, "normalizedFitCoef_3TermDistanceDensity_Semi_RED.csv")

writematrix(fit_goodness_3terms_distance_semi, "fit_goodness_normalizedFitCoef_3TermDistance_Semi_RED.csv")
writematrix(fit_goodness_3terms_distancedensity_semi, "fit_goodness_normalizedFitCoef_3TermDistanceDensity_Semi_RED.csv")


