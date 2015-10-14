global ds theta g u d50 delta k_visc depth dt1t c0 chezy r_trans porosity im;

%% 
% disp(' ')
% disp('Please select which roughness formula do you want to choose:')
% disp('Type 1 for - EINSTEIN modified by ENGELUND & HANSEN')
% disp('Type 2 for - VAN RIJN')
% rough_formula = input('Select now: ');
rough_formula = 1;

%%
% disp(' ')
% disp('Please select bedform roughness:')
% disp('Type 1 for - FLAT BED (no correction)')
% disp('Type 2 for - DUNES')
% % disp('Type 3 for - RIPPLES')
% bed_rough = input('Select now: ');
bed_rough = 1;

% if choise is wrong, then ..
% include if you want to correct the roughness due to dunes or ripple
%%

%% Choice 1: EINSTEIN [1950] modified by ENGELUND & HANSEN [1967]
k_chezy = 1;
if(rough_formula == 1)
    if   (bed_rough == 1) % FLAT BED case
        depth_theta = 1;    % depth could be funtion of theta in case of presence of bedforms
    elseif(bed_rough == 2) % DUNES
        depth_theta = 0.06 / theta + 0.3 * sqrt(theta); % ENGELUND & FREDSOE [1982] for dunes scaled with unperturbed theta
        if(depth_theta >= 1)
            depth_theta = 1;
        end
    end
    chezy = 6 + 2.5*log(depth_theta/(k_chezy*2.5*ds));
    c0 = 1 / (chezy^2 * depth_theta); % Friction coefficient unperturbed
    theta1 = theta * depth_theta;     % ENGELUND & FREDSOE [1982] for dunes
    dt1t = 1;
    if(depth_theta < 1)
        dt1t = 0.8*theta;
    end
    % Cd and Ct are defined in 18.a,b of COLOMBINI et al. [1987]
    c_d = 1 / c0 * (-5 * chezy^-3/depth_theta); %in the bracket there is the partial derivative of Cf over depth_theta
    c_t = (1 - dt1t/depth_theta) * (1 + 5 / chezy);
    % Extra results for Ashida & Michiue [1971]
    depth_AM = 1;
    chezy_AM = 6 + 2.5*log(depth_AM/(k_chezy*2.5*ds));
    c0_AM = 1 / (chezy_AM^2 * depth_AM); % Friction coefficient unperturbed
    theta1_AM = theta * depth_AM;     % ENGELUND & FREDSOE [1982] for dunes
    c_d_AM = 1 / c0_AM * (-5 * chezy_AM^-3/depth_AM); %in the bracket there is the partial derivative of Cf over depth_theta
    c_t_AM = (1 - 1/depth_AM) * (1 + 5 / chezy_AM);
end

%% Choice 2: VAN RIJN [1984]

if(rough_formula == 2)
    if(bed_rough == 1) % FLAT BED case
        depth_theta = 1;
        chezy = 18 * log(12*depth_theta/(3*ds));
    end
    if(bed_rough == 2) % DUNES
        depth_theta = 1;
        d_star = d50*(delta*g/k_visc^2)^(1/3);    % dimensionless particle diameter
        if(d_star <= 4)
            theta_cr = 0.240 / d_star;
        elseif and(d_star >  4, d_star <= 10)
            theta_cr = 0.140 / d_star^0.64;
        elseif and(d_star > 10, d_star <= 20)
            theta_cr = 0.040 / d_star^0.10;
        elseif and(d_star > 20, d_star <= 150)
            theta_cr = 0.013 / d_star^0.29;
        else
            theta_cr = 0.055;
        end
        u_star_cr = sqrt(theta_cr*delta*g*d50);
        chezy_grain = 18 * log(12*depth_theta/(3*2.11*ds));     % 3*d90 is 3* 2.11*d50
        u_star = (g^0.5 / chezy_grain) * u;
        tsp = (u_star^2 - u_star_cr^2)/u_star_cr^2;             % transport stage parameter
        if(tsp< 0)
            disp(' ')
            disp('Case in which there is no sediment transport, using Van Rijn [1984]');
            tsp = 0;
        end
        psi = 0.015 * (ds)^0.3 * (1 - exp(-0.5*tsp)) * (25 - tsp);
        dune_length = 7.3 * depth;
        dune_height = psi * dune_length;
        ds_eff = 3* 2.11 *d50 + 1.1 * dune_height * (1 - exp(-25*psi));
        chezy = 18 * log(12/(2.11*ds_eff));
    end
    c0 = 1 / (chezy^2 * depth_theta); % Friction coefficient unperturbed
    theta1 = theta * depth_theta;     % ENGELUND & FREDSOE [1982] for dunes
    dt1t = 1;
    if(depth_theta < 1)
        dt1t = 0.8*theta;
    end
    % Cd and Ct are defined in 18.a,b of COLOMBINI et al. [1987]
    c_d = 1 / c0 * (-36 * chezy^-3/depth_theta); % in the bracket there is the partial derivative of Cf over depth_theta
    c_t = (1 - dt1t/depth_theta) * (1 + 5 / chezy); % Ask GUIDO !!
        % Extra results for Ashida & Michiue [1971]
    depth_AM = 1;
    chezy_AM = 18 * log(12*depth_AM/(3*ds));
    c0_AM = 1 / (chezy_AM^2 * depth_AM); % Friction coefficient unperturbed
    theta1_AM = theta * depth_AM;     % ENGELUND & FREDSOE [1982] for dunes
    c_d_AM = 1 / c0_AM * (-36 * chezy_AM^-3/depth_AM); % in the bracket there is the partial derivative of Cf over depth_theta
    c_t_AM = (1 - 1/depth_AM) * (1 + 5 / chezy_AM);

end
