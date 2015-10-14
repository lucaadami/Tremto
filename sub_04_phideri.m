global ds theta g u d50 delta k_visc depth dt1t c0 chezy r_trans porosity im;

%%
cont = 0;
% if cont = 0 it is calculated only the formula selected
% if cont = 1 all the formulations are calculated
if(cont == 0)
    bed_max = 1;
    %     disp(' ')
    %     disp('Please select which bedload formula do you want to choose:')
    %     disp('Type  1 for - MEYER PETER & MULLER [1948]')
    %     disp('Type  2 for - WONG & PARKER [2006]     (MPM modified)')
    %     disp('Type  3 for - HUNZIKER & JAEGGI [2002] (MPM modified)')
    %     disp('Type  4 for - HUANG [2010]             (MPM modified)')
    %     disp('Type  5 for - ENGELUD & HANSEN [1972]')
    %     disp('Type  6 for - VAN RIJN [1984]')
    %     disp('Type  7 for - PARKER [1990]')
    %     disp('Type  8 for - ASHIDA & MICHUE [1971]')
    %     disp('Type  9 for - BAGNOLD [1980]')
    %     bedload_formula = input('Select now: ');
    bedload_formula = 2;
elseif(cont == 1)
    bed_max = 9;    % this value has to be the total number of the bedload formula implemented
    bedload_formula = 1;
end

%% common

phi_0 = 0;
phi_t = 0;
phi_d = 0;


for i = 1:bed_max
    %% Choice 1: MEYER PETER & MULLER [1948]
    
    if(bedload_formula == 1)
        theta_cr_1 = 0.047;
        theta_thres = theta - theta_cr_1;
        if(theta_thres > 0)
            phi_0 = 8 * theta_thres ^ 1.5;
            phi_t = 1.5 * theta * dt1t / theta_thres;
        end
        theta_crit = theta_cr_1;
    end
    
    %% Choice 2: WONG & PARKER [2006]     (MPM modified)
    
    if(bedload_formula == 2)
        %         theta_cr_2 = 0.047;
        theta_cr_2 = 0.0495;
        if (theta <= theta_cr_2)
            theta_thres = 0;
            phi_0 = 0;
            phi_t = 0;
        else
            theta_thres = theta - theta_cr_2;
            % Basement MPM modified for WP
            phi_0 = 3.97 * theta_thres ^ 1.5;
            phi_t = 1.5 * theta * dt1t / theta_thres;
            %             phi_0 = 4.93 * theta_thres ^1.6;
            %             phi_t = 1.6 * theta * dt1t / theta_thres;
        end
        theta_crit = theta_cr_2;
    end
    %% Choice 3: HUNZIKER & JAEGGI [2002] (MPM modified)
    
    if(bedload_formula == 3)
        theta_cr_3 = 0.05;
        theta_thres = theta - theta_cr_3;
        if(theta_thres > 0)
            phi_0 = 5 * theta_thres ^ 1.5;
            phi_t = 1.5 * theta * dt1t / theta_thres;
        end
        theta_crit = theta_cr_3;
    end
    %% Choice 4: HUANG [2010]
    
    if(bedload_formula == 4)
        theta_cr_4 = 0.047;
        theta_thres = theta - theta_cr_4;
        if(theta_thres > 0)
            phi_0 = 6 * theta_thres ^ (5/3);
            phi_t = (5/3) * theta * dt1t / theta_thres;
        end
        theta_crit = theta_cr_4;
    end
    %% Choice 5: ENGELUD & HANSEN [1972]
    
    if(bedload_formula == 5)
        theta_cr_5 = 0.04;          % Should be zero but it doesn't converge
        theta_thres = theta - theta_cr_5;
        phi_0 = 0.05 * theta ^ (5/2) / c0;
        phi_d =     - c_t;
        phi_t = 2.5 - c_t;
        theta_crit = theta_cr_5;
    end
    %% Choice 6: VAN RIJN [1984]
    
    if(bedload_formula == 6)
        d_star = d50*(delta*g/k_visc^2)^(1/3);          % dimensionless particle diameter (by Van Rijn)
        if(d_star <= 4)
            theta_cr_6 = 0.240 / d_star;
        elseif and(d_star >  4, d_star <= 10)
            theta_cr_6 = 0.140 / d_star^0.64;
        elseif and(d_star > 10, d_star <= 20)
            theta_cr_6 = 0.040 / d_star^0.10;
        elseif and(d_star > 20, d_star <= 150)
            theta_cr_6 = 0.013 / d_star^0.29;
        else
            theta_cr_6 = 0.055;
        end
        u_star_cr = sqrt(theta_cr_6*delta*g*ds);                % d50 instead of ds?
        chezy_grain = 18 * log(12*depth_theta/(3*2.11*ds));     % 3*d90 is 3* 2.11*d50
        u_star = (g^0.5 / chezy_grain) * u;
        tsp = (u_star^2 - u_star_cr^2)/u_star_cr^2;             % transport stage parameter
        if(tsp< 0)
            disp(' ')
            disp('Case in which there is no sediment transport, using Van Rijn [1984]');
            tsp = 0;
        end
        phi_0 = 0.053 * tsp ^ 2.1 / d_star ^ 0.3;
        phi_t = 2.1 * theta / (tsp * theta_cr_6);
        theta_thres = tsp;
        theta_crit = theta_cr_6;
    end
    %% Choice 7: PARKER [1990]
    
    if(bedload_formula == 7)
        theta_ref = 0.0386; % Parker 1990 reference Shields
        phi_parker = theta / theta_ref;
        phi_0 = 0.00218 * theta^(3/2);
        if(phi_parker > 1.588)
            phi_0 = phi_0 * 5474 * (1 - 0.853 / phi_parker)^(9/2);
            phi_t = 1.5 + 4.5 *  0.853 / phi_parker / (1 - 0.853/phi_parker);
        end
        if and(phi_parker > 1, phi_parker <= 1.588)
            phi_0 = phi_0 * exp(14.2 * (phi_parker - 1) - 9.28 * (phi_parker - 1) ^ 2);
            phi_t = 1.5 + phi_parker * (14.2 - 18.56 * (phi_parker - 1));
        end
        if(phi_parker <= 1)
            phi_0 = phi_0 * phi_parker ^ 14.2;
            phi_t = 15.7;
        end
        theta_thres = theta - theta_ref; % limit value of convergence
        theta_crit = theta_cr_7;
    end
    %% Choice 8: ASHIDA & MICHUE [1971]
    
    if(bedload_formula == 8)
        if and(d50>=0.0003,d50<0.007)
            disp('You are into the limits of the bedload formula of ASHIDA & MICHUE')
        else
            disp(' ')
            disp('ALERT!!!')
            disp('ALERT: d50 is out of the limits of validity of the bedload formula of ASHIDA & MICHUE')
        end
        chezy_8 = chezy;  %
        theta_8 = u^2 /(chezy_8^2 * delta * d50);         % theta is the mobility parameter calculated referring only to grain roughness
        d_star = d50*(delta*g/k_visc^2)^(1/3);            % dimensionless particle diameter (by Van Rijn)
        if(d_star <= 4)
            theta_cr_8 = 0.240 / d_star;
        elseif and(d_star >  4, d_star <= 10)
            theta_cr_8 = 0.140 / d_star^0.64;
        elseif and(d_star > 10, d_star <= 20)
            theta_cr_8 = 0.040 / d_star^0.10;
        elseif and(d_star > 20, d_star <= 150)
            theta_cr_8 = 0.013 / d_star^0.29;
        else
            theta_cr_8 = 0.055;
        end
        theta_thres = theta - theta_cr_8;
        chezy_grain = 18 * log(12*depth_theta/(3*2.11*ds));     % 3*d90 is 3* 2.11*d50
        u_star_prime = (g^0.5 / chezy_grain) * u;
        u_star       = (g^0.5 / chezy_8    ) * u;
        % Now we need to recall cderi in order to get c (function of bedforms)
        % and c' (flatbed case)
        phi_0 = 17 * (u_star_prime / u_star * theta_8) ^ (3/2) * (1 - theta_cr_8/theta_8) * (1 - sqrt(theta_cr_8/theta_8));
        phi_d = 17 *                          theta_8  ^ (3/2) * (1 - theta_cr_8/theta_8) * (1 - sqrt(theta_cr_8/theta_8)) ...
            * 1/chezy_grain * c_d * c0 - chezy_8 / chezy_grain ^ 2 * c_d_AM * c0;
        phi_t = (chezy_8/chezy_grain)^1.5 * (8.5 * (3 * sqrt(theta_8)- 2 * sqrt(theta_cr_8) - theta_cr_8/sqrt(theta_8))) ...
            + 1.5 * 17 * theta_8  ^ (3/2) * (1 - theta_cr_8/theta_8) * (1 - sqrt(theta_cr_8/theta_8)) *  sqrt(chezy_8) / chezy_grain ...
            * c_t * c0 / theta_8;
        %     phi_t = - (0.5 * (3 * theta * sqrt(theta_cr/theta)-2.*theta_cr-1*theta_cr ...
        %         *sqrt(theta_cr/theta)))/((-1+sqrt(theta_cr/theta))*(theta-1.*theta_cr) ...
        %         *sqrt(theta_cr/theta));
        theta_crit = theta_cr_8;
    end
    %% Choice 9: BAGNOLD [1980]
    
    % if(bedload_formula == 9)
    %     e_b = 0.15;      % efficiency factor e_b = [0.1,0.2]
    %     tan_psi = 0.6;    % dynamic coefficient of friction, 0.6 for naturally shaped sediment
    %     chezy_9 = chezy;
    %     phi_0 = % u ^ 2 / chezy_9 ^2 *
    %     phi_t = 0; % to complete
    %     theta_crit = theta_cr_9;
    % end
    
    %% Case of calculation of all the bedload formula
    if(cont==1)
        phi0(i) = phi_0;
        phid(i) = phi_d;
        phit(i) = phi_t;
        bedload_formula = bedload_formula +1;
    end
end
%% display results
% disp(' ')
% disp(['phi_0 = ', num2str(phi_0)]);
% disp(['phi_d = ', num2str(phi_d)]);
% disp(['phi_t = ', num2str(phi_t)]);
