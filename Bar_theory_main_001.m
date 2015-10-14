clear all; close all; clc

global ds theta g u d50 delta k_visc depth dt1t c0 chezy r_trans porosity im;

% Program that solves linear theory of bars
% using the scheme of equations and perturbations
% proposed by Colombini et al. [1987]

disp('Welcome to the program: TREMTO (v0.1)!')
disp(' ')
disp('Please chose the parameters [defaults]')
choice = input('Use dimensional parameters? Y/N [Y]: ','s');
if isempty(choice)
    choice = 'Y';
end
choice=upper(choice);
 
if choice=='N'
    beta_num = input('Aspect ratio (beta)       [20]: ');
    theta    = input('Shields stress (theta)   [0.1]: ');
    ds       = input('Relative roughness (ds) [0.01]: ');
   
    %Default parameters
    if isempty(beta_num)
        beta_num = 20;
    end
    if isempty(theta)
        theta = 0.1;
    end
    if isempty(ds)
        ds = 0.01;
    end
else
    % Define the dimensional parameters of the real case river
    discharge = input('Discharge (m^3/s) [700]: ');
    width = input('Width (m) [100]: ');
    d50 = input('D50 of the sediment (m) [0.04]: ');
    slopex = input('Longitudinal slope (m/m) [0.002]: ');
  
    %Default parameters
    if isempty(discharge)
        discharge = 700;
    end
    if isempty(width)
        width = 100;
    end
    if isempty(d50)
        d50 = 0.04;
    end
    if isempty(slopex)
        slopex = 0.002;
    end
    
    g = 9.81;

end
r_trans = input('Transversal bedslope correction coefficient [0.3]: ');
if isempty(r_trans)
    r_trans = 0.3;
end

% Declaration of global variables
delta = 1.65;
porosity = 0.4;
im = sqrt(-1);  % imaginery part

% Parameters of figures
beta_m=1;
beta_M=80;
Nbeta=100;

lambda_m = 0;
lambda_M = 2;
Nlambda=100;

sub_01_input
sub_03_cderi
sub_04_phideri


if theta<=theta_crit
    % case: theta < theta critical
    % Bedload transport is zero
    % Linear theory is not applyable
    disp(' ');
    disp(['theta    = ' num2str(theta)]);
    disp(['ds       = ' num2str(ds)]);
    disp(['beta     = ' num2str(beta_num)]);
    disp(['theta < theta_critical --> ' num2str(theta) ' < ' num2str(theta_crit)])
    disp(' ');
    disp('Linear bar theory is then not applyable')
    disp('Please modify your parameters in order to have theta>theta_cr')
else
    sub_05_computations
    disp(' ');
    disp(['theta    = ' num2str(theta)]);
    disp(['ds       = ' num2str(ds)]);
    disp(['beta     = ' num2str(beta_num)]);
    disp(['beta_cr  = ' num2str(beta_cr)]);
    disp(['beta_res = ' num2str(beta_res)]);
    
end
