% Definition of scaled parameters and choice of equations that describes
% bedload transport
global ds theta g u d50 delta k_visc depth dt1t c0 chezy r_trans porosity im;

% Uniform flow
if choice=='N'
    % Non dimensional case
    chezy=6.+2.5*log(1/(2.5*ds));
    fr=sqrt(theta*ds*delta*chezy^2);
    
    Q_0 = ds/(1-porosity)*1/(chezy*sqrt(theta));
else 
    % Dimensional case
    depth = uniform(discharge,width,d50,slopex);
    u = discharge/width/depth;
    fr = u/(9.81*depth)^0.5;
    % real_lambda = pi*width/bar_lenght;
    beta_num = width/2/depth;
    ds = d50/depth;
    theta = slopex/(delta*ds); 
    Q_0 = sqrt(delta * g * d50^3) / ((1-porosity) * depth * u);
    
   
end
