function d = uniform(q,w,d50,s)
fd = 1;
eps = 10 ^-10;
ks = 21/(d50^(1/6));
d = (q/(ks*s^0.5*w))^(3/5);
delta = 20;  % delta >=1
d_min = d*1/delta;
d_max = d*delta;
k_un=0;
while(abs(fd) > eps)
    k_un =k_un + 1;
    chez = 6+2.5*log(d/(2.5*d50));
    rh = (w*d)/(w+2*d);
    fd = w*d*chez*(9.81*s*rh)^0.5-q;
    if fd>0
        d_old =d;
        d = d - abs(fd/q);
        if d<=d_min
            %             d = (d_max+d_min)/2;
            d = (d_max+d_min)/2;
        end
        d_max = d_old;
    else
        d_old =d;
        d = d + abs(fd/q);
        if d>=d_max
            %             d = (d_min+d_max)/2;
            d = (d_min+d_max)/2;
        end
        d_min = d_old;
    end
%     disp(['fd         : ' num2str(fd,'%15.10f')    ' [m]']);
%     disp(['d_old      : ' num2str(d_old,'%15.10f') ' [m]']);
%     disp(['d_new      : ' num2str(d,'%15.10f')     ' [m]']);
%     disp(['d_min      : ' num2str(d_min,'%15.10f') ' [m]']);
%     disp(['d_max      : ' num2str(d_max,'%15.10f') ' [m]']);
%     disp(['k          : ' num2str(k_un,'%15.0f')   ' [-]']);
%     disp(' ');
%     pause() %pause(0.3)
    if(k_un == 20)
        d = (d_min+d_max)/2;
        k_un = 0;
    end
end
% chez = @(z) 6+2.5*log(z/(2.5*d50));
% rh = @(z) w*z/(w+2*z);
% fd = @(z) q-w*z*chez(z)*(9.81*s*rh(z))^0.5;
% ks = 21/(d50^(1/6));
% f0 = (q/(ks*s^0.5*w))^(3/5);
% d = fzero(fd,f0);
