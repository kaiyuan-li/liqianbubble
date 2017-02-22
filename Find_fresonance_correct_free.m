% resonance freq by sweeping freq using marmottant model, fixed pressure
% Qian Li
% benchmark with paper"nonlinear shell behavior of phospholipid coated microbubbles"


%? resonance? same pressure, get most oscillation energy?
%close all;
clear all;

global  R0 P0 GAMMA MU RHO c  kexi0  Pa f t_end tspan_excitation window_length t_interval hannwin


% experimental fres-diameter data
linestyle = {'k-*','k--+','k:<','k-.>','r-o','r:o','r--o','r-.o','b-','b--',...
    'b:','b-.','m-','m:','m--','m-.'};
co = [  0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
set(groot,'defaultAxesColorOrder',co,...
    'defaultAxesLineStyleOrder','-|--|:|o')

%% Constant Parameters
P0 = 1.01e5;                    % hydrostatic pressure (ambient pressure), (unit: Pa), (1.01e5 Pa)
GAMMA = 1.07;                   % polytropic exponent of the gas in the bubble, (1.4 for air, 1.07 for perfluorobutane (c4F10) gas)
MU = 0.001;                     % shear (dynamic) liquid viscosity, (unit: Pa), (0.001 Pa for water)
kexi0 = 0.0728;%725;%0.0725;                  % surface tension, (unit: N/m), (0.072 N/m for water/air)
RHO = 1000;                      % density of the liquid, (unit: kg/m3), (998 kg/m3 for water)
c = 1500;                       % speed of the sound in the liquid, (unit: m/s), (1500 m/s for water)



%% controls
Max_fun = 1; % 1= using fundamental response; 0= using maxium raidus excurstion as resonance reference
             %2 energy absorbtion.
win  = 0; % 1 = tukey win; 0 =no win driving
% use win = 0 and fixed cycles and Max_fun =1;
% or use win =1 and fixed time



%% Variable Parameters


    
R0_range = 2.8e-6;%[1.7 1.95 2.2 2.6 2.91 3.74 4 5]*1e-6;

for R0ii = 1:length(R0_range)
    R0 = R0_range(R0ii);
    display(strcat('R0 = ',num2str(R0*1e6),'um'));
    
       
    f0_uncoated(R0ii) = 1/(2*pi*R0)*sqrt(1/RHO*(3*GAMMA*P0+2*(3*GAMMA-1)*kexi0/R0));
    
    fdrive_range = [0.005]*1e6;
    cycle = 20;
    
    Pa_range = [120]*1e3;
    
    for Paii = 1:length(Pa_range)
        Pa = Pa_range(Paii);
             
        
             
        % sweep frequency, build maxium amplitude matrix of each oscillation to find the resonant
            % sweep frequency  found resonance.
            
            for ii = 1:length(fdrive_range);
                f = fdrive_range(ii);                                 % acoutic driving frequency, (unit: Hz)
                
                % Numerical integration
                x0 = [R0;0];                          % intial conditions of [x(0),x(1)], radius and velocity of bubble wall
                t_end = cycle*1/f;                               % integer * cycles. drivewave end time
                fs = 50e6;
                t_interval = 1/fs;
                tspan_excitation = [0:t_interval:t_end];     % tspan of acoustic excitation
                window_length = length(tspan_excitation);     % define hanning window length over tspan of acoustic excitation
                n = tspan_excitation/t_interval+1;
                if win == 1; % hanning wind
                    tukey = tukeywin(window_length,0.15);
                    hannwin = tukey;
                    %hannwin = 0.5*(1-cos(2*pi*n/window_length))'; % hann window for modulating driving pulse
                else % no window
                    hannwin = n./n;
                end
                
                tspan = [0:t_interval:1.5*t_end];   % interval of integration
                [t,x] = ode15s(@modFreebubble,tspan,x0);
                
                %Energy_p = sum((Pa*sin(2*pi*f*tspan_excitation).*hannwin(floor(tspan_excitation/t_interval+1))').^2)/1e20;
                
                x1 = x(:,1);  % R - t
                x2 = x(:,2);  % v - t
                
                %Energy_b = sum((x1-R0).^2);
                
                    % spectra of R-R0 versus t
                    L = length(x1);
                    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
                    Y = fft(x1-R0,NFFT)/L; % R0 would induce DC, large amplitude at freq=0
                    F = fs/2*linspace(0,1,NFFT/2+1);
                    
                    amplitude = 2*abs(Y(1:NFFT/2+1));
                    
                if Max_fun == 1;
                    
                    % fundamental response
                    spectra_resolution = F(2)-F(1);
                    f_position = floor([0.9*f/spectra_resolution 1.1*f/spectra_resolution]);
                    fundamentalResponse =  max(amplitude(f_position(1):f_position(2)));
                    
                    maxA(ii) = fundamentalResponse;
                elseif Max_fun == 0;
                    maxA(ii) = max(x1)-min(x1);
                elseif Max_fun == 2;
                    maxA(ii) = Energy_b/Energy_p;% bubble energy / excitation energy
                end
                
                maxE(ii) = max(x1);
                minE(ii) = min(x1);
                 figure(2);plot(t,x1); hold on;
                 legendinfo{ii} = [num2str(f/1e6),'MHz'];
                 hold on;
          
                
            end
            % legend(legendinfo)
            % resonance freq for each R0
            [A,B] = max(maxA);
            f_resonance(R0ii,Paii) = fdrive_range(B);
            display(strcat('resonance freq = ',num2str(f_resonance(R0ii,Paii)/1e6),'Mhz'));

        
        
%         figure(2);
%         plot(fdrive_range,maxA*1e6)
%         xlabel('freq')
%         ylabel('amplitude (um)')
%         title('resonance curve')
%         hold on;
%         
        %% evaluate if bubble in the linear oscillating region for all freq sweep
        figure(4);
        plot(fdrive_range,maxE-R0,fdrive_range,minE-R0)% maxium and mininum excursion
        hold on;
      
     
        title([' R0= ', num2str(R0*1e6), 'um, Pressure = ', num2str(Pa/1e3),'kPa'])
    end
%   disp('paused')

 % close all
end
%% simulation done alert
%sound(y,sound_fs); % sound alert


figure(100)
hold on;
plot(Pa_range/1e3,f_resonance/1e6,'LineWidth',2)

figure(2);
legend(legendinfo)
figformat;


