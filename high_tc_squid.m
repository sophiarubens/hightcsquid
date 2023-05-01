%% 1. preprocessing
gainR=10;
tol=0.1*0.01;
Rmtol=gainR-tol*gainR;
Rptol=gainR+tol*gainR;

% a. read data from files into arrays
vi_data=readmatrix('thurs6p_v_i.csv');
vphi_data=readmatrix('thurs6p_v_phi.csv');
shap_data_acc=readmatrix('shapiro12_25apr.csv');
shap_data_comp=readmatrix('shapiro6_25apr.csv');

% b. use this data for determining I_C and R_N:
vi_curr=vi_data(:,1)'; % raw values in volts --> manual says this indicates 100 microamps = 10^2 *10^-6 amps = 10^-4 amps so 
                       % use this prefactor for physically meaningful values [[use in axis labels, later on]]
vi_volt=vi_data(:,2)'; % raw values in volts --> manual says this indicates 100 microvolts = 10^-4 volts so use this prefactor 
                       % for physically meaningful values [[use in axis labels, later on]]
del_vi_curr=0.001; % since the oscilloscope-produced CSVs have currents specified to the nearest 0.01
del_vi_volt=0.001; % since the oscilloscope-produced CSVs have voltages specified to the nearest 0.01
delI_add_VI_m=abs(vi_curr./Rptol)';
delI_add_VI_p=abs(vi_curr./Rmtol)';
errors_vi_curr_m=(del_vi_curr)*ones(length(vi_curr),1)+delI_add_VI_m; %FOR FIGURES 1 AND 2
errors_vi_curr_p=(del_vi_curr)*ones(length(vi_curr),1)+delI_add_VI_p;
errors_vi_voltages=del_vi_volt*ones(length(vi_volt),1);

% c. use this data for determining deltaV:
vphi_flux=vphi_data(:,1)'; % use the same ch1 correction factor as above
vphi_volt=vphi_data(:,2)'; % use the same ch2 correction factor as above
del_vphi_flux=0.001; % since the oscilloscope-produced CSVs have fluxes specified to the nearest 0.0001
del_vphi_volt=0.00001; % since the oscilloscope-produced CSVs have voltages specified to the nearest 0.01
% see below for error bars --> they only need to be applied to the useful part of the data set --> those will be for FIGURE 3

% L isn't actually calculated, and calculating beta_L relies on the parameters determined using the above data sets 
% --> there is no more data to import

% d. use this data for determining e/h in the Shapiro steps section:
shap_curr_acc=shap_data_acc(:,1)';
shap_volt_acc=shap_data_acc(:,2)';
shap_curr_comp=shap_data_comp(:,1)';
shap_volt_comp=shap_data_comp(:,2)';

del_shap_curr=0.001; % applies to both Shapiro data sets considered here- just reflects the uncertainty in the oscilloscope data
del_shap_volt=0.001; % see below for error est. on the size (distance along voltage axis) of the Shapiro steps
delI_add_Shap6_m=abs(shap_curr_acc./Rptol)';
delI_add_Shap6_p=abs(shap_curr_acc./Rmtol)';
delI_add_Shap12_m=abs(shap_curr_comp./Rptol)';
delI_add_Shap12_p=abs(shap_curr_comp./Rmtol)';
errors_shap_curr6_m=del_shap_curr*ones(length(shap_curr_acc),1)+delI_add_Shap6_m; %FOR FIGURE 4
errors_shap_curr6_p=del_shap_curr*ones(length(shap_curr_acc),1)+delI_add_Shap6_p;
errors_shap_curr12_m=del_shap_curr*ones(length(shap_curr_acc),1)+delI_add_Shap12_m;
errors_shap_curr12_p=del_shap_curr*ones(length(shap_curr_acc),1)+delI_add_Shap12_p;
errors_shap_voltages=del_shap_volt*ones(length(shap_volt_acc),1);

%% 2. critical current (I_C)
% a. figure 1 -- plot the entire V-I curve to get a feel for the data
figure
plot(vi_curr,vi_volt,'.','MarkerFaceColor','b','LineStyle','none')
errorbar(vi_curr,vi_volt,-errors_vi_voltages,errors_vi_voltages,-errors_vi_curr_m,errors_vi_curr_p,'.')
title('all critical current data')
xlabel('current; A e-4')
ylabel('voltage; V e-4')

% b. filter so only the nonzero-slope part of the V-I curve just after the knee remains
vi_curr_linear=NaN;
vi_volt_linear=NaN;
errors_vi_lin_curr_m=NaN;
errors_vi_lin_curr_p=NaN;
j=1;
for i=1:2500
    if(vi_curr(i)>0.6 && vi_curr(i)<1.32)
        % fprintf('for i=%d, entered loop to preserve a value b/c CH1(i)=%f\n',i,CH1(i))
        vi_curr_linear(j)=vi_curr(i); % I know it's inefficient to resize the arrays which hold linearized data with every 
                                      % iteration, but this is a case where pre-declaring an array of zeroes won't work, because 
                                      % the nonlinear spacing and non-monotonic traversal of values makes it quite difficult to 
                                      % predict the size of the necessary linearized array in advance
        vi_volt_linear(j)=vi_volt(i);
        errors_vi_lin_curr_m(j)=errors_vi_curr_m(i); %must set inside loop to get the proper value--again, can't preallocate size
        errors_vi_lin_curr_p(j)=errors_vi_curr_p(i);
        j=j+1; % increment the index used to store values in the arrays of just the linear part of the data set
    end
end
vi_curr_linear=vi_curr_linear';
vi_volt_linear=vi_volt_linear';
errors_vi_lin_voltages=del_vi_volt*ones(length(vi_volt_linear),1); % okay to set voltage error bars manually b/c the voltages 
                                                                   % are not subject to the additional contribution to 
                                                                   % uncertainty from the gain resistor

% d. linear regression on the isolated part of the data set + extract coefficients
icmdl=fitlm(vi_curr_linear,vi_volt_linear);
icint=icmdl.Coefficients{1,1};
icintunc=icmdl.Coefficients{1,2};
icslope=icmdl.Coefficients{2,1};
icslunc=icmdl.Coefficients{2,2};
icrsq=icmdl.Rsquared.Ordinary;
ic_eqn=['y=',num2str(icslope,4),'x+',num2str(icint,4),'; R-squared=',num2str(icrsq,4)];

% e. calculate the intercepts
intercept=-icint/icslope; % x-int when y=0, so y=m*x+b becomes 0=m*xint*b --> xint=-b/m
I_C=(intercept/2)*(1e-4); % critical current is half the value of the intercept b/c there are two Josephson junctions

% f. prepare vectors that show this fit for more currents
curr=linspace(0,4,500);
volt=icslope*curr+icint;

% g. figure 2 -- plot the linear regression model (data, fit, and confidence bounds)
figure
hold on
plot(curr,volt); % although in the region of the actual data and fitlm() this is redundant, it helps visualize the extrapolation
                 % to the horizontal intercept
handle_ic=plot(icmdl,'MarkerSize',3,'Marker','o','MarkerFaceColor','blue');
fitHandle = findobj(handle_ic,'DisplayName','Fit');
cbHandles = findobj(handle_ic,'DisplayName','Confidence bounds');
cbHandles = findobj(handle_ic,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
fitHandle.Color = 'blue';
set(cbHandles, 'Color', 'r', 'LineWidth', 1)
text(vi_curr_linear(3),vi_volt_linear(3),ic_eqn,'Position',[0.2,0.46,0])
axis([0,1.5,0,1])
errorbar(vi_curr_linear,vi_volt_linear,-errors_vi_lin_voltages,errors_vi_lin_voltages,...
         -errors_vi_lin_curr_m,errors_vi_lin_curr_p,'.')
title('positive-slope portion of critical current data just above positive knee')
xlabel('current; A e-4')
ylabel('voltage; V e-4')
legend('linear regression with extrapolation','data','',' confidence bounds','','Location','northwest')
hold off

% h. error propagation: uncertainty in I_C
del_I_C=(((-icint*icslunc)/icslope)+icintunc)/((2*10^4)*icslope);

%% 3. normal state resistance (R_N)
% a. establish endpoints of the V-I curve -- again, do this manually, b/c the points aren't all stored in order
mincurr=0;
mincurridx=0;
maxcurr=0;
maxcurridx=0;
npts=length(vphi_data);
for j=1:1:npts %ich1 and ich2 are the data of interest here
    if(vi_curr(j)<mincurr)
        mincurr=vi_curr(j);
        mincurridx=j;
    end
    if(vi_curr(j)>maxcurr)
        maxcurr=vi_curr(j);
        maxcurridx=j;
    end
end
deltaI_width=vi_curr(maxcurridx)-vi_curr(mincurridx);
deltaV_height=vi_volt(maxcurridx)-vi_volt(mincurridx);
R_N=2*deltaV_height/deltaI_width; %no need for scale factors here (both voltage and current should be scaled by 1e-4 --> like 
                                  % scale factors cancel out)

% b. error propagation in normal state resistance
del_vi_curr_tot=del_vi_curr+0.5*mean(errors_vi_curr_m+errors_vi_curr_p); %use mean to get a simple estimate for the error on a current in this data
del_R_N=(4/deltaI_width)*(del_vi_volt+((deltaV_height*del_vi_curr_tot)/(deltaI_width)));

%% 4. voltage modulation depth (deltaV)
% a. exclude the noise at either extreme (+/-5.08) because only the sinusoidal part of the data is of interest
distilled_vphi_volt=NaN;
distilled_vphi_flux=NaN;
n_vphi_pts=length(vphi_flux);
% the below loop (once again) uses the slower-than-optimal tactic of resizing arrays after certain loop trips, but since I don't
% know the size of the "distilled" version of these arrays in advance, I and I don't want a bunch of placeholding zeroes to be
% plotted, this is the best option.
for q=1:1:n_vphi_pts
    if(vphi_flux(q)~=-5.08 && vphi_flux(q)~=5.08)
        distilled_vphi_flux(q)=vphi_volt(q);
        distilled_vphi_volt(q)=vphi_flux(q);
    end
end

errors_distilled_vphi_voltages=del_vphi_volt*ones(length(distilled_vphi_volt),1);
errors_distilled_vphi_fluxes=del_vphi_flux*ones(length(distilled_vphi_flux),1);

% b. figure 3 -- plot the meaningful part of the data
figure
plot(distilled_vphi_volt,distilled_vphi_flux,'o','MarkerSize',5,'Color','b','MarkerFaceColor','b')
errorbar(distilled_vphi_volt,distilled_vphi_flux,-errors_distilled_vphi_fluxes,errors_distilled_vphi_fluxes,...
         -errors_distilled_vphi_voltages,errors_distilled_vphi_voltages,'.')
xlabel('flux; A e-4')
ylabel('voltage; V e-4')
title('voltage as a function of flux')
axis([-5.08,5.08,-0.04,0.04])

% c. from inspection of the plot prepared in 4a., populate a matrix of minima and the subsequent maxima for as many complete 
% cycles appear in the data set, and use these matrices to obtain several values of deltaV
v_minmag=[0.034 ,0.0336,0.0332,0.0332,0.034, 0.0336,0.0348,0.0336,0.0316,0.0348]; %using the thursday 6pm data
v_maxmag=[0.0308,0.0312,0.0308,0.0312,0.0304,0.0316,0.0316,0.0316,0.0328,0.0312];
deltaV_values=v_minmag+v_maxmag;

% d. take the final value for deltaV to be the average of the values taken from different periods
deltaV=mean(deltaV_values)/(1e4); %correct the units, as specified above

% e. error propagation in deltaV
del_deltaV=del_vphi_volt/5000;

%% 5. inductance (L)
% a. doesn't require any new calculations: these values are given + there are no errors to propagate

%% 6. device parameter (beta_L)
% a. additional useful constants
L=73e-12; %H
Phi_0=2.07e-15; %Wb (in MKS units)
k_B=1.381e-23; %J/K
T=77; %K
prefac=(4*I_C*R_N)/(pi*deltaV);
correction=3.57*sqrt(k_B*T*L)/Phi_0;

% b-ii, c-ii, and d-ii don't need gain resistor tolerance updates because del_I_C isn't impacted by including the gain resistor tolerance
% b. i. way #1: eq. 3.1 (definition of the parameter)
beta1_L=(2*I_C*L)/(Phi_0);

% b. ii. error propagation for way 1
del_beta1_L=2*L*del_I_C/Phi_0;

% c. i. way #2: eq. 3.2 (correct only when thermal noise is negligible)
beta2_L=prefac-1;

% c. ii. error propagation for way #2
del_beta2_L=(4/(pi*deltaV))*(R_N*del_I_C+I_C*del_R_N+((I_C*R_N*del_deltaV)/(deltaV))); %account for gain resistors

% d. i. way #3: eq. 3.5 (eq. 3.2 corrected for thermal noise; this value should agree more closely with that from eq. 3.1)
beta3_L=prefac*(1-correction)-1;

% d. ii. error propagation for way #2
alpha=1-correction;
del_beta3_L=((4*alpha)/(pi*deltaV))*(R_N*del_I_C+I_C*del_R_N+((I_C*R_N*del_deltaV)/(deltaV)));

%% 7. e/h from Shapiro steps
% a. prepare some arrays to plot horizontal lines to reflect the approximate location of the Shapiro steps
stepx=linspace(-1.5,1.5,5e3);
step0=linspace(-0.264,-0.264,5e3);
step1=linspace(-0.016,-0.016,5e3);
step2=linspace(0.248,0.248,5e3);
step00=linspace(-0.424,-0.424,5e3);
step11=linspace(-0.216,-0.216,5e3);
step22=linspace(-0.008,-0.008,5e3);
step33=linspace(0.208,0.208,5e3);
step44=linspace(0.408,0.408,5e3);

% b. figure 4 -- plot the V-I data

figure
subplot(1,2,1)
plot(shap_curr_acc,shap_volt_acc,'.','MarkerFaceColor','b','LineStyle','none')
hold on
errorbar(shap_curr_acc,shap_volt_acc,-errors_shap_voltages,errors_shap_voltages,...
         -errors_shap_curr12_m,errors_shap_curr12_p,'.','Color','b')
plot(stepx,step0,'-g')
text(shap_curr_acc(3),shap_volt_acc(3),'step at -0.264e-4 V','Position',[0.4,-0.2,0])
plot(stepx,step1,'-g')
text(shap_curr_acc(3),shap_volt_acc(3),'step at -0.016e-4 V','Position',[0.7,0.05,0])
plot(stepx,step2,'-g')
text(shap_curr_acc(3),shap_volt_acc(3),'step at 0.248e-4 V','Position',[-0.8,0.3,0])
xlabel('current; A e-4')
ylabel('voltage; V e-4')
title('Shapiro steps data with -18.2 dBm of 12.312 GHz MW signal')
axis([-1.5,1.5,-1.5,1.5]);
hold off

subplot(1,2,2)
plot(shap_curr_comp,shap_volt_comp,'.','MarkerFaceColor','b','LineStyle','none')
hold on
errorbar(shap_curr_comp,shap_volt_comp,-errors_shap_voltages,errors_shap_voltages,...
         -errors_shap_curr6_m,errors_shap_curr6_p,'.','Color','b')
plot(stepx,step00,'-g')
text(shap_curr_comp(3),shap_volt_comp(3),'step at -0.424e-4 V','Position',[-0.1,-0.38,0])
plot(stepx,step11,'-g')
text(shap_curr_comp(3),shap_volt_comp(3),'step at -0.216e-4 V','Position',[0.2,-0.15,0])
plot(stepx,step22,'-g')
text(shap_curr_comp(3),shap_volt_comp(3),'step at -0.008e-4 V','Position',[-1.4,0.05,0])
plot(stepx,step33,'-g')
text(shap_curr_comp(3),shap_volt_comp(3),'step at 0.208e-4 V','Position',[-0.7,0.25,0])
plot(stepx,step44,'-g')
text(shap_curr_comp(3),shap_volt_comp(3),'step at 0.408e-4 V','Position',[-0.25,0.45,0])
xlabel('current; A e-4')
ylabel('voltage; V e-4')
title('Shapiro steps data with -18.2 dBm of 9.990 GHz MW signal')
axis([-1.5,1.5,-1.5,1.5]);
hold off

% c. calculate V_0 based on averaging the amplitude of the visible steps
V_0=(1/2)*(1e-4)*(0.248+0.264);
delV_0=((1e-4)+0.008)*(1e-4); % summed contributions from the resolution of oscilloscope CSV data and +/- one vertical point 
                              % (from estimating the center of each step / point of inflection from the plot) --> this second 
                              % kind of error is taken from the plot

% d. calculate e/h
nu_a=12.312*10^9; %Hz
eoverh=(nu_a)/(2*V_0);

% e. error propagation for e/h
deleoverh= @(V_0,delV_0,nu_a,delnu_a) ((1)/(2*V_0))*((delnu_a)+((nu_a*delV_0)/(V_0)));
delnu_a=1e-5; %b/c 10 decimal places are known in MHz (uncertainty = 1e-11), 4 decimal places are known in Hz (uncertainty= 1e-5)
deleoverh_val=deleoverh(V_0,delV_0,nu_a,delnu_a);