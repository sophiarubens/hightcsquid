%% formulas
pcterr= @(exp,theo) abs(exp-theo)/theo*100;
del_eoverh= @(V_0,delV_0,nu_a,delnu_a) ((1)/(2*V_0))*((delnu_a)+((nu_a*delV_0)/(V_0)));

%% set up receptacles for results and their associated errors
ntrials=14;
delnu_a=1e-5; %b/c 10 decimal places are known in MHz (uncertainty = 1e-11), 4 decimal places are known in Hz (uncertainty = 1e-5)
eoverh_values=zeros(ntrials,1);
deleoverh_values=zeros(ntrials,1);

%% analyze each data set in turn
figure
for k=1:1:ntrials
    if(k==1) %friday
        shap_data=readmatrix('shapiro_21apr.csv');
        nu_app=9.990;
        pwr='-17.0';
        V_0=(1/4)*(1e-4)*(0.44+0.4);
        delV_0=0.001+0.02; %each deltaV has summed contributions from the resolution of oscilloscope CSV data and +/- one vertical point (from estimating the center of each step / point of inflection from the plot) --> this second kind of error is taken from the plot
    elseif(k==2) %tuesday 1 (number represents order)
        shap_data=readmatrix('shapiro1_am25apr.csv');
        nu_app=10.109;
        pwr='-18.0';
        V_0=(1/4)*(1e-4)*(0.392+0.424);
        delV_0=0.0001+0.008;
    elseif(k==3) %tuesday 2 (number represents order)
        shap_data=readmatrix('shapiro2_am25apr.csv');
        nu_app=10.120;
        pwr='-18.0';
        V_0=(1/4)*(1e-4)*(0.4+0.424);
        delV_0=0.0001+0.008;
    elseif(k==4) %tuesday 4 (from here onward, the number represents the file number)
        shap_data=readmatrix('shapiro4_25apr.csv');
        nu_app=10.000;
        pwr='-18.2';
        V_0=(1/4)*(1e-4)*(0.416+0.456);
        delV_0=(1e-4)+0.008;
    elseif(k==5) %tuesday 5
        shap_data=readmatrix('shapiro5_25apr.csv');
        nu_app=10.009;
        pwr='-18.2';
        V_0=(1/4)*(1e-4)*(0.4+0.416);
        delV_0=(1e-4)+0.008;
    elseif(k==6) %tuesday 6
        shap_data=readmatrix('shapiro6_25apr.csv');
        nu_app=9.990;
        pwr='-18.2';
        V_0=(1/4)*(1e-4)*(0.408+0.424);
        delV_0=(1e-4)+0.008;
    elseif(k==7) %tuesday 7
        shap_data=readmatrix('shapiro7_25apr.csv');
        nu_app=10.013;
        pwr='-18.2';
        V_0=(1/4)*(1e-4)*(0.408+0.432);
        delV_0=(1e-4)+0.008;
    elseif(k==8) %tuesday 8
        shap_data=readmatrix('shapiro8_25apr.csv');
        nu_app=11.271;
        pwr='-18.2';
        V_0=(1/4)*(1e-4)*(0.344+0.248);
        delV_0=(1e-4)+0.008;
    elseif(k==9) %tuesday 9
        shap_data=readmatrix('shapiro9_25apr.csv');
        nu_app=11.691;
        pwr='-18.2';
        V_0=(1/4)*(1e-4)*(0.48+0.496);
        delV_0=(1e-4)+0.008;
    elseif(k==10) %tuesday 10
        shap_data=readmatrix('shapiro10_25apr.csv');
        nu_app=11.702;
        pwr='-18.2';
        V_0=(1/2)*(1e-4)*(0.24+0.24);
        delV_0=(1e-4)+0.008;
    elseif(k==11) %tuesday 11
        shap_data=readmatrix('shapiro11_25apr.csv');
        nu_app=11.835;
        pwr='-18.2';
        V_0=(1/2)*(1e-4)*(0.232+0.248);
        delV_0=(1e-4)+0.008;
    elseif(k==12) %tuesday 12  !! this trial has the lowest associated error ... although only two steps were visible here, they were far clearer than any of
                  %            !! the data sets where four steps were visible
        shap_data=readmatrix('shapiro12_25apr.csv');
        nu_app=12.312;
        pwr='-18.2';
        V_0=(1/2)*(1e-4)*(0.248+0.264);
        delV_0=(1e-4)+0.008;
    elseif(k==13) %tuesday 13
        shap_data=readmatrix('shapiro13_25apr.csv');
        nu_app=12.472;
        pwr='-18.2';
        V_0=(1/2)*(1e-4)*(0.272+0.272);
        delV_0=(1e-4)+0.008;
    elseif(k==14) %tuesday 14
        shap_data=readmatrix('shapiro14_25apr.csv');
        nu_app=12.512;
        pwr='-18.2';
        V_0=(1/2)*(1e-4)*(0.272+0.272);
        delV_0=(1e-4)+0.008;
    end
    delV_0=delV_0*1e-4;
    
    %% calculate e/h
    eoverh=(nu_app*10^9)/(2*V_0); %since applied frequencies were stored in GHz

    %% propagate error
    deleoverh=del_eoverh(V_0,delV_0,nu_app*10^9,delnu_a);

    %% store in value-holders
    eoverh_values(k)=eoverh;
    deleoverh_values(k)=deleoverh;

    %% add to plot
    titlestring=[pwr,' dBm of ',num2str(nu_app),' GHz'];
    shap_curr=shap_data(:,1)';
    shap_volt=shap_data(:,2)';
    subplot(3,5,k)
    plot(shap_curr,shap_volt,'.','MarkerFaceColor','b','LineStyle','none')
    xlabel('current; A e-4')
    ylabel('voltage; V e-4')
    title(titlestring)
end
sgtitle('several sets of Shapiro steps data')

%% percent error
eoverh_theo=2.4179671e14;
pcterr_alltrials=pcterr(eoverh_values,eoverh_theo);