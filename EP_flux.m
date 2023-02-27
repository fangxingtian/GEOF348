%Script computes EP fluxes from NCEP/NCAR Reanalysis 1 u,v and T
%fields from: "https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.pressure.html"
%Created 20.11.2019 Ashneel Chandra (AC)
%Modified 24.02.2021 AC
%Modified for computational efficiency 28.02.2022 AC
%Equations used to compute EP Fluxes are eq 3.1 and 3.2 from [Edmon et al., 1980]

clear all

%limits the number of cores to 4 if you are using cyclone 
%you don't need to run this if you are using your personal computer
maxNumCompThreads(4);

addpath(genpath('/Data/gfi/work/GEOF348/EP_flux/'));

%Constants
a = 6.37122e6; %radius of earth
omega = 7.2921e-05; % earth angular velocity
R = 287; %gas constant for dry air
c_p = 1004; %specific heat of dry air at constant pressure 

%This script is made to work with data for Years 2008 - 2018 daily data


%Predefine variables
T(1:144,1:73,1:17,1:4018) = NaN; 
u(1:144,1:73,1:17,1:4018) = NaN;
v(1:144,1:73,1:17,1:4018) = NaN;
theta(1:144,1:73,1:17,1:4018) = NaN;
du_dt(1:73,1:17,1:4018) = NaN;


%Read u, v and T
%Read T,u,v

path_to_data = '/Data/gfi/work/GEOF348/EP_flux';
%path_to_data = '/Users/ashneelchandra/OneDrive/Documents/PhD Research Fellow in Climate Dynamics (159124)/Data/NCEP Reanalysis 1/Pressure level variables';
T_file_path = [path_to_data,'/Air_Temperature/'];
T_list_1 = dir([path_to_data,'/Air_Temperature/air_*.nc']);
T_list_2 = char(T_list_1.name);
u_file_path = [path_to_data,'/U-Wind/'];
u_list_1 = dir([path_to_data,'/U-Wind/uwnd_*.nc']);
u_list_2 = char(u_list_1.name);
v_file_path = [path_to_data,'/V-Wind/'];
v_list_1 = dir([path_to_data,'/V-Wind/vwnd_*.nc']);
v_list_2 = char(v_list_1.name);
k=1;

Tlim = size(T_list_2,1);
for i = 1:Tlim
    i
    T_list_3 = strtrim(T_list_2(i,:));  
    Tfilename = [T_file_path T_list_3];
    u_list_3 = strtrim(u_list_2(i,:));  
    ufilename = [u_file_path u_list_3];
    v_list_3 = strtrim(v_list_2(i,:));  
    vfilename = [v_file_path v_list_3];
    
time_ = double(ncread(Tfilename,'time')/24 + datenum('01-Jan-1800'));
lat = ncread(Tfilename,'lat');
lon = ncread(Tfilename,'lon');
level = ncread(Tfilename,'level'); % millibar = hPa
air = ncread(Tfilename,'air'); %lon,lat,level,time , degK 
uwnd = ncread(ufilename,'uwnd'); %lon,lat,level,time , m/s
vwnd = ncread(vfilename,'vwnd'); %lon,lat,level,time , m/s


for j = 1:size(air,4)
    i
    j
    time(k,1)= time_(j);
    T(:,:,:,k) = air(:,:,:,j); 
    u(:,:,:,k) = uwnd(:,:,:,j);
    v(:,:,:,k) = vwnd(:,:,:,j); 
    k=k+1;
end

clear air uwnd vwnd time_   
end

time_vec=datevec(time);

%compute tendency of zonal mean wind m s^{-2}
[~,~,du_dt] = gradient(squeeze(mean(u,1,'omitnan')),1,1,(time*(24*60*60))); 

clear i k j 


%Boreal Winter (DJF)
ind_win = find(time_vec(:,2)==12 | time_vec(:,2)==1 | time_vec(:,2)==2);
%Boreal Summer (JJA)
ind_sum = find(time_vec(:,2)==6 | time_vec(:,2)==7 | time_vec(:,2)==8);

%Compute potential temperature (theta)
tic
for s = 1:size(level,1)            
    s
    theta(:,:,s,:) = squeeze(T(:,:,s,:))*(level(s)/1000)^(-R/c_p);
end
toc

%Convert from mbar/hPa to Pa
level_Pa = level*100;
lnP = log(level_Pa); %natural logarithm of pressure

%change 4th index of u, v, theta and 1st index of time here to 'ind_win' (for boreal winter), 'ind_sum' (for boreal
%summer) or as 'time' (for long term mean)

ind_season = ind_win;
%ind_season = ind_sum;
%ind_season = [1:size(time,1)]';


u_ = u(:,:,:,ind_season);
v_ = v(:,:,:,ind_season);
theta_ = theta(:,:,:,ind_season);
u_t_ = du_dt(:,:,ind_season);
time_ = time(ind_season,1);

clear u v theta time du_dt

u = u_;
v = v_;
theta = theta_;
du_dt = u_t_; 
time = time_;
clear u_ v_ theta_ time_ u_t_

%time mean of zonal mean wind tendency
du_dt_avg = squeeze(mean(du_dt,3,'omitnan'));


%Compute zonal means
tic
u_zm = squeeze(mean(u,1,'omitnan'));
v_zm = squeeze(mean(v,1,'omitnan'));
theta_zm = squeeze(mean(theta,1,'omitnan'));
toc

%Compute time mean of theta_zm (for stationary component)
theta_zm_stat = squeeze(mean(theta_zm,3,'omitnan'));

    
%Compute time anomaly of theta_zm (for transient component)
tic
for r = 1:size(time,1)
    r
    theta_zm_trans(:,:,r) = squeeze(theta_zm(:,:,r)) - theta_zm_stat;
end
toc
   

%Compute zonal mean anomalies
tic
for p = 1:size(lon,1)
    p
    u_zma(p,:,:,:) = squeeze(u(p,:,:,:)) - u_zm;
    v_zma(p,:,:,:) = squeeze(v(p,:,:,:)) - v_zm;
    theta_zma(p,:,:,:) = squeeze(theta(p,:,:,:)) - theta_zm;         
end
toc

clear v theta

%Compute time mean for zonal mean anomalies (stationary component)
 u_zma_stat = squeeze(mean(u_zma,4,'omitnan'));
 v_zma_stat = squeeze(mean(v_zma,4,'omitnan'));
 theta_zma_stat = squeeze(mean(theta_zma,4,'omitnan'));
           

    
%Compute time anomalies for zonal mean anomalies (transient component)
 tic   
for s = 1:size(time,1)
    s
    u_zma_trans(:,:,:,s) = squeeze(u_zma(:,:,:,s)) - u_zma_stat;
    v_zma_trans(:,:,:,s) = squeeze(v_zma(:,:,:,s)) - v_zma_stat;
    theta_zma_trans(:,:,:,s) = squeeze(theta_zma(:,:,:,s)) - theta_zma_stat;      
end
toc



%Compute products of zonal mean anomalies
tic
uv = u_zma.*v_zma;
v_theta = v_zma.*theta_zma;
uv_trans = u_zma_trans.*v_zma;
v_theta_trans = v_zma_trans.*theta_zma;
toc

%stationary
tic
uv_stat = u_zma_stat.*v_zma_stat;
v_theta_stat = v_zma_stat.*theta_zma_stat;
toc


%Compute zonal mean of products of anomalies
tic
uv_zm = squeeze(mean(uv,1,'omitnan'));
v_theta_zm = squeeze(mean(v_theta,1,'omitnan'));
uv_zm_trans = squeeze(mean(uv_trans,1,'omitnan'));
v_theta_zm_trans = squeeze(mean(v_theta_trans,1,'omitnan'));     
toc

%stationary
uv_zm_stat = squeeze(mean(uv_stat,1,'omitnan'));
v_theta_zm_stat = squeeze(mean(v_theta_stat,1,'omitnan'));


%Differentiate zonal mean of theta wrt lnP using centered difference method
%P or level is in Pascal
for s = 1:size(time,1)
    for p = 1:size(lat,1)
        s
        p
    k = 2;
    theta_zm_p(p,1,s) = (1/level_Pa(1))*((theta_zm(p,2,s) - theta_zm(p,1,s))/(lnP(2)-lnP(1)));
    theta_zm_p_trans(p,1,s) = (1/level_Pa(1))*((theta_zm_trans(p,2,s) - theta_zm_trans(p,1,s))/(lnP(2)-lnP(1)));
    
    for q = 1:size(level,1)-2
             
        theta_zm_p(p,k,s) = (1/level_Pa(k))*((theta_zm(p,q+2,s) - theta_zm(p,q,s))/(lnP(q+2)-lnP(q)));
        theta_zm_p_trans(p,k,s) = (1/level_Pa(k))*((theta_zm_trans(p,q+2,s) - theta_zm_trans(p,q,s))/(lnP(q+2)-lnP(q)));
        k = k+1;
    end
    theta_zm_p(p,k,s) = (1/level_Pa(k))*((theta_zm(p,k,s) - theta_zm(p,k-1,s))/(lnP(k)-lnP(k-1)));
    theta_zm_p_trans(p,k,s) = (1/level_Pa(k))*((theta_zm_trans(p,k,s) - theta_zm_trans(p,k-1,s))/(lnP(k)-lnP(k-1)));
    clear k
    end
end

%stationary
for p = 1:size(lat,1)
    p
    k = 2;
    theta_zm_p_stat(p,1) = (1/level_Pa(1))*((theta_zm_stat(p,2) - theta_zm_stat(p,1))/(lnP(2)-lnP(1)));
   
    
    for q = 1:size(level,1)-2
        p
        q     
        theta_zm_p_stat(p,k) = (1/level_Pa(k))*((theta_zm_stat(p,q+2) - theta_zm_stat(p,q))/(lnP(q+2)-lnP(q)));
       
        k = k+1;
    end
    theta_zm_p_stat(p,k) = (1/level_Pa(k))*((theta_zm_stat(p,k) - theta_zm_stat(p,k-1))/(lnP(k)-lnP(k-1)));
  
    clear k
end

  

%Compute F_phi and F_p
tic
for p = 1:size(lat,1)
    p
    F_phi(p,:,:) = -a*cosd(lat(p))*uv_zm(p,:,:) ;
    F_p(p,:,:) = (2*omega*sind(lat(p))*a*v_theta_zm(p,:,:))./theta_zm_p(p,:,:);
    F_phi_trans(p,:,:) = -a*cosd(lat(p))*uv_zm_trans(p,:,:) ;
    F_p_trans(p,:,:) = (2*omega*sind(lat(p))*a*v_theta_zm_trans(p,:,:))./theta_zm_p(p,:,:);  
end
toc

%stationary
for p = 1:size(lat,1)
    p
    F_phi_stat(p,:) = -a*cosd(lat(p))*uv_zm_stat(p,:) ;
    F_p_stat(p,:) = (2*omega*sind(lat(p))*a*v_theta_zm_stat(p,:))./theta_zm_p_stat(p,:);
end


%Compute time  mean of the EP flux components
F_phi_avg = mean(F_phi,3,'omitnan');
F_p_avg = mean(F_p,3,'omitnan');
F_phi_avg_trans = mean(F_phi_trans,3,'omitnan');
F_p_avg_trans = mean(F_p_trans,3,'omitnan');


%Compute div of F_phi 
for p = 1:size(level,1)
   
    k = 2;
    div_F_phi(1,p) = (1/(a*cosd(lat(1))))*((F_phi_avg(2,p)*cosd(lat(1)) - F_phi_avg(1,p)*cosd(lat(1)))/(lat(2)-lat(1)));
    div_F_phi_trans(1,p) = (1/(a*cosd(lat(1))))*((F_phi_avg_trans(2,p)*cosd(lat(1)) - F_phi_avg_trans(1,p)*cosd(lat(1)))/(lat(2)-lat(1)));
    div_F_phi_stat(1,p) = (1/(a*cosd(lat(1))))*((F_phi_stat(2,p)*cosd(lat(1)) - F_phi_stat(1,p)*cosd(lat(1)))/(lat(2)-lat(1)));
    
    for q = 1:size(lat,1)-2
        p
        q    
        div_F_phi(k,p) = (1/(a*cosd(lat(k))))*((F_phi_avg(q+2,p)*cosd(lat(k)) - F_phi_avg(q,p)*cosd(lat(k)))/(lat(q+2)-lat(q)));
        div_F_phi_trans(k,p) = (1/(a*cosd(lat(k))))*((F_phi_avg_trans(q+2,p)*cosd(lat(k)) - F_phi_avg_trans(q,p)*cosd(lat(k)))/(lat(q+2)-lat(q)));
        div_F_phi_stat(k,p) = (1/(a*cosd(lat(k))))*((F_phi_stat(q+2,p)*cosd(lat(k)) - F_phi_stat(q,p)*cosd(lat(k)))/(lat(q+2)-lat(q)));
        k = k+1;
    end
     div_F_phi(k,p) = (1/(a*cosd(lat(k))))*((F_phi_avg(k,p)*cosd(lat(k)) - F_phi_avg(k-1,p)*cosd(lat(k)))/(lat(k)-lat(k-1)));
     div_F_phi_trans(k,p) = (1/(a*cosd(lat(k))))*((F_phi_avg_trans(k,p)*cosd(lat(k)) - F_phi_avg_trans(k-1,p)*cosd(lat(k)))/(lat(k)-lat(k-1)));
     div_F_phi_stat(k,p) = (1/(a*cosd(lat(k))))*((F_phi_stat(k,p)*cosd(lat(k)) - F_phi_stat(k-1,p)*cosd(lat(k)))/(lat(k)-lat(k-1)));
    clear k
end


%Compute div of F_p
for p = 1:size(lat,1)
   
    k = 2;
    div_F_p(p,1) = (F_p_avg(p,2) - F_p_avg(p,1))/(level_Pa(2)-level_Pa(1));
    div_F_p_trans(p,1) = (F_p_avg_trans(p,2) - F_p_avg_trans(p,1))/(level_Pa(2)-level_Pa(1));
    div_F_p_stat(p,1) = (F_p_stat(p,2) - F_p_stat(p,1))/(level_Pa(2)-level_Pa(1));
    
    for q = 1:size(level,1)-2
        p
        q    
        div_F_p(p,k) = (F_p_avg(p,q+2) - F_p_avg(p,q))/(level_Pa(q+2)-level_Pa(q));
        div_F_p_trans(p,k) = (F_p_avg_trans(p,q+2) - F_p_avg_trans(p,q))/(level_Pa(q+2)-level_Pa(q));
        div_F_p_stat(p,k) = (F_p_stat(p,q+2) - F_p_stat(p,q))/(level_Pa(q+2)-level_Pa(q));
        k = k+1;
    end
    div_F_p(p,k) = (F_p_avg(p,k) - F_p_avg(p,k-1))/(level_Pa(k)-level_Pa(k-1));
    div_F_p_trans(p,k) = (F_p_avg_trans(p,k) - F_p_avg_trans(p,k-1))/(level_Pa(k)-level_Pa(k-1));
    div_F_p_stat(p,k) = (F_p_stat(p,q+2) - F_p_stat(p,q))/(level_Pa(q+2)-level_Pa(q));
    clear k
end


%Compute div of F
div_F = div_F_phi + div_F_p;
div_F_trans = div_F_phi_trans + div_F_p_trans;
div_F_stat = div_F_phi_stat + div_F_p_stat;


%scale averaged EP flux components for plotting
for p = 1:size(lat,1)
    for q = 1:size(level,1)
        p
        q
        scaled_F_phi(p,q) = cosd(lat(p))*sqrt(1000/level_Pa(q))*F_phi_avg(p,q)/(a*pi); 
        scaled_F_p(p,q) = cosd(lat(p))*sqrt(1000/level_Pa(q))*F_p_avg(p,q)/100000;
        scaled_F_phi_trans(p,q) = cosd(lat(p))*sqrt(1000/level_Pa(q))*F_phi_avg_trans(p,q)/(a*pi); 
        scaled_F_p_trans(p,q) = cosd(lat(p))*sqrt(1000/level_Pa(q))*F_p_avg_trans(p,q)/100000;
        scaled_F_phi_stat(p,q) = cosd(lat(p))*sqrt(1000/level_Pa(q))*F_phi_stat(p,q)/(a*pi); 
        scaled_F_p_stat(p,q) = cosd(lat(p))*sqrt(1000/level_Pa(q))*F_p_stat(p,q)/100000;
    end
end



%Plot EP flux vectors as arrows and divergence as contours
vneg = [-800,-700,-600,-500,-400,-300,-200,-100];
v_o = [0,0];
vpos = [100,200,300,400,500,600,700,800];
dd = 1e-3;
cc=[-5:0.5:5]*dd;
scale = [cc(1,1) cc(1,size(cc,2))];
[X,Y] = meshgrid(lat,level);
X = X';
Y = Y';
factor=1e3;

hold on
subplot(3,3,1)
hold on
imagescn(lat,level,du_dt_avg'*factor);
contour(X,Y,div_F,vneg,'--k','LineWidth',1.5,'ShowText','on') %you can change _stat to _trans to plot transient component
contour(X,Y,div_F,v_o,'-k','LineWidth',2,'ShowText','on')
contour(X,Y,div_F,vpos,'-k','LineWidth',1.5,'ShowText','on') %you can change _stat to _trans to plot transient component
%%contour(X,Y,div_F,'LineWidth',3,'ShowText','on')
hold on
%pcolor(X,Y,du_dt_avg*1e3)
quiversc(X,Y,scaled_F_phi,scaled_F_p,'y','density',50,'MaxHeadSize',1.2,'AutoScaleFactor',1,'LineWidth',1.5) %you can change _stat to _trans to plot transient component
%%qq=quiver(X,Y,scaled_F_phi,scaled_F_p,'Color','y','LineWidth',1.5,'MaxHeadSize',0.5,'AutoScale','on','AutoScaleFactor',1.8);
hold on
cb=colorbar('southoutside','Ticks',cc,'FontSize',14);
cb.Label.String = '\bfm s^{-2}';
caxis(scale)
set(gca, 'ytick',[50 100 200 300 500 700 1000], 'ylim',[50 1000], 'YScale', 'log', 'Ydir', 'reverse', 'Fontsize', 14,'FontWeight','bold')
set(gca, 'xlim',[-90 10]); %you can change longitude limits here to show northern or southern hemishpere
cmocean('balance',10*numel(cc)-10)
colorbar('off')
title({['South'];['Total']}, 'Fontsize', 18,'FontWeight','bold')
%xlabel('Latitude', 'Fontsize', 15,'FontWeight','bold')
ylabel('Level (hPa)',  'Fontsize', 15,'FontWeight','bold')
box on
hold on

subplot(3,3,4)
hold on
imagescn(lat,level,du_dt_avg'*factor);
contour(X,Y,div_F_trans,vneg,'--k','LineWidth',1.5,'ShowText','on') %you can change _stat to _trans to plot transient component
contour(X,Y,div_F_trans,v_o,'-k','LineWidth',2,'ShowText','on')
contour(X,Y,div_F_trans,vpos,'-k','LineWidth',1.5,'ShowText','on') %you can change _stat to _trans to plot transient component
%%contour(X,Y,div_F_trans,'LineWidth',3,'ShowText','on')
hold on
quiversc(X,Y,scaled_F_phi_trans,scaled_F_p_trans,'y','density',50,'MaxHeadSize',1.2,'AutoScaleFactor',1,'LineWidth',1.5) %you can change _stat to _trans to plot transient component
%%qq=quiver(X,Y,scaled_F_phi_trans,scaled_F_p_trans,'Color','y','LineWidth',1.5,'MaxHeadSize',0.5,'AutoScale','on','AutoScaleFactor',1.8);
hold on
cb=colorbar('southoutside','Ticks',cc,'FontSize',14);
cb.Label.String = '\bfm s^{-2}';
caxis(scale)
set(gca, 'ytick',[50 100 200 300 500 700 1000], 'ylim',[50 1000], 'YScale', 'log', 'Ydir', 'reverse', 'Fontsize', 14,'FontWeight','bold')
set(gca, 'xlim',[-90 10]); %you can change longitude limits here to show northern or southern hemishpere
cmocean('balance',10*numel(cc)-10)
colorbar('off')
title('Transient', 'Fontsize', 18,'FontWeight','bold')
%xlabel('Latitude', 'Fontsize', 15,'FontWeight','bold')
ylabel('Level (hPa)',  'Fontsize', 15,'FontWeight','bold')
box on
hold on

subplot(3,3,7)
hold on
imagescn(lat,level,du_dt_avg'*factor);
contour(X,Y,div_F_stat,vneg,'--k','LineWidth',1.5,'ShowText','on') %you can change _stat to _trans to plot transient component
contour(X,Y,div_F_stat,v_o,'-k','LineWidth',2,'ShowText','on')
contour(X,Y,div_F_stat,vpos,'-k','LineWidth',1.5,'ShowText','on') %you can change _stat to _trans to plot transient component
%%contour(X,Y,div_F_stat,'LineWidth',3,'ShowText','on')
hold on
quiversc(X,Y,scaled_F_phi_stat,scaled_F_p_stat,'y','density',50,'MaxHeadSize',1.2,'AutoScaleFactor',1,'LineWidth',1.5) %you can change _stat to _trans to plot transient component
%%qq=quiver(X,Y,scaled_F_phi_stat,scaled_F_p_stat,'Color','y','LineWidth',1.5,'MaxHeadSize',0.5,'AutoScale','on','AutoScaleFactor',1.8);
hold on
cb=colorbar('southoutside','Ticks',cc,'FontSize',14);
cb.Label.String = '\bfm s^{-2}';
caxis(scale)
set(gca, 'ytick',[50 100 200 300 500 700 1000], 'ylim',[50 1000], 'YScale', 'log', 'Ydir', 'reverse', 'Fontsize', 14,'FontWeight','bold')
set(gca, 'xlim',[-90 10]); %you can change longitude limits here to show northern or southern hemishpere
cmocean('balance',10*numel(cc)-10)
title('Stationary', 'Fontsize', 18,'FontWeight','bold')
xlabel('Latitude', 'Fontsize', 15,'FontWeight','bold')
ylabel('Level (hPa)',  'Fontsize', 15,'FontWeight','bold')
box on
hold on

hold on
subplot(3,3,2)
hold on
imagescn(lat,level,du_dt_avg'*factor);
contour(X,Y,div_F,vneg,'--k','LineWidth',1.5,'ShowText','on') %you can change _stat to _trans to plot transient component
contour(X,Y,div_F,v_o,'-k','LineWidth',2,'ShowText','on')
contour(X,Y,div_F,vpos,'-k','LineWidth',1.5,'ShowText','on') %you can change _stat to _trans to plot transient component
%%contour(X,Y,div_F,'LineWidth',3,'ShowText','on')
hold on
quiversc(X,Y,scaled_F_phi,scaled_F_p,'y','density',50,'MaxHeadSize',1.2,'AutoScaleFactor',1,'LineWidth',1.5) %you can change _stat to _trans to plot transient component
%%qq=quiver(X,Y,scaled_F_phi,scaled_F_p,'Color','y','LineWidth',1.5,'MaxHeadSize',0.5,'AutoScale','on','AutoScaleFactor',1.8);
hold on
cb=colorbar('southoutside','Ticks',cc,'FontSize',14);
cb.Label.String = '\bfm s^{-2}';
caxis(scale)
set(gca, 'ytick',[50 100 200 300 500 700 1000], 'ylim',[50 1000], 'YScale', 'log', 'Ydir', 'reverse', 'Fontsize', 14,'FontWeight','bold')
set(gca, 'xlim',[-10 90]); %you can change longitude limits here to show northern or southern hemishpere
cmocean('balance',10*numel(cc)-10)
colorbar('off')
title({['North'];['Total']}, 'Fontsize', 18,'FontWeight','bold')
%xlabel('Latitude', 'Fontsize', 15,'FontWeight','bold')
ylabel('Level (hPa)',  'Fontsize', 15,'FontWeight','bold')
box on
hold on

subplot(3,3,5)
hold on
imagescn(lat,level,du_dt_avg'*factor);
contour(X,Y,div_F_trans,vneg,'--k','LineWidth',1.5,'ShowText','on') %you can change _stat to _trans to plot transient component
contour(X,Y,div_F_trans,v_o,'-k','LineWidth',2,'ShowText','on')
contour(X,Y,div_F_trans,vpos,'-k','LineWidth',1.5,'ShowText','on') %you can change _stat to _trans to plot transient component
%%contour(X,Y,div_F_trans,'LineWidth',3,'ShowText','on')
hold on
quiversc(X,Y,scaled_F_phi_trans,scaled_F_p_trans,'y','density',50,'MaxHeadSize',1.2,'AutoScaleFactor',1,'LineWidth',1.5) %you can change _stat to _trans to plot transient component
%%qq=quiver(X,Y,scaled_F_phi_trans,scaled_F_p_trans,'Color','y','LineWidth',1.5,'MaxHeadSize',0.5,'AutoScale','on','AutoScaleFactor',1.8);
hold on
cb=colorbar('southoutside','Ticks',cc,'FontSize',14);
cb.Label.String = '\bfm s^{-2}';
caxis(scale)
set(gca, 'ytick',[50 100 200 300 500 700 1000], 'ylim',[50 1000], 'YScale', 'log', 'Ydir', 'reverse', 'Fontsize', 14,'FontWeight','bold')
set(gca, 'xlim',[-10 90]); %you can change longitude limits here to show northern or southern hemishpere
cmocean('balance',10*numel(cc)-10)
colorbar('off')
title('Transient', 'Fontsize', 18,'FontWeight','bold')
%xlabel('Latitude', 'Fontsize', 15,'FontWeight','bold')
ylabel('Level (hPa)',  'Fontsize', 15,'FontWeight','bold')
box on
hold on

subplot(3,3,8)
hold on
imagescn(lat,level,du_dt_avg'*factor);
contour(X,Y,div_F_stat,vneg,'--k','LineWidth',1.5,'ShowText','on') %you can change _stat to _trans to plot transient component
contour(X,Y,div_F_stat,v_o,'-k','LineWidth',2,'ShowText','on')
contour(X,Y,div_F_stat,vpos,'-k','LineWidth',1.5,'ShowText','on') %you can change _stat to _trans to plot transient component
%%contour(X,Y,div_F_stat,'LineWidth',3,'ShowText','on')
hold on
quiversc(X,Y,scaled_F_phi_stat,scaled_F_p_stat,'y','density',50,'MaxHeadSize',1.2,'AutoScaleFactor',1,'LineWidth',1.5) %you can change _stat to _trans to plot transient component
%%qq=quiver(X,Y,scaled_F_phi_stat,scaled_F_p_stat,'Color','y','LineWidth',1.5,'MaxHeadSize',0.5,'AutoScale','on','AutoScaleFactor',1.8);
hold on
cb=colorbar('southoutside','Ticks',cc,'FontSize',14);
cb.Label.String = '\bfm s^{-2}';
caxis(scale)
set(gca, 'ytick',[50 100 200 300 500 700 1000], 'ylim',[50 1000], 'YScale', 'log', 'Ydir', 'reverse', 'Fontsize', 14,'FontWeight','bold')
set(gca, 'xlim',[-10 90]); %you can change longitude limits here to show northern or southern hemishpere
cmocean('balance',10*numel(cc)-10)
title('Stationary', 'Fontsize', 18,'FontWeight','bold')
xlabel('Latitude', 'Fontsize', 15,'FontWeight','bold')
ylabel('Level (hPa)',  'Fontsize', 15,'FontWeight','bold')
box on
hold on

%export_fig('EP_Flux_DJF_2008-2018.jpg','-transparent','-r600');
