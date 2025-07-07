% Here, user should load measured wind data
% It is assumed that 'dateme' was constructed using datetime function
load('wind_measured.mat','dateme','vme','dme');

% Here, user should load the reanalysis wind data
% It is assumed that 'datere' was constructed using datetime function
load('wind_reanalysis.mat','datere','vre','dre');

% Interpolate reanalysis to measured time points
vrea = interp1(datere,vre,dateme);
drea = interp1(datere,dre,dateme);

% Bias correction BC_02 in Solari & Alonso (2025) (following Lemos et al.
% 2020)
% wind speed
pmin = .01;
pmax = .9999; 
prob = exp(-exp(-linspace(-log(-log(pmin)),-log(-log(pmax)),20)));
qx   = quantile(vrea,prob);
qy   = quantile(vme,prob);
vre  = interp1(qx,qy,vre,'linear','extrap');
% wind direction
prob = .01:.01:.99;
qux  = quantile(cosd(drea),prob);
quy  = quantile(cosd(dme),prob);
qvx  = quantile(sind(drea),prob);
qvy  = quantile(sind(dme),prob);
u2   = interp1(qux,quy,cosd(dre),'linear','extrap');
v2   = interp1(qvx,qvy,sind(dre),'linear','extrap');
dre  = wrapTo360(atan2d(v2,u2));

% Saves bias corected reanalysis
save wind_reanalysis_bc02.mat fechare vre dre;
