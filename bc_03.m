% Here, user should load measured wind data
% It is assumed that 'dateme' was constructed using datetime function
load('wind_measured.mat','dateme','vme','dme');

% Here, user should load the reanalysis wind data
% It is assumed that 'datere' was constructed using datetime function
load('wind_reanalysis.mat','datere','vre','dre');

% Interpolate reanalysis to measured time points
vrea = interp1(datere,vre,dateme);
drea = interp1(datere,dre,dateme);

% Bias correction following Solari & Alonso 2025 (BC_03)
% Iterates 500 times in the IRPM. User should check if this is enough or
% requires more iteratuins <--- IMPORTANT
xs = vrea.*[cosd(drea) sind(drea)];
xt = vme.*[cosd(dme) sind(dme)];
xc = xs;
for iter = 1:500
    [Q,R] = qr(randn(2));
    Q     = Q*diag(sign(diag(R)));
    xca   = xc*Q;
    xta   = xt*Q;
    xca   = qqmap(xca,xta);
    xc    = xca/Q;
end
[dc,vc] = cart2pol(xc(:,1),xc(:,2));
dc      = wrapTo360(rad2deg(dc));
% Fits a GP regression model and uses it to correct the full reanalysis
mod_dv = fitrgp([vrea drea],vc-vrea,'Standardize',true);
mod_dd = fitrgp([vrea drea],wrapTo180(dc-drea),'Standardize',true);
dvre   = predict(mod_dv,[vre dre]);
ddre   = predict(mod_dd,[vre dre]);
vre    = vre + dvre;
dre    = wrapTo360(dre + ddre);

% Saves bias corected reanalysis
save wind_reanalysis_bc02.mat fechare vre dre;