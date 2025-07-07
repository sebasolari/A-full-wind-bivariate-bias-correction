% Here, user should load measured wind data
% It is assumed that 'dateme' was constructed using datetime function
load('wind_measured.mat','dateme','vme','dme');

% Here, user should load the reanalysis wind data
% It is assumed that 'datere' was constructed using datetime function
load('wind_reanalysis.mat','datere','vre','dre');

% Interpolate reanalysis to measured time points
vrea = interp1(datere,vre,dateme);
drea = interp1(datere,dre,dateme);

% Bias correction BC-01 in Solari & Alonso (2025)
dcen = 0:5:360;
dven = 30;
prob = [0.01 0.05 0.1:.1:0.9 0.95 0.99]';
fun  = fittype('A*x^B+C');
a    = zeros(numel(dcen),1);
b    = zeros(numel(dcen),1);
c    = zeros(numel(dcen),1);
for id = 1:numel(dcen)
    daux = wrapTo180(drea-dcen(id));
    idau = find(daux>-dven/2 & daux<dven/2);
    qx   = quantile(vrea(idau),prob);
    qy   = quantile(vme(idau),prob);
    out  = fit(qx,qy,fun,fitoptions('METHOD','NonlinearLeastSquares','StartPoint',[1 1 0]));
    a(id) = out.A;
    b(id) = out.B;
    c(id) = out.C;
end
ft              = fittype('fourier6');
opts            = fitoptions(ft);
opts.Display    = 'Off';
opts.Lower      = [-Inf.*ones(1,13) 2*pi/362];
opts.StartPoint = [zeros(1,13) 2*pi];
opts.Upper      = [Inf.*ones(1,13) 2*pi/358];
fitres_A        = fit(dcen(:),a,ft,opts);
fitres_B        = fit(dcen(:),b,ft,opts);
fitres_C        = fit(dcen(:),c,ft,opts);
a               = fitres_A(dre);
b               = fitres_B(dre);
c               = fitres_C(dre);
vre             = a.*vre.^b+c;
vre             = max(0,vre);
figure

% Saves bias corected reanalysis
save wind_reanalysis_bc01.mat fechare vre dre;
