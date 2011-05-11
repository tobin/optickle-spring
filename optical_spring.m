%% Modeling an Optical Spring
% This script builds a simple two-mirror cavity in Optickle and then finds
% the (opto)mechanical transfer functions at various cavity detunings.
%
% Tobin Fricke
% February 10, 2011

%% Parameters
% Here we specify some basic parameters of the optical setup, such as the
% reflectivities of the optics.  These will be used later, when the
% Optickle model is created.

vFrf     = 0;           % carrier only
P_laser  = 12;          % Laser power [Watts]

par.IX.T = 0.02930;     % Power transmission coefficient
par.EX.T = 10e-6;

par.w = 2 * pi * 0.74;  % Suspension resonant frequency [rad/s]
par.mass = 10;          % Test mass mass [kg]


%% Build the Optickle model

% Create the initial Optickle object by giving Optickle a vector of the RF
% frequencies that will be in the model.  In our case, we have only the
% carrier, so vFrf = [0].  (These frequencies are all relative to the
% carrier.)
opt = Optickle(vFrf);

c = opt.c;              % speed of light [m/2]
lambda = opt.lambda;    % wavelength [m]

% Add the optics to the model
opt = addSource(opt, 'Laser', sqrt(P_laser));
%[opt, sn] = addMirror(opt, name, aio, Chr, Thr, Lhr, Rar, Lmd, Nmd)
opt = addMirror(opt, 'ITM', 0, 0, par.IX.T, 0, 0, 0);
opt = addMirror(opt, 'ETM', 0, 0, par.EX.T, 0, 0, 0);

% Link them together
opt = addLink(opt, 'Laser', 'out', 'ITM', 'bk', 0);
opt = addLink(opt, 'ITM', 'fr', 'ETM', 'fr', 3995);
opt = addLink(opt, 'ETM', 'fr', 'ITM', 'fr', 3995);

% We don't actually need any probes for our purposes here, but Optickle
% complains if we don't add any.  (Probes are like photodiodes, except that
% they can observe the fields without disturbing them.)
opt = addProbeOut(opt, 'CAVITY DC',  'ITM', 'fr', 0, 0);

% Frequencies of interest
f = logspace(log10(0.1), log10(10), 101);

% Tell Optickle the mechanical transfer functions of the suspensions
dampRes = [0.01 + 1i, 0.01 - 1i];
mechTFobj = zpk([], -par.w * dampRes, 1 / par.mass);
opt = setMechTF(opt, 'ITM', mechTFobj);
opt = setMechTF(opt, 'ETM', mechTFobj);
mechTF = squeeze(freqresp(mechTFobj, 2*pi*f));

% for convenience, get the mirror amplitude reflectivities:
t1 = sqrt(par.IX.T);
t2 = sqrt(par.EX.T);
r1 = sqrt(1 - par.IX.T);
r2 = sqrt(1 - par.EX.T);

% Compute the finesse and cavity gain:
F = 4 * r1 * r2 / (1 - r1*r2)^2;
finesse = (pi/2)*sqrt(F);
g = t1 / (1 - r1*r2);

%% sweepLinear
% To validate the analytical model, check that it agrees with the DC fields
% computed by Optickle over various detunings.

detuning = 20e-12;
Nsweep = 31;
pos = zeros(opt.Ndrive, 1);
pos(getDriveNum(opt, 'ITM')) = detuning/2;
pos(getDriveNum(opt, 'ETM')) = detuning/2;

dx = linspace(-detuning, detuning, 101);
k = (2*pi)/lambda;

[xPos, sigDC, fDC] = sweepLinear(opt, -pos, pos, Nsweep);

plot(1e12 * dx, g^2./(1 + F*sin(k*dx).^2), '-', 'LineWidth', 2.5);
hold all
plot(1e12 * (xPos(getDriveNum(opt, 'ITM'), :) + xPos(getDriveNum(opt, 'ETM'), :)), ...
     sigDC(getProbeNum(opt, 'CAVITY DC'), :)/P_laser, 'o', ...
     'LineWidth', 2.5, 'MarkerSize', 5);
hold off
xlabel('detuning (picometers)');
ylabel('cavity power buildup');
legend('Analytic', 'Optickle', 'Location', 'Best');

%% Plot the derivative

x = (xPos(getDriveNum(opt, 'ITM'), :) + xPos(getDriveNum(opt, 'ETM'), :));

dx = x;
phi = dx * k;
plot(1e12*dx,  - P_laser * (2*F*g^2) * cos(phi) .* sin(phi) * k ./ (1 + F*sin(phi).^2).^2, '-', ...
     'LineWidth', 2.5);
hold all

P = sigDC(getProbeNum(opt, 'CAVITY DC'), :);
plot( 1e12 * (x(2:end)+(x(1:(end-1))))/2,  diff(P)./diff(x), 'o', 'LineWidth', 2.5,'MarkerSize', 5);
hold off
grid on;

%% Tickle
% Now the Optickle model is constructed, and all that's left is to call
% the "tickle" function to have it compute fields and transfer functions.

results = struct('pos', [],'mMech', [], 'mMechDOF', []);

detunings = linspace(-10, 10, 11)*1e-12;
nETM = getDriveNum(opt, 'ETM');
nITM = getDriveNum(opt, 'ITM');
% Loop over various a set of cavity detunings (in [pico]meters)
for ii=1:length(detunings),    
    detuning = detunings(ii);
    
    pos = zeros(opt.Ndrive, 1);    
    pos(nITM) = detuning/2;
    pos(nETM) = detuning/2;
      
    % Call tickle -- this is where all computations happen
    [fDC, sigDC, sigAC, mMech, noiseAC, noiseMech] = tickle(opt, pos, f);
    
    result.pos = pos;
    result.mMech = mMech;

    % Change to the common/differential basis.   Note that the coordinate 
    % for the ETM is opposite that of the ITM, because they are facing 
    % opposite directions.  So ETM+ITM gives the differential mode.
    S = [1 -1; 1 1];    
    mMechDOF = zeros(size(mMech));
    jj = [nITM nETM];
    for kk=1:size(mMech,3)
        mMechDOF(jj,jj,kk) = S * mMech(jj,jj,kk) / S;
    end
    
    result.mMechDOF = mMechDOF;
    results(ii) = result;
end

%% Plot the transfer functions in the mirror basis
close all;
clim = [-10 10];
caxis(clim);
cm   = colormap();

for ii=1:length(results);
    pos = results(ii).pos;
    mMech = results(ii).mMech;
    
    detuning = pos(nETM) + pos(nITM);
    
    % pick a pretty color
    color = interp1(linspace(clim(1), clim(2), size(cm,1)), cm, detuning*1e12);
         
    dofnames = {'ITM', 'ETM'};
    ax = [0 0 0 0];
    for nFrom=1:2,
        for nTo = 1:2,
            plotnum = (nTo-1)*2 + nFrom;
            ax(plotnum) = subplot(2,2, plotnum);
                        
            % Plot the Optickle results
            semilogx(f, db(getTF(mMech,nTo,nFrom) .* mechTF), ...
                  '-', 'Color', color, 'MarkerSize', 5, 'LineWidth', 2.5);
            hold all
            title(sprintf('%s to %s', dofnames{nFrom}, dofnames{nTo}));
        end
    end
end
linkaxes(ax, 'xy');
%% Plot the transfer functions in the common/differential basis
close all;
clim = [-10 10];
caxis(clim);
cm   = colormap();

for ii=1:length(results);
    pos = results(ii).pos;
    mMechDOF = results(ii).mMechDOF;
    
    detuning = pos(nETM) + pos(nITM);
    
    % pick a pretty color
    color = interp1(linspace(clim(1), clim(2), size(cm,1)), cm, detuning*1e12);
         
    dofnames = {'common', 'differential'};
    ax = [0 0 0 0];
    for nFrom=1:2,
        for nTo = 1:2,
            plotnum = (nTo-1)*2 + nFrom;
            ax(plotnum) = subplot(2,2, plotnum);
                        
            % Plot the Optickle results
            semilogx(f, db(getTF(mMechDOF,nTo,nFrom) .* mechTF), ...
                  '-', 'Color', color, 'MarkerSize', 5, 'LineWidth', 2.5);
            hold all
            title(sprintf('%s to %s', dofnames{nFrom}, dofnames{nTo}));
        end
    end
end
linkaxes(ax, 'xy');
%%  Plot the differential mode along with the prediction
%colorbar;
%xlabel('frequency [Hz]');
%ylabel('differential mode transfer function');
%title('mechanical transfer function versus cavity detuning [picometers]');
%legend('Analytic model','Optickle simulation', 'location', 'best');
% Make the fonts bigger
%set(findall(gcf, 'Type','text'), 'FontSize', 12)

close all

% Pick a pretty color
clim = [-10 10];
caxis(clim);
cm   = colormap();
for ii=1:length(results);
    pos = results(ii).pos;
    mMechDOF = results(ii).mMechDOF;
    detuning = pos(nETM) + pos(nITM);
    
    color = interp1(linspace(clim(1), clim(2), size(cm,1)), cm, detuning*1e12);
    
    % Calculate the optical spring constant
    phi = detuning * (2*pi/lambda);
    k_opt1 = -(2*P_laser/c) * (2*F*g^2) * (2*pi/lambda)^2 * detuning;
    k_opt2 = -(2*P_laser/c) * (2*F*g^2) * cos(phi) * sin(phi) * (2*pi/lambda) / (1 + F*sin(phi)^2)^2;
    k_opt = k_opt2;

    % Calculate the expected mechanical transfer function in the presense
    % of radiation pressure.    We assume that the damping coefficient is
    % unmodified, and the only action of the optical spring is to change
    % the effective spring constant.
%      poles = -par.w * dampRes;    
%      poles = real(poles) + 1i*[-1 1].*sqrt(imag(poles).^2 + k_opt/par.mass);     
    poles = -sqrt(par.w^2 + k_opt/par.mass) * dampRes;
    sys = zpk([], poles, 1 / par.mass);
    
    % Plot using more points, for increased prettiness
    f2 = logspace(log10(min(f)), log10(max(f)), 1001);
    sysTF = squeeze(freqresp(sys, 2*pi*f2));
    semilogx(f2, db(sysTF), '-', 'Color', (color + [1 1 1])/2, 'LineWidth', 2.5);
    hold all
    
    % Plot the Optickle results
    nTo = 2;
    nFrom = 2;
    semilogx(f, db(getTF(mMechDOF,nTo,nFrom).*mechTF), ...
        'o', 'Color', color, 'MarkerSize', 5, 'LineWidth', 2.5);
    title(sprintf('%s to %s', dofnames{nFrom}, dofnames{nTo}));
end

colorbar;
xlabel('frequency [Hz]');
ylabel('differential mode transfer function');
title('mechanical transfer function versus cavity detuning [picometers]');
legend('Analytic model','Optickle simulation', 'location', 'best');
% Make the fonts bigger
set(findall(gcf, 'Type','text'), 'FontSize', 12)
hold off