%% Radiation pressure check
% To check whether Optickle computes the radiation pressure force in the
% way I expect, this script constructs a simple model consisting of only a
% laser, an AM modulator, and a perfectly-reflective mirror.   We ask
% Optickle for the transfer function from power fluctuations to mirror
% position, multiply by (2/c) to get the radiation pressure to position
% transfer function, and check whether this agrees with the mechanical
% transfer function we supplied.
%
% Tobin Fricke - March 15, 2011

%% Laser, Modulator, Mirror
% Construct the Optickle model

vFrf     = 0;               % carrier only
P_laser  = 12;              % Laser power [Watts]
par.w     = 2 * pi * 0.74;  % Suspension resonant frequency [rad/s]
par.mass = 10;              % Test mass mass [kg]

% Build the Optickle model
opt = Optickle(vFrf);

% Add the optics to the model
opt = addSource(opt, 'Laser', sqrt(P_laser));
opt = addModulator(opt, 'AM', 1);
opt = addMirror(opt, 'Mirror', 0, 0, 0, 0, 0, 0);  % perfect mirror

% Link them together
opt = addLink(opt, 'Laser', 'out', 'AM', 'in', 0);
opt = addLink(opt, 'AM', 'out', 'Mirror', 'fr', 0);

% Probes
opt = addProbeIn(opt, 'Mirror DC',  'Mirror', 'fr', 0, 0);

% Frequencies of interest
f = logspace(log10(0.01), log10(10), 101);

% Tell Optickle the mechanical transfer function of the suspension
dampRes = [0.01 + 1i, 0.01 - 1i];
mechTFobj = zpk([], -par.w * dampRes, 1 / par.mass);
opt = setMechTF(opt, 'Mirror', mechTFobj);
mechTF = squeeze(freqresp(mechTFobj, 2*pi*f));

%% Tickle
[fDC, sigDC, sigAC, mMech, noiseAC, noiseMech] = tickle(opt, [], f);

%% Look at the results
% The reaction matrix is something like
%
%  [ 1  0
%    0  H(s) ]  where H(s) is the mechanical transfer function
%
% i.e. AM drive causes AM with a flat transfer function, and mechanical
% drive causes motion per whatever transfer function we supplied.
%
% The modified transfer function will couple AM to mirror position,
% so mMech should look like:
%
%  [1    0
%   G(s) 1]
%
% Multiplying these together gives us
%
%  [1 0  [1 0      [1 0
%   G 1]  0 H]  =   G H] 
%
% So all we need to do is read-out mMech(2,1,:) to get the transfer
% function of (AM drive) to (Mirror position).  We can further use our
% probe after the modulator to get the transfer function from (AM drive) to
% (Watts modulation) (which will be flat).  Divide that by (2/c) to get (AM
% drive) to (Newtons of radiation pressure).  Then divide those two
% transfer functions to get the mechanical transfer function of the
% suspension.  Should agree with what we set up in the first place.

% get (Radiation pressure Newtons) / (AM drive)
tf1 = (2/opt.c) * ...
    getTF(sigAC, getProbeNum(opt, 'Mirror DC'), getDriveNum(opt, 'AM'));

% get (Mirror meters) / (AM drive)
tf2 = getTF(mMech, getDriveNum(opt, 'Mirror'), getDriveNum(opt, 'AM'));

% get (Mirror meters) / (Radiation pressure Newtons)
mechTF_inferred = (tf2./tf1);

semilogx(f, db(mechTF_inferred), 'r--', 'linewidth', 2.5);
hold on
plot(f, db(mechTF), 'k', 'linewidth', 2.5);
hold off
legend('Inferred mech TF', 'Supplied mech TF', 'Location', 'Best');
xlabel('frequency [Hz]');
ylabel(' dB(meters per Newton) ');
title('radiation pressure check');
grid on
set([gca;findall(gca, 'Type','text')], 'FontSize', 16)

%% Plot the residual
plot(f, abs(mechTF_inferred ./ mechTF), 'o-', 'LineWidth', 2.5);
ylim([0 1]);
grid on;
title('residual');
legend('| inferred / supplied |');
ylabel('ratio');
xlabel('frequency');
set([gca;findall(gca, 'Type','text')], 'FontSize', 16)

%%
% There's a 2X disagreement