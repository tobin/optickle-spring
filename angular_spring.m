%% Angular spring

%% Parameters
% Here we specify some basic parameters of the optical setup, such as the
% reflectivities of the optics.  These will be used later, when the
% Optickle model is created.

vFrf     = 0;           % carrier only
P_laser  = 3500;          % Laser power [Watts]

par.IX.T = 0.02930;     % Power transmission coefficient
par.EX.T = 10e-6;

par.w = 2 * pi * 0.74;  % Suspension resonant frequency [rad/s]
par.mass = 10;          % Test mass mass [kg]


par.w_pit = 2 * pi * 0.6;   % pitch mode resonance frequency

% Mirror dimensions
par.rTM = 0.25/2;           % test-mass radius
par.tTM = 0.1;              % test-mass thickness
par.iTM = (3 * par.rTM^2 + par.tTM^2) / 12;  % TM moment / mass

par.iI = par.mass * par.iTM;  % moment of mirrors

par.IX.ROC = 14600;
par.EX.ROC = 7400;
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
opt = addMirror(opt, 'ITM', 0, 1/par.IX.ROC, par.IX.T, 0, 0, 0);
opt = addMirror(opt, 'ETM', 0, 1/par.EX.ROC, par.EX.T, 0, 0, 0);

% Link them together
opt = addLink(opt, 'Laser', 'out', 'ITM', 'bk', 0);
opt = addLink(opt, 'ITM', 'fr', 'ETM', 'fr', 3995);
opt = addLink(opt, 'ETM', 'fr', 'ITM', 'fr', 3995);

% We don't actually need any probes for our purposes here, but Optickle
% complains if we don't add any.  (Probes are like photodiodes, except that
% they can observe the fields without disturbing them.)
opt = addProbeOut(opt, 'CAVITY DC',  'ITM', 'fr', 0, 0);



% Tell Optickle the mechanical transfer functions of the suspensions
dampRes = [0.01 + 1i, 0.01 - 1i];
mechTFobj = zpk([], -par.w_pit * dampRes, 1 / par.iI);

opt = setMechTF(opt, 'ITM', mechTFobj, 2);
opt = setMechTF(opt, 'ETM', mechTFobj, 2);

opt = setCavityBasis(opt, 'ITM', 'ETM');
%%

f = logspace(log10(0.1), log10(10), 301);

mechTF = squeeze(freqresp(mechTFobj, 2*pi*f));

[fDC, sigDC]   = tickle(opt, [], f);
[sigAC, mMech] = tickle01(opt, [], f);

%% Plot the transfer functions in the mirror basis
clf
for ii=1:2
    for jj=1:2
        H = getTF(mMech, ii, jj) .* mechTF;
        subplot(2,1,1);
        loglog(f, abs(H), '.-');
        hold all
        subplot(2,1,2);
        semilogx(f, unwrap(angle(H))*180/pi, '.-');
        hold all
    end
end

%%

% Change to the common/differential basis.   Note that the coordinate
% for the ETM is opposite that of the ITM, because they are facing
% opposite directions.  So ETM+ITM gives the differential mode.

%S = [1 -1; 1 1];
[S, ~] = eig(mMech(:,:,1));
S = inv(real(S));

mMechDOF = zeros(size(mMech));
jj = [1 2];
for kk=1:size(mMech,3)
    mMechDOF(jj,jj,kk) = S * mMech(jj,jj,kk) / S;
end

clf
drive_names = {'STABLE', 'UNSTABLE'};
for ii=1:2
    for jj=1:2
        H = getTF(mMechDOF, ii, jj) .* mechTF;
        
        [B,A] = invfreqs(H, 2*pi*f, 1, 2);
        [Z,P,K] = tf2zpk(B,A);
        fprintf('%s --> %s:\n', drive_names{jj}, drive_names{ii});
        zpk(Z,P,K)
        subplot(2,1,1);
        loglog(f, abs(H), '.-');
        hold all
        subplot(2,1,2);
        semilogx(f, unwrap(angle(H))*180/pi, '.-');
        hold all
    end
end
