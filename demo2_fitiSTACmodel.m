% demo2_fitiSTACmodel.m
%
% Second demo script that illustrates fitting iSTAC model (including
% filters and exponentiated-quadratic) to data from a simulated neuron with
% three stimulus filters, stimulated with Gaussian white noise.

%% 1a. Make filters for simulated LNP neuron

nkt = 30;  % number of time bins in filter

dtBin = .01; % bin size for representing time (here .01s => 100 Hz stimulus frame rate)
RefreshRate = 1/dtBin; % Refresh rate
tk = (-nkt+1:0)'; % vector of time indices (in units of stim frames)

% Make some stimulus filters for simulated neuron
filt1 = exp(-((tk+nkt/5)/(nkt/12)).^2)-.25*exp(-((tk+nkt/2)/(nkt/5)).^2); % biphasic filter
filt1 = filt1./norm(filt1);  %normalize to make a unit vector

filt2 = [filt1(2:end); 0];  % make 2nd filter (shift filter 1 by 2)
filt2 = filt2-filt1*(filt1'*filt2); % orthogonalize w.r.t. filter 1
filt2 = filt2./norm(filt2); 

filt3 = [filt2(2:end); 0];  % make 3rd filter (shift filter 2 by 1)
filt3 = filt3-[filt1 filt2]*([filt1 filt2]'*filt3); % orthogonalize w.r.t. filter 1 & 2
filt3 = filt3./norm(filt3); % normalize to make a unit vector

filts_true = [filt1 filt2 filt3];

% Plot these filters
plot(tk*dtBin, filts_true, 'o-');
title('filters for simulation');
xlabel('time before spike (s)'); ylabel('filter coefficient');
axis tight;


%% 1b.  Simulate spike responses from LNP neuron

% Create stimulus 
slen = 10000;   % Stimulus length (make longer to see better recovery of true filters)
Stim = randn(slen,1); % 1D Gaussian white noise (temporal) stimulus
%Stim = conv2(Stim,normpdf(-3:3,0,1)','same'); % smooth stimulus

% Convolve stimulus with filters
f1 = sameconv(Stim,filt1);
f2 = sameconv(Stim,filt2);
f3 = sameconv(Stim,filt3);

% Compute output of nonlinearity
softrect = @(x)(log(1+exp(x))); % soft-rectification function
fnlin = @(x1,x2,x3)(softrect(40*x1 + 12*x2.^2 + 10*(x3).^2));
lam = fnlin(f1,f2,f3);

%  Simulate spike train
sps = poissrnd(lam/RefreshRate); % generate spikes

% Report number of spikes and spike rate
nsp = sum(sps); % number of spikes
fprintf('LNP simulation: %d time bins, %d spikes  (%.2f sp/s)\n', slen,nsp, nsp/slen*RefreshRate);

%% 2. Fit iSTAC filters

nktfilt = 30; % number of time bins to use for filter 
% This is important: normally would want to vary this to find optimal filter length

nfilts = 3; % number of filters to estimate

% Compute STA and STC
[sta,stc,rawmu,rawcov] = simpleSTC(Stim,sps,nktfilt);  % compute STA and STC

% Compute iSTAC filter estimates
fprintf('\nComputing iSTAC estimate\n');
nFilts = 3; % number of filters to compute
[istacFilts,vals,DD] = compiSTAC(sta(:),stc,rawmu,rawcov,nFilts); % find iSTAC filters

% Plot filters against true filter
tt = -nkt+1:0;  ttf = -nktfilt+1:0;
aa = sign(diag((filts_true'*istacFilts))); % sign of dot products (to flip filter estimates to visually match, if necessary)
subplot(231); plot(tt, filts_true(:,1), 'k--', ttf, aa(1)*istacFilts(:,1), 'linewidth', 2);
title('filter 1'); xlabel('time before spike (bins)');
subplot(232); plot(tt, filts_true(:,2), 'k--', ttf, aa(2)*istacFilts(:,2), 'linewidth', 2);
title('filter 2'); 
subplot(233); plot(tt, filts_true(:,3), 'k--', ttf, aa(3)*istacFilts(:,3), 'linewidth', 2);
title('filter 3'); 
%% 3. Fit "exponentiated quadratic" nonlinearity for two 1-filter models

% Fit and plot nonlinearity along each filter
pp1 = fitNlin_expquad_ML(Stim,sps,istacFilts(:,1),RefreshRate); % fit nonlinearity
[x1,nl1] = compNlin_1D(istacFilts(:,1),pp1.nlfun,Stim); % gridded nonlinearity for plotting
pp2 = fitNlin_expquad_ML(Stim,sps,istacFilts(:,2),RefreshRate); % fit nonlinearity
[x2,nl2] = compNlin_1D(istacFilts(:,1),pp2.nlfun,Stim); % gridded nonlinearity for plotting
pp3 = fitNlin_expquad_ML(Stim,sps,istacFilts(:,3),RefreshRate); % fit nonlinearity
[x3,nl3] = compNlin_1D(istacFilts(:,1),pp3.nlfun,Stim); % gridded nonlinearity for plotting

subplot(234); 
plot(x1,nl1,x2,nl2,x3,nl3,'linewidth',2);
xlabel('normalized filter output');
ylabel('firing rate (Hz)'); 
legend('filter 1', 'filter 2', 'filter3', 'location', 'northwest');
title('1D nonlinearites iSTAC filters');

%% 4. Fit nonlinearity for two 2-filter models 

% Fit iSTAC model nonlinearity using filters 1 and 2
pp12 = fitNlin_expquad_ML(Stim,sps,istacFilts(:,1:2),RefreshRate); % LNP model struct
[xgrd,ygrd,nlvals] = compNlin_2D(pp12.k,pp12.nlfun,Stim); % compute gridded nonlinearity

% Plot gridded 2D nonlinearity
subplot(235); mesh(xgrd,ygrd,nlvals); 
zlm = [0 max(nlvals(:))*1.01];
axis tight; set(gca,'zlim',zlm);
xlabel('filter 1 axis');ylabel('filter 2 axis'); 
zlabel('spike rate (sps/sec)'); title('2D iSTAC nonlinearity (filters 1 & 2)')

% Fit iSTAC model nonlinearity using filters 2 and 3
pp23 = fitNlin_expquad_ML(Stim,sps,istacFilts(:,2:3),RefreshRate); % LNP model struct
[xgrd,ygrd,nlvals] = compNlin_2D(pp23.k,pp23.nlfun,Stim);

% Plot gridded 2D nonlinearity
subplot(236); mesh(xgrd,ygrd,nlvals); 
axis tight; set(gca,'zlim',zlm);
xlabel('filter 2 axis');ylabel('filter 3 axis'); 
zlabel('spike rate (sps/sec)'); title('2D iSTAC nonlinearity (filters 2 & 3)')

%% 5. Fit exp-quadratic nonlinearity with all 3 filters

pp123 = fitNlin_expquad_ML(Stim,sps,istacFilts,RefreshRate); % LNP model struct


%% 6. Now compare firing rates and log-likelihood of different models

% Compute istac filter outputs
ss1 = sameconv(Stim,istacFilts(:,1));  % filter 1 output
ss2 = sameconv(Stim,istacFilts(:,2));  % filter 2 output
ss3 = sameconv(Stim,istacFilts(:,3));  % filter 3 output

% Compute firing rates
rr1 = pp1.nlfun(ss1); % model with filter 1 only
rr2 = pp1.nlfun(ss2); % model with filter 2 only
rr12 = pp12.nlfun([ss1 ss2]); % model with filters 1 and 2
rr23 = pp23.nlfun([ss2 ss3]); % model with filters 2 and 3
rr123 = pp123.nlfun([ss1 ss2 ss3]); % model with all 3 filters

% Plot rates
tt = 1:100;
subplot(211); stem(tt,sps(tt),'k'); title('spike counts per bin');
subplot(212); 
plot(tt,rr1(tt),tt,rr2(tt),tt,rr12(tt),tt,rr23(tt),tt,rr123(tt),'linewidth', 2); 
xlabel('time bin'); 
ylabel('spike rate (sps/s)'); 
title('rates of fitted models');
legend('istac-1', 'istac-2', 'istac-12', 'istac-23', 'istac-123');

% Compute log-likelihoods
flogli = @(rr)(log(rr)'*sps - sum(rr)/RefreshRate);

loglis = [flogli(rr1), flogli(rr2), flogli(rr12), flogli(rr23), flogli(rr123)];
fprintf('\nLog-likelihood of fitted iSTAC models\n-----------------\n');
fprintf(' istac1   = %.0f\n istac2   = %.0f\n istac12  = %.0f\n istac23  = %.0f\n istac123 = %.0f\n',loglis);

% Note that these are TRAINING log-likelihoods. To select a model, one
% should divide data into "training" and "test" sets, fit the iSTAC filters
% and nonlinearities on the training data, and then compute log-likelihood
% on the held out test data.
