% demo1_iSTAC.m
%
% Simple script for testing out the iSTAC code

% 0. Set up parameters for an example Linear-Nonlinear-Poisson (LNP) neuron 
nt = 20;        % number of temporal elements of filter
tvec = (-nt+1:0)'; % vector of time indices (in units of stim frames)

% Make some temporal filters
filt1 = exp(-((tvec+4.5)/1.5).^2/2) -.2*exp(-((tvec+nt/2)/3).^2/2); %1st filter
filt1 = filt1./norm(filt1);  %normalize

filt2 = [diff(filt1); 0];  % 2nd filter
filt2 = filt2- filt1*(filt1'*filt2); %orthogonalize to 1st filter
filt2 = filt2./norm(filt2); % normalize

% Make plots
plot(tvec, [filt1 filt2])  
title('filters for simulation');
xlabel('time before spike'); ylabel('filter coeff');


%% 1.  1st example:  rectified linear LNP neuron
% ===================================================

% Create stimulus ------------------
slen = 2000;   % Stimulus length
Stim = randn(slen,1);
RefreshRate = 100; % refresh rate

linresp = sameconv(Stim, filt1);  % filter output
r = max(linresp,0)*50; % instantaneous spike rate (rectified linear nonlinearity)
spikes = poissrnd(r/RefreshRate); % generate spikes

[sta,stc,rawmu,rawcov] = simpleSTC(Stim,spikes, nt); % Compute STA and STC
[u,~,~] = svd(stc); % Compute eigenvectors of STC matrix

% Compute iSTAC estimator
nfilts = 1;
istacfilts = compiSTAC(sta, stc, rawmu, rawcov, nfilts);
plot(tvec, filt1, 'k--', tvec, sta./norm(sta), tvec, u(:,end), tvec, istacfilts);
legend('true k', 'STA', 'STC', 'iSTAC', 'location', 'northwest');
xlabel('time before spike (bins)'); ylabel('filter coeff');
title('filter estimates');

% Compute error (using subspace angle between true filter and estimates)
Errs = [subspace(filt1,sta) subspace(filt1, u(:,end)) subspace(filt1, istacfilts)];
fprintf(1, 'Errors: STA=%.3f, STC=%.3f, iSTAC=%.3f\n', Errs(1), Errs(2), Errs(3));


%% 2. 2nd example:  2-filter LNP model with quadratic nonlinearity, correlated stimuli
% ===================================================
%
% In this example, the STA also lies in the space spanned by the two
% filters, but by using iSTAC we can pull out the 2 relevant axes much more
% accurately

% % Create stimulus ------------------
slen = 10000;   % Stimulus length (Better convergence w/ longer stimulus)
Stim = randn(slen,1);
Stim = conv2(Stim,normpdf(-3:3,0,1)','same'); % smooth stimulus
RefreshRate = 100; % refresh rate

DC = [0.75 .5];  % Setting to zero means expected STA is zero
linresp = [sameconv(Stim, filt1)+DC(1) sameconv(Stim,filt2)+DC(2)];  % filter output
r = 10*linresp(:,1).^2 + 8*linresp(:,2).^2; % instantaneous spike rate
spikes = poissrnd(r/RefreshRate); % generate spikes

[sta,stc,rawmu,rawcov] = simpleSTC(Stim,spikes,nt);
[u,s,v] = svd(stc);

nfilts = 10;  % (Only need 2, but compute 10 for demonstration purposes)
eigvalthresh = 0.05; % eigenvalue cutoff threshold (for pruning dims from raw stimulus)
[istacfilts,vals] = compiSTAC(sta, stc, rawmu, rawcov, nfilts,eigvalthresh);
KLcontributed = [vals(1); diff(vals)];
nfilts = length(vals);

subplot(221);  plot(1:nfilts, KLcontributed, 'o');
title('KL contribution');
xlabel('subspace dimensionality');

subplot(221);
plot(tvec, filt1, 'k--', tvec, u(:,1:2)*u(:,1:2)'*filt1, ...
    tvec, istacfilts(:,1:2)*istacfilts(:,1:2)'*filt1, 'r');
title('Reconstruction of 1st filter'); ylabel('filter coeff');
legend('true k', 'STC', 'iSTAC', 'location', 'northwest');

subplot(223);
plot(tvec, filt2, 'k--', tvec, u(:,1:2)*u(:,1:2)'*filt2, ...
    tvec, istacfilts(:,1:2)*istacfilts(:,1:2)'*filt2, 'r');
title('Reconstruction of 2nd filter');
xlabel('time before spike'); ylabel('filter coeff');

subplot(222);  
plot(1:nfilts, KLcontributed, 'o');
title('KL contribution');
xlabel('subspace dimensionality');

subplot(224);
plot(tvec, istacfilts(:,1:2));
title('iSTAC filters');
legend('1st', '2nd', 'location', 'northwest');

Errs = [subspace([filt1 filt2], u(:,1:2)) subspace([filt1 filt2], istacfilts(:,1:2))];
fprintf(1, 'Errors: STC=%.3f, iSTAC=%.3f\n', Errs(1), Errs(2));


%% 3. Illustrate convergence behavior with 1D linear-rectified model
% ===================================================

% Run several experiments of different length
slens = [1000 4000 16000 64000]; % number of samples for different expts
nexpts = length(slens);  % number of experiments
niters = 20;  % number of times to repeat each experiment
nfilts = 1; % number of filters to estimate 

% Initialize matrices for storing errors
ErrsSTA = zeros(niters, nexpts);
ErrsiSTAC = zeros(niters, nexpts);

% Run simulated experiments
for i = 1:nexpts
    slen = slens(i);  % length of this block of experiments
    fprintf('\nStimulus length: %d\n------\n', slen);
    for j = 1:niters

        % Generate stimulus and LNP response
        Stim = randn(slen,1);
        linresp = sameconv(Stim, filt1);  % filter output
        r = max(linresp,0)*50; % instantaneous spike rate
        spikes = poissrnd(r/RefreshRate); % generate spikes

        [sta,stc,rawmu,rawcov] = simpleSTC(Stim,spikes, nt); % Compute STA and STC
        istacFilt =  compiSTAC(sta, stc, rawmu, rawcov, nfilts);  % compute iSTAC estimate

        ErrsSTA(j,i) = subspace(filt1,sta); 
        ErrsiSTAC(j,i)= subspace(filt1, istacFilt);
    end
end

% Make plot
clf;
plot(slens, mean(ErrsSTA), 'bo-', slens, mean(ErrsiSTAC), 'ro-');
set(gca, 'xscale', 'log', 'yscale', 'linear');
legend('STA', 'iSTAC', 'location', 'northeast');
xlabel('stim length (samples)');
ylabel('error (radians)');
title('comparison of STA and iSTAC');
