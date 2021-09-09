function [mest, psig, pest, mse, res, ln_like1, ln_like2, F_Bayes, F_Bayes2Ln] = geostatistical_inversion2(A,d,Q,R,Xb)
%% perform geostatistical inversion
% inputs:
%     A == design matrix
%     R == data covariance matrix
%     Q == prior model covariance matrix
%     d == data
%     Xb = X matrix for beta 

% outputs:
%     mest         == Mode of estimated model parameter values
%     psig         == diagonal elments of posterior model covariance
%     pest         == posterior estimate??
%     mse          == mean squared error (square of residuals)
    
% References

% Reinisch, E. C., M. Cardiff, and K. L. Feigl (2018), Characterizing
% Volumetric Strain at Brady Hot Springs, Nevada, USA Using Geodetic Data,
% Numerical Models, and Prior Information, Geophysical Journal
% International, 1501?1513. http://dx.doi.org/10.1093/gji/ggy347 

% Kitanidis, P. K. (2007), On stochastic inverse modeling, Geophysical
% Monograph-American Geophysical Union, 171, 19.

% Cardiff, M. (2020), Bayesian Approaches to Inverse Problems, edited.

%% verify sizes
[ndata,mparam] = size(A)

[nR,mR] = size(R)
if nR ~= ndata || mR ~= ndata
    nR
    mR
    ndata
    error(sprintf('Data covariance matrix R must be %d by %d\n',ndata,ndata));
end

[nQ,mQ] = size(Q)
if nQ ~= mparam || mQ ~= mparam
    nQ
    mQ
    mparam
    error(sprintf('Model covariance matrix Q must be %d by %d\n',mparam,mparam));
end

% multiply 
AR = A'*R;

fprintf(1,'Start defining Ggeo and Dgeo\n')
Ggeo = [(A'*A*Q)*A'+AR+Xb*(A*Xb)', A'*(A*Xb)+Xb];
%clear Q
disp('Ggeo defined')

% 2021/09/07 Kurt simplify 
Dgeo = A'*d + Xb;
%Dgeo = A'*d;
%clear b
disp('Dgeo defined')

%clearvars -except A AR ARfac Ggeo Dgeo include_AR V0 vfactor nf depthsz centx centy dx dy datx daty nrowsA ncolsA datx_vec daty_vec grdx grdy dz dt sig2 mu bi ncells lx ly irefcount zcol topz

% what is pest??
 % from regular ls with full dVdt3=unit(-23107,'m^3/year'), from 4Okada -3.1624e+04
 % from regular ls with < 500, DV =  -37107.3 m^3
% from Ali et al. (2016) result for In20130513_20140511: 2013.3616 2014.3562 -19244.912 3863.12 TSX T53_32785_38296 0.898067 0.906401 0.9946
%pest = pinv(Ggeo'*Ggeo)*Ggeo'*Dgeo; %doulbe version: DV = -376305 m^3 %pinv(Ggeo'*Ggeo,  1e-18)*Ggeo'*Dgeo; % DV =  -37236.9 m^3
%load('Q.mat')
disp('Start inversion')
%pest = pinv((Ggeo))*(Dgeo);
pest = pinv(Ggeo'*Ggeo)*Ggeo'*Dgeo; 

mpest=size(pest)

disp('inversion completed')


disp('finding posterior model covariance')

% store only diagonal elements
psig = diag(Q - [A*Q; Xb']'*pinv(Ggeo)*[A', Xb]*[A*Q; Xb']);

disp('defining beta and xi')
%clear Ggeo
%load('b.mat')
beta = pest(end)
xi = pest(1:end-1);
% beta = pest(end-1:end);
% xi = pest(1:end-2);

% mode 
disp('Mode of estimated parameter values')
mest = Xb*beta + Q*A'*xi;

% compare to true model parameters
%mse0 = (b-A*mest)'*(b-A*mest)
res = (d-A*mest);
%load('AR-50.mat')
[Arows, Acols] = size(A);
%degfree = Arows - Acols;
%mse_pre = ((res)'*(res))/((.005/dt)^2*ARfac)/degfree
%mse = (b-A*mest)'*pinv(AR)*A'*(b-A*mest)*ARfac
mse = (d-A*mest)'*pinv(AR)*A'*(d-A*mest)  % simplifies because ARfac = 1
%clear AR

%% find logs of likelihoods
%
% if bi == 2
%     llikelihood_dP = -1/2*mse
%     lprior_dP = -1/2*(mest - Xb*beta)'*pinv(Q)*(mest-Xb*beta)
%     lpost_dP = lprior_dP*llikelihood_dP
% end
% end

ln_like1 =  -0.5*(mest - Xb*beta)'*pinv(Q)*(mest-Xb*beta)

ln_like2  = -0.5*mse*ln_like1

% Bayes Factor
%  Reinisch et al. [2018] state
% "the Bayes Factor F_B which is the ratio of the likelihood of the model
% under the alternative hypothesis to the likelihood of the model under the
% null hypothesis (e.g. Kass & Raftery 1995; Gelman et al. 2013, p. 183).
% To interpret the mechanism driving the volume changes at depth, we set
% the null hypothesis H0 to be decreasing pore fluid pressure and the
% alternative hypothesis H1 to be thermal contraction. We find the value of
% FB = 211.7 leading to (2lnFB) = 10.7. According to the scale by Kass &
% Raftery (1995), a value of (2ln FB ) greater than 10 gives ?very strong
% evidence? against the null hypothesis."

% Kass, R., and A. Raftery (1995), Bayes Factors, J. Am. Stat. Assoc, 90,
% 773-795. https://www.jstor.org/stable/2291091 
% In a 1935 paper and in his
% book Theory of Probability, Jeffreys developed a methodology for
% quantifying the evidence in favor of a scientific theory. The centerpiece
% was a number, now called the Bayes factor, which is the posterior odds of
% the null hypothesis when the prior probability on the null is one-half.
% Although there has been much discussion of Bayesian hypothesis testing in
% the context of criticism of P-values, less attention has been given to
% the Bayes factor as a practical tool of applied statistics. In this
% article we review and discuss the uses of Bayes factors in the context of
% five scientific applications in genetics, sports, ecology, sociology, and
% psychology. We emphasize the following points:
% * From Jeffreys' Bayesian viewpoint, the purpose of hypothesis testing is
% to evaluate the evidence in favor of a scientific theory. * Bayes factors
% offer a way of evaluating evidence in favor of a null hypothesis.
% * Bayes factors provide a way of incorporating external information into
% the evaluation of evidence about a hypothesis.
% * Bayes factors are very general and do not require alternative models to
% be nested.
% * Several techniques are available for computing Bayes factors, including
% asymptotic approximations that are easy to compute using the output from
% standard packages that maximize likelihoods.
% * In "nonstandard" statistical models that do not satisfy common
% regularity conditions, it can be technically simpler to calculate Bayes
% factors than to derive non-Bayesian significance tests.
% * The Schwarz criterion (or BIC) gives a rough approximation to the
% logarithm of the Bayes factor, which is easy to use and does not require
% evaluation of prior distributions.
% * When one is interested in estimation or prediction, Bayes factors may
% be converted to weights to be attached to various models so that a
% composite estimate or prediction may be obtained that takes account of
% structural or model uncertainty.
% * Algorithms have been proposed that allow model uncertainty to be taken
% into account when the class of models initially considered is very large.
% * Bayes factors are useful for guiding an evolutionary model-building
% process.
% * It is important, and feasible, to assess the sensitivity of conclusions
% to the prior distributions used.

F_Bayes = ln_like2 - ln_like1

% (2ln(FB)) This value is used for assessment
F_Bayes2Ln = 2.0 * log(F_Bayes)

return
end



