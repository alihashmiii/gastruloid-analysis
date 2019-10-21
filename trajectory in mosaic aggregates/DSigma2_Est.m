function [ varargout ] = DSigma2_Est( xMatrix , R , dT , varargin )
%DSigma2_Est
%
%Estimate diffusion coefficient (D) and/or variance of
%localization noise (Sigma2) from single-particle tracking
%data. 
%
%If you use this code for your work, please reference the appropriate
%literature [1-3]. Details of the calculations for maximum likelihood
%estimation are found in Ref. [1] and [3]. Details of the Optimal Least
%Squares Fit (OLSF) are found in Refs. [2] and [3].
%
%Andrew Berglund, NIST Center for Nanoscale Science and Technology
%Created: March 17, 2010
%Latest version: July 6, 2011
%Please send comments, questions, and problems to ajberglund@gmail.com
%
%SYNTAX
%   [EstVals1 , EstVals2 , ...] = 
%       DSigma2_Est(xMatrix , R , dT , Est1 , FixedVals1 , Est2, FixedVals2 ,...)
%
%       xMatrix =   Matrix of d-dimensional particle tracking data.  
%                   size(xMatrix) = [N d] where N is the number of data
%                   points and d is the dimension.
%       R       =   Motion blur coefficient R between 0 and 1/4. 
%                   See reference [1] for details
%       dT      =   Frame interval
%
%       Choice of estimator is specificied by the strings Est1,Est2,... and
%       corresponding fixed parameter values FixedVals1, FixedVals2,...
%
%       Est1 =   Estimator name chosen from the following
%                   'DSigma2_OLSF'  =   OLSF Estimate of D and Sigma2
%                   'DSigma2_MLE'   =   MLE of D and Sigma2
%                   'D_MLE'         =   MLE of D (Sigma2 fixed)
%                   'D_Sigma2'      =   MLE of Sigma2 (D fixed)
%
%       FixedVals provides fixed values of D or Sigma2 for some estimators:
%
%                     FixedVals    = []   
%                       for DSigma2_* estimators this value is ignored
%                     FixedVals    = [Sigma2]     
%                       for D_* estimators with Sigma2 fixed
%                     FixedVals    = [D]    
%                       for Sigma2_* estimators with D fixed
%
%       EstVals is a Matlab STRUCTURE containing the estimated (and fixed)
%       values:
%                   EstVals.D       = Estimated D value
%                   EstVals.Sigma2  = Estimated Sigma^2 value
%                   EstVals.v       = Estimated drift velocity *vector*
%                   EstVals.flag    = Error flag set to 1 if estimator
%                                       attained convergence criteria. 
%EXAMPLES:
%   EstVals = DSigma2_Est(xM,R,dT,'DSigma2_MLE',[])
%       returns estimates of D and Sigma^2 in the structure EstVals
%
%   [evMLE,evOLSF] = DSigma2_Est(xM,R,dT,'DSigma2_MLE',[],'DSigma2_OLSF',[])
%       returns estimates of D and Sigma^2 using both MLE and OLSF
%       estimators
%
%   EstVals = DSigma2_Est(xM,R,dT,'D_MLE',Sigma2)
%       returns MLE estimate of D with Sigma2 fixed at specified value
%
%HISTORY:
%   Last Update: 7/27/2011
%
%REFERENCES:
%   [1] "Statistics of Single-Particle Tracking," A.J. Berglund 
%   Physical Review E, 82, 011917 (2010).
%   http://pre.aps.org/abstract/PRE/v82/i1/e011917
%
%   [2] "Mean square displacement analysis of single-particle trajectories
%   with localization error: Brownian motion in an isotropic medium," X.
%   Michalet, Physical Review E, 82, 041914 (2010).
%   http://pre.aps.org/abstract/PRE/v82/i4/e041914
%
%   [3] "Optimal Diffusion Coefficient Estimation in Single-Particle
%   Tracking," X. Michalet and A.J. Berglund , Physical Review E (2012).
%
%DISCLAIMERS: This software was developed at the National Institute of
%Standards and Technology by employees of the Federal Government in the
%course of their official duties. Pursuant to title 17 Section 105 of the
%United States Code this software is not subject to copyright protection
%and is in the public domain. DSigma2_Est is an experimental system. NIST
%assumes no responsibility whatsoever for its use by other parties, and
%makes no guarantees, expressed or implied, about its quality, reliability,
%or any other characteristic. We would appreciate acknowledgement if the
%software is used. Certain commercial software is identified in this system
%in order to specify the procedure adequately. Such identification is not
%intended to imply recommendation or endorsement by the National Institute
%of Standards and Technology, nor is it intended to imply that the
%materials or equipment identified are necessarily the best available for
%the purpose.

    d = size(xMatrix,2); % dimension
    numEstimators = numel(varargin)/2; % I'm not checking that valid number of arguments is passed
    EstVals = cell(numEstimators,1);
    
    xInt = xMatrix;%(:,kk); % X is the kk_th trajectory

    for pp = 1:numEstimators;
        thisEst = varargin{2*pp-1};
        InitVals = varargin{2*pp};
        switch thisEst
            case 'DSigma2_MLE'
                InitVals = [];
                if isempty(InitVals)
                    D_i = mean(var(diff(xInt,1))/2/dT/2);
                    Sigma2_i = mean(var(diff(xInt,1))/2);
                else
                    D_i = InitVals(1);
                    Sigma2_i = InitVals(2);
                end

                dX=diff(xInt,1)-repmat(mean(diff(xInt,1)),size(xInt,1)-1,1);
                
                options = optimset('MaxFunEvals',100000,'MaxIter',100000);
                [b,~,errflag] = fminsearch(@(b)-Likelihood_subfunction(...
                    dst(dX).^2,b(1),b(2),dT,R ), [D_i Sigma2_i],options);

                EstVals{pp}.D = b(1);
                EstVals{pp}.Sigma2 = b(2);
                EstVals{pp}.v = mean(diff(xInt))/dT;
                EstVals{pp}.flag = errflag;
            case 'D_MLE'
                if isempty(InitVals)
                    flag = 0;
                    error('D_MLE estimator option requires a non-empty InitVals parameters.');
                else 
                    D_i = mean(var(diff(xInt,1))/2/dT/2);
                    Sigma2 = InitVals;
                end

                dX=diff(xInt,1)-repmat(mean(diff(xInt,1)),size(xInt,1)-1,1);

                options = optimset('MaxFunEvals',100000,'MaxIter',100000);
                [b,~,errflag] = fminsearch(@(D)-Likelihood_subfunction(dst(dX).^2,D,Sigma2,dT,R) , [D_i],options);

                EstVals{pp}.D = b(1);
                EstVals{pp}.Sigma2 = Sigma2;
                EstVals{pp}.v = mean(diff(xInt))/dT;
                EstVals{pp}.flag = errflag;
            case 'Sigma2_MLE'
                if isempty(InitVals)
                    flag = 0;
                    error('Sigma2_MLE estimator option requires a non-empty InitVals parameters.');
                else
                    Sigma2_i = mean(var(diff(xInt,1))/2);
                    D = InitVals;
                end

                dX=diff(xInt,1)-repmat(mean(diff(xInt,1)),size(xInt,1)-1,1);

                options = optimset('MaxFunEvals',100000,'MaxIter',100000);
                [Sigma2,~,errflag] = fminsearch(@(b)-Likelihood_subfunction(dst(dX).^2,D,b,dT,R) , [Sigma2_i],options);

                EstVals{pp}.D = D;
                EstVals{pp}.Sigma2 = Sigma2;
                EstVals{pp}.v = mean(diff(xInt))/dT;
                EstVals{pp}.flag = errflag;
            case 'DSigma2_OLSF' 
                % Xavier Michalet's algorithm, see Ref. [2]

                NXM = numel(xInt); % XM uses N = number of trajectory points
                % Calculate MSD at each iteration, since in
                % general pmin will be << N

                donea = 0;
                doneb = 0;

                pa = floor(NXM/10);
                pb = floor(NXM/10);

                iter = 1;
                maxiter = 100;
                flag = 1;

%                     [pam,pbm] = PMin_XM(Inf,NXM);
%                     rho = MSDoverlap(xInt);%,max(pam,pbm));
%                     rho = rho(:);
                rho = [];
                while ~donea && ~doneb
                    if iter == maxiter;
                        donea = 1;
                        doneb = 1;

                        D=NaN;
                        sig2=NaN;

                        flag=0;
                        disp(['OLSF did not converge in ' num2str(maxiter) 'iterations']);
                    else
%                             rho = MSDoverlap(xInt , max(pa(end),pb(end)));
                        if max(pa(end),pb(end)) > numel(rho)
                            % calculate more rho points if necessary
                            rho = MSDoverlap(xInt , max(pa(end),pb(end)));
%                             rho = rho(:);
                        end
                        rho = rho(:);

                        % one line least squares fit to MSD = a+b*rho ;
                        if ~donea
                            A = [ones(pa(end),1) (1:pa(end)).' ];
                            B = rho(1:pa(end));
                            X = A\B;
                            aa = X(1); % a estimate for pa
                            ba = X(2); % b estimate for pa

                            if aa<0
                                xa = 0;
                            elseif ba < 0
                                xa = Inf;
                            else
                                xa = aa/ba;
                            end

                            newpa = PMin_XM(xa , NXM);
                            if ismember(newpa,pa);
                                Da = ba/2/dT/d; %this is not the D estimate, but used below
                                sig2 = (aa+4*Da*R*dT)/2/d;%***check the /d's here***
                                donea = 1;
                            end

                            pa = [pa(:); newpa;];                            
                        end

                        if ~doneb
                            A = [ones(pb(end),1) (1:pb(end)).' ];
                            B = rho(1:pb(end));
                            X = A\B;
                            ab = X(1); % a estimate for pb
                            bb = X(2); % b estimate for pb

%                             pb = [pb(:); PMin_XM(ab/bb , NXM);];

                            if ab<0
                                xb = 0;
                            elseif bb < 0
                                xb = Inf;
                            else
                                xb = ab/bb;
                            end

                            newpb = PMin_XM(xb , NXM);
                            if ismember(newpb,pb);
                                D = bb/2/dT/d;
%                                 sig2 = (ab+4*D*R*dT)/2;
                                doneb = 1;
                            end

                            pb = [pb(:); newpb;];   
                        end
                        iter = iter+1;
                    end
                end                                        
                EstVals{pp}.D = D;
                EstVals{pp}.Sigma2 = sig2; % need to return sigma for consistency
                EstVals{pp}.v = mean(diff(xInt,1))/dT;
                EstVals{pp}.flag=flag;

            otherwise
                EstVals{pp}.flag=-1;
                warning(['Unrecognized estimator "' thisEst '"']);
        end
    end;

    for pp=1:numEstimators;
        varargout{pp} = EstVals{pp};
    end;
end


function L = Likelihood_subfunction1D( dXX , D, sig2, dT, R )
% dXX is square dst of diff of trajectory with mean subtracted and
% k_th element multiplied by (-1)^k 

    if D<0 || sig2 < 0;
        L = -Inf;
    else
        N = numel(dXX);

        alpha = 2*D*dT*(1-2*R) + 2*sig2;
        beta = 2*R*D*dT - sig2;

        eigvec = alpha + 2*beta*cos(pi*(1:N)./(N+1));
        eigvec = reshape(eigvec, size(dXX));

        L = -.5*sum( log(eigvec) + (2/(N+1)).*dXX./eigvec );
    end;
end

function L = Likelihood_subfunction( dXX , varargin )
    L = 0;
    for kk=1:size(dXX,2);
        L = L+Likelihood_subfunction1D(dXX(:,kk),varargin{:});
    end
end

function [rho] = MSDoverlap(X , n)
% return mean-square displacement vector rho for positions X, out to sample
% interval n. size(X) = [N d] where d is the dimensionality
    [N,~] = size(X);
    if nargin<2;
        n = N-1;
    end;

    rho = zeros(n,1);
    for kk=1:n;
        ix =1;
        num = 0;
        while (ix+kk)<=N;
            rho(kk) = rho(kk) + sum((X(ix+kk,:)-X(ix,:)).^2,2);
            ix = ix+1;
            num = num+1;
        end;
        rho(kk) = rho(kk)/num;
    end;
end

function [pa , pb] = PMin_XM(x,N)
    % Calculate moptimal fit point from XM formulas.
    % These formulas use XM conventions for N. These differ from AB.

    fa = 2+1.6*x.^0.51;
    La = 3+(4.5.*N.^.4-8.5).^1.2;

    fb = 2+1.35*x.^0.6;
    Lb = 0.8 + 0.564*N;

    if isinf(x);
        pa = floor(La);
        pb = floor(Lb);
    else
        pa = floor(fa.*La./(fa.^3+La.^3).^.33);
        pb = min(floor(Lb), floor(fb.*Lb./(fb.^3+Lb.^3).^.33));
    end

    % make sure nothing is zero
    pa = max(2,pa);
    pb = max(2,pb);
end