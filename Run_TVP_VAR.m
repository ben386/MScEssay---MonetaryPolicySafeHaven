%% This File runs the Time-Varying Parameter VARs %%

% When using this code or parts of it, please cite Paul, Pascal. 2019. 
% "The Time-Varying Effect of Monetary Policy on Asset Prices". 
% Review of Economics & Statistics.

% The code largely builds on the one provided on the website of Gary Koop,
% but is adapted to include the possibility of an exogenous variable.

clear all;
clc;

%% Choose Functions to Run

choose_set_figures      = 1; %1;
colored_figures         = 1;
        % 1: TVP-VARX in the main text, produces figures 3,4,5,14,16-19
        % 2: TVP-VAR Cholesky identification, produces figure 15
        % 3: TVP-VARX sensitivity "Priors", produces figures 20-22
        % 4: TVP-VARX sensitivity "ELB Episode", produces figures 23-25
        % 5: TVP-VARX sensitivity "Timing of Actions", produces figures 26-28

figure_settings

% VAR Settings
p               = 2;     % Number of lags
IRF_T           = 20;%40;    % Impulse response horizon
T               = 20;%40;   % Computed IRF Horizon

% Settings Gibbs Sampler
nrep            = 3000; % Number of replications
nburn           = 2000;  % Number of burn-in-draws

%% Load Data

%Set working directory
cd 'C:\Users\osull\Dropbox\Economics\Postgrad\Dissertation\MATLAB\Paul (2019) Code\Replication Files'

%cd '..\Data\Macro Time Series'
addpath([cd '\Data\Macro Time Series'])

% Load Macro series (1960 M1 - 2015 M10)
macro           = xlsread('macrodata.xlsx','Sheet 1','A2:E575');

% Convert stock prices, dividends, and house prices to real units
%macro(:,2)      = 100*macro(:,2)./macro(:,5);
%macro(:,3)      = 100*macro(:,3)./macro(:,5);
%macro(:,4)      = 100*macro(:,4)./macro(:,5);

%cd '..\Monetary Policy Surprises'
addpath([cd '\Data\Monetary Policy Surprises'])

% Load Monetary Policy Surprises
surprises   = xlsread('Monetary_Surprises.xls','Scheduled_Meetings','C2:D405');

%cd '..\..\Time-Varying Parameter VAR'
addpath([cd '\Time-Varying Parameter VAR'])

%% Define VAR

if cholesky == 0
    ydata           = [macro(:,1), macro(:,2), macro(:,3), macro(:,4), macro(:,5)];
    labels          = {'Federal Funds Rate','Gold Price','Yield Spread','Industrial Production','Consumer Price Index'};
    %labels          = {'Federal Funds Rate','Real Interest Rate','Stock Prices','Dividends','House Prices','Consumer Price Index','Industrial Production','Fundamental Stock Prices','Stock Prices - Fundamental Stock Prices'};
    differencing    = [0, 1, 0, 1, 1];
else
    ydata           = [macro(:,1), macro(:,2), macro(:,3), macro(:,4), macro(:,5)];
    labels          = {'Federal Funds Rate','Gold Price','Yield Spread','Industrial Production','Consumer Price Index'};
    differencing    = [0, 1, 0, 1, 1];
    %ydata           = [macro(:,3), macro(:,5), macro(:,6), macro(:,1), macro(:,4), macro(:,2)];
    %labels          = {'Dividends','Real Interest Rate','Consumer Price Index','Industrial Production','Federal Funds Rate','House Prices','Stock Prices','Fundamental Stock Prices','Stock Prices - Fundamental Stock Prices'};
    %differencing    = [1, 1, 1, 0, 1, 1];
end

%% Pre-processing

% Take Logs and First-Difference Data
for jj=1:size(ydata,2)
    if differencing(1,jj) == 1 
        ydata(:,jj)     = 100*log(ydata(:,jj));
        ydata(:,jj)     = [0;diff(ydata(:,jj))];
    end
end
ydata           = ydata(2:end,:);

% Define Yearlab
yearlab         = zeros(size(ydata,1),1);
months          = [1/12;2/12;3/12;4/12;5/12;0.5;7/12;8/12;9/12;10/12;11/12;1];
yearlab(1,1)    = 1960;
count           = 0;
k               = 0;
for j=1:size(yearlab,1)
    count = 1 + count;
    yearlab(j,1) = 1960 + k + months(count,1);
    if count == 12
        count = 0;
        k     = k + 1;
    end
end

% Restrict Sample
start_sample    = 12*(start_year-1970)+start_month-1;
end_sample      = size(ydata,1)-12*(2017-end_year)-(9-end_month);
Y               = ydata(start_sample:end_sample,:);
yearlab         = yearlab(start_sample:end_sample,:);

% Pick Instrument and Restrict Sample
end_instr       = size(surprises,1)-12*(2017-end_year)-(9-end_month);
instr           = surprises(1:end_instr,instr_use);
instr           = [zeros(size(Y,1)-size(instr,1),1);instr];

% Orthogonalizing Instrument against Lags of Y
X_clean         = lagmatrix(ydata(start_sample-p:end_sample,:),1:p);
X_clean         = X_clean(p+1:end,:);
beta_instr      = X_clean\instr;
instr           = instr-X_clean*beta_instr;

% Number of observations and dimension of X and Y
t               = size(Y,1); % t is the time-series observations of Y
M               = size(Y,2); % M is the dimensionality of Y

 % Size of the training sample
tau             = 12*(start_year_TVP-start_year)-(start_month-start_month_TVP)-p;
%tau             = 12*(1990-start_year+1)-(start_month-start_month_TVP)-p;

% ===================================| VAR EQUATION |==============================
% Generate lagged Y matrix. This will be part of the X matrix
ylag = mlag2(Y,p); % Y is [T x M]. ylag is [T x (Mp)]
% Form RHS matrix X_t = [1 y_t-1 y_t-2 ... y_t-k] for t=1:T
ylag = ylag(p+tau+1:t,:);

if cholesky == 1
    K = M + p*(M^2); % K is the number of elements in the state vector
else
    K         = M + p*(M^2) + M;
    contemp   = zeros(M,M);
    for mm = 1:M
        contemp(mm,mm)  = 1;
    end
end

% Create Z_t matrix: This creates the regressor matrix in the right format
% and includes constant and instrument
Z = zeros((t-tau-p)*M,K);

for i = 1:t-tau-p
    if cholesky == 0
        ztemp = [eye(M) contemp*instr(p+tau+i,1)];
    else
        ztemp = eye(M);
    end
    for j = 1:p        
        xtemp = ylag(i,(j-1)*M+1:j*M);
        xtemp = kron(eye(M),xtemp);
        ztemp = [ztemp xtemp];
    end
    Z((i-1)*M+1:i*M,:) = ztemp;
end

% Redefine FAVAR variables y
y       = Y(tau+p+1:t,:)';
yearlab = yearlab(tau+p+1:t);
% Time series observations
t       = size(y,2);   % t is now length_sample - p - tau

%% Preliminary TVP-VAR Settings

it_print    = 100;  % Print in the screen every "it_print"-th iteration

%========= PRIORS:
% To set up training sample prior as Primiceri, use the following subroutine
if cholesky == 1
    [B_OLS,VB_OLS,A_OLS,sigma_OLS,VA_OLS] = ts_prior(Y,tau,M,p);
else
    [B_OLS,VB_OLS,A_OLS,sigma_OLS,VA_OLS,hbar_OLS] = ts_prior_X(Y,tau,M,p,instr,contemp);
end

% Or use uninformative values
% B_OLS     = zeros(K,1);
% VB_OLS    = eye(K);

%OLS estimation
% Mdl = varm(5,2);
% EstMdl = estimate(Mdl, ydata(1:132,:))

%-------- Now set prior means and variances (_prmean / _prvar)
% This is the Kalman filter initial condition for the time-varying
% parameters B(t)
% B_0 ~ N(B_OLS, 4Var(B_OLS))
B_0_prmean = B_OLS;
B_0_prvar  = 4*VB_OLS;

% Note that for IW distribution I keep the _prmean/_prvar notation...
% Q is the covariance of B(t)
% Q ~ IW(k2_Q*size(subsample)*Var(B_OLS),size(subsample))
Q_prmean    = ((kappa_Q).^2)*tau*VB_OLS;
Q_prvar     = tau;

% Sigma is the covariance of the VAR covariance, SIGMA
% Sigma ~ IW(I,M+1)
Sigma_prmean = eye(M);
Sigma_prvar  = M+1;

%========= INITIALIZE MATRICES:
% Specify covariance matrices for measurement and state equations
consQ       = 0.0001;
Qdraw       = consQ*eye(K);
Qchol       = sqrt(consQ)*eye(K);
Btdraw      = zeros(K,t);
Sigmadraw   = 0.1*eye(M);

% Storage matrices for posteriors and stuff
Bt_postmean = zeros(K,t);
Qmean       = zeros(K,K);
Sigmamean   = zeros(M,M);

% Storage for trace plots and autocorrelation functions

%========= IMPULSE RESPONSES:
% Note that impulse response and related stuff involves a lot of storage
% and, hence, put istore=0 if you do not want them

bigj  = zeros(M,M*p);
bigj(1:M,1:M) = eye(M);

imp_save_all = zeros(nrep,M,T,t);

%% Gibbs Sampler

tic;
disp('Number of iterations');

for irep = 1:nrep + nburn    % GIBBS iterations starts here
    % Print iterations
    
    if mod(irep,it_print) == 0
        disp(irep);toc;
    end
    
    % -----------------------------------------------------------------------------------------
    %   STEP I: Sample B_t from p(B_t|y,Sigma) (Drawing coefficient states, pp. 844-845)
    % -----------------------------------------------------------------------------------------
    [Btdraw,log_lik] = carter_kohn_hom(y,Z,Sigmadraw,Qdraw,K,M,t,B_0_prmean,B_0_prvar);
    % This follows Appendix A.6 in Primiceri (2005)
    
    Btemp = Btdraw(:,2:t)' - Btdraw(:,1:t-1)';
    sse_2Q = zeros(K,K);
    for i = 1:t-1
        sse_2Q = sse_2Q + Btemp(i,:)'*Btemp(i,:);
    end
    
    Qinv        = inv(sse_2Q + Q_prmean);
    %Qinvdraw    = wish(Qinv,t+Q_prvar);
    Qinvdraw    = wish((Qinv+1e-3*eye(size(Qinv))), t+Q_prvar);
    Qdraw       = inv(Qinvdraw);
    Qchol       = chol(Qdraw);
    
    % -----------------------------------------------------------------------------------------
    %   STEP I: Sample Sigma from p(Sigma|y,B_t) which is i-Wishart
    % ----------------------------------------------------------------------------------------
    
    yhat        = zeros(M,t);
    Btdraw_2    = [Btdraw(1:6,:);zeros(6,size(Btdraw,2));Btdraw(13:end,:)];
    for i = 1:t
        yhat(:,i)   = y(:,i) - Z((i-1)*M+1:i*M,:)*Btdraw(:,i);
    end
    
    sse_2S      = zeros(M,M);
    for i = 1:t
        sse_2S      = sse_2S    + yhat(:,i)*yhat(:,i)';
    end
    
    Sigmainv        = inv(sse_2S + Sigma_prmean);
    Sigmainvdraw    = wish(Sigmainv,t+Sigma_prvar);
    Sigmadraw       = inv(Sigmainvdraw);
    Sigmachol       = chol(Sigmadraw);
    
    %----------------------------SAVE AFTER-BURN-IN DRAWS AND IMPULSE RESPONSES -----------------
    if irep > nburn
        % Save only the means of B(t), Q and SIGMA
        Bt_postmean = Bt_postmean + Btdraw;
        Qmean       = Qmean + Qdraw;
        Sigmamean   = Sigmamean + Sigmadraw;
        
        if istore == 1
            % Impulse response analysis

            for i = 1:t % Get impulses recursively for each time period
                    
                biga = zeros(M*p,M*p);
                for j = 1:p-1
                    biga(j*M+1:M*(j+1),M*(j-1)+1:j*M) = eye(M);
                end

                bbtemp = Btdraw(M+M+1:K,i);  % Get the draw of B(t) at time i=1,...,T  (exclude intercept and exogenous variables)

                splace = 0;

                for ii = 1:p
                    for iii = 1:M
                        biga(iii,(ii-1)*M+1:ii*M) = bbtemp(splace+1:splace+M,1)';
                        splace = splace + M;
                    end
                end

                % Normalization to 20bp in 1991 M1 (0.56... is the
                % posterior mean estimate of the FFR for 1991 M1)
                shock = Btdraw(M+1:2*M,i)*0.2/0.567939496315211;

                % Now get impulse responses for 1 through IRF_T future periods
                impresp = zeros(M,T);  % matrix to store initial response at each period
                impresp(:,1) = shock;
                bigai = biga;
                for j = 2:T
                    impresp(:,j) = bigj*bigai*bigj'*shock;
                    bigai = bigai*biga;
                end
                imp_save_all(irep-nburn,:,:,i) = impresp;                               
            end
        end
    end     
end

clc;
toc; % Stop timer and print total time

%% Save File

%cd '..\Output\TVP'

save('TVP_Solution.mat');

%cd '..\..\Time-Varying Parameter VAR'

%% Compute Impulse Responses

Bt_postmean = Bt_postmean./nrep;   % Posterior mean of B(t) (VAR regression coeff.)
Qmean       = Qmean./nrep;         % Posterior mean of Q (covariance of B(t))
Sigmamean   = Sigmamean./nrep;     % Posterior mean of SIGMA (VAR covariance matrix)

imp_all_norm    = zeros(t,M,T);
imp_all_chol    = zeros(t,M,T);

if cholesky == 0
    for i = 1:t

        % Normalization to 20bp shock in 1991 M1
        shock   = Bt_postmean(M+1:2*M,i)*0.2/Bt_postmean(M+1,1);

        % Impulse responses with constant posterior mean of B for each t
        bbtemp  = Bt_postmean(M+M+1:K,i);
        splace  = 0;

        biga = zeros(M*p,M*p);

        for j = 1:p-1
            biga(j*M+1:M*(j+1),M*(j-1)+1:j*M) = eye(M);
        end
        
        for ii = 1:p
            for iii = 1:M
                biga(iii,(ii-1)*M+1:ii*M) = bbtemp(splace+1:splace+M,1)';
                splace = splace + M;
            end
        end

        impresp      = zeros(M,T);
        impresp(:,1) = shock;
        bigai        = biga;

        for j = 2:T
            impresp(:,j) = bigj*bigai*bigj'*shock;
            bigai = bigai*biga;
        end

        imp_all_norm(i,:,:) = impresp;

    end
end

if cholesky == 1
    for i = 1:t    
        bbtemp = Btdraw(M+1:K,i);
        splace = 0;

        biga = zeros(M*p,M*p);

        for j = 1:p-1
            biga(j*M+1:M*(j+1),M*(j-1)+1:j*M) = eye(M);
        end

        for ii = 1:p
            for iii = 1:M
                biga(iii,(ii-1)*M+1:ii*M) = bbtemp(splace+1:splace+M,1)';
                splace = splace + M;
            end
        end

        impresp      = zeros(M,T);
        
        Sigmachol   = chol(Sigmamean);
        shock       = Sigmachol';   % First shock is the Cholesky of the VAR covariance
        diagonal    = diag(diag(shock));
        shock       = inv(diagonal)*shock;    % Unit initial shock 

        % Now get impulse responses for 1 through IRF_T future periods
        
        impresp(:,1) = shock(:,4);
        bigai        = biga;

        for j = 2:T
            impresp(:,j) = bigj*bigai*bigj'*shock(:,4);
            bigai = bigai*biga;
        end
        
        imp_all_chol(i,:,:) = impresp;
    end
end

for vv=1:M
	if differencing(1,vv) == 1
        imp_all_norm(:,vv,:) =cumsum(imp_all_norm(:,vv,:),3);
        imp_all_chol(:,vv,:) =cumsum(imp_all_chol(:,vv,:),3);
	end
end

if cholesky == 1
    imp_all     = imp_all_chol;
else
    imp_all     = imp_all_norm;
end

imp_1       = zeros(size(imp_all,1),size(imp_all,3));
imp_1(:,:)  = imp_all(:,1,:);
imp_1_surf  = rot90(imp_1(:,1:IRF_T),2);

imp_2       = zeros(size(imp_all,1),size(imp_all,3));
imp_2(:,:)  = imp_all(:,2,:);
imp_2_surf  = rot90(imp_2(:,1:IRF_T),2);

imp_3       = zeros(size(imp_all,1),size(imp_all,3));
imp_3(:,:)  = imp_all(:,3,:);
imp_3_surf  = rot90(imp_3(:,1:IRF_T),2);
    
imp_4       = zeros(size(imp_all,1),size(imp_all,3));
imp_4(:,:)  = imp_all(:,4,:);
imp_4_surf  = rot90(imp_4(:,1:IRF_T),2);

imp_5       = zeros(size(imp_all,1),size(imp_all,3));
imp_5(:,:)  = imp_all(:,5,:);
imp_5_surf  = rot90(imp_5(:,1:IRF_T),2);

%imp_6       = zeros(size(imp_all,1),size(imp_all,3));
%imp_6(:,:)  = imp_all(:,6,:);
%imp_6_surf  = rot90(imp_6(:,1:IRF_T),2);

%% Compute Real Interest Rate and Fundamental Stock Price Response

% Computing Real Interest Rate Response
%imp_7       = zeros(size(imp_all,1),size(imp_all,3));
%for jj=1:IRF_T
    %imp_7(:,jj)     =imp_1(:,jj)-(imp_5(:,jj+12)-imp_5(:,jj));
%end
%imp_7_surf = rot90(imp_7(:,1:IRF_T),2);

% Computing Fundamental Stock Price Response
%Lambda      = 1/(1+0.04/12);
%imp_8       = zeros(size(imp_all,1),size(imp_all,3));
%J           = 250;

%for k=0:IRF_T

     %Funda=zeros(t,J+1);

     %for j=0:J
          %Funda(:,j+1)=Lambda^j*[(1-Lambda)*imp_3(:,k+j+2)-imp_7(:,k+j+1)./12];
     %end

     %imp_8(:,k+1) = sum(Funda,2);

%end
%imp_8_surf  = rot90(imp_8(:,1:IRF_T),2);

%imp_9       = imp_2-imp_8;    
%imp_9_surf  = rot90(imp_9(:,1:IRF_T),2);

%% Figure: 3D Impulse Responses

cd 'C:\Users\osull\Dropbox\Economics\Postgrad\Dissertation\MATLAB\Paul (2019) Code\Replication Files\Output\TVP'

figure(1)%FED FUNDS
surf(rot90(imp_1_surf,3))
if colored_figures == 0
    colormap('gray')
end
xlim([1 t])
ylim([1 IRF_T])
set(gca,'YTick',1:10:IRF_T)
set(gca,'YTickLabel',{'20','10'})%set(gca,'YTickLabel',{'40','30','20','10'})
if start_year_TVP == 1990 && end_year == 2007
    set(gca,'XTick',13:24:t)
    set(gca,'XTickLabel',{'1991','1993','1995','1997','1999','2001','2003','2005','2007'})
elseif start_year_TVP == 1991 && end_year >= 2016
    set(gca,'XTick',13:36:t)
    set(gca,'XTickLabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016'})
end
axis tight
set(gca,'Box','off');
set(gca,'Fontsize',20);
set(gcf, 'Position', get(0,'Screensize'));
y_pos = get(get(gca, 'YLabel'), 'Position');
set(get(gca, 'YLabel'), 'Position', y_pos + [5 0 0]); % adjust position [5 0 0]
title_fig = title(labels(1,1));
set(title_fig,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
ylabel('Month','Interpreter','Latex','Fontsize', 24,'FontName','Palatino')
set(get(gca,'ylabel'),'rotation',0)
zlabel('Percent','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
h = gcf;
set(h,'PaperPositionMode','auto')
set(h,'PaperOrientation','portrait')
set(h,'Position',[0 0 1500 1000])
tightfig
view(20,40) %Remove to get original plots
if choose_set_figures == 1
    print(h,'-depsc','-r300')
elseif choose_set_figures == 2
    print(h,'-djpeg','..\Output\TVP\Figure_15a','-r75')
elseif choose_set_figures == 3
    print(h,'-djpeg','..\Output\TVP\Figure_20a','-r75')
elseif choose_set_figures == 4
    print(h,'-djpeg','..\Output\TVP\Figure_23a','-r75')
elseif choose_set_figures == 5
    print(h,'-djpeg','..\Output\TVP\Figure_26a','-r75')
end

figure(2)%GOLD PRICE
surf(rot90(imp_2_surf, 3))%surf(rot90(imp_7_surf,3))
if colored_figures == 0
    colormap('gray')
end
xlim([1 t])
ylim([1 IRF_T])
set(gca,'YTick',1:10:IRF_T)
set(gca,'YTickLabel',{'20','10'})%set(gca,'YTickLabel',{'40','30','20','10'})
if start_year_TVP == 1990 && end_year == 2007
    set(gca,'XTick',13:24:t)
    set(gca,'XTickLabel',{'1991','1993','1995','1997','1999','2001','2003','2005','2007'})
elseif start_year_TVP == 1991 && end_year >= 2016
    set(gca,'XTick',13:36:t)
    set(gca,'XTickLabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016'})
end
axis tight
set(gca,'Box','off');
set(gca,'Fontsize',20);
set(gcf, 'Position', get(0,'Screensize'));
y_pos = get(get(gca, 'YLabel'), 'Position');
set(get(gca, 'YLabel'), 'Position', y_pos + [5 0 0]); % adjust position
title_fig = title(labels(1,2));
set(title_fig,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
ylabel('Month','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
set(get(gca,'ylabel'),'rotation',0)
zlabel('Percent','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
h = gcf;
set(h,'PaperPositionMode','auto')
set(h,'PaperOrientation','portrait')
set(h,'Position',[0 0 1500 1000])
tightfig
view(20,40)
if choose_set_figures == 1
    print(h,'-depsc','-r300')
elseif choose_set_figures == 2
    print(h,'-djpeg','..\Output\TVP\Figure_15b','-r75')
elseif choose_set_figures == 3
    print(h,'-djpeg','..\Output\TVP\Figure_20b','-r75')
elseif choose_set_figures == 4
    print(h,'-djpeg','..\Output\TVP\Figure_23b','-r75')
elseif choose_set_figures == 5
    print(h,'-djpeg','..\Output\TVP\Figure_26b','-r75')
end

figure(3)%YIELD SPREAD
surf(rot90(imp_3_surf,3))%surf(rot90(imp_2_surf,3))
if colored_figures == 0
    colormap('gray')
end
xlim([1 t])
ylim([1 IRF_T])
set(gca,'YTick',1:10:IRF_T)
set(gca,'YTickLabel',{'20','10'})%set(gca,'YTickLabel',{'40','30','20','10'})
if start_year_TVP == 1990 && end_year == 2007
    set(gca,'XTick',13:24:t)
    set(gca,'XTickLabel',{'1991','1993','1995','1997','1999','2001','2003','2005','2007'})
elseif start_year_TVP == 1991 && end_year >= 2016
    set(gca,'XTick',13:36:t)
    set(gca,'XTickLabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016'})
end
axis tight
set(gca,'Box','off');
set(gca,'Fontsize',20);
set(gcf, 'Position', get(0,'Screensize'));
y_pos = get(get(gca, 'YLabel'), 'Position');
set(get(gca, 'YLabel'), 'Position', y_pos + [5 0 0]); % adjust position
title_fig = title(labels(1,3));
set(title_fig,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
ylabel('Month','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
set(get(gca,'ylabel'),'rotation',0)
zlabel('Percent','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
h = gcf;
set(h,'PaperPositionMode','auto')
set(h,'PaperOrientation','portrait')
set(h,'Position',[0 0 1500 1000])
tightfig
view(60,65)
%if choose_set_figures == 1
    print(h,'-depsc','-r300')
%elseif choose_set_figures == 2
%    print(h,'-djpeg','..\Output\TVP\Figure_15c','-r75')
%elseif choose_set_figures == 3
%    print(h,'-djpeg','..\Output\TVP\Figure_20c','-r75')
%elseif choose_set_figures == 4
%    print(h,'-djpeg','..\Output\TVP\Figure_23c','-r75')
%elseif choose_set_figures == 5
%    print(h,'-djpeg','..\Output\TVP\Figure_26c','-r75')
%end

figure(4)%IndProd
surf(rot90(imp_4_surf,3))%surf(rot90(imp_3_surf,3))
if colored_figures == 0
    colormap('gray')
end
xlim([1 t])
ylim([1 IRF_T])
set(gca,'YTick',1:10:IRF_T)
set(gca,'YTickLabel',{'20','10'})%set(gca,'YTickLabel',{'40','30','20','10'})
if start_year_TVP == 1990 && end_year == 2007
    set(gca,'XTick',13:24:t)
    set(gca,'XTickLabel',{'1991','1993','1995','1997','1999','2001','2003','2005','2007'})
elseif start_year_TVP == 1991 && end_year >= 2016
    set(gca,'XTick',13:36:t)
    set(gca,'XTickLabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016'})
end
axis tight
set(gca,'Box','off');
set(gca,'Fontsize',20);
set(gcf, 'Position', get(0,'Screensize'));
y_pos = get(get(gca, 'YLabel'), 'Position');
set(get(gca, 'YLabel'), 'Position', y_pos + [5 0 0]); % adjust position
title_fig = title(labels(1,4));
set(title_fig,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
ylabel('Month','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
set(get(gca,'ylabel'),'rotation',0)
zlabel('Percent','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
h = gcf;
set(h,'PaperPositionMode','auto')
set(h,'PaperOrientation','portrait')
set(h,'Position',[0 0 1500 1000])
view(20, 40)
tightfig
%if choose_set_figures == 1
    print(h,'-depsc','-r300')
%elseif choose_set_figures == 2
%    print(h,'-djpeg','..\Output\TVP\Figure_15d','-r75')
%elseif choose_set_figures == 3
%    print(h,'-djpeg','..\Output\TVP\Figure_20d','-r75')
%elseif choose_set_figures == 4
%    print(h,'-djpeg','..\Output\TVP\Figure_23d','-r75')
%elseif choose_set_figures == 5
%    print(h,'-djpeg','..\Output\TVP\Figure_26d','-r75')
%end

figure(5)%CPI
surf(rot90(imp_5_surf,3))%surf(rot90(imp_4_surf,3))
if colored_figures == 0
    colormap('gray')
end
xlim([1 t])
ylim([1 IRF_T])
set(gca,'YTick',1:10:IRF_T)
set(gca,'YTickLabel',{'20','10'})%set(gca,'YTickLabel',{'40','30','20','10'})
if start_year_TVP == 1990 && end_year == 2007
    set(gca,'XTick',13:24:t)
    set(gca,'XTickLabel',{'1986','1989','1991','1993','1995','1997','1999','2001','2003','2005','2007'})
elseif start_year_TVP == 1991 && end_year >= 2016
    set(gca,'XTick',13:36:t)
    set(gca,'XTickLabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016'})
end
axis tight
set(gca,'Box','off');
set(gca,'Fontsize',20);
set(gcf, 'Position', get(0,'Screensize'));
y_pos = get(get(gca, 'YLabel'), 'Position');
set(get(gca, 'YLabel'), 'Position', y_pos + [5 0 0]); % adjust position
title_fig = title(labels(1,5));
set(title_fig,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
ylabel('Month','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
set(get(gca,'ylabel'),'rotation',0)
zlabel('Percent','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
h = gcf;
set(h,'PaperPositionMode','auto')
set(h,'PaperOrientation','portrait')
set(h,'Position',[0 0 1500 1000])
tightfig
view(35,60)
%if choose_set_figures == 1
    print(h,'-depsc','-r300')
%elseif choose_set_figures == 2
%    print(h,'-djpeg','..\Output\TVP\Figure_15e','-r75')
%elseif choose_set_figures == 3
%    print(h,'-djpeg','..\Output\TVP\Figure_20e','-r75')
%elseif choose_set_figures == 4
%    print(h,'-djpeg','..\Output\TVP\Figure_23e','-r75')
%elseif choose_set_figures == 5
%    print(h,'-djpeg','..\Output\TVP\Figure_26e','-r75')
%end

figure(6)
surf(rot90(imp_5_surf,3))
if colored_figures == 0
    colormap('gray')
end
xlim([1 t])
ylim([1 IRF_T])
set(gca,'YTick',1:10:IRF_T)
set(gca,'YTickLabel',{'40','30','20','10'})
if start_year_TVP == 1990 && end_year == 2007
    set(gca,'XTick',13:24:t)
    set(gca,'XTickLabel',{'1991','1993','1995','1997','1999','2001','2003','2005','2007'})
elseif start_year_TVP == 1991 && end_year >= 2016
    set(gca,'XTick',13:36:t)
    set(gca,'XTickLabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016'})
end
axis tight
set(gca,'Box','off');
set(gca,'Fontsize',20);
set(gcf, 'Position', get(0,'Screensize'));
y_pos = get(get(gca, 'YLabel'), 'Position');
set(get(gca, 'YLabel'), 'Position', y_pos + [5 0 0]); % adjust position
title_fig = title(labels(1,5));%title_fig = title(labels(1,6));
set(title_fig,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
ylabel('Months','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
set(get(gca,'ylabel'),'rotation',0)
zlabel('Percent','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
h = gcf;
set(h,'PaperPositionMode','auto')
set(h,'PaperOrientation','portrait')
set(h,'Position',[0 0 1500 1000])
tightfig
view(35,45)
%if choose_set_figures == 1
    %print(h,'-djpeg','-r75')
%elseif choose_set_figures == 2
%    print(h,'-djpeg','..\Output\TVP\Figure_15f','-r75')
%elseif choose_set_figures == 3
%    print(h,'-djpeg','..\Output\TVP\Figure_20f','-r75')
%elseif choose_set_figures == 4
%    print(h,'-djpeg','..\Output\TVP\Figure_23f','-r75')
%elseif choose_set_figures == 5
%    print(h,'-djpeg','..\Output\TVP\Figure_26f','-r75')
%end

%figure(7)
%surf(rot90(imp_6_surf,3))
%if colored_figures == 0
%    colormap('gray')
%end
%xlim([1 t])
%ylim([1 IRF_T])
%set(gca,'YTick',1:10:IRF_T)
%set(gca,'YTickLabel',{'40','30','20','10'})
%if start_year_TVP == 1990 && end_year == 2007
%    set(gca,'XTick',13:24:t)
%    set(gca,'XTickLabel',{'1991','1993','1995','1997','1999','2001','2003','2005','2007'})
%elseif start_year_TVP == 1991 && end_year >= 2016
%    set(gca,'XTick',13:36:t)
%    set(gca,'XTickLabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016'})
%end
%axis tight
%set(gca,'Box','off');
%set(gca,'Fontsize',20);
%set(gcf, 'Position', get(0,'Screensize'));
%y_pos = get(get(gca, 'YLabel'), 'Position');
%set(get(gca, 'YLabel'), 'Position', y_pos + [5 0 0]); % adjust position
%title_fig = title(labels(1,7));
%set(title_fig,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
%ylabel('Months','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
%set(get(gca,'ylabel'),'rotation',0)
%zlabel('Percent','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
%h = gcf;
%set(h,'PaperPositionMode','auto')
%set(h,'PaperOrientation','portrait')
%set(h,'Position',[0 0 1500 1000])
%tightfig
%if choose_set_figures == 1
%    print(h,'-djpeg','-r75')
%elseif choose_set_figures == 2
%    print(h,'-djpeg','..\Output\TVP\Figure_15g','-r75')
%elseif choose_set_figures == 3
%    print(h,'-djpeg','..\Output\TVP\Figure_20g','-r75')
%elseif choose_set_figures == 4
%    print(h,'-djpeg','..\Output\TVP\Figure_23g','-r75')
%elseif choose_set_figures == 5
%    print(h,'-djpeg','..\Output\TVP\Figure_26g','-r75')
%end

%figure(8)
%surf(rot90(imp_8_surf,3))
%if colored_figures == 0
%    colormap('gray')
%end
%xlim([1 t])
%ylim([1 IRF_T])
%set(gca,'YTick',1:10:IRF_T)
%set(gca,'YTickLabel',{'40','30','20','10'})
%if start_year_TVP == 1990 && end_year == 2007
%    set(gca,'XTick',13:24:t)
 %   set(gca,'XTickLabel',{'1991','1993','1995','1997','1999','2001','2003','2005','2007'})
%elseif start_year_TVP == 1991 && end_year >= 2016
%    set(gca,'XTick',13:36:t)
%    set(gca,'XTickLabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016'})
%end
%axis tight
%set(gca,'Box','off');
%set(gca,'Fontsize',20);
%set(gcf, 'Position', get(0,'Screensize'));
%y_pos = get(get(gca, 'YLabel'), 'Position');
%set(get(gca, 'YLabel'), 'Position', y_pos + [5 0 0]); % adjust position
%title_fig = title(labels(1,8));
%set(title_fig,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
%ylabel('Months','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
%set(get(gca,'ylabel'),'rotation',0)
%zlabel('Percent','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
%h = gcf;
%set(h,'PaperPositionMode','auto')
%set(h,'PaperOrientation','portrait')
%set(h,'Position',[0 0 1500 1000])
%tightfig
%if choose_set_figures == 1
%    print(h,'-djpeg','-r75')
%elseif choose_set_figures == 2
%    print(h,'-djpeg','..\Output\TVP\Figure_15h','-r75')
%elseif choose_set_figures == 3
%    print(h,'-djpeg','..\Output\TVP\Figure_20h','-r75')
%elseif choose_set_figures == 4
%    print(h,'-djpeg','..\Output\TVP\Figure_23h','-r75')
%elseif choose_set_figures == 5
%    print(h,'-djpeg','..\Output\TVP\Figure_26h','-r75')
%end

%figure(9)
%surf(rot90(imp_9_surf,3))
%if colored_figures == 0
%    colormap('gray')
%end
%xlim([1 t])
%ylim([1 IRF_T])
%set(gca,'YTick',1:10:IRF_T)
%set(gca,'YTickLabel',{'40','30','20','10'})
%if start_year_TVP == 1990 && end_year == 2007
%    set(gca,'XTick',13:24:t)
%    set(gca,'XTickLabel',{'1991','1993','1995','1997','1999','2001','2003','2005','2007'})
%elseif start_year_TVP == 1991 && end_year >= 2016
%    set(gca,'XTick',13:36:t)
%    set(gca,'XTickLabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016'})
%end
%axis tight
%set(gca,'Box','off');
%set(gca,'Fontsize',20);
%set(gcf, 'Position', get(0,'Screensize'));
%y_pos = get(get(gca, 'YLabel'), 'Position');
%set(get(gca, 'YLabel'), 'Position', y_pos + [5 0 0]); % adjust position
%title_fig = title(labels(1,9));
%set(title_fig,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
%ylabel('Months','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
%set(get(gca,'ylabel'),'rotation',0)
%zlabel('Percent','Interpreter','Latex','Fontsize',24,'FontName','Palatino')
%h = gcf;
%set(h,'PaperPositionMode','auto')
%set(h,'PaperOrientation','portrait')
%set(h,'Position',[0 0 1500 1000])
%tightfig
%if choose_set_figures == 1
%    print(h,'-djpeg','-r75')
%elseif choose_set_figures == 2
%    print(h,'-djpeg','..\Output\TVP\Figure_15i','-r75')
%elseif choose_set_figures == 3
%    print(h,'-djpeg','..\Output\TVP\Figure_20i','-r75')
%elseif choose_set_figures == 4
%    print(h,'-djpeg','..\Output\TVP\Figure_23i','-r75')
%elseif choose_set_figures == 5
%    print(h,'-djpeg','..\Output\TVP\Figure_26i','-r75')
%end

%% Figure: 2D Impulse Responses
%% 

% Data in log-levels
start_sample_yoy    = 12*(start_year_TVP-1960)+start_month_TVP;
end_sample_yoy      = size(macro,1)-12*(2017-end_year)-(9-end_month);
D_data              = [macro(:,1) macro(:,2) macro(:,3) macro(:,4) macro(:,5)];
for jj=1:size(ydata,2)
    if differencing(1,jj)==1
        D_data(:,jj) = 100*log(D_data(:,jj));
    end
end
data = D_data(start_sample_yoy:end_sample_yoy,:);

%% Start in 1990 M1

if choose_set_figures == 4

    stock_min   = min(imp_2_surf(:,28));
    stock_max   = max(imp_2_surf(:,28));

    stock_lb_1  = stock_min+(stock_min-stock_max)/20;
    stock_ub_1  = stock_max-(stock_min-stock_max)/20;

    house_min   = min(imp_4_surf(:,28));
    house_max   = max(imp_4_surf(:,28));

    house_lb_1  = house_min+(house_min-house_max)/20;
    house_ub_1  = house_max-(house_min-house_max)/20;

    stock_index_max = max(data(:,2)./100)+(max(data(:,2)./100)-min(data(:,2)./100))/20;
    stock_index_min = min(data(:,2)./100)-(max(data(:,2)./100)-min(data(:,2)./100))/20;   

    house_index_max = max(data(:,4)./100)+(max(data(:,2)./100)-min(data(:,2)./100))/20;
    house_index_min = min(data(:,4)./100)-(max(data(:,2)./100)-min(data(:,2)./100))/20;

    fig = figure;
    if colored_figures == 0
        left_color  = [0 0 0];
        right_color = [0.41 0.41 0.41];
        set(fig,'defaultAxesColorOrder',[left_color;right_color]);
    end
    subplot(2,1,1)
    yyaxis left
    h_area_1 = area([8 15],[stock_lb_1 stock_lb_1],'FaceColor',[.9 .9 .9]);
    hold on
    h_area_2 = area([135 140],[stock_lb_1 stock_lb_1],'FaceColor',[.9 .9 .9]);
    h_area_3 = area([216 234],[stock_lb_1 stock_lb_1],'FaceColor',[.9 .9 .9]);
    imp_3 = plot(flipud(imp_2_surf(:,28)),'Linewidth',3,'Color',rgb(irf_color),'Linestyle','-.');
    ylabel('Percent','Fontsize',20,'FontName','Palatino','Interpreter','Latex')
    ylim([stock_lb_1 stock_ub_1])
    yyaxis right
    plot(data(:,2)./100,'Linewidth',3,'Color',rgb(dash_color),'Linestyle','-') ;
    axis tight
    set(gca,'XTick',13:36:t)
    set(gca,'XTickLabel',{'1991','1994','1997','2000','2003','2006','2009','2012','2015'})
    tttt = title('Stock Prices');
    set(gca,'Box','on','FontSize', 20)
    ylabel('log(real S\&P 500)','Fontsize',20,'FontName','Palatino','Interpreter','Latex')
    set(tttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    legend(imp_3,{'1 year'},'Location','Northwest')
    ylim([stock_index_min stock_index_max])
    xlim([1 t])
    subplot(2,1,2)
    yyaxis left
    h_area_1 = area([8 15],[house_lb_1 house_lb_1], 'FaceColor', [.9 .9 .9]);
    hold on
    h_area_2 = area([135 140],[house_lb_1 house_lb_1], 'FaceColor', [.9 .9 .9]);
    h_area_3 = area([216 234],[house_lb_1 house_lb_1], 'FaceColor', [.9 .9 .9]);
    imp_6 = plot(flipud(imp_4_surf(:,28)),'Linewidth',3,'Color',rgb(irf_color),'Linestyle','-.');
    ylabel('Percent','Fontsize',20,'FontName','Palatino','Interpreter','Latex')
    ylim([house_lb_1 house_ub_1])
    xlim([1 t])
    yyaxis right
    plot(data(:,4)./100,'Linewidth',3,'Color',rgb(dash_color),'Linestyle','-') 
    hold on
    axis tight
    set(gca,'XTick',13:36:t)
    set(gca,'XTickLabel',{'1991','1994','1997','2000','2003','2006','2009','2012','2015'})
    tttt = title('House Prices');
    set(gca,'Box','on','FontSize', 20)
    ylabel('log(real Case-Shiller Index)','Fontsize',20,'FontName','Palatino','Interpreter','Latex')
    set(tttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    legend(imp_6,{'1 year'},'Location','Northwest')
    ylim([house_index_min house_index_max])
    xlim([1 t])
    h = gcf;
    set(h,'PaperPositionMode','auto')
    set(h,'PaperOrientation','portrait')
    set(h,'Position',[0 0 1500 1000])
    tightfig
    print(h,'-depsc','..\Output\TVP\Figure_24')
    
    stock_min   = min(imp_6_surf(:,28)./imp_2_surf(:,28));
    stock_max   = max(imp_6_surf(:,28)./imp_2_surf(:,28));
    
    stock_lb_2  = stock_min+(stock_min-stock_max)/20;
    stock_ub_2  = stock_max-(stock_min-stock_max)/20;
    
    house_min   = min(imp_6_surf(:,28)./imp_4_surf(:,28));
    house_max   = max(imp_6_surf(:,28)./imp_4_surf(:,28));
    
    house_lb_2  = house_min+(house_min-house_max)/20;
    house_ub_2  = house_max-(house_min-house_max)/20;
    
    fig = figure;
    if colored_figures == 0
        left_color  = [0 0 0];
        right_color = [0.41 0.41 0.41];
        set(fig,'defaultAxesColorOrder',[left_color;right_color]);
    end
    subplot(2,1,1)
    yyaxis left
    h_area_1 = area([8 15], [stock_ub_2 stock_ub_2],'FaceColor', [.9 .9 .9]);
    hold on
    h_area_2 = area([135 140], [stock_ub_2 stock_ub_2],'FaceColor', [.9 .9 .9]);
    h_area_3 = area([216 234], [stock_ub_2 stock_ub_2],'FaceColor', [.9 .9 .9]);
    imp_9 = plot(flipud(imp_6_surf(:,28)./imp_2_surf(:,28)),'Linewidth',3,'Color',rgb(irf_color),'Linestyle','-.');
    ylim([stock_lb_2 stock_ub_2])
    yyaxis right
    plot(data(:,2)./100,'Linewidth',3,'Color',rgb(dash_color),'Linestyle','-') 
    axis tight
    set(gca,'XTick',13:36:t)
    set(gca,'XTickLabel',{'1991','1994','1997','2000','2003','2006','2009','2012','2015'})
    tttt = title('Sacrifice Ratio - $\frac{\mathrm{Output}}{\mathrm{Stock \thinspace Prices}}$');
    set(gca,'Box','on','FontSize', 20)
    ylabel('log(real S\&P 500)','Fontsize',20,'FontName','Palatino','Interpreter','Latex')
    set(tttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    legend(imp_9,{'1 year'},'Location','Northwest')
    ylim([stock_index_min stock_index_max])
    xlim([1 t])
    subplot(2,1,2)
    yyaxis left
    h_area_1 = area([8 15], [house_ub_2 house_ub_2], 'FaceColor', [.9 .9 .9]);
    hold on
    h_area_2 = area([135 140], [house_ub_2 house_ub_2], 'FaceColor', [.9 .9 .9]);
    h_area_3 = area([216 234], [house_ub_2 house_ub_2], 'FaceColor', [.9 .9 .9]);
    imp_12 = plot(flipud(imp_6_surf(:,28)./imp_4_surf(:,28)),'Linewidth',3,'Color',rgb(irf_color),'Linestyle','-.');
    ylim([house_lb_2 house_ub_2])
    yyaxis right
    plot(data(:,4)./100,'Linewidth',3,'Color',rgb(dash_color),'Linestyle','-') 
    axis tight
    set(gca,'XTick',13:36:t)
    set(gca,'XTickLabel',{'1991','1994','1997','2000','2003','2006','2009','2012','2015'})
    tttt = title('Sacrifice Ratio - $\frac{\mathrm{Output}}{\mathrm{House \thinspace Prices}}$');
    set(gca,'Box','on','FontSize', 20)
    ylabel('log(real Case-Shiller Index)','Fontsize',20,'FontName','Palatino','Interpreter','Latex')
    set(tttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    legend(imp_12,{'1 year'},'Location','Northwest')
    ylim([house_index_min house_index_max])
    xlim([1 t])
    h = gcf;
    set(h,'PaperPositionMode','auto')
    set(h,'PaperOrientation','portrait')
    set(h,'Position',[0 0 1500 1000])
    tightfig
    print(h,'-depsc','..\Output\TVP\Figure_25')
    
end

%% Start in 1991 M1

if start_year_TVP == 1991 && choose_set_figures ~= 2
    
    stock_min   = min([imp_2_surf(:,4);imp_2_surf(:,28)]);
    stock_max   = max([imp_2_surf(:,4);imp_2_surf(:,28)]);
    
    stock_lb_1  = stock_min+(stock_min-stock_max)/20;
    stock_ub_1  = stock_max-(stock_min-stock_max)/20;
    
    house_min   = min([imp_4_surf(:,4);imp_4_surf(:,28)]);
    house_max   = max([imp_4_surf(:,4);imp_4_surf(:,28)]);
    
    house_lb_1  = house_min+(house_min-house_max)/20;
    house_ub_1  = house_max-(house_min-house_max)/20;

    stock_index_max = max(data(:,2)./100)+(max(data(:,2)./100)-min(data(:,2)./100))/20;
    stock_index_min = min(data(:,2)./100)-(max(data(:,2)./100)-min(data(:,2)./100))/20;   
    
    house_index_max = max(data(:,4)./100)+(max(data(:,2)./100)-min(data(:,2)./100))/20;
    house_index_min = min(data(:,4)./100)-(max(data(:,2)./100)-min(data(:,2)./100))/20;
    
    fig = figure;
    if colored_figures == 0
        left_color  = [0 0 0];
        right_color = [0.41 0.41 0.41];
        set(fig,'defaultAxesColorOrder',[left_color;right_color]);
    end
    subplot(2,1,1)
    yyaxis left
    h_area_1 = area([1 3], [stock_lb_1 stock_lb_1], 'FaceColor', [.9 .9 .9]);
    hold on
    h_area_2 = area([123 128], [stock_lb_1 stock_lb_1], 'FaceColor', [.9 .9 .9]);
    h_area_3 = area([204 222], [stock_lb_1 stock_lb_1], 'FaceColor', [.9 .9 .9]);
    imp_1 = plot(flipud(imp_2_surf(:,4)),'Linewidth',3,'Color',rgb(irf_color),'Linestyle',':');
    imp_3 = plot(flipud(imp_2_surf(:,28)),'Linewidth',3,'Color',rgb(irf_color),'Linestyle','-.');
    hold on
    ylabel('Percent','Fontsize',20,'FontName','Palatino','Interpreter','Latex')
    ylim([stock_lb_1 stock_ub_1])
    yyaxis right
    plot(data(:,2)./100,'Linewidth',3,'Color',rgb(dash_color),'Linestyle','-') ;
    axis tight
    set(gca,'XTick',13:36:t)
    set(gca,'XTickLabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016'})
    tttt = title('Stock Prices');
    set(gca,'Box','on','FontSize', 20)
    ylabel('log(real S\&P 500)','Fontsize',20,'FontName','Palatino','Interpreter','Latex')
    set(tttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    legend([imp_1 imp_3],{'3 years','1 year'},'Location','Northwest')
    ylim([stock_index_min stock_index_max])
    xlim([1 t])
    subplot(2,1,2)
    yyaxis left
    h_area_1 = area([1 3], [house_lb_1 house_lb_1], 'FaceColor', [.9 .9 .9]);
    hold on
    h_area_2 = area([123 128], [house_lb_1 house_lb_1], 'FaceColor', [.9 .9 .9]);
    h_area_3 = area([204 222], [house_lb_1 house_lb_1], 'FaceColor', [.9 .9 .9]);
    imp_4 = plot(flipud(imp_4_surf(:,4)),'Linewidth',3,'Color',rgb(irf_color),'Linestyle',':');
    imp_6 = plot(flipud(imp_4_surf(:,28)),'Linewidth',3,'Color',rgb(irf_color),'Linestyle','-.');
    ylabel('Percent','Fontsize',20,'FontName','Palatino','Interpreter','Latex')
    ylim([house_lb_1 house_ub_1])
    xlim([1 t])
    yyaxis right
    plot(data(:,4)./100,'Linewidth',3,'Color',rgb(dash_color),'Linestyle','-') 
    hold on
    axis tight
    set(gca,'XTick',13:36:t)
    set(gca,'XTickLabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016'})
    tttt = title('House Prices');
    set(gca,'Box','on','FontSize', 20)
    ylabel('log(real Case-Shiller Index)','Fontsize',20,'FontName','Palatino','Interpreter','Latex')
    set(tttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    legend([imp_4 imp_6],{'3 years','1 year'},'Location','Northwest')    
    ylim([house_index_min house_index_max])
    xlim([1 t])
    h = gcf;
    set(h,'PaperPositionMode','auto')
    set(h,'PaperOrientation','portrait')
    set(h,'Position',[0 0 1500 1000])
    tightfig
    if choose_set_figures == 1
        print(h,'-depsc','..\Output\TVP\Figure_4')
    elseif choose_set_figures == 3
        print(h,'-depsc','..\Output\TVP\Figure_21')
    elseif choose_set_figures == 5
        print(h,'-depsc','..\Output\TVP\Figure_27')
    end

    stock_min   = min([imp_6_surf(:,4)./imp_2_surf(:,4);imp_6_surf(:,28)./imp_2_surf(:,28)]);
    stock_max   = max([imp_6_surf(:,4)./imp_2_surf(:,4);imp_6_surf(:,28)./imp_2_surf(:,28)]);

    stock_lb_2  = stock_min+(stock_min-stock_max)/20;
    stock_ub_2  = stock_max-(stock_min-stock_max)/20;

    house_min   = min([imp_6_surf(:,4)./imp_4_surf(:,4);imp_6_surf(:,28)./imp_4_surf(:,28)]);
    house_max   = max([imp_6_surf(:,4)./imp_4_surf(:,4);imp_6_surf(:,28)./imp_4_surf(:,28)]);

    house_lb_2  = house_min+(house_min-house_max)/20;
    house_ub_2  = house_max-(house_min-house_max)/20;

    fig = figure;
    if colored_figures == 0
        left_color  = [0 0 0];
        right_color = [0.41 0.41 0.41];
        set(fig,'defaultAxesColorOrder',[left_color;right_color]);
    end
    subplot(2,1,1)
    yyaxis left
    h_area_1 = area([1 3], [stock_ub_2 stock_ub_2], 'FaceColor', [.9 .9 .9]);
    hold on
    h_area_2 = area([123 128], [stock_ub_2 stock_ub_2], 'FaceColor', [.9 .9 .9]);
    h_area_3 = area([204 222], [stock_ub_2 stock_ub_2], 'FaceColor', [.9 .9 .9]);
    imp_7 = plot(flipud(imp_6_surf(:,4)./imp_2_surf(:,4)),'Linewidth',3,'Color',rgb(irf_color),'Linestyle',':');
    imp_9 = plot(flipud(imp_6_surf(:,28)./imp_2_surf(:,28)),'Linewidth',3,'Color',rgb(irf_color),'Linestyle','-.');
    ylim([stock_lb_2 stock_ub_2])
    yyaxis right
    plot(data(:,2)./100,'Linewidth',3,'Color',rgb(dash_color),'Linestyle','-') 
    axis tight
    set(gca,'XTick',13:36:t)
    set(gca,'XTickLabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016'})
    tttt = title('Sacrifice Ratio - $\frac{\mathrm{Output}}{\mathrm{Stock \thinspace Prices}}$');
    set(gca,'Box','on','FontSize', 20)
    ylabel('log(real S\&P 500)','Fontsize',20,'FontName','Palatino','Interpreter','Latex')
    set(tttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    legend([imp_7 imp_9],{'3 years','1 year'},'Location','Northwest')
    ylim([stock_index_min stock_index_max])
    xlim([1 t])
    subplot(2,1,2)
    yyaxis left
    h_area_1 = area([1 3], [house_ub_2 house_ub_2], 'FaceColor', [.9 .9 .9]);
    hold on
    h_area_2 = area([123 128], [house_ub_2 house_ub_2], 'FaceColor', [.9 .9 .9]);
    h_area_3 = area([204 222], [house_ub_2 house_ub_2], 'FaceColor', [.9 .9 .9]);
    imp_10  = plot(flipud(imp_6_surf(:,4)./imp_4_surf(:,4)),'Linewidth',3,'Color',rgb(irf_color),'Linestyle',':');
    imp_12  = plot(flipud(imp_6_surf(:,28)./imp_4_surf(:,28)),'Linewidth',3,'Color',rgb(irf_color),'Linestyle','-.');
    ylim([house_lb_2 house_ub_2])
    yyaxis right
    plot(data(:,4)./100,'Linewidth',3,'Color',rgb(dash_color),'Linestyle','-') 
    axis tight
    set(gca,'XTick',13:36:t)
    set(gca,'XTickLabel',{'1992','1995','1998','2001','2004','2007','2010','2013','2016'})
    tttt = title('Sacrifice Ratio - $\frac{\mathrm{Output}}{\mathrm{House \thinspace Prices}}$');
    set(gca,'Box','on','FontSize', 20)
    ylabel('log(real Case-Shiller Index)','Fontsize',20,'FontName','Palatino','Interpreter','Latex')
    set(tttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    legend([imp_10 imp_12],{'3 years','1 year'},'Location','Northwest')    
    ylim([house_index_min house_index_max])
    xlim([1 t])
    h = gcf;
    set(h,'PaperPositionMode','auto')
    set(h,'PaperOrientation','portrait')
    set(h,'Position',[0 0 1500 1000])
    tightfig
    if choose_set_figures == 1
        print(h,'-depsc','..\Output\TVP\Figure_5')
    elseif choose_set_figures == 3
        print(h,'-depsc','..\Output\TVP\Figure_22')
    elseif choose_set_figures == 5
        print(h,'-depsc','..\Output\TVP\Figure_28')
    end

end

%% Figure 16: Impulse Responses for selected period with confidence bands

if istore == 1
   
    %start_confidence    = 12*(year_confidence-start_year_TVP)+month_confidence-start_month_TVP+1;
    start_confidence    = 12*(2016-start_year_TVP)+month_confidence-start_month_TVP+1;
    imp_confidence      = zeros(nrep,M,T);
    labels_confidence   = {'Federal Funds Rate', 'Gold Price', 'Yield Spread', 'Industrial Production', 'Consumer Price Index'};
    %labels_confidence   = {'Federal Funds Rate','Stock Prices','Dividends','House Prices','Consumer Price Index','Industrial Production'};
    
    for vv=1:M
        if differencing(1,vv) == 1
            imp_confidence(:,vv,:) = cumsum(imp_save_all(:,vv,:,start_confidence),3);
        else
            imp_confidence(:,vv,:) = imp_save_all(:,vv,:,start_confidence);
        end
    end
    
    zero_line = zeros(1,IRF_T);
    
    IRF_chosen_1    = squeeze(prctile(imp_confidence,2.5));
    IRF_chosen_2    = squeeze(prctile(imp_confidence,16));
    IRF_chosen_3    = squeeze(prctile(imp_confidence,50));
    IRF_chosen_4    = squeeze(prctile(imp_confidence,84));
    IRF_chosen_5    = squeeze(prctile(imp_confidence,97.5));
    
    % Computing Real Interest Rate Response
    imp_RR       = zeros(1,size(IRF_chosen_2,2));
    for jj=1:IRF_T
        imp_RR(1,jj)     =IRF_chosen_2(1,jj)-(IRF_chosen_2(5,jj+12)-IRF_chosen_2(5,jj));
    end

    % Computing Fundamental Stock Price Response
    Lambda       = 1/(1+0.04/12);
    imp_FS       = zeros(1,size(IRF_chosen_2,2));

    for k=0:IRF_T

         Funda=zeros(1,J+1);

         for j=0:J
              Funda(1,j+1)=Lambda^j*[(1-Lambda)*IRF_chosen_2(3,k+j+2)-imp_RR(1,k+j+1)./12];
         end

         imp_FS(1,k+1) = sum(Funda);

    end
    
    figure(16)
    for mm = 1:5
        subplot(3,2,mm)
        plot(IRF_chosen_1(mm,:),'Color',rgb(irf_color),'LineWidth',1,'LineStyle','-')
        hold on
        plot(IRF_chosen_2(mm,:),'Color',rgb(irf_color),'LineWidth',1,'LineStyle','-')
        plot(IRF_chosen_3(mm,:),'Color',rgb(irf_color),'LineWidth',3,'LineStyle','-')
        plot(IRF_chosen_4(mm,:),'Color',rgb(irf_color),'LineWidth',1,'LineStyle','-')
        plot(IRF_chosen_5(mm,:),'Color',rgb(irf_color),'LineWidth',1,'LineStyle','-')
        %if mm == 1
        %    plot(imp_RR(1,:),'Linewidth',3,'Linestyle','-.','Color',rgb(dash_color))
        %elseif mm == 2
        %    plot(imp_FS(1,:),'Linewidth',3,'Linestyle','-.','Color',rgb(dash_color))
        %end
        ciplot(IRF_chosen_1(mm,:),IRF_chosen_5(mm,:),band_color);
        ciplot(IRF_chosen_2(mm,:),IRF_chosen_4(mm,:),band_color);
        plot(zero_line,'Color',rgb('Black'))
        axis tight
        xlim([1 IRF_T])
        ttt= title(labels_confidence(1,mm));
        set(gca,'Box','on','FontSize', 16)
        ylabel('Percent','Fontsize',20,'FontName','Palatino','Interpreter','Latex')
        if mm == 5 || mm == 6
            xlabel('Months','Fontsize',20,'FontName','Palatino','Interpreter','Latex')
        end
        set(ttt,'fontsize',24,'FontName','Palatino','Interpreter','Latex')
        grid on
        set(gca,'GridLineStyle','--')
        ax = gca;
        ax.GridAlpha = .5;
    end
    h = gcf;
    set(h,'PaperPositionMode','auto')
    set(h,'PaperOrientation','landscape')
    set(h,'Position',[0 0 1500 1000])
    tightfig
    print(h,'-depsc','..\Output\TVP\Figure_16')

end

%% Figures 17-19: Differences in Impulse Responses

if istore == 1
    
    %% Preliminaries Differences Asset Prices
    
    date_stock_1  = 12*(year_compare_stock_1-start_year_TVP)+month_compare_stock_1-start_month_TVP+1;
    date_stock_2  = 12*(year_compare_stock_2-start_year_TVP)+month_compare_stock_2-start_month_TVP+1;
    date_stock_3  = 12*(year_compare_stock_3-start_year_TVP)+month_compare_stock_3-start_month_TVP+1;
    
    date_house_1  = 12*(year_compare_house_1-start_year_TVP)+month_compare_house_1-start_month_TVP+1;
    date_house_2  = 12*(year_compare_house_2-start_year_TVP)+month_compare_house_2-start_month_TVP+1;
    date_house_3  = 12*(year_compare_house_3-start_year_TVP)+month_compare_house_3-start_month_TVP+1;
    
    imp_stock_1(:,:) = cumsum(imp_save_all(:,2,1:40,date_stock_1),3);
    imp_stock_2(:,:) = cumsum(imp_save_all(:,2,1:40,date_stock_2),3);
    imp_stock_3(:,:) = cumsum(imp_save_all(:,2,1:40,date_stock_3),3);
    
    imp_house_1(:,:) = cumsum(imp_save_all(:,4,1:40,date_house_1),3);
    imp_house_2(:,:) = cumsum(imp_save_all(:,4,1:40,date_house_2),3);
    imp_house_3(:,:) = cumsum(imp_save_all(:,4,1:40,date_house_3),3);

    %% Figure 17: Differences Stock Prices
    
    figure(17)
    subplot(2,2,1)
    plot(median(imp_stock_1),'Linewidth',1,'Color',rgb(irf_color),'Marker','o')
    hold on
    plot(median(imp_stock_2),'Linewidth',1,'Color',rgb(irf_color),'Marker','x')
    plot(median(imp_stock_3),'Linewidth',1,'Color',rgb(irf_color))
    legend('1991 M1','2002 M1','2007 M12','location','Northeast')
    set(gca,'XTick',0:10:IRF_T)
    xlim([1 IRF_T])
    set(gca,'Box','on','FontSize',24)
    ylabel('Percent','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    ttt= title('Stock Prices');
    set(ttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    grid on
    set(gca,'GridLineStyle','--')
    ax = gca;
    ax.GridAlpha = .5;
    
    subplot(2,2,2)
    plot(squeeze(prctile(imp_stock_1-imp_stock_2,16)),'Linewidth',1,'Color',rgb(irf_color));
    hold on
    plot(squeeze(prctile(imp_stock_1-imp_stock_2,31)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_stock_1-imp_stock_2,50)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_stock_1-imp_stock_2,69)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_stock_1-imp_stock_2,84)),'Linewidth',1,'Color',rgb(irf_color));
    ciplot(squeeze(prctile(imp_stock_1-imp_stock_2,16)),squeeze(prctile(imp_stock_1-imp_stock_2,84)),band_color);
    ciplot(squeeze(prctile(imp_stock_1-imp_stock_2,31)),squeeze(prctile(imp_stock_1-imp_stock_2,69)),band_color);
    set(gca,'XTick',0:10:IRF_T)
    axis tight
    xlim([1 IRF_T])
    ylim([min(min(squeeze(prctile(imp_stock_1-imp_stock_2,16))),min(squeeze(prctile(imp_stock_1-imp_stock_2,84))))-1 ...
          max(max(squeeze(prctile(imp_stock_1-imp_stock_2,16))),max(squeeze(prctile(imp_stock_1-imp_stock_2,84))))+1])
    set(gca,'Box','on','FontSize',24)
    ylabel('Percent','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    ttt= title('1991 M1-2002 M1');
    set(ttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    grid on
    set(gca,'GridLineStyle','--')
    ax = gca;
    ax.GridAlpha = .5;
    
    subplot(2,2,3)
    plot(squeeze(prctile(imp_stock_1-imp_stock_3,16)),'Linewidth',1,'Color',rgb(irf_color));
    hold on
    plot(squeeze(prctile(imp_stock_1-imp_stock_3,31)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_stock_1-imp_stock_3,50)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_stock_1-imp_stock_3,69)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_stock_1-imp_stock_3,84)),'Linewidth',1,'Color',rgb(irf_color));
    ciplot(squeeze(prctile(imp_stock_1-imp_stock_3,16)),squeeze(prctile(imp_stock_1-imp_stock_3,84)),band_color);
    ciplot(squeeze(prctile(imp_stock_1-imp_stock_3,31)),squeeze(prctile(imp_stock_1-imp_stock_3,69)),band_color);
    set(gca,'XTick',0:10:IRF_T)
    axis tight
    xlim([1 IRF_T])
    ylim([min(min(squeeze(prctile(imp_stock_1-imp_stock_3,16))),min(squeeze(prctile(imp_stock_1-imp_stock_3,84))))-1 ...
          max(max(squeeze(prctile(imp_stock_1-imp_stock_3,16))),max(squeeze(prctile(imp_stock_1-imp_stock_3,84))))+1])
    set(gca,'Box','on','FontSize',24)
    ylabel('Percent','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    xlabel('Months','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    ttt= title('1991 M1-2007 M12');
    set(ttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    grid on
    set(gca,'GridLineStyle','--')
    ax = gca;
    ax.GridAlpha = .5;
    
    subplot(2,2,4)
    plot(squeeze(prctile(imp_stock_2-imp_stock_3,16)),'Linewidth',1,'Color',rgb(irf_color));
    hold on
    plot(squeeze(prctile(imp_stock_2-imp_stock_3,31)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_stock_2-imp_stock_3,50)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_stock_2-imp_stock_3,69)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_stock_2-imp_stock_3,84)),'Linewidth',1,'Color',rgb(irf_color));
    ciplot(squeeze(prctile(imp_stock_2-imp_stock_3,16)),squeeze(prctile(imp_stock_2-imp_stock_3,84)),band_color);
    ciplot(squeeze(prctile(imp_stock_2-imp_stock_3,31)),squeeze(prctile(imp_stock_2-imp_stock_3,69)),band_color);
    set(gca,'XTick',0:10:IRF_T)
    axis tight
    xlim([1 IRF_T])
    ylim([min(min(squeeze(prctile(imp_stock_2-imp_stock_3,16))),min(squeeze(prctile(imp_stock_2-imp_stock_3,84))))-1 ...
          max(max(squeeze(prctile(imp_stock_2-imp_stock_3,16))),max(squeeze(prctile(imp_stock_2-imp_stock_3,84))))+1])
    set(gca,'Box','on','FontSize',24)
    ylabel('Percent','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    xlabel('Months','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    ttt= title('2002 M1-2007 M12');
    set(ttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    grid on
    set(gca,'GridLineStyle','--')
    ax = gca;
    ax.GridAlpha = .5;
    h = gcf;
    set(h, 'PaperPositionMode','auto')
    set(h,'PaperOrientation','landscape')
    set(h,'Position',[0 0 1500 1000])
    tightfig
    print(h,'-depsc','..\Output\TVP\Figure_17')

    %% Figure 18: Differences House Prices
    
    figure(18)
    subplot(2,2,1)
    plot(median(imp_house_1),'Linewidth',1,'Color',rgb(irf_color),'Marker','o')
    hold on
    plot(median(imp_house_2),'Linewidth',1,'Color',rgb(irf_color),'Marker','x')
    plot(median(imp_house_3),'Linewidth',1,'Color',rgb(irf_color))
    legend('1991 M1','1995 M1','2007 M12','location','Northeast')
    set(gca,'XTick',0:10:IRF_T)
    xlim([1 IRF_T])
    set(gca,'Box','on','FontSize',24)
    ylabel('Percent','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    xlabel('Months','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    ttt= title('House Prices');
    set(ttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    grid on
    set(gca,'GridLineStyle','--')
    ax = gca;
    ax.GridAlpha = .5;
    
    subplot(2,2,2)
    plot(squeeze(prctile(imp_house_1-imp_house_2,16)),'Linewidth',1,'Color',rgb(irf_color));
    hold on
    plot(squeeze(prctile(imp_house_1-imp_house_2,31)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_house_1-imp_house_2,50)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_house_1-imp_house_2,69)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_house_1-imp_house_2,84)),'Linewidth',1,'Color',rgb(irf_color));
    ciplot(squeeze(prctile(imp_house_1-imp_house_2,16)),squeeze(prctile(imp_house_1-imp_house_2,84)),band_color);
    ciplot(squeeze(prctile(imp_house_1-imp_house_2,31)),squeeze(prctile(imp_house_1-imp_house_2,69)),band_color);
    set(gca,'XTick',0:10:IRF_T)
    axis tight
    xlim([1 IRF_T])
    ylim([min(min(squeeze(prctile(imp_house_1-imp_house_2,16))),min(squeeze(prctile(imp_house_1-imp_house_2,84))))-1 ...
          max(max(squeeze(prctile(imp_house_1-imp_house_2,16))),max(squeeze(prctile(imp_house_1-imp_house_2,84))))+1])
    set(gca,'Box','on','FontSize',24)
    ylabel('Percent','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    xlabel('Months','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    ttt= title('1991 M1-1995 M1');
    set(ttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    grid on
    set(gca,'GridLineStyle','--')
    ax = gca;
    ax.GridAlpha = .5;
    
    subplot(2,2,3)
    plot(squeeze(prctile(imp_house_1-imp_house_3,16)),'Linewidth',1,'Color',rgb(irf_color));
    hold on
    plot(squeeze(prctile(imp_house_1-imp_house_3,31)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_house_1-imp_house_3,50)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_house_1-imp_house_3,69)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_house_1-imp_house_3,84)),'Linewidth',1,'Color',rgb(irf_color));
    ciplot(squeeze(prctile(imp_house_1-imp_house_3,16)),squeeze(prctile(imp_house_1-imp_house_3,84)),band_color);
    ciplot(squeeze(prctile(imp_house_1-imp_house_3,31)),squeeze(prctile(imp_house_1-imp_house_3,69)),band_color);
    set(gca,'XTick',0:10:IRF_T)
    axis tight
    xlim([1 IRF_T])
    ylim([min(min(squeeze(prctile(imp_house_1-imp_house_3,16))),min(squeeze(prctile(imp_house_1-imp_house_3,84))))-1 ...
          max(max(squeeze(prctile(imp_house_1-imp_house_3,16))),max(squeeze(prctile(imp_house_1-imp_house_3,84))))+1])
    set(gca,'Box','on','FontSize',24)
    ylabel('Percent','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    xlabel('Months','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    ttt= title('1991 M1-2007 M12');
    set(ttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    grid on
    set(gca,'GridLineStyle','--')
    ax = gca;
    ax.GridAlpha = .5;
    
    subplot(2,2,4)
    plot(squeeze(prctile(imp_house_2-imp_house_3,16)),'Linewidth',1,'Color',rgb(irf_color));
    hold on
    plot(squeeze(prctile(imp_house_2-imp_house_3,31)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_house_2-imp_house_3,50)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_house_2-imp_house_3,69)),'Linewidth',1,'Color',rgb(irf_color));
    plot(squeeze(prctile(imp_house_2-imp_house_3,84)),'Linewidth',1,'Color',rgb(irf_color));
    ciplot(squeeze(prctile(imp_house_2-imp_house_3,16)),squeeze(prctile(imp_house_2-imp_house_3,84)),band_color);
    ciplot(squeeze(prctile(imp_house_2-imp_house_3,31)),squeeze(prctile(imp_house_2-imp_house_3,69)),band_color);
    set(gca,'XTick',0:10:IRF_T)
    axis tight
    xlim([1 IRF_T])
    ylim([min(min(squeeze(prctile(imp_house_2-imp_house_3,16))),min(squeeze(prctile(imp_house_2-imp_house_3,84))))-1 ...
          max(max(squeeze(prctile(imp_house_2-imp_house_3,16))),max(squeeze(prctile(imp_house_2-imp_house_3,84))))+1])
    set(gca,'Box','on','FontSize',24)
    ylabel('Percent','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    xlabel('Months','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    ttt= title('1995 M1-2007 M12');
    set(ttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    grid on
    set(gca,'GridLineStyle','--')
    ax = gca;
    ax.GridAlpha = .5;
    h = gcf;
    set(h, 'PaperPositionMode','auto')
    set(h,'PaperOrientation','landscape')
    set(h,'Position',[0 0 1500 1000])
    tightfig
    print(h,'-depsc','..\Output\TVP\Figure_18')

    %% Preliminaries Sacrifice Ratios
    
    date_sacr_1  = 12*(year_compare_sacr_1-start_year_TVP)+month_compare_sacr_1-start_month_TVP+1;
    date_sacr_2  = 12*(year_compare_sacr_2-start_year_TVP)+month_compare_sacr_2-start_month_TVP+1;
    
    imp_sacr_stock_1(:,:) = cumsum(imp_save_all(:,2,1:40,date_sacr_1),3);
    imp_sacr_stock_2(:,:) = cumsum(imp_save_all(:,2,1:40,date_sacr_2),3);
    imp_sacr_house_1(:,:) = cumsum(imp_save_all(:,4,1:40,date_sacr_1),3);
    imp_sacr_house_2(:,:) = cumsum(imp_save_all(:,4,1:40,date_sacr_2),3);
    
    imp_indus_s_1(:,:) = cumsum(imp_save_all(:,6,1:40,date_sacr_1),3);
    imp_indus_s_2(:,:) = cumsum(imp_save_all(:,6,1:40,date_sacr_2),3);
    imp_indus_h_1(:,:) = cumsum(imp_save_all(:,6,1:40,date_sacr_1),3);
    imp_indus_h_2(:,:) = cumsum(imp_save_all(:,6,1:40,date_sacr_2),3);
    
    imp_stock_1_pos = [];
    imp_stock_2_pos = [];
    imp_house_1_pos = [];
    imp_house_2_pos = [];
    
    imp_indus_s_1_pos = [];
    imp_indus_s_2_pos = [];
    imp_indus_h_1_pos = [];
    imp_indus_h_2_pos = [];
    
    for jj=1:3000
        
        check_s = 0;
        check_h = 0;
        
        for kk=1:40
                
            if imp_sacr_stock_1(jj,kk)>=0 || imp_indus_s_1(jj,kk)>=0 || imp_sacr_stock_2(jj,kk)>=0 || imp_indus_s_2(jj,kk)>=0
                check_s = 1;
            end
            if imp_sacr_house_1(jj,kk)>=0 || imp_indus_h_1(jj,kk)>=0 || imp_sacr_house_2(jj,kk)>=0 || imp_indus_h_2(jj,kk)>=0
                check_h = 1;
            end
                
        end
        
        if check_s == 0
        
            imp_stock_1_pos     = [imp_stock_1_pos;imp_sacr_stock_1(jj,:)];
            imp_stock_2_pos     = [imp_stock_2_pos;imp_sacr_stock_2(jj,:)];
            imp_indus_s_1_pos   = [imp_indus_s_1_pos;imp_indus_s_1(jj,:)];
            imp_indus_s_2_pos   = [imp_indus_s_2_pos;imp_indus_s_2(jj,:)]; 
            
        end
        
        if check_h == 0
            
            imp_house_1_pos     = [imp_house_1_pos;imp_sacr_house_1(jj,:)];
            imp_house_2_pos     = [imp_house_2_pos;imp_sacr_house_2(jj,:)];
            imp_indus_h_1_pos   = [imp_indus_h_1_pos;imp_indus_h_1(jj,:)];
            imp_indus_h_2_pos   = [imp_indus_h_2_pos;imp_indus_h_2(jj,:)];
            
        end
    end
    
    imp_sacr_stock_1 = imp_indus_s_1_pos./imp_stock_1_pos;
    imp_sacr_stock_2 = imp_indus_s_2_pos./imp_stock_2_pos;
    
    imp_sacr_house_1 = imp_indus_h_1_pos./imp_house_1_pos;
    imp_sacr_house_2 = imp_indus_h_2_pos./imp_house_2_pos;
    
    diff_sacr_stock = imp_sacr_stock_2-imp_sacr_stock_1;
    diff_sacr_house = imp_sacr_house_2-imp_sacr_house_1;
    
    diff_sacr_stock_llb   = squeeze(prctile(diff_sacr_stock,16));
    diff_sacr_stock_lb    = squeeze(prctile(diff_sacr_stock,31));
    diff_sacr_stock_me    = squeeze(prctile(diff_sacr_stock,50));
    diff_sacr_stock_ub    = squeeze(prctile(diff_sacr_stock,69));
    diff_sacr_stock_uub   = squeeze(prctile(diff_sacr_stock,84));
    
    min_sacr_stock        = min(min(diff_sacr_stock_uub,diff_sacr_stock_llb));
    max_sacr_stock        = max(max(diff_sacr_stock_uub,diff_sacr_stock_llb));
    
    diff_sacr_house_llb   = squeeze(prctile(diff_sacr_house,16));
    diff_sacr_house_lb    = squeeze(prctile(diff_sacr_house,31));
    diff_sacr_house_me    = squeeze(prctile(diff_sacr_house,50));
    diff_sacr_house_ub    = squeeze(prctile(diff_sacr_house,69));
    diff_sacr_house_uub   = squeeze(prctile(diff_sacr_house,84));

    min_sacr_house        = min(min(diff_sacr_house_uub,diff_sacr_house_llb));
    max_sacr_house        = max(max(diff_sacr_house_uub,diff_sacr_house_llb));
    
    %% Figure 19: Differences Sacrifice Ratios
    
    figure(15)
    subplot(2,2,1)
    plot(median(imp_sacr_stock_1),'Linewidth',1,'Color',rgb(irf_color),'Marker','o')
    hold on
    plot(median(imp_sacr_stock_2),'Linewidth',1,'Color',rgb(irf_color),'Marker','x')
    legend('1995 M1','2007 M12','location','Northeast')
    set(gca,'XTick',0:10:IRF_T)
    xlim([1 IRF_T])
    set(gca,'Box','on','FontSize',24)
    ylabel('Percent','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    xlabel('Months','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    ttt= title('Sacrifice Ratio - $\frac{\mathrm{Output}}{\mathrm{Stock \thinspace Prices}}$');
    set(ttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    grid on
    set(gca,'GridLineStyle','--')
    ax = gca;
    ax.GridAlpha = .5;
    
    subplot(2,2,2)
    plot(median(imp_sacr_house_1),'Linewidth',1,'Color',rgb(irf_color),'Marker','o')
    hold on
    plot(median(imp_sacr_house_2),'Linewidth',1,'Color',rgb(irf_color),'Marker','x')
    legend('1995 M1','2007 M12','location','Northeast')
    set(gca,'XTick',0:10:IRF_T)
    xlim([1 IRF_T])
    set(gca,'Box','on','FontSize',24)
    ylabel('Percent','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    xlabel('Months','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    ttt= title('Sacrifice Ratio - $\frac{\mathrm{Output}}{\mathrm{House \thinspace Prices}}$');
    set(ttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    grid on
    set(gca,'GridLineStyle','--')
    ax = gca;
    ax.GridAlpha = .5;
    
    subplot(2,2,3)
    plot(diff_sacr_stock_llb,'Linewidth',1,'Color',rgb(irf_color));
    hold on
    plot(diff_sacr_stock_lb,'Linewidth',1,'Color',rgb(irf_color));
    plot(diff_sacr_stock_me,'Linewidth',3,'Color',rgb(irf_color));
    plot(diff_sacr_stock_ub,'Linewidth',1,'Color',rgb(irf_color));
    plot(diff_sacr_stock_uub,'Linewidth',1,'Color',rgb(irf_color));
    ciplot(diff_sacr_stock_llb,diff_sacr_stock_uub,band_color);
    ciplot(diff_sacr_stock_lb,diff_sacr_stock_ub,band_color);
    set(gca,'XTick',0:10:IRF_T)
    axis tight
    xlim([1 IRF_T])
    ylim([min_sacr_stock-(max_sacr_stock-min_sacr_stock)/5 max_sacr_stock+(max_sacr_stock-min_sacr_stock)/5])
    ttt= title('2007 M12-1995 M1');
    set(gca,'Box','on','FontSize',24)
    ylabel('Percent','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    xlabel('Months','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    set(ttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    grid on
    set(gca,'GridLineStyle','--')
    ax = gca;
    ax.GridAlpha = .5;
    
    subplot(2,2,4)
    plot(diff_sacr_house_llb,'Linewidth',1,'Color',rgb(irf_color));
    hold on
    plot(diff_sacr_house_lb,'Linewidth',1,'Color',rgb(irf_color));
    plot(diff_sacr_house_me,'Linewidth',3,'Color',rgb(irf_color));
    plot(diff_sacr_house_ub,'Linewidth',1,'Color',rgb(irf_color));
    plot(diff_sacr_house_uub,'Linewidth',1,'Color',rgb(irf_color));
    ciplot(diff_sacr_house_llb,diff_sacr_house_uub,band_color);
    ciplot(diff_sacr_house_lb,diff_sacr_house_ub,band_color);
    set(gca,'XTick',0:10:IRF_T)
    axis tight
    xlim([1 IRF_T])
    ylim([min_sacr_house-(max_sacr_house-min_sacr_house)/5 max_sacr_house+(max_sacr_house-min_sacr_house)/5])
    ttt= title('2007 M12-1995 M1');
    set(gca,'Box','on','FontSize',24)
    ylabel('Percent','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    xlabel('Months','Fontsize',24,'FontName','Palatino','Interpreter','Latex')
    set(ttt,'fontsize',30,'FontName','Palatino','Interpreter','Latex')
    grid on
    set(gca,'GridLineStyle','--')
    ax = gca;
    ax.GridAlpha = .5;
    
    h = gcf;
    set(h, 'PaperPositionMode','auto')
    set(h,'PaperOrientation','landscape')
    set(h,'Position',[0 0 1500 1000])
    tightfig
    print(h,'-depsc','..\Output\TVP\Figure_19')
    
end
