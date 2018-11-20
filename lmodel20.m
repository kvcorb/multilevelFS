%% 
load carsmall

y = MPG;
X = [ones(size(Weight)), Weight];
Z = ones(size(y));
lme = fitlmematrix(X,y,Z,Model_Year)

Z = double([Model_Year==70, Model_Year==76, Model_Year==82]);
lme = fitlmematrix(X,y,Z,[],'Covariancepattern','Isotropic')

% porting del codice R

xdata=readtable("simBRAmatlab.txt");
xdata.Properties.VariableNames = {'id', 'logP',	'country',	't',	'm1',	'm2',...
    'm3',	'm4',	'm5',	'm6',	'm7',	'm8',	'm9',	'm10',	'm11',	'logQ',	'logV'};

% componenti di \beta
p=14;
% dimensione (numero dei passi) della FS
dimFS=181;
% nsimul
nsimul=1;
% bsb offset to avoid rank deficient error
bsboffset=10;


%%%%%%%%%%%%%% MC matrici 3D

BV3D=NaN(p,476,nsimul);
S23D=NaN(1,476,nsimul);
mdrB3D=NaN(1,475,nsimul);
mdrF3D=NaN(1,484,nsimul);
mdrB3Di=NaN(1,484,nsimul);
mdrF3Di=NaN(1,484,nsimul);

% \beta
BV3D=NaN(p,dimFS,nsimul);
% random effects \u_i
U3D=NaN(5,dimFS,nsimul);
% squared residuals
RESSQ3D=NaN(300,dimFS,nsimul);




%%%%%% start MONTE CARLO LOOP
for ind = 1:nsimul

%coffebra.cont=fx_simBRA(xdata);
coffebra.cont=xdata;
% nn: numbero of units
nn=size(coffebra.cont,1);
countrynames=unique(coffebra.cont.country);
% number of units belonging to a group
grpunt=nn/length(countrynames);

% crea matrice d'appoggio
tabcountry=zeros(5,length(countrynames));
%colnames(tabcountry)=countrynames

% numero di subsamples
nsubsamp=50;


% contamination
% logV[5]=logV[5]+1.5
% logV[6]=logV[6]+1.5
% logV[7]=logV[7]+1.5
% logV[8]=logV[8]+1.5

%%%%% plot(coffebra.cont$logQ, coffebra.cont$logV, cex = 1.5, pch = 20, xlab = "W", ylab = "V")
lme = fitlme(xdata,'logV ~ logQ + t+m1+m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+(1 | country)','FitMethod','ML')


y=xdata.logV;
X=[xdata.logQ, xdata.t, xdata.m1, ];

Z = double([strcmp(xdata.country,'G'), strcmp(xdata.country,'F'), strcmp(xdata.country,'I'), strcmp(xdata.country,'S'), strcmp(xdata.country,'UK')]);

%Z = double([Model_Year==70, Model_Year==76, Model_Year==82]);

lme2 = fitlmematrix(X,y,Z,G,'FixedEffectPredictors', ...
{'logQ',  't', 'm1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11'}, ...
'RandomEffectPredictors',{{'Intercept'}},'RandomEffectGroups',{'Country'});


y = coffebra.cont.logV;
%ngroups = max(fit.coffebra.cont@Gp)
%X = fit.coffebra.cont@pp$X
%Z = t(fit.coffebra.cont@pp$Zt)

yhat=fitted(lme);
res=y-yhat;

% extract covariance parameters
[psi,mse,stats] = covarianceParameters(lme);
sigma2eps=mse;
sigma2U=cell2mat(psi);

% predict y from subset
ypred = predict(lme,tblnewsub);

end
