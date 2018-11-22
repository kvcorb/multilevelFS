function out=fx_simBRA(xdata)

% ricordarsi di chiamare mtR per fissare il seed e lo state!



ncountries = length(unique(xdata.country));
N =size(xdata,1)/ncountries;

%a0 = randn(ncountries, 0.8, sd = 0.01);
% a=mean, b=sd that is r = b.*randn(n,1) + a
a0 = 0.01.*randn(ncountries,1) + 0.8;


sdG = 2.58/500;
sdF = 2.65/500;
sdI = 3.33/500;
sdS = 3.29/500;
sdUK = 3.17/500;

WG=sdG .* randn(N,1);
WF=sdF .* randn(N,1);
WI=sdI .* randn(N,1);
WS= sdS .* randn(N,1);
WUK=sdUK .* randn(N,1);

betaf=[0.01, 0.02, -0.02, -0.03, -0.05, 0.06, 0.02, -0.05, -0.01, -0.06, -0.01, 0.04, 0.03];

betaG = [a0(1),betaf];
betaF = [a0(2),betaf];
betaI = [a0(3),betaf];
betaS = [a0(4),betaf];
betaUK = [a0(5),betaf];

% xa=table2array(x(:,[2 4:17]))


xID=1:(N*ncountries);
%xstate=(xdata(:,3));

% con la table
% IdUK=strcmp(xdata{:,3},'UK')


% con i char
% char(xdata{1:10,3})=='G'


% in MATLAB le matrici non possono contenere codici alfanumerici!
% Codifichiamo le nazioni
xstate=[repmat(1,60,1);repmat(2,60,1);repmat(3,60,1);repmat(4,60,1);repmat(5,60,1)];

% ok
intercept=ones(N,1);


% da creare una matrice con 1a col. logQ e le altre con le dummy vars.
IdG=strcmp(xdata{:,3},'G');
TlnPG=table2array([xdata(IdG,16) xdata(IdG,4:15)] );
lnPG=[intercept TlnPG]* betaG' + WG;
IdF=strcmp(xdata{:,3},'F');
TlnPF=table2array([xdata(IdF,16) xdata(IdF,4:15)] );
lnPF=[intercept TlnPF]* betaF' + WF;
IdI=strcmp(xdata{:,3},'I');
TlnPI=table2array([xdata(IdI,16) xdata(IdI,4:15)] );
lnPI=[intercept TlnPI]* betaI' + WI;
IdS=strcmp(xdata{:,3},'S');
TlnPS=table2array([xdata(IdS,16) xdata(IdS,4:15)] );
lnPS=[intercept TlnPS] * betaS' + WS;
IdUK=strcmp(xdata{:,3},'UK');
TlnPUK=table2array([xdata(IdUK,16) xdata(IdUK,4:15)] );
lnPUK=[intercept TlnPUK]* betaUK' + WUK;



x1 =[lnPG; lnPF; lnPI; lnPS; lnPUK];

x1=[x1 xstate table2array(xdata(:,4:16))];

% rowNames = {'a','b','c'};
% colNames = {'x','y','z'};
% sTable = array2table(x1,'RowNames',rowNames,'VariableNames',colNames)

% V=P*Q
logV= x1(:,1).* x1(:,15);

x1=[x1,logV];

%colnames(x1)[16] <- "logV"


scatter(x1(:,15), x1(:,16))


out=x1;
end

