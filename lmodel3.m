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
bsboffset=0;
% number of groups
ngroups=5;

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
%for ind = 1:nsimul

coffebracont=fx_simBRA(xdata);
%coffebracont=xdata;
% nn: numbero of units
nn=size(coffebracont,1);
countrynames=unique(coffebracont(:,2));
% number of units belonging to a group
grpunt=nn/length(countrynames);

% crea matrice d'appoggio
tabcountry=zeros(5,length(countrynames));
%colnames(tabcountry)=countrynames

% numero di subsamples
nsubsamp=50;


% contamination (to be added later!)
% logV[5]=logV[5]+1.5
% logV[6]=logV[6]+1.5
% logV[7]=logV[7]+1.5
% logV[8]=logV[8]+1.5

%%%%% plot(coffebra.cont$logQ, coffebra.cont$logV, cex = 1.5, pch = 20, xlab = "W", ylab = "V")
lme = fitlme(xdata,'logV ~ logQ + t+m1+m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+(1 | country)','FitMethod','ML');



% %%%%% regression in math format
% X = table2array(xdata(:,4:16));
% y =table2array(xdata(:,17));
% %Z = num2cell(ones(size(y)));
% X = [ones(size(y)), X];
% Z = ones(size(y));
% v = [1 2 3 4 5];
% u = repelem(v,60);
% G = u';
% lme2 = fitlmematrix(X,y,Z,G);
 y = xdata.logV;
% %ngroups = max(fit.coffebra.cont@Gp)
% %X = fit.coffebra.cont@pp$X
% %Z = t(fit.coffebra.cont@pp$Zt)

yhat=fitted(lme);
res=y-yhat;

% extract covariance parameters
[psi,mse,stats] = covarianceParameters(lme);
sigma2eps=mse;
sigma2U=cell2mat(psi);

% predict y from subset
%ypred = predict(lme,tblnewsub);

% or use Z instead from lme2!


positions=1:nn;
positions = reshape(positions, [grpunt, ngroups]);
temp = NaN(nsubsamp,p*ngroups);
for l = 1:nsubsamp
    subset=NaN(p+bsboffset, ngroups);
    % set.seed(l)
    sss = datasample(1:(grpunt - (p+bsboffset)+1), ngroups);
    for ii = 1:size(positions,2)
        subset(:,ii) = positions(sss(ii):(sss(ii)+(p+bsboffset)-1), ii);
    end
    temp(l,:) = reshape(subset, [1,p*ngroups]);
end

resorder =NaN(nn, size(temp,1));

for j = 1:size(temp,1)
    if (mod(j,5)==0)
        disp([num2str(j) '\t']);
    end
    mysub = temp(j,:);
    %fitcoffebracontsub = fitlmematrix(X(mysub,:),y(mysub,:),Z(mysub,:),G(mysub,:));
    fitcoffebracontsub = fitlme(xdata(mysub,:),'logV ~ logQ + t+m1+m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+(1 | country)','FitMethod','ML');
    %ngroups = max(fit.coffebra.cont.sub@Gp)
%     ysub = y(mysub,:);
%     Xsub = X(mysub,:);
%     Zsub = Z(mysub,:);
    beta.sub = fitcoffebracontsub.Coefficients;
    [psi,mse,stats] = covarianceParameters(fitcoffebracontsub);
    sigma2U_sub = cell2mat(psi);
    sigma2eps_sub = mse;
    %Ghat_sub = sigma2U_sub * diag(ngroups);
    
    %%% This is to check if beta is obtained
    %Vhat_sub = sigma2U_sub * Zsub * Zsub' + sigma2eps_sub * diag(length(ysub));
    
    %%% uhat.sub = Ghat.sub1%*%t(Z)%*%ginv(as.matrix(Vhat.sub1)) %*% (y - X%*%beta.sub)
    [uhat_sub, uhat_subnames] = randomEffects(fitcoffebracontsub);
    yhat_sub = predict(fitcoffebracontsub,xdata);
    res_sub = y - yhat_sub;
    resorder(:,j) = sort(res_sub.^2);
end

% median maybe best mean!
[B,index]=min(median(resorder, 1));
bsb = temp(:,index);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%        END INIT PART

%%% Storage of quantities in FS
res_fwd=NaN(nn, nn-length(bsb) + 1);
units_in=NaN(nn, nn-length(bsb) + 1);
beta_fwd = NaN(nn, nn-length(bsb) + 1);
%ttest_fwd = beta.fwd
sigma2U_fwd = NaN(1, nn-length(bsb)+ 1);
sigma2eps_fwd = NaN(1, nn-length(bsb)+ 1);
row_inside_fwd = NaN(nn, nn-length(bsb)+ 1);
uhat_fwd = NaN(ngroups, nn-length(bsb)+ 1);
%ores_fwd=res.fwd


fscycle=1:length(bsb)

% original:
first=1;
for m=1:length(bsb)
    
    % tabcountry[2:3,1:5]=0
    
    % cat("Subset size m of the FWD = ", m, "\n")
    if (m == length(bsb))
        %%%% first step of the forward search
        step = 1;
        ttt = bsb;
        
        fitcoffebracontsub = fitlme(xdata(ttt,:),'logV ~ logQ + t+m1+m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+(1 | country)','FitMethod','ML');
        %ngroups = max(fit.coffebra.cont.sub@Gp)
        %ysub = y(ttt,:);
        %Xsub = X(ttt,:);
        %Zsub = Z(ttt,:);
        beta_sub = fitcoffebracontsub.Coefficients(:,2);
        ttest_sub = fitcoffebracontsub.Coefficients(:,6);
        [psi,mse,~] = covarianceParameters(fitcoffebracontsub);
        sigma2U_sub = cell2mat(psi);
        sigma2eps_sub = mse;
        %Ghat_sub = sigma2U_sub * diag(ngroups);
        %loglikhat_sub = fitcoffebracontsub.LogLikelihood;
        [uhat_sub, ~] = randomEffects(fitcoffebracontsub);
        yhat_sub = predict(fitcoffebracontsub, xdata);
        res_sub = y - yhat_sub;
        % posm = which(order(res.sub^2)<= m)
        
        %%% save FS data
        
        res_fwd(:, step) = res_sub;
        beta_fwd(:, step) = beta_sub;
        %ttest_fwd(:, step) = ttest_sub;
        sigma2U_fwd(:, step) = sigma2U_sub;
        sigma2eps_fwd(:, step) = sigma2eps_sub;
        row_inside_fwd[1:length(bsb), step] = bsb;
        uhat_fwd(:, step) = uhat_sub;
    else
        % steps > 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ORDERING RESIDUALS
        
        % creiamo una matrice con 3 colonne: 1a indici delle units, 2a gruppo di apparteneneza, 3a SQRES
        tmpres=as.data.frame(cbind(1:nn,coffebra.cont$country,as.numeric(res.sub^2)))
        names(tmpres)=c("unit","country","res2")
        
        % ordiniamo gli SQRES per gruppo
        otmpres=tmpres[order(tmpres$country, tmpres$res2),]
        otmpres=cbind(otmpres,1:nn)
        names(otmpres)=c("unit","country","res2","position")
        
        % prendiamo gli indici delle units ordinati in senso crescente
        ordbsb=as.numeric(otmpres$unit)
        
        if (first==1)
            % al primo run:
            % dividiamo lo step [m] della FS in gruppi di uguale dimensione
            m1=floor(m/5);
            % ordind cont. indici del tipo: 1, 2, 3, ..., ngroup1, (ngroup1+1), (ngroup1+2), ..., (ngroup1+ngroup2), (ngroup2+1), ..., (ngroup2+ngroup3)
            ordind=[(1:m1)+0*grpunt (1:m1)+1*grpunt (1:m1)+2*grpunt (1:m1)+3*grpunt (1:m1)+4*grpunt];
            % tabcountry row2 contiene il numero delle units assegnate a ciascun gruppo
            tabcountry(2,1:5)=m1;
            % che vengono immesse nelle variabili nc1,  nc2 etc.
            nc1=tabcountry(2,1);
            nc2=tabcountry(2,2);
            nc3=tabcountry(2,3);
            nc4=tabcountry(2,4);
            nc5=tabcountry(2,5);
            
            % ordind viene popolato con gli indici trovati dalla FS
            ordind=[(1:nc1)+0*grpunt (1:nc2)+1*grpunt (1:nc3)+2*grpunt (1:nc4)+3*grpunt (1:nc5)+4*grpunt];
            
        else
            % generic 1+m step
            % nc1 e nc2 cont. il numero di units assegnate a ciascun gruppo
            nc1=tabcountry(2,1);
            nc2=tabcountry(2,2);
            nc3=tabcountry(2,3);
            nc4=tabcountry(2,4);
            nc5=tabcountry(2,5);
            
            % ordind (lungo come la somma delle units assegnate) cont. gli indici delle units assegnate
            ordind = [(1:nc1)+0*grpunt (1:nc2)+1*grpunt (1:nc3)+2*grpunt (1:nc4)+3*grpunt (1:nc5)+4*grpunt];
        end
        
        
        % prendiamo le prime 1:nc1 e 1:nc2 units con il minore SQRES
        % NB: ordind contiene sequenze concatenate tipo 1:5, (1+ng1:8+ng1), (1+(ng1+ng2):8+(ng1+ng2), ...
        
        % bsblast cont. gli ID delle unità (ordinate per SQRES) selezionate per il subset
        bsblast=ordbsb(ordind);
        % selinf cont. gli indici (non gli ID) delle prossime 2 units che competono per entrare nel subset
        
        nc1s=nc1;
        nc2s=nc2;
        nc3s=nc3;
        nc4s=nc4;
        nc5s=nc5;
        
        if (nc1==grpunt)
            nc1=grpunt-1;
            nc1s=NaN;
        end
        if (nc2==grpunt)
            nc2=grpunt-1;
            nc2s=NaN;
        end
        if (nc3==grpunt)
            nc3=grpunt-1;
            nc3s=NaN;
        end
        if (nc4==grpunt)
            nc4=grpunt-1;
            nc4s=NaN;
        end
        if (nc5==grpunt)
            nc5=grpunt-1;
            nc5s=NaN;
        end
        
        selind= [(nc1s+1)+0*grpunt (nc2s+1)+1*grpunt (nc3s+1)+2*grpunt (nc4s+1)+3*grpunt (nc5s+1)+4*grpunt];
        % ordindnew contiene gli indici delle units gia' incluse e delle prossime 2 units che competono per entrare nel subset
        % ordindnew=c((1:nc1+1)+0*grpunt, (1:nc2+1)+1*grpunt, (1:nc3+1)+2*grpunt, (1:nc4+1)+3*grpunt, (1:nc5+1)+4*grpunt)
        
        % mincountry=which.min(otmpres[selind,3]) # not needed?
        
        % idx cont. l'indice della unit da aggiungere al subset con il minore SQRES
        idx=selind[which.min(otmpres[selind,3])];
        % idy cont. ID della unit da aggiungere al subset con il minore SQRES
        % idy=as.numeric(as.vector(otmpres[idx,1]))
        
        nbsbind=sort(c(ordind,idx));
        
        %   id_ctry=which(match(tabcountry[1,],tmpres[ordindx[cv_idx],2])==1)
        id_ctry=which(match(countrynames, otmpres[idx,2])==1)
        cv_idx=otmpres[idx,1]
        
        if (mod(m,100) == 0)
            cat("Progression: ", step, "Monte Carlo loop: ", ind, "\n")
            cat("Subset size m of the FWD = ", m, "\n")
            cat("m:", m, " id_ctry:", countrynames[id_ctry], " cv_idx:", cv_idx,"\n")
        end
        
        % debug code ignore!
        % if (m>=100)
        %   print("break!")
        %   end
        
        % numero di units incluse per ciascun gruppo
        
        tabcountry(2,1)=sum(otmpres(nbsbind,2)==countrynames(1));
        tabcountry(2,2)=sum(otmpres(nbsbind,2)==countrynames(2));
        tabcountry(2,3)=sum(otmpres(nbsbind,2)==countrynames(3));
        tabcountry(2,4)=sum(otmpres(nbsbind,2)==countrynames(4));
        tabcountry(2,5)=sum(otmpres(nbsbind,2)==countrynames(5));
        
        
        % per aggiornamenti futuri del codice
        tabcountry(4,id_ctry)=1;
        % valori delle ultime units incluse per ciascun gruppo
        tabcountry(5,id_ctry)=otmpres(idx,3);
        % id_best=ordindx[cv_idx]
        
        
        % tt contiene gli ID delle units incluse nel subset
        ttt = otmpres(nbsbind,1);
        % units.in matrice che per ogni passo m della FS segnala le units presenti (nel passo corrente)
        units.in(ttt, step) = 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% init.fit = lmer(logV ~ logQ +(1 | country), data = coffebra.cont, REML=F, subset = ttt)
        init.fit = lmer(logV ~ logQ+t+m1+m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+(1 | country), data = coffebra.cont, REML=T, subset = ttt)
        %%%        init.fit = lmer(w ~ t +(1 | pigname), data = pig, REML=F, subset = ttt)
        check.ngroups = max(init.fit@Gp)
        
        % cat("\n", "Progression = ", step, "\n")
        
        y.sub = init.fit@resp$y
        X.sub = init.fit@pp$X
        Z.sub = t(init.fit@pp$Zt)
        beta.sub = init.fit@beta
        ttest.sub = beta.sub/diag(summary(init.fit)$vcov)^0.5
        sigma2U.sub = summary(init.fit)$varcor[[1]][1]
        sigma2eps.sub = summary(init.fit)$sigma^2
        loglikhat.sub = summary(init.fit)$logLik[1]
        Ghat.sub = sigma2U.sub * diag(ngroups)
        
        %%% Ste equal to zero rows unused
        
        not.mysub = setdiff(1:nn, ttt)
        %%%    Z.reduced = Z
        
        
        Vhat.sub1 = sigma2U.sub * Z %*% t(Z) + sigma2eps.sub * diag(nn)
        Ghat.sub1 = sigma2U.sub * diag(ngroups)
        
        
        %%% identical to ranef(ranef(fit.coffebra.cont.sub)$pigname[[1]])
        uhat.sub = Ghat.sub1%*%t(Z)%*%ginv(as.matrix(Vhat.sub1)) %*% (y - X%*%beta.sub)
        %%%        yhat.sub = X%*%beta.sub + Z.reduced%*%uhat.sub
        yhat.sub = X%*%beta.sub + Z%*%uhat.sub
        res.sub = y - yhat.sub
        res.fwd[, step] = as.vector(res.sub)
        beta.fwd[, step] = as.vector(beta.sub)
        ttest.fwd[, step] = as.vector(ttest.sub)
        sigma2U.fwd[1, step] = as.vector(sigma2U.sub)
        sigma2eps.fwd[1, step] = as.vector(sigma2eps.sub)
        %row.inside.fwd[1:length(nbsbind), step] = nbsbind
        row.inside.fwd[1:length(ttt), step] = ttt
        uhat.fwd[, step] = as.vector(uhat.sub)
        
        % verifichiamo l'ordinamento
        ores.fwd[, step]=otmpres$res2
        
        first=0;
    end
    
end


%}