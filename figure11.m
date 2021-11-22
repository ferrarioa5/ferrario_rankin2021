

load case1

p.a=2; p.b=2.8; p.c=5.5;
p.delta=0.015; p.TD=0.022; p.tau=0.01; p.taui=0.25; p.theta=0.5; p.m=6;
sigma=20;


%% utility functions
Si=@(x) 1./(1+exp(sigma*(-x)));
H=@(z) z>=0;

ddehist=@(t) [1,0,1,0];

numDF=12*5;
numPR=12*5;

rtol=1e-7; atol=1e-7;

dfs=linspace(0,1,numDF);
PRs=linspace(0.5,20,numPR);
Acr=zeros(numPR,numDF);
Bcr=zeros(numPR,numDF);

tic
parfor i=1:numDF
    i
    for j=1:numPR
        
        df=dfs(i)
        PR=PRs(j)
        TR=1/PR; delta=pi*PR;
        
        del=0.001;
        tstart=0; 
        tstop=70*TR-(TR-p.TD)/2;
        tspan=[0,tstop];
        
        FA=@(t,p,df)(Si(sin(delta*t)).*Si(sin(delta*(p.TD-t)))+(1-df.^(1/p.m)).*Si(-sin(delta*t)).*Si(-sin(delta*(p.TD-t)))).*H(t-tstart);
        FB=@(t,p,df)(Si(-sin(delta*t)).*Si(-sin(delta*(p.TD-t)))+(1-df.^(1/p.m)).*Si(sin(delta*t)).*Si(sin(delta*(p.TD-t)))).*H(t-tstart);
        
        dydt=@(t,y,Z,p,df) [...
            (-y(1)+Si(p.a*y(2)-p.b*Z(4,1)+p.c*FA(t,p,df)-p.theta))/p.tau;...
            (-y(2)+Si(p.a*y(1)-p.b*Z(3,1)+p.c*FB(t,p,df)-p.theta))/p.tau;...
            Si(y(1)-p.theta)*(1-y(3))/p.tau-y(3)/p.taui;...
            Si(y(2)-p.theta)*(1-y(4))/p.tau-y(4)/p.taui;...
            ];
        
        njumps=floor(PR*tstop);
        jumps=zeros(1,3*njumps);
        labels=strings(1,3*njumps);
        
        k=1;
        while(k<njumps+1)
            if (mod(k,2))
                jumps(3*k-2)=k/PR; labels(3*k-2)="0";
                jumps(3*k-1)=k/PR+p.delta; labels(3*k-1)="D";
                jumps(3*k)=k/PR+p.TD; labels(3*k)="TD";
            else
                jumps(3*k-2)=k/PR; labels(3*k-2)="TR";
                jumps(3*k-1)=k/PR+p.delta; labels(3*k-1)="TR+D";
                jumps(3*k)=k/PR+p.TD; labels(3*k)="TR+TD";
            end
            k=k+1;
        end
        
        [~,idx]=sort(jumps);
        jumps=jumps(idx);
        labels=labels(idx);
        
        lags=[p.delta];
        options = ddeset('RelTol',rtol,'AbsTol',atol,'Jumps',jumps);
        
        der=@(t,y,Z,p) dydt(t,y,Z,p,df);
        
        sol = dde23(@(t,y,Z) der(t,y,Z,p), lags, @(t) ddehist(t), tspan, options);
        t=sol.x; y=sol.y; % z=sol.yp;
        
        xlimits=[tstop-2*TR-2*TR,tstop-2*TR];
        
        idx_time1=find(t(1:end-1)<xlimits(1) & t(1,2:end)>=xlimits(1));
        idx_time2=find(t(1:end-1)<xlimits(2) & t(1,2:end)>=xlimits(2));
        
        t=t(idx_time1:idx_time2); y=y(:,idx_time1:idx_time2);
        
        idxA=find(y(1,1:end-1)<p.theta & y(1,2:end)>=p.theta);
        idxB=find(y(2,1:end-1)<p.theta & y(2,2:end)>=p.theta);
        Acr(j,i)=length(idxA);
        Bcr(j,i)=length(idxB);
    end
end
toc

Z=Acr+Bcr;
h=surf(PRs,dfs,Z');
xlim([PRs(1),PRs(end)])
ylim([dfs(1),dfs(end)])
set(h,'edgecolor','none')
xlabel('PR')
ylabel('df')
view(2);



