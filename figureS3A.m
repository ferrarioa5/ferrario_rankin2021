

% clear, close all

PR=17; p.b=2; p.taui=0.2; p.theta=0.5; p.delta=0.03; p.TD=0.025; p.tau=0.001; p.a=0.6;

TR=1/PR;
tstart=1/PR;
tstop=7/PR-0.005;

Si=@(z) z>=0;
ddehist=@(t) [1,0,1,0];

t=0:0.00001:tstop;

FA=@(t,p) Si(t-tstart).* ( p.c*Si(mod(t,2*TR)-TR).*Si(TR+p.TD-mod(t,2*TR))+p.d*Si(mod(t-TR,2*TR)-TR).*Si(TR+p.TD-mod(t-TR,2*TR)) );
FB=@(t,p) Si(t-tstart).* ( p.c*Si(mod(t-TR,2*TR)-TR).*Si(TR+p.TD-mod(t-TR,2*TR))+p.d*Si(mod(t,2*TR)-TR).*Si(TR+p.TD-mod(t,2*TR)) );

dydt=@(t,y,Z,p) [...
    (-y(1)+Si(p.a*y(2)-p.b*Z(4,1)+FA(t,p)-p.theta))/p.tau;...
    (-y(2)+Si(p.a*y(1)-p.b*Z(3,1)+FB(t,p)-p.theta))/p.tau;...
    Si(y(1)-p.theta)*(1-y(3))/p.tau-y(3)/p.taui;...
    Si(y(2)-p.theta)*(1-y(4))/p.tau-y(4)/p.taui;...
    ];


lags=[p.delta];
tspan=[0,tstop];

q=[...
    2.9, 0.44 ... % APcAS
    ];

njumps=floor(PR*tstop);
jumps=zeros(1,3*njumps);
labels=strings(1,3*njumps);

k=1;
while(k<njumps+1)
    if (mod(k,2))
        jumps(4*k-3)=p.taui*log((p.a+q(1)*q(2)-p.theta)/p.b)+(2/PR-p.delta)+k/PR; labels(4*k-3)="t*";
        jumps(4*k-2)=k/PR; labels(4*k-2)="0";
        jumps(4*k-1)=k/PR+p.delta; labels(4*k-1)="D";
        jumps(4*k)=k/PR+p.TD; labels(4*k)="TD";
    else
        jumps(4*k-3)=p.taui*log((p.a+q(1)*q(2)-p.theta)/p.b)+(2/PR-p.delta)+k/PR; labels(4*k-3)="TR+t*";
        jumps(4*k-2)=k/PR; labels(4*k-2)="TR";
        jumps(4*k-1)=k/PR+p.delta; labels(4*k-1)="TR+D";
        jumps(4*k)=k/PR+p.TD; labels(4*k)="TR+TD";
    end
    k=k+1;
end

[~,idx]=sort(jumps);
jumps=jumps(idx);
labels=labels(idx);

options = ddeset('RelTol',1e-7,'AbsTol',1e-7,'Jumps',jumps);


names={'APcAS'}; 
nnames=length(names);

out=cell(nnames,1);
for i=1:nnames
    p.c = q(i,1); p.d = p.c*q(i,2);    
    sol = dde23(@(t,y,Z) dydt(t,y,Z,p), lags, @(t) ddehist(t), tspan, options);
    t=sol.x; y=sol.y; % z=sol.yp;
    out{i}=[t;y];
end



%% plotting
make_colors
font_size=15;

% figure(1)
% tight_margin_plot
xlimits=[tstop-2/PR,tstop];
lw=2.5;

for i=1:nnames
%     subplot(5,2,i); 
    hold on
    tmp=out{i};
    t=tmp(1,:); y=tmp(2:end,:);
    plot(t,y(1,:),'-','Color',dark_blue,'LineWidth',lw)
    plot(t,y(2,:),'-','Color',inred,'LineWidth',lw)
    plot(t+p.delta,y(3,:),'-','Color',light_blue,'LineWidth',lw)
    plot(t+p.delta,y(4,:),'-','Color',orange,'LineWidth',lw)
    for element=jumps
        q=plot([element,element],[-0.1,1.1],'LineStyle','-','Color','k','LineWidth',1);
        q.Color(4) = 0.25;
    end
    
    if (i>0)
        xticks(jumps)
        xticklabels(labels)
%         legend('u_A','u_B','s_A','s_B','FontSize',font_size)
    else
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    end
    yticks([0,1])
    xlim(xlimits)
    title(names{i},'fontsize',font_size)
    ylim([-0.1,1.1])
    set(gca,'FontSize',font_size)
end












                
                
                
                
                
                
                
                
                
                
                
