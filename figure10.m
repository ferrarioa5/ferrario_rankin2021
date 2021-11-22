

clear, close all

addpath('..')

npoints=1e2;

% b=inh strength, TI=inh time scale, DE=inh delay, TD=tone duration, PR=presentation rate
TI=0.1; TH=0.5; DE=0.02; TD=0.03; PR=15; TR=1/PR; a=0.3; b=2;

% LIM=[2.999,3.001];
LIM=[b+TH,3];
YLIM=[0,1];

% p=+, m=-
N_m=@(TR)exp(-(TR-DE-TD)/TI); N_t=@(TR)exp(-(TR-2*DE)/TI); N_p=@(TR)exp(-(TR-DE)/TI);
M_m=@(TR)exp(-(2*TR-DE-TD)/TI); M_p=@(TR)exp(-(2*TR-TD)/TI);

states = struct;

%% MAIN states 
states.I=@(a,b,c,d,TR) (-b.*N_m(TR)+d-TH>=0).*(a-b+d-TH>=0);
states.IS=@(a,b,c,d,TR) (-b.*N_t(TR)+d-TH>=0).*(a-b+d-TH<0);
states.ID=@(a,b,c,d,TR) (-b.*N_m(TR)+d-TH<0).*(a-b+d-TH>=0);
states.IDS=@(a,b,c,d,TR) (-b.*N_t(TR)+d-TH<0).*(a-b.*N_t(TR)+d-TH>=0).*(a-b+d-TH<0);
states.AS=@(a,b,c,d,TR) (a-b.*N_p(TR)+d-TH<0).*(-b.*M_m(TR)+d-TH>=0);
states.ASD=@(a,b,c,d,TR) (a-b.*N_p(TR)+d-TH<0).*(-b.*M_m(TR)+d-TH<0).*(a-b.*M_m(TR)+d-TH>=0);
states.AP=@(a,b,c,d,TR) (a-b.*M_p(TR)+d-TH<0);

%% CONNECT states
states.AScINT=@(a,b,c,d,TR) (a-b.*N_p(TR)+d-TH>=0).*(a-b.*N_t(TR)+d-TH<0);
states.APcAS=@(a,b,c,d,TR) (a-b.*M_p(TR)+d-TH>=0).*(a-b.*M_m(TR)+d-TH<0);


%% vary c,DF -- plotting
font_size=14; kmax=4;
names=fieldnames(states);
nstates=length(names);

figure('Position', [50 50 900 700])
% subplot(5,3,[3,6]); hold on
%subplot(5,3,[10,13]); hold on

dLIM=(LIM(2)-LIM(1))/npoints;
cs=LIM(1):dLIM:LIM(2);
DFs=YLIM(1):dLIM:YLIM(2);
[X,Y]=meshgrid(cs,DFs); 
S=zeros(length(DFs),length(cs),nstates); 

cs_tmp=(LIM(1)-dLIM):dLIM:(LIM(2)+dLIM);
DFs_tmp=(YLIM(1)-dLIM):dLIM:(YLIM(2)+dLIM);
[X_tmp,Y_tmp]=meshgrid(cs_tmp,DFs_tmp);

text_xs=zeros(nstates,1); text_ys=zeros(nstates,1);

figure(2)
%% generate structure containing periods + changing names of the states (when necessary) + plot states
periods = struct('name',{'None'},'period',{0});
count=1;
for i=1:nstates
    name_tmp=names{i};
    f_tmp=states.(name_tmp);
    
    v_undscr=strfind(name_tmp,'_'); n_undscr=length(v_undscr);
    if n_undscr==2 % W/K state
        m=str2num(name_tmp(v_undscr(1)+1:v_undscr(2)-1));
        k=str2num(name_tmp(v_undscr(2)+1:end));
        if contains(name_tmp,'W')
            periods(end+1)=struct('name',{name_tmp},'period',{2*(m+2*k)*TR});
        elseif contains(name_tmp,'K')
            periods(end+1)=struct('name',{name_tmp},'period',{(m+2*k)*TR});
        end
        name_tmp=strcat(name_tmp(1:v_undscr(1)),'{',num2str(m),',',num2str(k),'}'); % in this case also change name to "_{...}"
    elseif n_undscr==1 % cycle skipping state
        k=str2num(name_tmp(v_undscr(1)+1:end));
        if contains(name_tmp,'INT')
            periods(end+1)=struct('name',{name_tmp},'period',{2*(2*k-1)*TR});
        elseif contains(name_tmp,'S')
            periods(end+1)=struct('name',{name_tmp},'period',{2*k*TR});
        end
        name_tmp=strcat(name_tmp(1:v_undscr(1)),'{',num2str(k),'}'); % in this case also change name to "_{...}"
        
        % rename states to make it clearer when plotting
        name_tmp=switch_name(name_tmp(1:v_undscr(1)-1),name_tmp(v_undscr(1):end));
        
        if k>kmax
            name_tmp="";
        end
        
    else % 2TR-periodic state   
        periods(end+1)=struct('name',{name_tmp},'period',{2*TR});
    end
    
    
%     if ~isempty(str2num(name_tmp(end))) && ~isempty(str2num(name_tmp(end-1)))
%         name_tmp=strcat(name_tmp(1:end-2),'{',name_tmp(end-1),',',name_tmp(end),'}');
%     end
    
    % rewrite K and W states c-d-g-cd-.... and convert them to ^1,^2,...
    lower_v=isstrprop(name_tmp,'lower');  % vector of lower cases 
    first_idx=find(isstrprop(name_tmp,'lower')==0,1); % first upper case index
    if first_idx>1
        lower_str=name_tmp(1:first_idx-1);
        name_tmp=strcat(name_tmp(first_idx:end),'^{',name_tmp(1:first_idx-1),'}');
    end
    
    z_value=periods(end).period/(2*TR);
%     S(:,:,i)=f_tmp(a,b,X,X.*Y,TR)*z_value;
    S(:,:,i)=f_tmp(a,b,X,X.*Y,TR);
    

    
    Z_tmp=zeros(length(DFs)+2,length(cs)+2);
    if sum(S(:,:,i))==0
        disp(strcat('Conditions for ', name_tmp, ' not met'))
    else
        disp(strcat('Conditions for ', name_tmp, ' are met.. plotting'))
        Z_tmp(2:end-1,2:end-1)=count*S(:,:,i);
%         Z_tmp(2:end-1,2:end-1)=S(:,:,i);
        surf(X_tmp,Y_tmp,Z_tmp)
        alpha(.5)
        C=contour(X_tmp,Y_tmp,Z_tmp,[0.5 0.5],'k');
        idx=Z_tmp>0.5; xval=X_tmp(idx); yval=Y_tmp(idx);
        
        xmin=min(xval); xmax=max(xval); dx=xmax-xmin;
        ymin=min(yval); ymax=max(yval); dy=ymax-ymin;
        
        if ~isempty(xmin+dx/2) && ~isempty(ymin+dy/2)
            text_xs(i)=xmin+dx/2; text_ys(i)=ymin+dy/2;
        end
        count=count+1;
    end
    
    if text_xs(i)~=0 && text_ys(i)~=0
        text(text_xs(i),text_ys(i),count+4,name_tmp,'HorizontalAlignment','center','Color','k','fontsize',font_size,'fontweight','bold')
%         text(text_xs(i),text_ys(i),z_value+5,name_tmp,'HorizontalAlignment','center','Color','k','fontsize',font_size,'fontweight','bold')
    end
end

% plot(LIM,[0.8,0.8],'k--','LineWidth',1.5)

colormap('jet')
% colbar=colorbar('horiz'); 
shading interp
view(2)
xlim(LIM)
ylim(YLIM)
xlabel('C','fontsize',font_size)
ylabel('DF','fontsize',font_size)

yticks([0,1])
xticks(LIM)
set(gca,'FontSize',font_size)







%% figure2b
p.b=b; p.taui=TI; p.theta=TH; p.delta=DE; p.TD=TD; p.tau=0.001;

tstart=1/PR;
tstop=7/PR-0.001;

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

options = ddeset('RelTol',1e-7,'AbsTol',1e-7,'Jumps',jumps);

names={'I','IS','ID','IDS','AS','ASD','AScINT','APcAS','AP'}; 
nnames=length(names);

q=[...
    3, 0.9; ... % I
    3, 0.715; ... % INTS
    3, 0.68; ... % INTD
    3, 0.625; ... % INTDS
    3, 0.47; ... % AS 
    3, 0.4; ... % ASD
    3, 0.53; ... % AScINT
    3, 0.33; ... % APcAS
    3, 0.2 ... % AP
    ];


out=cell(nnames,1);
for i=1:nnames
    p.c = q(i,1); p.d = p.c*q(i,2);
    if i==3
        p.a=0.6;
    else
        p.a=0.3;
    end
    
    if i==2 || i==9
        p.tau=0.0005;
    elseif i==4 || i==6
        p.tau=0.002;
    elseif i==3
        p.tau=0.004;
    else
        p.tau=0.0005;
    end
    
    sol = dde23(@(t,y,Z) dydt(t,y,Z,p), lags, @(t) ddehist(t), tspan, options);
    t=sol.x; y=sol.y; % z=sol.yp;
    out{i}=[t;y];
end



%% plotting
make_colors
font_size=15;

figure(1)
% tight_margin_plot
xlimits=[tstop-2/PR,tstop]-0.005;
lw=2;

for i=1:nnames
%     subplot(3,4,i+ceil(i/3)); hold on
    subplot(5,3,i+(ceil(i/2)-1)); hold on

%     subplot(3,4,i+ceil(i/3)-1); hold on
%     subplot(3,3,i); hold on
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
    
    if (i>=7)
        xticks(jumps)
        xticklabels(labels)
    else
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    end
    
    if(i==9)
        legend('u_A','u_B','s_A','s_B','FontSize',font_size,'Orientation','horizontal')
    end
    yticks([0,1])
    xlim(xlimits)
    title(names{i},'fontsize',font_size)
    ylim([-0.1,1.1])
    set(gca,'FontSize',font_size)
end





%% PR vs df
figure(3), hold on 
% TI=0.1; TH=0.5; DE=0.02; TD=0.03; PR=15; TR=1/PR; a=0.8; b=2;
c=5; m=6; a=1; TI=0.2; DE=0.01; TH=0.5; b=2;
lw=1;

% p=+, m=-
N_m=@(TR)exp(-(TR-DE-TD)/TI); N_t=@(TR)exp(-(TR-2*DE)/TI); N_p=@(TR)exp(-(TR-DE)/TI);
M_m=@(TR)exp(-(2*TR-DE-TD)/TI); M_p=@(TR)exp(-(2*TR-TD)/TI);

states = struct;

%% MAIN states 
states.INT=@(a,b,c,d,TR) (-b.*N_m(TR)+d-TH>=0).*(a-b+d-TH>=0);
states.INTS=@(a,b,c,d,TR) (-b.*N_t(TR)+d-TH>=0).*(a-b+d-TH<0);
states.INTD=@(a,b,c,d,TR) (-b.*N_m(TR)+d-TH<0).*(a-b+d-TH>=0);
states.INTDS=@(a,b,c,d,TR) (-b.*N_t(TR)+d-TH<0).*(a-b.*N_t(TR)+d-TH>=0).*(a-b+d-TH<0);
states.AS=@(a,b,c,d,TR) (a-b.*N_p(TR)+d-TH<0).*(-b.*M_m(TR)+d-TH>=0);
states.ASD=@(a,b,c,d,TR) (a-b.*N_p(TR)+d-TH<0).*(-b.*M_m(TR)+d-TH<0).*(a-b.*M_m(TR)+d-TH>=0);
states.AP=@(a,b,c,d,TR) (a-b.*M_p(TR)+d-TH<0);

%% CONNECT states
states.AScINT=@(a,b,c,d,TR) (a-b.*N_p(TR)+d-TH>=0).*(a-b.*N_t(TR)+d-TH<0);
states.APcAS=@(a,b,c,d,TR) (a-b.*M_p(TR)+d-TH>=0).*(a-b.*M_m(TR)+d-TH<0);


% TI=0.15;

% subplot(5,3,[9,12,15]); hold on
% subplot(5,3,[11,12,14,15]); hold on

font_size=14; kmax=4;
names=fieldnames(states);
nstates=length(names);

LIM=[3,25];
YLIM=[0,1];

dTR=0.1; TRs=LIM(1):dTR:LIM(2);
ddf=0.001; dfs=YLIM(1):ddf:YLIM(2);

[X,Y]=meshgrid(TRs,dfs);
S=zeros(length(dfs),length(TRs),nstates); 

count=1;
for i=1:nstates
    name_tmp=names{i};
    f_tmp=states.(name_tmp);
    
    S(:,:,i)=i*f_tmp(a,b,c,c*(1-Y.^(1/m)),1./X);

    if sum(S(:,:,i))==0
        disp(strcat('Conditions for ', name_tmp, ' not met'))
    else
        disp(strcat('Conditions for ', name_tmp, ' are met.. plotting'))
%         surf(X,Y,S(:,:,i))
        alpha(.5)
        C=contour(X,Y,S(:,:,i),(i-1)+[0.5,0.5],'k','LineWidth',lw);
        count=count+1;
    end
end

xlabel('PR','fontsize',font_size)
ylabel('df','fontsize',font_size)
% colormap('jet')
% colbar=colorbar('horiz'); 
shading interp
view(2)
xlim(LIM)
ylim(YLIM)
set(gca,'FontSize',font_size)

% yticks(YLIM(1):0.2:YLIM(2))
xticks(LIM(1):10:LIM(2))






















