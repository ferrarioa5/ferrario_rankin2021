

clear, close all

% b=inh strength, TI=inh time scale, DE=inh delay, TD=tone duration, PR=presentation rate
TI=0.4; TH=0.5; DE=0.015; TD=0.005; PR=5; TR=1/PR; a=0.4; b=3;
% LIM=[1.6,2.3];
% YLIM=[0.21,0.7];
LIM=[0,3];
YLIM=[0,1];
npoints=1e3;

% p=+, m=- 
N_p=@(TR)exp(-(TR-DE)/TI); N_m=@(TR)exp(-(TR-DE-TD)/TI); 
M_p=@(TR)exp(-(2*TR-DE)/TI); M_m=@(TR)exp(-(2*TR-DE-TD)/TI); 
L_p=@(TR,k)exp(-(k*TR-DE)/TI); L_m=@(TR,k)exp(-(k*TR-DE-TD)/TI); 
R_p=@(TR,k)exp(-(k*TR+TD-2*DE)/TI); R_m=@(TR,k)exp(-(k*TR-2*DE)/TI); 


% useful quantities
C1=@(a,b,c,d,TR) d-TH;
C2_p=@(a,b,c,d,TR) a-b.*M_p(TR)+d-TH; C2_m=@(a,b,c,d,TR) a-b.*M_m(TR)+d-TH;
C3_p=@(a,b,c,d,TR) -b.*N_p(TR)+c-TH; C3_m=@(a,b,c,d,TR) -b.*N_m(TR)+c-TH;
C4_p=@(a,b,c,d,TR) -b.*M_p(TR)+c-TH; C4_m=@(a,b,c,d,TR) -b.*M_m(TR)+c-TH;
C5_p=@(a,b,c,d,TR) a-b.*N_p(TR)+d-TH; C5_m=@(a,b,c,d,TR) a-b.*N_m(TR)+d-TH;
C6_p=@(a,b,c,d,TR) a-b.*N_p(TR)+c-TH; C6_m=@(a,b,c,d,TR) a-b.*N_m(TR)+c-TH;
C7_p=@(a,b,c,d,TR) -b.*N_p(TR)+d-TH; C7_m=@(a,b,c,d,TR) -b.*N_m(TR)+d-TH;
C8_p=@(a,b,c,d,TR) -b.*M_p(TR)+d-TH; C8_m=@(a,b,c,d,TR) -b.*M_m(TR)+d-TH;
C9=@(a,b,c,d,TR) a-b.*M_p(TR)-TH;
C10=@(a,b,c,d,TR) a-b.*N_p(TR)-TH;

states = struct;


%% 2TR-periodic states
% % MAIN states 
states.SB=@(a,b,c,d,TR) (C3_p(a,b,c,d,TR)<0).*(C8_m(a,b,c,d,TR)>=0).*(C9(a,b,c,d,TR)<0); % SAB (A/B)
states.SD=@(a,b,c,d,TR) (C4_m(a,b,c,d,TR)>=0).*(C2_m(a,b,c,d,TR)>=0).*(C3_p(a,b,c,d,TR)<0).*(C8_m(a,b,c,d,TR)<0).*(C9(a,b,c,d,TR)<0); % SABD (A/B)
states.I=@(a,b,c,d,TR) (C1(a,b,c,d,TR)>=0).*(C6_p(a,b,c,d,TR)<0); % INTA/B

% states.S=@(a,b,c,d,TR) (C1(a,b,c,d,TR)<0).*(C2_p(a,b,c,d,TR)<0).*(C3_p(a,b,c,d,TR)<0); % SA/B
% states.AP=@(a,b,c,d,TR) (C3_m(a,b,c,d,TR)>=0).*(C2_p(a,b,c,d,TR)<0); % AP
% states.AS=@(a,b,c,d,TR) (C3_m(a,b,c,d,TR)>=0).*(C8_m(a,b,c,d,TR)>=0).*(C5_p(a,b,c,d,TR)<0).*(C10(a,b,c,d,TR)<0); % AS (A/B)
% states.ASD=@(a,b,c,d,TR) (C2_m(a,b,c,d,TR)>=0).*(C3_m(a,b,c,d,TR)>=0).*(C5_p(a,b,c,d,TR)<0).*(C8_m(a,b,c,d,TR)<0).*(C10(a,b,c,d,TR)<0); % ASD (A/B)
% states.ID=@(a,b,c,d,TR) (C3_m(a,b,c,d,TR)>=0).*(C5_m(a,b,c,d,TR)>=0).*(C7_m(a,b,c,d,TR)<0).*(C10(a,b,c,d,TR)<0); % INTD
% states.IB=@(a,b,c,d,TR) (C7_m(a,b,c,d,TR)>=0).*(C10(a,b,c,d,TR)<0); % INT



%% vary c,DF -- plotting
font_size=18;
figure('Position', [50 50 900 700])

names=fieldnames(states);
nstates=length(names);

dLIM=(LIM(2)-LIM(1))/npoints;
cs=LIM(1):dLIM:LIM(2);
DFs=YLIM(1):dLIM:YLIM(2);
[X,Y]=meshgrid(cs,DFs); 
S=zeros(length(DFs),length(cs),nstates); 

cs_tmp=(LIM(1)-dLIM):dLIM:(LIM(2)+dLIM);
DFs_tmp=(YLIM(1)-dLIM):dLIM:(YLIM(2)+dLIM);
[X_tmp,Y_tmp]=meshgrid(cs_tmp,DFs_tmp);

figure(1), hold on 

text_xs=zeros(nstates,1); text_ys=zeros(nstates,1);

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
        text_xs(i)=xmin+dx/2; text_ys(i)=ymin+dy/2;
        text(text_xs(i),text_ys(i),nstates+1,name_tmp,'HorizontalAlignment','center','Color','w','fontsize',font_size,'fontweight','bold')
        %         text(text_xs(i),text_ys(i),z_value+2,name_tmp,'HorizontalAlignment','center','Color','k','fontsize',font_size,'fontweight','bold')
    end
end

colormap('hot')
caxis([0,5])
shading interp
view(2)
xlim(LIM)
ylim(YLIM)
xlabel('C','fontsize',font_size)
ylabel('DF','fontsize',font_size)

% xticks([0.5,1,1.5,2,2.5,3])
% yticks([0,1])
set(gca,'FontSize',font_size)































