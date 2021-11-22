
clear 
close all

addpath('../')
make_colors;
font_size=14;

s1=0.7;
s2=0.2;
a=1;
idx=1;
maxhaed=0;
s=0.05;

norma=15;
n=11;

uA=linspace(-0.03,1.03,n);
uB=linspace(-0.03,1.03,n);


u=linspace(min(uA),max(uA),17);

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
tstop=15; dt=0.01; tspan=0:dt:tstop;


%% sigmoid 
hh1=subplot(1,2,1);
hold on
H=@(x) x>=0;
rhs=@(t,y)[-y(1)+H(a*(y(2)-s2)); -y(2)+H(a*(y(1)-s1))];

for x=u
    [~,y] = ode45(rhs, tspan, [x,min(uB)], options);
    plot(y(:,1),y(:,2),'Color',return_color(y))
    v=rhs(0,[y(idx,1),y(idx,2)]); norma=sqrt(v(1)^2+v(2)^2);



    
    [~,y] = ode45(rhs, tspan, [x,max(uB)], options);
    plot(y(:,1),y(:,2),'Color',return_color(y))
    v=rhs(0,[y(idx,1),y(idx,2)]); norma=sqrt(v(1)^2+v(2)^2);



    
    [~,y] = ode45(rhs, tspan, [min(uA),x], options);
    plot(y(:,1),y(:,2),'Color',return_color(y))
    v=rhs(0,[y(idx,1),y(idx,2)]); norma=sqrt(v(1)^2+v(2)^2);



    
    [~,y] = ode45(rhs, tspan, [max(uA),x], options);
    plot(y(:,1),y(:,2),'Color',return_color(y))
    v=rhs(0,[y(idx,1),y(idx,2)]); norma=sqrt(v(1)^2+v(2)^2);



end

plot([s1,s1],[min(uA),max(uA)],'r','LineWidth',1) % uB nullcline
plot([min(uB),max(uB)],[s2,s2],'b','LineWidth',1) % uA nullcline
plot([0,0],[0,0],'.','Color','k','MarkerSize',25)
plot([1,1],[1,1],'.','Color','k','MarkerSize',25)

xticks([min(uA),s1,max(uA)])
xticklabels({num2str(min(uA)),'s_1',num2str(max(uA))})
yticks([min(uB),s2,max(uB)])
yticklabels({num2str(min(uB)),'s_2',num2str(max(uB))})

sep_left=@(x) (x-1)*s2/(s1-1); % left separatrix
xx=[min(uA),s1];
plot(xx,[sep_left(xx(1)),sep_left(xx(2))],'Color',yellow,'LineWidth',2)
right_left=@(x) x*(s2-1)/s1+1; % right separatrix
xx=[s1,max(uA)];
plot(xx,[right_left(xx(1)),right_left(xx(2))],'Color',orange,'LineWidth',2)

plot([s1,s1],[s2,s2],'.','Color','r','MarkerSize',25)

axis([min(uA),max(uA),min(uB),max(uB)])
title('Heaviside')


%% sigmoid
hh2=subplot(1,2,2);
hold on
lambda=20;
H=@(x) 1./(1+exp(-lambda*x));
rhs=@(t,y)[-y(1)+H(a*(y(2)-s2)); -y(2)+H(a*(y(1)-s1))];

for x=u
    [~,y] = ode45(rhs, tspan, [x,min(uB)], options);
    plot(y(:,1),y(:,2),'Color',return_color(y))
    v=rhs(0,[y(idx,1),y(idx,2)]); norma=sqrt(v(1)^2+v(2)^2);



    
    [~,y] = ode45(rhs, tspan, [x,max(uB)], options);
    plot(y(:,1),y(:,2),'Color',return_color(y))
    v=rhs(0,[y(idx,1),y(idx,2)]); norma=sqrt(v(1)^2+v(2)^2);



    
    [~,y] = ode45(rhs, tspan, [min(uA),x], options);
    plot(y(:,1),y(:,2),'Color',return_color(y))
    v=rhs(0,[y(idx,1),y(idx,2)]); norma=sqrt(v(1)^2+v(2)^2);



    
    [~,y] = ode45(rhs, tspan, [max(uA),x], options);
    plot(y(:,1),y(:,2),'Color',return_color(y))
    v=rhs(0,[y(idx,1),y(idx,2)]); norma=sqrt(v(1)^2+v(2)^2);



end

u=linspace(min(uA),max(uA),100);
plot(H(a*(u-s2)),u,'b','LineWidth',1) % uA-nullcline
plot(u,H(a*(u-s1)),'r','LineWidth',1) % uB-nullcline

axis([min(uA),max(uA),min(uB),max(uB)])

xticks([min(uA),s1,max(uA)])
xticklabels({num2str(min(uA)),'s_1',num2str(max(uA))})
yticks([min(uB),s2,max(uB)])
yticklabels({num2str(min(uB)),'s_2',num2str(max(uB))})
title('Sigmoid')


newtol=1e-6;
newmaxit=100;
hjac=1e-4;
f=@(y)[-y(1)+H(a*(y(2)-s2)); -y(2)+H(a*(y(1)-s1))];
df=@(y)MyJacobian(f,y,hjac);

% 0 and 1 fixed points 
[fp0,conv,~]=MySolve(f,[0;0],df,newtol,newmaxit);
if conv
    plot(fp0(1),fp0(2),'.','Color','k','MarkerSize',25)
end
[fp1,conv,~]=MySolve(f,[1;1],df,newtol,newmaxit);
if conv
    plot(fp1(1),fp1(2),'.','Color','k','MarkerSize',25)
end

%
[fp,conv,J]=MySolve(f,[s1;s2],df,newtol,newmaxit);

[V,D]=eig(df(fp));
h=0.00001;
vv=V(:,2);

point_sepL=fp+h*vv;
point_sepR=fp-h*vv;

[~,y] = ode45(rhs, -tspan, point_sepL, options);
plot(y(:,1),y(:,2),'Color',yellow,'LineWidth',2)

[~,y] = ode45(rhs, -tspan, point_sepR, options);
plot(y(:,1),y(:,2),'Color',orange,'LineWidth',2)

plot(fp(1),fp(2),'.','Color','r','MarkerSize',25)

figure(1)
subplot(1,2,1)
set(gca,'FontSize',font_size)
subplot(1,2,2)
set(gca,'FontSize',font_size)



p = get(hh1,'pos');
p(1) = p(1)-0.06;
p(3)=p(3)+0.07;
set(hh1,'pos',p);


p = get(hh2,'pos');
p(1)=p(1)-0.02;
p(3)=p(3)+0.07;
set(hh2,'pos',p);





