
%Starting point
x0=0;
y0=0;
fi0=0;
u0=0;
v0=0;
r0=0;
z10=cos(fi0)*x0+sin(fi0)*y0;
z20=-sin(fi0)*x0+cos(fi0)*y0;
z30=fi0;

[t,sol] = ode45('modelsim7ehv', [0,60], [z10, z20, z30, u0, v0, r0]);
z1 = sol(:,1);
z2 = sol(:,2);
z3 = sol(:,3);
u = sol(:,4);
v = sol(:,5);
w = sol(:,6);

x=[];
y=[];
fi=[];

for hihi=1:length(sol(:,1))
    trans=inv([cos(z3(hihi,1)) sin(z3(hihi,1)) ;-sin(z3(hihi, 1)) cos(z3(hihi, 1))]) ;
    x=[x trans(1, :)*[z1(hihi, 1) z2(hihi,1)]'];
    y=[y trans(2,:)*[z1(hihi,1) z2(hihi,1)]' ];
    fi=[fi z3(hihi,1)];
end

figure(1) ,plot(t,u) ,title('surge') ;
xlabel('t [ s ] '), ylabel('ve1ocity [ m/s ]')
hold on

plot([0 40] ,[2.25 2.25],'k',[40 40],[2.25 1.125],'k',[40 50] ,[1.125 1.125],'k',[50 50],[1.125 0],'k',[50 60],[0 0],'k')
plot([0 5],[0 0],'r',[5 5] ,[0 1.5] ,'r',[5 60] ,[1.5 1.5],'r')
figure(2) ,plot(t ,v),title('sway');
xlabel('t [ s ]') ,ylabel('velocity [m/s]')
figure(3) ,plot(t,w) ,title('yaw') ;
xlabel('t [ s ]') ,ylabel('angu1ar velocity [ rad/s ] ')

figure(4)
okj=find(t<5*5&t>0) ;
plot( [x(okj) x(okj (length(okj))+1)], [y(okj) y(okj (length(okj))+1)] ,'b')
hold on

kleurtjes=['r' 'r' 'r' 'r' 'r' 'r' 'r' 'k' 'k' 'k' 'k'];

for ay=1:11
    okj=find(t<(ay+1)*5&t>ay*5);
    plot ([x(okj) x(okj(length(okj))+1)], [y(okj) y(okj(length(okj))+1)] ,kleurtjes(ay))
hold on
end

title('Movement of the hovercraft in earth coordinates');
hold on
x2=x(find(t>=45));
y2=y(find(t>=45));

[i,j]=size(sol);

lll=0.3;%length of the ship

n=[30 34 82 99 105 109 112 114 116 119 124 131 138 150 190 227 245 259 578 4601];

colortjes=['b' 'b' 'r' 'r' 'r' 'r' 'r' 'r' 'r' 'r' 'r' 'r' 'r' 'r' 'r' 'r' 'k' 'k' 'k' 'k' 'k'];

for nij=1:1:length(n)
    Shipplot(x(1,n(nij)),y(1,n(nij)),fi(1,n(nij)),lll,colortjes(nij));
end

xlabel('X [m]'), ylabel('Y [m]')

figure(5),subplot(2,42,1:16) ;
plot(t,u),hold on,title('surge');
xlabel('t [ s ]') ,ylabel('velocity [ m/s ]'),
plot([0 40],[2.25 2.25],'k', [40 40],[2.25 1.125],'k',[40 50], [1.125 1.125],'k',[50 50],[1.125 0],'k',[50 60],[0 0],'k')
plot([0 5],[0 0],'r',[5 5],[0 1.5],'r',[5 60],[1.5 1.5],'r')

subplot(2,42,18), title('F1'), subplot(2,42,20), title('Rudder');

subplot(2,42,23:38)
plot (t,v),hold on,
title('sway');
xlabel('t [ s ]') ,ylabel('velocity [ m/s ]'),
plot([0 40] ,[2.25 2.25],'k',[40 40],[2.25 1.125],'k',[40 50], [1.125 1.125],'k',[50 50],[1.125 0],'k',[50 60],[0 0],'k')
plot([0 5],[0 0],'r',[5 5],[0 1.5],'r',[5 60],[1.5 1.5],'r')

subplot(2,42,40), title('F1'), subplot(2,42,42), title('Rudder');

subplot(2,42,43:58),plot(t,w) ,hold on,title('yaw'),
plot([0 40],[2.25 2.25],'k', [40 40],[2.25 1.125],'k', [40 50], [1.125 1.125],'k',[50 50],[1.125 0],'k',[50 60],[0 0],'k')
plot([0 5],[0 0],'r',[5 5],[0 1.5],'r',[5 60],[1.5 1.5],'r')
xlabel('t [ s ] ' ) ,ylabel('angular velocity [ rad/s ] '),

subplot(2,42,60),title('F1'),subplot(2,42,62),title('Rudder');

subplot(2,42,65:84)
okj=find(t<5*5&t>0);
plot([x(okj) x(okj(length(okj))+1)], [y(okj) y(okj(length(okj))+1)], 'b')
hold on
for ay=1:11
    okj=find(t<(ay+1)*5&t>ay*5) ;
    plot([x(okj) x(okj(length(okj))+1)],[y(okj) y(okj(length(okj))+1)] ,kleurtjes(ay))
    hold on
end

title('Movement of the hovercraft in earth coordinates');
hold on
[i,j]=size(sol) ;

xlabel('X [m]'), ylabel('Y [m]')

set(4, 'Position', [0 50 530-30 530-50]);
set(5, 'Position', [500 50 530-30 530-50]);

