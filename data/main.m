clc;clear;close all;
q0v = xlsread('data.xlsx');
q0 = q0v(:,1);
qv = q0v(:,2:4);
qv0 = [qv q0];

n = size(qv0,1);
Euler = zeros(n,3);
qNorm = zeros(n,1);
axis_spe = [0 0 1]';
AXIS = zeros(3,n);
for i = 1:n
    qv0(i,:) = qv0(i,:)/norm(qv0(i,:),2);
    Euler(i,:) = QtoEulerAngle(qv0(i,:));
    C_I2B = QtoC(qv0(i,:));
    AXIS(:,i) = C_I2B*axis_spe;
    qNorm(i,:) = norm(AXIS(:,i),2);
end

XTime = 0.1*(1:1:n)';
% figure(1)
% plot(XTime,Euler(:,1),'LineWidth',2);hold on
% plot(XTime,Euler(:,2),'g-.','LineWidth',2);hold on
% plot(XTime,Euler(:,3),'r--','LineWidth',2);hold on

figure(2)%四元数误差
plot(XTime,qv0(:,1),'LineWidth',2);hold on
plot(XTime,qv0(:,2),'LineWidth',2);hold on
plot(XTime,qv0(:,3),'g-.','LineWidth',2);hold on
plot(XTime,qv0(:,4),'r--','LineWidth',2);hold on

theta = linspace(0,pi,200);
phi = linspace(0,2*pi,300);
[Theta,Phi] = meshgrid(theta,phi);
R = 1;
X = R*sin(Theta).*cos(Phi);
Y = R*sin(Theta).*sin(Phi);
Z = R*cos(Theta);

figure(3)
plot3(AXIS(1,:),AXIS(2,:),AXIS(3,:),'LineWidth',2)
text(AXIS(1,1),AXIS(2,1),AXIS(3,1),'*','color','r','Fontsize',16);
text(AXIS(1,end),AXIS(2,end),AXIS(3,end),'o','color','r','Fontsize',10);
hold on
mesh(X,Y,Z)
axis equal



