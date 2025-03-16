function [rf,grad] = design_pTxSSWE(ID, SlabThickness, dwell_time)

%% Load data
addpath('data/');
gamma= 2.675e8;

path = ['data/s',num2str(ID),'/'];
filename = [path,'calibdata_axi'];
load(filename);
Nc = 8;
xx = 50;
yy = 40;

%%%Downsample calibration
US = 1;
UA = [102*US,xx,yy];
b0map_axi = imresize3(b0map_axi,UA);
mask_axi = imresize3(mask_axi,UA);
 temp2 = [];
for i = 1:Nc
    temp = squeeze(rfmap_axi(:,:,:,i));
    temp = imresize3(temp,UA,'box');
    temp2(:,:,:,i) = temp;
end
rfmap_axi = temp2;

%% Gradient 
%%% Const definition
SR = 160;
mg = 8/1000;
dt = 2e-6;
up = SR*dt:SR*dt:mg;
up = [up;zeros(2,length(up))];
down = mg-SR*dt:-SR*dt:0;
down = [down;zeros(2,length(down))];
ts = floor(0.5e-3/dt-size(up,2)-size(down,2));
step = [mg*ones(1,ts);zeros(2,ts)];
grad1 = [up,step,down];
grad = [grad1,-grad1,grad1,-grad1,grad1,-grad1];
lengthKP = size(grad,2);

%%% Add rewinder
SS = sum(grad1(1,:))/2;
rt = round(sqrt(SS/(SR*dt)));
up = SR*dt:SR*dt:SR*dt*rt;
up = [up;zeros(2,length(up))];
down = SR*dt*rt-SR*dt:-SR*dt:0;
down = [down;zeros(2,length(down))];
Grewinder = [up,down];
SS1 = sum(Grewinder(1,:))/2;
grad = [grad,Grewinder];

%%% Add Spiral
T = lengthKP*dt;
len = lengthKP;
kpr = zeros(3,len);
for i = 1:len
    ttt = i*dt;
    Kp = 20/(1+exp(10*(ttt/T-0.5)));
    Ktheta = 8/2.88*pi*1000*ttt;
    Kphi = 2/2.88*pi*1000*ttt;
    Kp_temp = [Kp*sin(Ktheta)*cos(Kphi);Kp*sin(Ktheta)*sin(Kphi);Kp*cos(Ktheta)];
    kpr(:,i) = Kp_temp;
end
dkval = diff(kpr,1,2);
gradr = dkval/gamma/dt;
gradr = [gradr,zeros(3,size(Grewinder,2)+1)];
grad = [grad(1,:);gradr(2:3,:)];

%% Mask location
Maskloc = double(mask_axi);
Maskloc = sum(Maskloc,2);
Maskloc = sum(Maskloc,3);

for i = length(Maskloc):-1:1
    if Maskloc(i)>0
        break
    end
end
corx = 52-round(SlabThickness/6);
cory = 51+round(SlabThickness/6);

%% Draw K-space Trajectory
kp3 = cumsum(grad(1,size(grad,2):-1:1))*dt*gamma;
kp2 = cumsum(grad(2,size(grad,2):-1:1))*dt*gamma;
kp1 = cumsum(grad(3,size(grad,2):-1:1))*dt*gamma;
kp = [kp1;kp2;kp3];
figure;
for i = 1:floor(size(kp,2)/3)
    g = 3;
    index = (i-1)*g;
    if i == 1
        plot3(kp(1,index+1:index+g),kp(2,index+1:index+g),kp(3,index+1:index+g),'-','linewidth',5,'MarkerIndices',[1,g],'Color',[105,55/96+100,155]/255);
    else
        plot3(kp(1,index:index+g),kp(2,index:index+g),kp(3,index:index+g),'-','linewidth',5,'MarkerIndices',[1,g+1],'Color',[105,55/96+100,155]/255);
    end
    hold on;
end
xlim([-50,50]);
ylim([-50,50]);
zlim([-600,600]);
grid on
set(gca,'Fontname','Times New Roman','fontsize',32);
ax = gca;
ax.LineWidth = 4; % 坐标轴线条磅数
title('k-space trajectory');

%% Calculation


%%%%---------------------Passband: Water inside the Slab-------------------------%%%%
kp1 = cumsum(grad(1,1:lengthKP))*dt*gamma;
kp2 = cumsum(grad(2,1:lengthKP))*dt*gamma;
kp3 = cumsum(grad(3,1:lengthKP))*dt*gamma;
kp = [kp1;kp2;kp3];


TFA = 13;

foxkt = [0.306,0.3,0.24];
poffset = [0 0 0]; %实际FOV与场图之间的偏移
soi = 1:102*US;
b1mapsMSn = rfmap_axi(soi,:,:,:);
b0mapMS = b0map_axi(soi,:,:)/1e6;
maskMS = mask_axi(soi,:,:);


list = 1:(corx-1)*US;
maskMS(list,:,:) = 0*maskMS(list,:,:);
list = US*cory+1:US*102;
maskMS(list,:,:) = 0*maskMS(list,:,:);

[b1arr,posarr] = create_array(b1mapsMSn,maskMS,foxkt,1e-3*poffset);



phasetrack = 1:lengthKP;
gObj = gradPulse(ones(1,phasetrack(end)),dt,'unitgrad','Tx');
kb0 = gObj.KspaceTrajectory;
nKT = size(kp,2);

A = [];b=[];
lB0 = 0.3;
kernalmat = 1i.* gamma.* dt.* exp(1i* ( posarr*kp + b0mapMS(maskMS)*kb0(phasetrack) - 0/(gamma/2/pi)*kb0(phasetrack) ) );
wB0 = diag(lB0*abs(b0mapMS(maskMS)*1e6)+1);
[nspa,nchs] = size(b1arr);
sysmat = complex(zeros(nspa,nchs*nKT)); % this prealloc results in much faster construct
for idx = 1:nchs
    iidx0 = (idx-1)*nKT + 1;
    sysmat(:,iidx0:(iidx0+nKT-1)) = wB0*diag(b1arr(:,idx)) * kernalmat;
end
A = [A;sysmat];

tar = wB0*(TFA/180)*pi*ones(size(sysmat,1),1);
b = [b;tar];
%%%%---------------------End of the Stopband Calculation-------------------------%%%%

%%%%---------------------Stopband1: Outside of the Slab-------------------------%%%%
filename = [path,'calibdata_axi'];
load(filename);
Nc = 8;
xx = 10;
yy = 8;
UA = [102,xx,yy];
b0map_axi = imresize3(b0map_axi,UA);

mask_axi = imresize3(mask_axi,UA);

temp2 = [];
for i = 1:Nc
    temp = squeeze(rfmap_axi(:,:,:,i));
    temp = imresize3(temp,UA,'box');
    temp2(:,:,:,i) = temp;
end
rfmap_axi = temp2;

kp1 = cumsum(grad(1,1:lengthKP))*dt*gamma;
kp2 = cumsum(grad(2,1:lengthKP))*dt*gamma;
kp3 = cumsum(grad(3,1:lengthKP))*dt*gamma;
kp = [kp1;kp2;kp3];

TFA = 0;

foxkt = [0.306,0.3,0.24];

soi = 1:102;
b1mapsMSn = rfmap_axi(soi,:,:,:);
b0mapMS = b0map_axi(soi,:,:)/1e6;
maskMS = mask_axi(soi,:,:);


list = corx-2:cory+2;
maskMS(list,:,:) = 0*maskMS(list,:,:);

[b1arr,posarr] = create_array(b1mapsMSn,maskMS,foxkt,1e-3*poffset);

phasetrack = 1:lengthKP;
gObj = gradPulse(ones(1,phasetrack(end)),dt,'unitgrad','Tx');
kb0 = gObj.KspaceTrajectory;
nKT = size(kp,2);

kernalmat = 1i.* gamma.* dt.* exp(1i* ( posarr*kp + b0mapMS(maskMS)*kb0(phasetrack) - 0/(gamma/2/pi)*kb0(phasetrack) ) );
[nspa,nchs] = size(b1arr);
sysmat = complex(zeros(nspa,nchs*nKT)); % this prealloc results in much faster construct
for idx = 1:nchs
    iidx0 = (idx-1)*nKT + 1;
    sysmat(:,iidx0:(iidx0+nKT-1)) = diag(b1arr(:,idx)) * kernalmat;
end
wR = 10;
A = [A;wR*sysmat];

tar = (TFA/180)*pi*ones(size(sysmat,1),1);
b = [b;tar];


%%%%---------------------Stopband2: Fat in the whole brain-------------------------%%%%
maskMS = mask_axi(soi,:,:);
list = corx-2:corx-1;
maskMS(list,:,:) = 0*maskMS(list,:,:);
list = cory+1:cory+2;
maskMS(list,:,:) = 0*maskMS(list,:,:);

[b1arr,posarr] = create_array(b1mapsMSn,maskMS,foxkt,1e-3*poffset);
phasetrack = 1:lengthKP;
gObj = gradPulse(ones(1,phasetrack(end)),dt,'unitgrad','Tx');
kb0 = gObj.KspaceTrajectory;
nKT = size(kp,2);
wtS = 8;
kernalmat = 1i.* gamma.* dt.* exp(1i* ( posarr*kp + b0mapMS(maskMS)*kb0(phasetrack) + 1000/(gamma/2/pi)*kb0(phasetrack) ) );
[nspa,nchs] = size(b1arr);
sysmat = complex(zeros(nspa,nchs*nKT)); % this prealloc results in much faster construct
for idx = 1:nchs
    iidx0 = (idx-1)*nKT + 1;
    sysmat(:,iidx0:(iidx0+nKT-1)) = diag(b1arr(:,idx)) * kernalmat;
end
A = [A;wtS*sysmat];

tar = 0*(TFA/180)*pi*ones(size(sysmat,1),1);
b = [b;tar];

[b1arr,posarr] = create_array(b1mapsMSn,maskMS,foxkt,1e-3*poffset);
phasetrack = 1:lengthKP;
gObj = gradPulse(ones(1,phasetrack(end)),dt,'unitgrad','Tx');
kb0 = gObj.KspaceTrajectory;
nKT = size(kp,2);
kernalmat = 1i.* gamma.* dt.* exp(1i* ( posarr*kp + b0mapMS(maskMS)*kb0(phasetrack) + 900/(gamma/2/pi)*kb0(phasetrack) ) );
[nspa,nchs] = size(b1arr);
sysmat = complex(zeros(nspa,nchs*nKT)); % this prealloc results in much faster construct
for idx = 1:nchs
    iidx0 = (idx-1)*nKT + 1;
    sysmat(:,iidx0:(iidx0+nKT-1)) = diag(b1arr(:,idx)) * kernalmat;
end
A = [A;wtS*sysmat];

tar = 0*(TFA/180)*pi*ones(size(sysmat,1),1);
b = [b;tar];


[b1arr,posarr] = create_array(b1mapsMSn,maskMS,foxkt,1e-3*poffset);
phasetrack = 1:lengthKP;
gObj = gradPulse(ones(1,phasetrack(end)),dt,'unitgrad','Tx');
kb0 = gObj.KspaceTrajectory;
nKT = size(kp,2);
kernalmat = 1i.* gamma.* dt.* exp(1i* ( posarr*kp + b0mapMS(maskMS)*kb0(phasetrack) + 1100/(gamma/2/pi)*kb0(phasetrack) ) );
[nspa,nchs] = size(b1arr);
sysmat = complex(zeros(nspa,nchs*nKT)); % this prealloc results in much faster construct
for idx = 1:nchs
    iidx0 = (idx-1)*nKT + 1;
    sysmat(:,iidx0:(iidx0+nKT-1)) = diag(b1arr(:,idx)) * kernalmat;
end
A = [A;wtS*sysmat];

tar = 0*(TFA/180)*pi*ones(size(sysmat,1),1);
b = [b;tar];


[b1arr,posarr] = create_array(b1mapsMSn,maskMS,foxkt,1e-3*poffset);
phasetrack = 1:lengthKP;
gObj = gradPulse(ones(1,phasetrack(end)),dt,'unitgrad','Tx');
kb0 = gObj.KspaceTrajectory;
nKT = size(kp,2);
kernalmat = 1i.* gamma.* dt.* exp(1i* ( posarr*kp + b0mapMS(maskMS)*kb0(phasetrack) + 1200/(gamma/2/pi)*kb0(phasetrack) ) );
[nspa,nchs] = size(b1arr);
sysmat = complex(zeros(nspa,nchs*nKT)); % this prealloc results in much faster construct
for idx = 1:nchs
    iidx0 = (idx-1)*nKT + 1;
    sysmat(:,iidx0:(iidx0+nKT-1)) = diag(b1arr(:,idx)) * kernalmat;
end
A = [A;wtS*sysmat];

tar = 0*(TFA/180)*pi*ones(size(sysmat,1),1);
b = [b;tar];

[b1arr,posarr] = create_array(b1mapsMSn,maskMS,foxkt,1e-3*poffset);
phasetrack = 1:lengthKP;
gObj = gradPulse(ones(1,phasetrack(end)),dt,'unitgrad','Tx');
kb0 = gObj.KspaceTrajectory;
nKT = size(kp,2);
kernalmat = 1i.* gamma.* dt.* exp(1i* ( posarr*kp + b0mapMS(maskMS)*kb0(phasetrack) + 800/(gamma/2/pi)*kb0(phasetrack) ) );
[nspa,nchs] = size(b1arr);
sysmat = complex(zeros(nspa,nchs*nKT)); % this prealloc results in much faster construct
for idx = 1:nchs
    iidx0 = (idx-1)*nKT + 1;
    sysmat(:,iidx0:(iidx0+nKT-1)) = diag(b1arr(:,idx)) * kernalmat;
end
A = [A;wtS*sysmat];

tar = 0*(TFA/180)*pi*ones(size(sysmat,1),1);
b = [b;tar];
%%%%---------------------End of the Stopband Calculation-------------------------%%%%


%%% MLS:RF calculation
[wt,rho,eta]= solve_mlstr(A,b,2000,1e-5);

wt11 = reshape(wt,[],nchs);
wt11 = wt11.';
myrf = [];
rf = {};
nt = 1;
for ind=1:lengthKP/nt
    rf{ind} = [];
    for j  =1:nt
        rf{ind} = [rf{ind},wt11(:,ind)];
    end
end
for ind=1:lengthKP/nt
    myrf = [myrf rf{ind}];
end
rf = myrf;

%% RF and Grad shown
% save rf rf
% save grad grad
figure;
fs = 24;
ax = gca;
ax.LineWidth = 3; % 坐标轴线条磅数
ax.FontName = 'Times New Roman'; % 字体格式
ax.FontSize = 14; % 字体大小
hold on
for i = 1:3
plot((1:size(grad,2))*dt*1000,(grad(i,:)*1000),'linewidth',3);
end
ylim([-20,20]);
xlim([0,3.21]);
xlabel('Time/ms');
ylabel('Grad/mT·m^{-1}');
box on
grid on
set(gca,'Fontname','Times New Roman','FontSize',fs);
legend('Gx','Gy','Gz','Location','Northeast');
title('Gradient');

rf_plot = 1*[rf,zeros(8,size(grad,2)-size(rf,2))];
figure;
fs = 24;
ax = gca;
ax.LineWidth = 3; % 坐标轴线条磅数
ax.FontName = 'Times New Roman'; % 字体格式
ax.FontSize = 14; % 字体大小
hold on
for i = 1:8
plot((1:size(rf_plot,2))*dt*1000,abs(rf_plot(i,:))*1e6,'linewidth',3);
end
xlim([0,3.21]);
xlabel('Time/ms');
ylabel('|RF|/V');
box on
grid on
legend('ch1','ch2','ch3','ch4','ch5','ch6','ch7','ch8','Location','Northeast');
set(gca,'Fontname','Times New Roman','FontSize',fs);
title('RF');

%%%%-----------------------------------------------%%%%
%%%%-----------------------------------------------%%%%
%%%%-----------------------------------------------%%%%
%%%%-----------------------------------------------%%%%
%%%%-----------------------------------------------%%%%
poffset = [0,0,0];

%% Bloch simulation for water and show
odt = dt;
dt = dwell_time;
CC = odt/dt;
grad = imresize(grad,[3,round(size(grad,2)*CC)]);
rf = imresize(rf,[8,round(size(rf,2)*CC)]);
% rf = [zeros(8,2),rf];
Nslab = 8;
Ci = Nslab+1;
% Wb = zeros(3,102*UA,100*UA,80*UA);
ID = 1;
filename = [path,'calibdata_axi'];
load(filename);

soi = 1:1:80;
maskMS = mask_axi(:,:,soi);
b1mapsMSn = rfmap_axi(:,:,soi,:);
b0mapMS = b0map_axi(:,:,soi)/1e6;
[Mxypat,Mzpat] = run_bloch_sim (rf,grad,b1mapsMSn,maskMS,[0.306,0.3,0.24],b0mapMS,...
        0,[],dt,1e-3*poffset);



%%%% ----------------------Show Sagittal View---------------------%%%%
mxypat = zeros(102,100,80);
mzpat = zeros(102,100,80);
mxypat = imresize3(Mxypat,1,'box');
mzpat = imresize3(Mzpat,1,'box');
mask_axi1 = imresize3(mask_axi,1,'box');

mxy = mxypat(corx:cory,:,:);
masks = mask_axi1(corx:cory,:,:);
vec = asin(abs(mxy(masks)))*180/pi;
RMSE_water = sqrt(sum((vec-13).^2)/length(vec));
Cov_water = std(vec)/mean(vec);


% save mxypatw4 mxypat
slice = 50;
Img = squeeze(mxypat(:,slice,:));
Img1 = squeeze(mzpat(:,slice,:));

Img2 = rot90(asin( abs(Img))/pi*180); 
Img2(Img2==90)=0;

fs = 16;
figure;
imagesc(Img2);
t = 0:0:0;
set(gca,'xtick',t);
set(gca,'ytick',t);
h=colorbar;
colormap hot
set(get(h,'title'),'string','[°]');
set(gca,'Fontname','Times New Roman','FontSize',fs);
caxis([0,20]);
axis image;
title('Water-Sagittal (M)');

Img2 = rot90(angle(Img)); 

fs = 16;
figure;
imagesc(Img2);
t = 0:0:0;
set(gca,'xtick',t);
set(gca,'ytick',t);
h=colorbar;
set(get(h,'title'),'string','[°]');
set(gca,'Fontname','Times New Roman','FontSize',fs);
caxis([-3.2,3.2]);
axis image;
colorcet('C2');
title('Water-Sagittal (P)');

%%%% ----------------------Show Coronal View---------------------%%%%
% mxypat = zeros(102,100,80);
% mzpat = zeros(102,100,80);
% mxypat = imresize3(Mxypat,1,'box');
% mzpat = imresize3(Mzpat,1,'box');
% mask_axi1 = imresize3(mask_axi,1,'box');
% 
% 
% for slice = corx-2:cory
%     Img = squeeze(mxypat(slice,:,:));
%     Img1 = squeeze(mzpat(slice,:,:));
% 
%     Img2 = rot90(asin( abs(Img))/pi*180); 
%     Img2(Img2==90)=0;
% 
%     fs = 16;
%     figure;
%     imagesc(Img2);
%     t = 0:0:0;
%     set(gca,'xtick',t);
%     set(gca,'ytick',t);
%     h=colorbar;
%     set(get(h,'title'),'string','[°]');
%     set(gca,'Fontname','Times New Roman','FontSize',fs);
%     caxis([0,15]);
%     axis image;
% 
%     Img2 = rot90(angle(Img)); 
%     fs = 16;
%     figure;
%     imagesc(Img2);
%     t = 0:0:0;
%     set(gca,'xtick',t);
%     set(gca,'ytick',t);
%     h=colorbar;
%     set(get(h,'title'),'string','[°]');
%     set(gca,'Fontname','Times New Roman','FontSize',fs);
%     caxis([-3.2,3.2]);
%     axis image;
%     colorcet('C2');
% end

%% Bloch simulation for fat and show
Nslab = 8;
Ci = Nslab+1;
% Wb = zeros(3,102*UA,100*UA,80*UA);
ID = 1;

filename = [path,'calibdata_axi'];
load(filename);

soi = 1:1:80;
maskMS = mask_axi(:,:,soi);
b1mapsMSn = rfmap_axi(:,:,soi,:);
b0mapMS = b0map_axi(:,:,soi)/1e6;
[Mxypat,Mzpat] = run_bloch_sim (rf,grad,b1mapsMSn,maskMS,[0.306,0.3,0.24],b0mapMS+1000/(gamma/2/pi),...
        0,[],dt,1e-3*poffset);



%%%% ----------------------Show Sagittal View---------------------%%%%
mxypat = zeros(102,100,80);
mzpat = zeros(102,100,80);
mxypat = imresize3(Mxypat,1,'box');
mzpat = imresize3(Mzpat,1,'box');
mask_axi1 = imresize3(mask_axi,1,'box');

masks = mask_axi1(:,:,:);
mxy = mxypat(:,:,:);
vec = asin(abs(mxy(masks)))*180/pi;
RMSE_fat = sqrt(sum((vec-0).^2)/length(vec));

% save mxypatf4 mxypat
slice = 50;
Img = squeeze(mxypat(:,slice,:));
Img1 = squeeze(mzpat(:,slice,:));

Img2 = rot90(asin( abs(Img))/pi*180); 
Img2(Img2==90)=0;

fs = 16;
figure;
imagesc(Img2);
t = 0:0:0;
set(gca,'xtick',t);
set(gca,'ytick',t);
h=colorbar;
set(get(h,'title'),'string','[°]');
set(gca,'Fontname','Times New Roman','FontSize',fs);
colormap hot
colorbar
caxis([0,6])
axis image;
title('Fat-Sagittal (M)');


% %% Bloch simulation for 1mm resolution and show
% Nslab = 8;
% Ci = Nslab+1;
% % Wb = zeros(3,102*UA,100*UA,80*UA);
% ID = 1;
% 
% filename = [path,'calibdata_axi'];
% load(filename);
% 
% %----------------------Upsample--------------------%
% UA = 3;
% b0map_axi = imresize3(b0map_axi,UA);
% 
% mask_axi = imresize3(mask_axi,UA);
% temp = false(size(mask_axi));
% ori = 'sag';
% if ori == 'sag'
%     ss = 50;
%     Rs = (ss-1)*UA+1:ss*UA;
%     temp(:,Rs,:) = mask_axi(:,Rs,:);
%     temp(:,Rs,:) = mask_axi(:,Rs,:);
% else
%     as = 51;
%     Rs = (as-1)*UA+1:as*UA;
%     temp(Rs,:,:) = mask_axi(Rs,:,:);
%     temp(Rs,:,:) = mask_axi(Rs,:,:);
% end
% 
% mask_axi = temp;
% temp2 = [];
% for i = 1:Nc
%     temp = squeeze(rfmap_axi(:,:,:,i));
%     temp = imresize3(temp,UA,'box');
%     temp2(:,:,:,i) = temp;
% end
% rfmap_axi = temp2;
% %----------------------Upsample--------------------%
% 
% for i = 4:4 
%     i
% %     mrf1 = modulate_rfphase (rf, grad(1,:), (i-Ci)*1.2+1.5, dt*1e6, 0);
% %     mrf2 = modulate_rfphase (rf, grad(1,:), (i-1)*1.2+1.5, dt*1e6, 0);
%     soi = 1:1:80*UA;
%     maskMS = mask_axi(:,:,soi);
%     b1mapsMSn = rfmap_axi(:,:,soi,:);
%     b0mapMS = b0map_axi(:,:,soi)/1e6;
%     [Mxypat,Mzpat] = run_bloch_sim (rf,grad,b1mapsMSn,maskMS,[0.306,0.3,0.24],b0mapMS+0/(gamma/2/pi),...
%             0,[],dt,1e-3*poffset);
% %     Wb(1,:,:,:) = real(mxypat);
% %     Wb(2,:,:,:) = imag(mxypat);
% %     Wb(3,:,:,:) = mzpat;
% %     [mxypat,mzpat] = run_bloch_sim (2*exp(1i*phase)*rf,grad,b1mapsMSn,maskMS,[0.3,0.3,0.24],b0mapMS/(gamma/2/pi),...
% %             0,[],dt,poffset,Wb);
% end
% 
% 
% %%%% Show sagittal view
% mxypat = zeros(102,100,80);
% mzpat = zeros(102,100,80);
% mxypat = imresize3(Mxypat,1,'box');
% mzpat = imresize3(Mzpat,1,'box');
% mask_axi1 = imresize3(mask_axi,1,'box');
% 
% 
% slice = 150;
% Img = squeeze(mxypat(:,slice,:));
% Img1 = squeeze(mzpat(:,slice,:));
% 
% Img2 = rot90(asin( abs(Img))/pi*180); 
% Img2(Img2==90)=0;
% 
% fs = 16;
% figure;
% imagesc(Img2);
% t = 0:0:0;
% set(gca,'xtick',t);
% set(gca,'ytick',t);
% h=colorbar;
% set(get(h,'title'),'string','[°]');
% set(gca,'Fontname','Times New Roman','FontSize',fs);
% caxis([0,15]);
% axis image;

%% Output
myrf = [rf,zeros(8,size(grad,2)-size(rf,2))];
mygrad = grad;
[Nc,Nt] = size(myrf);
rf_out = 1e6*reshape(myrf',[Nc*Nt,1]);
temp = mygrad;
temp(3,:) = -mygrad(1,:);
temp(1,:) = mygrad(3,:);
temp(2,:) = -mygrad(2,:);
grad_out = 1e3*temp';
% mygrad = 1e3*temp';
% 
save_pTXRFPulse_toINI(grad_out,rf_out,[]);
