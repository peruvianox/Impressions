COPnorms = imread('GreenvilleCOPnorms.jpg');

imshow(COPnorms); 

SD_1.med = ginput(10);
SD_1.lat = ginput(10);

SD_2.med = ginput(10);
SD_2.lat = ginput(10);


figure; % hamburger style
hold on; 
plot(SD_1.med(:,1),359-SD_1.med(:,2), 'b' ); 
plot(SD_1.lat(:,1),359-SD_1.lat(:,2), 'b' ); 
plot(SD_2.med(:,1),359-SD_2.med(:,2), 'r' ); 
plot(SD_2.lat(:,1),359-SD_2.lat(:,2), 'r' ); 
xlim([0 182]);
ylim([0 359]);

figure; % hot dog style
hold on; 
plot(359-SD_1.med(:,2),SD_1.med(:,1), 'b' ); 
plot(359-SD_1.lat(:,2),SD_1.lat(:,1), 'b' ); 
plot(359-SD_2.med(:,2),SD_2.med(:,1), 'r' ); 
plot(359-SD_2.lat(:,2),SD_2.lat(:,1), 'r' ); 
ylim([0 182]);
xlim([0 359]);

% data interpolated and extracted using matlab curve fitting toolbox

SD_1.x_med = fx_SD_1med/189; 
SD_1.y_med = x_SD_1med/359; 
SD_1.x_lat = fx_SD_1lat/189; 
SD_1.y_lat = x_SD_1lat/359; 

SD_2.x_med = fx_SD_2med/189; 
SD_2.y_med = x_SD_2med/359; 
SD_2.x_lat = fx_SD_2lat/189; 
SD_2.y_lat = x_SD_2lat/359; 

figure; 
imshow(COPnorms); 
hold on; 
plot(189*SD_1.x_med, 359 - (359*SD_1.y_med), 'b' ); 
plot(189*SD_1.x_lat(:,1), 359-(359*SD_1.y_lat), 'b' ); 
plot(189*SD_2.x_med, 359-(359*SD_2.y_med), 'r' ); 
plot(189*SD_2.x_lat, 359-(359*SD_2.y_lat), 'r' ); 

COPnormsData.Image = COPnorms; 
COPnormsData.SD_1 = SD_1; 
COPnormsData.SD_2 = SD_2; 

save('COPnormsData.mat', 'COPnormsData');  






