% Temporal Spatial Norms Organization 

[NUM, TXT, RAW] = xlsread('TempSpatNorms.xlsx');

%% Time Distance Parameters
% Columns are: Cadence, Stride Length, Walking Speed, Step Lengt
% Units are: steps/min, m, m/min, and cm

TempSpatNorms.Age2.TimeDist = NUM(1, 2:7);
TempSpatNorms.Age3.TimeDist = NUM(2, 2:7);
TempSpatNorms.Age4.TimeDist = NUM(3, 2:7);
TempSpatNorms.Age5.TimeDist = NUM(4, 2:7);
TempSpatNorms.Age6.TimeDist = NUM(5, 2:7);
TempSpatNorms.Age7.TimeDist = NUM(6, 2:7);
TempSpatNorms.Age8.TimeDist = NUM(7, 2:7);
TempSpatNorms.Age9.TimeDist = NUM(8, 2:7);
TempSpatNorms.Age10.TimeDist = NUM(9, 2:7);
TempSpatNorms.Age11.TimeDist = NUM(10, 2:7);
TempSpatNorms.Age12.TimeDist = NUM(11, 2:7);
TempSpatNorms.Age13.TimeDist = NUM(12, 2:7);
TempSpatNorms.Age14.TimeDist = NUM(13, 2:7);
TempSpatNorms.Age15.TimeDist = NUM(14, 2:7);
TempSpatNorms.Age16.TimeDist = NUM(15, 2:7);
TempSpatNorms.AgeAdult.TimeDist = NUM(16, 2:7);

%% Stance-Swing Parameters
% Columns are: Stance period, swing period, Initial double support, single
% support, and final double support. Units are all % of gait cycle

TempSpatNorms.Age2.StanceSwing = NUM(1, 8:12);
TempSpatNorms.Age3.StanceSwing = NUM(2, 8:12);
TempSpatNorms.Age4.StanceSwing = NUM(3, 8:12);
TempSpatNorms.Age5.StanceSwing = NUM(4, 8:12);
TempSpatNorms.Age6.StanceSwing = NUM(5, 8:12);
TempSpatNorms.Age7.StanceSwing = NUM(6, 8:12);
TempSpatNorms.Age8.StanceSwing = NUM(7, 8:12);
TempSpatNorms.Age9.StanceSwing = NUM(8, 8:12);
TempSpatNorms.Age10.StanceSwing = NUM(9, 8:12);
TempSpatNorms.Age11.StanceSwing = NUM(10, 8:12);
TempSpatNorms.Age12.StanceSwing = NUM(11, 8:12);
TempSpatNorms.Age13.StanceSwing = NUM(12, 8:12);
TempSpatNorms.Age14.StanceSwing = NUM(13, 8:12);
TempSpatNorms.Age15.StanceSwing = NUM(14, 8:12);
TempSpatNorms.Age16.StanceSwing = NUM(15, 8:12);
TempSpatNorms.AgeAdult.StanceSwing = NUM(16, 8:12);


%% Save data in .mat file
save('TempSpatNorms.mat', 'TempSpatNorms'); 
