function[Output] = GetTempSpatNorms(age)

%% Ask for age if not present
a = exist('age');
if a == 0
prompt = 'What is the age of the subject?';
age = inputdlg(prompt);
end

if iscell(age)
    age = cell2mat(age);
end

if isstr(age)
    age = str2double(age);
end

%% Load data
load ('TempSpatNorms.mat'); 

%% Get Temp Spat data depending on age
if age == 2 % if they are 2
    TSN = TempSpatNorms.Age2;
end

if age == 3 % if they are 3
    TSN = TempSpatNorms.Age3;
end

if age == 4 % if they are 4
    TSN = TempSpatNorms.Age4;
end

if age == 5 % if they are 5
    TSN = TempSpatNorms.Age5;
end

if age == 6 % if they are 6
    TSN = TempSpatNorms.Age6;
end

if age == 7 % if they are 7
    TSN = TempSpatNorms.Age7;
end

if age == 8 % if they are 8
    TSN = TempSpatNorms.Age8;
end

if age == 9 % if they are 9
    TSN = TempSpatNorms.Age9;
end

if age == 10 % if they are 10
    TSN = TempSpatNorms.Age10;
end

if age == 11 % if they are 11
    TSN = TempSpatNorms.Age11;
end

if age == 12 % if they are 12
    TSN = TempSpatNorms.Age12;
end

if age == 13 % if they are 13
    TSN = TempSpatNorms.Age13;
end

if age == 14 % if they are 14
    TSN = TempSpatNorms.Age14;
end

if age == 15 % if they are 15
    TSN = TempSpatNorms.Age16;
end

if age == 16 % if they are 16
    TSN = TempSpatNorms.Age16;
end

if age == 0 || 16 < age  % if they are older than 16
    TSN = TempSpatNorms.AgeAdult;
end


Output = TSN;

end
    
