
clear all; close all;clc;


Subjects = [1 2 3 4 5 8 9];

nEpisodes = 9;

for i=Subjects
    for j=2:nEpisodes
        eval(['Train1_ang_cor{',int2str(i),',',int2str(j),'} = csvread(''Segmented Movements\Vicon\Angles\m07_s0',int2str(i),'_e0',int2str(j),'_angles.txt'');']);
    end
end

for i=Subjects
    for j=nEpisodes+1
        eval(['Train1_ang_cor{',int2str(i),',',int2str(j),'} = csvread(''Segmented Movements\Vicon\Angles\m07_s0',int2str(i),'_e',int2str(j),'_angles.txt'');']);
    end
end

for i=Subjects
    for j=2:nEpisodes
        eval(['Test1_ang_inc{',int2str(i),',',int2str(j),'} = csvread(''Incorrect Segmented Movements\Vicon\Angles/m07_s0',int2str(i),'_e0',int2str(j),'_angles_inc.txt'');']);
    end
end

for i=Subjects
    for j=nEpisodes+1
        eval(['Test1_ang_inc{',int2str(i),',',int2str(j),'} = csvread(''Incorrect Segmented Movements\Vicon\Angles/m07_s0',int2str(i),'_e',int2str(j),'_angles_inc.txt'');']);
    end
end


k = 1;
for i=Subjects
    for j=2:nEpisodes+1
        Correct_Ini{k} = Train1_ang_cor{i,j};
        k = k+1;
    end
end

N = length(Correct_Ini);

nDim = size(Correct_Ini{1},2);

k = 1;
for i=Subjects
    for j=2:nEpisodes+1
        Incorrect_Ini{k} = Test1_ang_inc{i,j};
        k = k+1;
    end
end

nSubjects = length(Subjects);

for i = 1:nSubjects
    Len_Tr(i) = size(Correct_Ini{i},1);
end
Length_mean = 229;

for i=1:N
    for j=1:nDim
        Correct_Aligned{i}(:,j)=interp1([1:size(Correct_Ini{i},1)],Correct_Ini{i}(:,j),linspace(1,size(Correct_Ini{i},1),Length_mean));
    end
end

for i=1:N
    for j=1:nDim
        Incorrect_Aligned{i}(:,j)=interp1([1:size(Incorrect_Ini{i},1)],Incorrect_Ini{i}(:,j),linspace(1,size(Incorrect_Ini{i},1),Length_mean));
    end
end


Correct_Xm = [Correct_Aligned{1}'];
for i = 2:N
    Correct_Xm = [Correct_Xm; Correct_Aligned{i}'];
end

Incorrect_Xm = [Incorrect_Aligned{1}'];
for i = 2:N
    Incorrect_Xm = [Incorrect_Xm; Incorrect_Aligned{i}'];
end


Data_mean = repmat(mean(Correct_Xm,2), 1, size(Correct_Xm,2));
centered_Correct_Data = Correct_Xm - Data_mean;

Data_mean_inc = repmat(mean(Incorrect_Xm,2), 1, size(Incorrect_Xm,2));
centered_Incorrect_Data = Incorrect_Xm - Data_mean_inc;

scaling_value = ceil(max(max(max(centered_Correct_Data)),abs(min(min(centered_Correct_Data)))));
Data_Correct = centered_Correct_Data/scaling_value;
Data_Incorrect = centered_Incorrect_Data/scaling_value;

csvwrite('Data_Correct.csv',Data_Correct)
csvwrite('Data_Incorrect.csv',Data_Incorrect)

