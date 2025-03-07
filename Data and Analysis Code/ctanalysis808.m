clear all

%%
data = xlsread("concentrationdata808.xlsx");

adaptivedata = data(3:42,2:23);
rotordata = data(3:42, 24:31);
basecycle = 1:1:40;

cycle = 1:0.01:40;
adaptivedata(isnan(adaptivedata)) = 0;
rotorinterp = interp1(basecycle,rotordata,cycle,"spline");
adaptiveinterp = interp1(basecycle,adaptivedata,cycle,"spline");


baseline = [adaptivedata(1:15,:),rotordata(1:15,:)];
threshold = 10*std(baseline(:));
rotorgeneindices = zeros(1,8);
adaptiveindices = zeros(1,22);

trash = 0;
for i = 1:8
    if(max(rotorinterp(:,i)) > threshold)
        rotorgeneindices(i) = find(rotorinterp(:,i) > threshold,1);
    else
        rotorgeneindices(i) = length(rotorinterp(:,i));
    end

end

for i = 1:22
    if(max(adaptiveinterp(:,i)) > threshold)
        adaptiveindices(i) = find(adaptiveinterp(:,i) > threshold,1);
    else
        adaptiveindices(i) = length(adaptiveinterp(:,i));
    end
end

rotorgeneCts = cycle(rotorgeneindices);
adaptiveCts = cycle(adaptiveindices);

averagesadaptive = [mean(adaptiveCts(1:6)) mean(adaptiveCts(7:12)) mean(adaptiveCts(13:18)) mean(adaptiveCts(19:22))];
averagesrotor = [mean(rotorgeneCts(1:2)),mean(rotorgeneCts(3:4)),mean(rotorgeneCts(5:6)),mean(rotorgeneCts(7:8))];

adaptiveaverage = mean([adaptiveCts(1:12),adaptiveCts(15)]);
rotoraverage = mean(rotorgeneCts(1:6));

%% For overall consistency Cq data
data = xlsread("Consistency for 808 PCR Additional Experiments.xlsx");

%%
adaptivedata = data(:,2:20);
basecycle = 1:1:40;

cycle = 1:0.01:40;
adaptivedata(isnan(adaptivedata)) = 0;
adaptiveinterp = interp1(basecycle,adaptivedata,cycle,"spline");


baseline = [adaptivedata(1:15,:)];
threshold = 10*std(baseline(:));
adaptiveindices = zeros(1,19);
%%
trash = 0;

for i = 1:19
    if(max(adaptiveinterp(:,i)) > threshold)
        adaptiveindices(i) = find(adaptiveinterp(:,i) > threshold,1);
    else
        adaptiveindices(i) = length(adaptiveinterp(:,i));
    end
end

adaptiveCts = cycle(adaptiveindices);


adaptiveaverage = mean(adaptiveCts);
