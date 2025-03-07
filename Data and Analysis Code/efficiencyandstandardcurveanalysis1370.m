clear all 

%% Pull in data from Excel Files
rotordata = xlsread("Normalized 1370 Dilution Curve.xlsx",1);
adaptivedata = xlsread("Normalized 1370 Dilution Curve.xlsx",2);
timescurve = xlsread("Normalized 1370 Dilution Curve.xlsx",3);
timescurve = timescurve(1:41,35:66);

cyclerotor = 1:0.01:40;
cycleadaptive = 1:0.01:40;

rotorinterp = interp1(rotordata(:,1),rotordata(:,2:33),cyclerotor,"spline");
adaptiveinterp = interp1(adaptivedata(:,1),adaptivedata(:,2:33),cycleadaptive,"spline");
standardtimeinterp = interp1(adaptivedata(:,1),timescurve(2:41,:),cycleadaptive,"spline");
%% Cq and Efficiency Analysis
baseline = [adaptivedata(1:15,2:33),rotordata(1:15,2:33)];
threshold = 10*std(baseline(:));
rotorgeneindices = zeros(1,32);
adaptiveindices = zeros(1,32);
trash = 0;
for i = 1:32
    if(max(rotorinterp(:,i)) > threshold)
        rotorgeneindices(i) = find(rotorinterp(:,i) > threshold,1);
    else
        rotorgeneindices(i) = length(rotorinterp(:,i));
    end
    if(max(adaptiveinterp(:,i)) > threshold)
        adaptiveindices(i) = find(adaptiveinterp(:,i) > threshold,1);
    else
        adaptiveindices(i) = length(adaptiveinterp(:,i));
    end
end

rotorgeneCts = cyclerotor(rotorgeneindices);
adaptiveCts = cycleadaptive(adaptiveindices);
adaptiveCtime = zeros(1,32);

adaptiveexponentialphase = zeros(5,24);
rotorexponentialphase = zeros(5,24);
expphaserange = 4;
displacedx = 0:1:expphaserange;
adaptivewindowtest = zeros(24,5);
rotorwindowtest = zeros(24,5);
for i = 1:24
    target = round(adaptiveCts(i));
    adaptiveexponentialphase(:,i) = adaptivedata(target-2:target+2,i+1);
    target = round(rotorgeneCts(i));
    rotorexponentialphase(:,i) = rotordata(target-2:target+2,i+1);
end

adaptivelogexpphase = log(adaptiveexponentialphase);
rotorlogexpphase = log(rotorexponentialphase);
onesforfit = ones(length(displacedx),1);
fittest = [displacedx',onesforfit];
adaptivefits = zeros(2,24);
rotorfits = zeros(2,24);

for i = 1:24
    adaptivefits(:,i) = fittest\adaptivelogexpphase(:,i);
    rotorfits(:,i) = fittest\rotorlogexpphase(:,i);
end
adaptotherefficiencies = exp(adaptivefits(1,:));
avgadaptefficiency = mean(adaptotherefficiencies);
rotorotherefficiencies = exp(rotorfits(1,:));
avgrotorefficiency = mean(rotorotherefficiencies);
ctconca = exp(adaptivefits(2,:));
ctconca = ctconca./min(ctconca);
ctconcr = exp(rotorfits(2,:));
ctconcr = ctconcr./max(ctconcr);

for i = 1:32
    adaptiveCtime(i) = standardtimeinterp(adaptiveindices(i),i);
end

rotor105eff = mean(rotorotherefficiencies(1:4));
adaptive105eff = mean(adaptotherefficiencies(1:4));
rotor105std = std(rotorotherefficiencies(1:4));
adaptive105std = std(adaptotherefficiencies(1:4));

concentrations = [10^5,10^4,10^3,10^2,10^1,1];
averagerotor = [mean(rotorgeneCts(1:4)), mean(rotorgeneCts(5:8)), mean(rotorgeneCts(9:12)), mean(rotorgeneCts(13:16)), mean(rotorgeneCts(17:20)), mean(rotorgeneCts(21:24))];
stdrotor = [std(rotorgeneCts(1:4)), std(rotorgeneCts(5:8)), std(rotorgeneCts(9:12)), std(rotorgeneCts(13:16)), std(rotorgeneCts(17:20)), std(rotorgeneCts(21:24))];
averageadaptive = [mean(adaptiveCts(1:4)), mean(adaptiveCts(5:8)), mean(adaptiveCts(9:12)), mean(adaptiveCts(13:16)), mean(adaptiveCts(17:20)), mean(adaptiveCts(21:24))];
stdadaptive = [std(adaptiveCts(1:4)), std(adaptiveCts(5:8)), std(adaptiveCts(9:12)), std(adaptiveCts(13:16)), std(adaptiveCts(17:20)), std(adaptiveCts(21:24))];
rngadaptive = [range(adaptiveCts(1:4)), range(adaptiveCts(5:8)), range(adaptiveCts(9:12)), range(adaptiveCts(13:16)), range(adaptiveCts(17:20)), range(adaptiveCts(21:24))];

averageadaptivetime = [mean(adaptiveCtime(1:4)), mean(adaptiveCtime(5:8)), mean(adaptiveCtime(9:12)), mean(adaptiveCtime(13:16)), mean(adaptiveCtime(17:20)), mean(adaptiveCtime(21:24))];
stdadaptivetime = [std(adaptiveCtime(1:4)), std(adaptiveCtime(5:8)), std(adaptiveCtime(9:12)), std(adaptiveCtime(13:16)), std(adaptiveCtime(17:20)), std(adaptiveCtime(21:24))];


%% Linear efficiency fits and PCR plotting
figure
hold on
mdlrotor = fitlm(log10(concentrations),averagerotor);
mdladaptive = fitlm(log10(concentrations),averageadaptive);
plot([6,5,4,3,2,-1],mdlrotor.Coefficients.Estimate(1)+mdlrotor.Coefficients.Estimate(2)*[6,5,4,3,2,-1],LineWidth=2,Color="#0072BD",LineStyle='--')
plot([6,5,4,3,2,-1],mdladaptive.Coefficients.Estimate(1)+mdladaptive.Coefficients.Estimate(2)*[6,5,4,3,2,-1],LineWidth=2,Color="#A2142F",LineStyle='--')
eb = errorbar(log10(concentrations), averagerotor, stdrotor, 'ob',"LineWidth",1);
eb2 = errorbar(log10(concentrations), averageadaptive, stdadaptive, 'or',"LineWidth",1);
ylabel("Cq Value")
xlabel("Log_{10} [Copies/\muL]")
set(gca,'FontName' , 'Arial')
set(gca,'FontSize',20)
set(gca,'FontWeight','bold')
xlim([-1,6])

rotorefficiency = 10^(-1/mdlrotor.Coefficients.Estimate(2))
adaptiveefficiency = 10^(-1/mdladaptive.Coefficients.Estimate(2))


cycle40time = timescurve(41,:)/60;
cycle40avgs = mean(cycle40time);
std40time = std(cycle40time);
figure
hold on
mdltime = fitlm(log10(concentrations),averageadaptivetime/60);
plot([6,5,4,3,2,-1],mdltime.Coefficients.Estimate(1)+mdltime.Coefficients.Estimate(2)*[6,5,4,3,2,-1],LineWidth=2,Color="#A2142F",LineStyle='--')
eb3 = errorbar(log10(concentrations), averageadaptivetime/60, stdadaptivetime/60, 'ob',"LineWidth",1,Color='red');
ylabel("Quantification Time (min)")
xlabel("Log_{10} [Copies/\muL]")
set(gca,'FontName' , 'Arial')
set(gca,'FontSize',20)
set(gca,'FontWeight','bold')
xlim([-1,6])


figure
hold on
plot(cyclerotor,rotorinterp(:,1:4),Color="#0072BD",LineWidth=2)
plot(cyclerotor,rotorinterp(:,5:8),Color="#D95319",LineWidth=2)
plot(cyclerotor,rotorinterp(:,9:12),Color="#EDB120",LineWidth=2)
plot(cyclerotor,rotorinterp(:,13:16),Color="#7E2F8E",LineWidth=2)
plot(cyclerotor,rotorinterp(:,17:20),Color="#77AC30",LineWidth=2)
plot(cyclerotor,rotorinterp(:,21:24),Color="#4DBEEE",LineWidth=2)
plot(cyclerotor,rotorinterp(:,25:28),Color="#A2142F",LineWidth=2)
plot(cyclerotor,rotorinterp(:,29:32),Color="#FFC0CB",LineWidth=2)
ylabel("Fluorescence (Normalized)")
xlabel("Cycle")
set(gca,'FontName' , 'Arial')
set(gca,'FontSize',20)
set(gca,'FontWeight','bold')
xlim([0,40])
ylim([-1,100])
yline(threshold,LineStyle="--",LineWidth=2, color="black")

figure
hold on
plot(cycleadaptive,adaptiveinterp(:,1:4),Color="#0072BD",LineWidth=2)
plot(cycleadaptive,adaptiveinterp(:,5:8),Color="#D95319",LineWidth=2)
plot(cycleadaptive,adaptiveinterp(:,9:12),Color="#EDB120",LineWidth=2)
plot(cycleadaptive,adaptiveinterp(:,13:16),Color="#7E2F8E",LineWidth=2)
plot(cycleadaptive,adaptiveinterp(:,17:20),Color="#77AC30",LineWidth=2)
plot(cycleadaptive,adaptiveinterp(:,21:24),Color="#4DBEEE",LineWidth=2)
plot(cycleadaptive,adaptiveinterp(:,25:28),Color="#A2142F",LineWidth=2)
plot(cycleadaptive,adaptiveinterp(:,29:32),Color="#FFC0CB",LineWidth=2)
ylabel("Fluorescence (Normalized)")
xlabel("Cycle")
set(gca,'FontName' , 'Arial')
set(gca,'FontSize',20)
set(gca,'FontWeight','bold')
xlim([0,40])
ylim([-1,100])
yline(threshold,LineStyle="--",LineWidth=2, color="black")

%% Speed testing figure
timedata = xlsread("Speedtest for 1370 normalized.xlsx",1);
speedfluor = xlsread("Speedtest for 1370 normalized.xlsx",2);

cyclespeed = 1:0.01:40;
speedinterp = interp1(speedfluor(:,1),speedfluor(:,2:8),cyclespeed,"spline");
timeinterp = interp1(timedata(2:41,1),timedata(2:41,10:16),cyclespeed,"spline");
baseline = speedfluor(1:15,2:8);
threshold = 10*std(baseline(:));
speedindices = zeros(1,7);

for i = 1:7
    if(max(speedinterp(:,i)) > threshold)
            speedindices(i) = find(speedinterp(:,i) > threshold,1);
        else
            speedindices(i) = length(speedinterp(:,i));
    end
end
speedCts = cyclespeed(speedindices);

figure
hold on
plot(timeinterp(:,1:4)/60,speedinterp(:,1:4),Color="#0072BD",LineWidth=2)
plot(timeinterp(:,5:7)/60,speedinterp(:,5:7),Color="#FFC0CB",LineWidth=2)
ylabel("Fluorescence (Normalized)")
xlabel("Time (min)")
set(gca,'FontName' , 'Arial')
set(gca,'FontSize',20)
set(gca,'FontWeight','bold')
xlim([0,16])
ylim([-1,80])


speedavg = mean(speedCts(1:4));
speedstd = std(speedCts(1:4));