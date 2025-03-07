clear all

%% Extract time and fluorescence data from binary.
filelist = dir('Consistency data 808\**\*.bin');
timelist = dir('Consistency data 808\**\coeff*');
numbersamples = length(filelist);
maximacompilation = zeros(15,numbersamples);
minimacompilation = zeros(15,numbersamples);
coolingcompilation = zeros(15,numbersamples);
heatingcompilation = zeros(15,numbersamples);
for i = 1:length(filelist)
    handler = fopen(strcat(filelist(i).folder, "\", filelist(i).name));
    A = fread(handler,'double');
    frewind(handler);
    B = fread(handler, 'int64');
    frewind(handler);
    fseek(handler, 8, "bof");
    C = fread(handler,'double',8);
    frewind(handler);
    D = fread(handler,'int64',8);
    D = D - min(D);

    times = D(1:2:length(D))/1000;
    fluors = C(1:2:length(C));
    data = readmatrix(strcat(timelist(i).folder, "\", timelist(i).name));
    coolingcompilation(:,i) = data(4:2:32,5);
    heatingcompilation(:,i) = data(3:2:32,5);

    maxima = islocalmax(fluors, 'MinSeparation',15,'SamplePoints',times,'MinProminence',50);
    fluormax = fluors(maxima);
    gaps = diff(times);
    idx = find(gaps > 0.5);
    idx(fluors(idx) > mean(fluors)) = [];
    fluormin = fluors(idx);


    maximacompilation(:,i) = fluormax(2:16);
    minimacompilation(:,i) = fluormin(2:16);
end


averagemaxima = mean(maximacompilation);
stddevmaxima = std(maximacompilation);
avgmax = mean(averagemaxima);

averageminima = mean(minimacompilation);
stddevminima = std(minimacompilation);
avgmin = mean(averageminima);
%% Grab cycling values from the fitting data.

heatingvals = (heatingcompilation(1:15,:)-mean(heatingcompilation(1:15,:)))/1000;
coolingvals = (coolingcompilation(1:15,:)-mean(coolingcompilation(1:15,:)))/1000;
avgheat = mean(mean(heatingcompilation(1:15,:)))/1000;
avgcool = mean(mean(coolingcompilation(1:15,:)))/1000;

%% Create the histograms.

minimadeviation = minimacompilation - averageminima;
maximadeviation = maximacompilation -averagemaxima;

binedges = -40:2:40;

figure
hold on
histogram(minimadeviation,binedges,'Normalization','percentage')
histogram(maximadeviation,binedges,'Normalization','percentage')
ylabel("Frequency (%)")
xlabel("Deviation (mV)")
ylim([0,25])
stddeviationannealing = std(minimadeviation(:));
stddeviationmelting = std(maximadeviation(:));
set(gca,'FontName' , 'Arial')
set(gca,'FontSize',20)
set(gca,'FontWeight','bold')


binedges2 = -2:0.1:2;

figure
hold on
histogram(coolingvals,binedges2,'Normalization','percentage')
histogram(heatingvals,binedges2,'Normalization','percentage')
ylabel("Frequency (%)")
xlabel("Deviation (s)")
stddeviationcooling = std(coolingvals(:));
stddeviationheating = std(heatingvals(:));
set(gca,'FontName' , 'Arial')
set(gca,'FontSize',20)
set(gca,'FontWeight','bold')

minfactored = stddeviationannealing/avgmin;
maxfactored = stddeviationmelting/avgmax;

outliersminima=minimadeviation(abs(minimadeviation)>40);
percentoutliermin = length(outliersminima)/length(minimadeviation(:))*100
outliersmaxima=maximadeviation(abs(maximadeviation)>40);
percentoutliermax = length(outliersmaxima)/length(minimadeviation(:))*100
outlierscoolingtime=coolingvals(abs(coolingvals)>2);
percentoutliercool = length(outlierscoolingtime)/length(minimadeviation(:))*100
outliersheatingtime=heatingvals(abs(heatingvals)>2);
percentoutlierheat = length(outliersheatingtime)/length(minimadeviation(:))*100