clear all

%%
filelist = dir('1370 nm raw data dilution curve\**\*.bin');
timelist = dir('1370 nm raw data dilution curve\**\coeff*');

maximacompilation = zeros(15,32);
minimacompilation = zeros(15,32);
coolingcompilation = zeros(45,32);
heatingcompilation = zeros(45,32);
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
    coolingcompilation(:,i) = data(2:2:length(data(:,6)),6);
    heatingcompilation(:,i) = data(1:2:length(data(:,6)),6);

    maxima = islocalmax(fluors, 'MinSeparation',15,'SamplePoints',times,'MinProminence',50);
    fluormax = fluors(maxima);
    gaps = diff(times);
    idx = find(gaps > 0.5);
    idx(fluors(idx) > mean(fluors)) = [];
    fluormin = fluors(idx);


    maximacompilation(:,i) = fluormax(2:16);
    minimacompilation(:,i) = fluormin(2:16);
end
%%
averagemaxima = mean(maximacompilation);
stddevmaxima = std(maximacompilation);
avgmax = mean(averagemaxima);

averageminima = mean(minimacompilation);
stddevminima = std(minimacompilation);
avgmin = mean(averageminima);

heatingvals = (heatingcompilation(2:16,:)-mean(heatingcompilation(2:16,:)))/1000;
coolingvals = (coolingcompilation(2:16,:)-mean(coolingcompilation(2:16,:)))/1000;
avgheat = mean(mean(heatingcompilation(2:16,:)))/1000;
avgcool = mean(mean(coolingcompilation(2:16,:)))/1000;

minimadeviation = minimacompilation - averageminima;
maximadeviation = maximacompilation -averagemaxima;


binedges = -40:2:40;

figure
hold on
histogram(minimadeviation,binedges,'Normalization','percentage')
histogram(maximadeviation,binedges,'Normalization','percentage')
ylabel("Frequency (%)")
xlabel("Deviation (mV)")
stddeviationannealing = std(minimadeviation(:));
stddeviationmelting = std(maximadeviation(:));
set(gca,'FontName' , 'Arial')
set(gca,'FontSize',20)
set(gca,'FontWeight','bold')


binedges = -2:0.1:2;
figure
hold on
histogram(coolingvals, binedges,'Normalization','percentage')
histogram(heatingvals,binedges,'Normalization','percentage')
ylabel("Frequency (%)")
xlabel("Deviation (s)")
ylim([0,40])
stddeviationcooling = std(coolingvals(:));
stddeviationheating = std(heatingvals(:));
set(gca,'FontName' , 'Arial')
set(gca,'FontSize',20)
set(gca,'FontWeight','bold')

