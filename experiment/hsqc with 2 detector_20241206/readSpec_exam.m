% read real spectrum from topspin exported txt
%
% Mengjia He, 2024.10.25

% topspin examdata/20241206, glycine at detector 1, glucose at detector 2
% normal                        spec-80               spec-14     
% GC                            spec-82               spec-16
% SL                            spec-84               spec-63

close all; clearvars; clc;

% specify sample, sequence and nuclei
sample = 'glucose';
seqName = 'HSQC';
nuclei = '13C';

% load FID
fileName = append('spec-',sample,'-',seqName,'-',nuclei,'.txt');
spec = readTopSpinTxt('spec-14-996.txt');
spec = [spec;readTopSpinTxt('spec-16-996.txt')];
spec = [spec;readTopSpinTxt('spec-63-996.txt')];


% define frequency axis
switch nuclei
    case '1H'
        
        SW = 10;
        switch sample
            case 'glucose'; freq_shift = 0; xlimts = [2,6]; freq_ref = 1;
                spec = spec/1e8;
            case 'glycine'; freq_shift = 0; xlimts = [2,6]; freq_ref = 1.5;
                spec = spec/3e8;
        end
    case '13C'
        
        SW = 150;
        switch sample
            case 'glucose'; freq_shift = 0; xlimts = [50,110]; freq_ref = -4;
                spec = spec/1e8;
            case 'glycine'; freq_shift = 0; xlimts = [20,60]; freq_ref = 65;
                spec = spec/3e8;

        end
end
freq_axis = linspace(0,SW,size(spec,2))-freq_ref;

% visulize data
close all;
fig1 = figure('Name',fileName);
fig1.Position = [100 100 400 180];
plot(freq_axis, spec(1,:),'LineWidth',0.1,'color','k'); hold on;
plot(freq_axis+freq_shift, spec(2,:),'LineWidth',0.1,'color',lineColor('lightorange')); hold on;
plot(freq_axis+freq_shift*2, spec(3,:),'LineWidth',0.1,'color',lineColor('lightblue')); 
axis tight;
xlim(xlimts);
ylim([-0.1,5]);
yticks([]);
set(gca, 'XDir','reverse');
xlabel(append(nuclei,' Chemical shift, ppm'));
set(gca,'FontSize',11,'FontName','Arial','linewidth',1);
box off;
set(get(gca,'YAxis'),'visible','off');
% lgd = legend('Normal', 'GC', 'SL');
% set(lgd, 'Box', 'off');        
% set(lgd, 'Color', 'none');   
