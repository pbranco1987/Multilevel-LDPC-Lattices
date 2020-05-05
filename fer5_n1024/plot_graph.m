% Authors: P. R. Branco da Silva and D. Silva
% Paper: Multilevel LDPC lattices with efficient encoding and decoding and a generalization of construction D'
% Journal: IEEE Transactions on Information Theory, vol. 65, no. 5, pp. 3246–3260, May 2019
% DOI: 10.1109/TIT.2018.2883119

% This script plots Fig. 2 of the paper (see caption for details)
% Design for P_e <= 10^-5 and n = 1024


addpath('../Library')

% line width
lw = 1.0;

%% Load data

% Generalized Construction D' LDPC (2-level, n=1000) with multistage decoding
% Level 0
results_0 = load('results_ldpc_cs_C0.mat');

% Level 1
results_1 = load('results_ldpc_cs_C1.mat');

% Level 2 (uncoded)
% WER is calculated theoretically
n = 1024;
m0 = 788;
m1 = 103;
R0 = 1-(m0/n); % R0 = 0.23
R1 = 1-(m1/n); % R1 = 0.90
R = R0+R1;
VNR_dB = results_1.VNR_dB;
fer2 = zeros(size(VNR_dB));
for i=1:length(VNR_dB)
    vnr = 10^(VNR_dB(i)/10);
    sigma = sqrt((2^(4-2*R))/(2*pi*exp(1)*vnr));
    fer2(i) = pf_unc(sigma/4,n);
end
results_2.VNR_dB = VNR_dB;
results_2.fer = fer2;

% Total WER
% Approximated by the union bound of the WERs of the individual levels
results.VNR_dB = results_1.VNR_dB;
fer0 = [results_0.fer zeros(1,length(results_1.fer)-length(results_0.fer))];
fer1 = results_1.fer;
results.fer = min(fer0+fer1+fer2,1);

% Original Construction D' LDPC (2-level, n=1024) with multistage decoding
results_ldpc_nested = load('results_ldpc_nested.mat');

% Polar (2-level, n=1024)
% Kindly provided by Cong Ling
data_polar = load('data_polar.mat');

% Original Construction D' LDPC (2-level, n=1000) with joint min-sum decoding
% Extracted from [9]
data_ldpc_minsum = load('data_ldpc_minsum.mat');

%% Plotting

hold off;

% LDPC Orig. Const. D' joint min-sum
semilogy(data_ldpc_minsum.VNR_dB,data_ldpc_minsum.fer,'-+','Color','m','linewidth',lw);
hold on;

% LDPC Orig. Const. D' MSD
semilogy(results_ldpc_nested.VNR_dB,results_ldpc_nested.fer,'-','Marker','s','Color',[0,0.5,0],'LineWidth',lw);

% Polar
semilogy(data_polar.VNR_dB,data_polar.fer,'-*','Color','r','LineWidth',lw);

% LDPC Gen. Const. D' MSD
semilogy(results.VNR_dB,results.fer,'-','Marker','o','Color',[0 0 1],'LineWidth',lw);
semilogy(results_0.VNR_dB,results_0.fer,'--','Color',[0 0.5 1],'LineWidth',lw);
semilogy(results_1.VNR_dB,results_1.fer,'-.','Color',[0,0.5,0.5],'LineWidth',lw);
semilogy(results_2.VNR_dB,results_2.fer,'-','Color',[0.5 0.5 0.5],'LineWidth',lw);

% Poltyrev Limit
y1 = get(gca,'ylim');
semilogy([0 0],y1,'k--','LineWidth',2);

% Graph specifications
xlabel('VNR (dB)','Interpreter', 'LaTex');
ylabel('Word Error Rate (WER)','Interpreter', 'LaTex');
legend({'LDPC / joint min-sum [9]',...
        'LDPC',...
        'Polar [13]',... 
        'LDPC (generalized D$''$)','Level 0 ($R_0 = 0.23$)','Level 1 ($R_1 = 0.90$)','Level 2 (uncoded)',...
        'Poltyrev limit'},...
        'units','pixels',...
        'Location','SouthWest',...
        'Interpreter', 'LaTex');
axis([0 3.5 -inf inf]);
ax = gca;
ax.GridColor = 'k';
ax.GridAlpha = 0.18;
grid on;
set(gca,'XMinorGrid','off')
hold off;
ylim([1e-7-eps,1e0])
figure(gcf)

% export
set(gcf,'Color','w');
% uses export_fig library to export to pdf
%addpath('../Library/export_fig')
%export_fig('fer5-n1024.pdf');
