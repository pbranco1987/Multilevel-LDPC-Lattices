% Authors: P. R. Branco da Silva and D. Silva
% Paper: Multilevel LDPC lattices with efficient encoding and decoding and a generalization of construction D'
% Journal: IEEE Transactions on Information Theory, vol. 65, no. 5, pp. 3246–3260, May 2019
% DOI: 10.1109/TIT.2018.2883119

% This script plots Fig. 4 of the paper (see caption for details)
% Design for P_e <= 10^-2 and n = 10000


addpath('../Library')

% line width
lw = 1.0;

%% Load data

% Generalized Construction D' LDPC (2-level, n=10000) with multistage decoding
% Level 0
results_0 = load('results_ldpc_cs_C0.mat');

% Level 1
results_1 = load('results_ldpc_cs_C1.mat');

% Level 2 (uncoded)
% WER is calculated theoretically
n = 10000;
R0 = 0.4094;
R1 = 0.9730;
R = R0+R1;
VNR_dB = results_1.VNR_dB;
vnr = 10.^(VNR_dB/10);
sigma = sqrt((2^(4-2*R))./(2*pi*exp(1).*vnr));

fer2 = [];
for sig = sigma
    fer2 = [fer2 pf_unc(sig/4, n)];    
end
results_2.VNR_dB = VNR_dB;
results_2.fer = fer2;

% Total WER
% Approximated by the union bound of the WERs of the individual levels
pad = zeros(1,(length(results_1.VNR_dB)-length(results_0.VNR_dB)));
fer0 = [results_0.fer pad];
fer1 = results_1.fer;
results.VNR_dB = results_1.VNR_dB;
results.fer = min(fer0+fer1+fer2,1);

% LDLC n=10000
% Kindly provided by Brian Kurkoski
data_ldlc = load('data_ldlc.mat');

% LDA n=10000
% Kindly provided by N. di Pietro and J. J. Boutros
data_lda = load('data_lda.mat');

%% Plotting

hold off;

% LDPC Gen. Const. D'
semilogy(results.VNR_dB,results.fer,'-','Marker','o','Color',[0,0,1],'linewidth',lw);
hold on;
semilogy(results_0.VNR_dB,results_0.fer,'--','Color',[0 0.5 1],'LineWidth',lw);
semilogy(results_1.VNR_dB,results_1.fer,'-.','Color',[0,0.5,0.5],'LineWidth',lw);
semilogy(results_2.VNR_dB,results_2.fer,'-','Color',[0.5 0.5 0.5],'LineWidth',lw);

% LDLC
semilogy(data_ldlc.VNR_dB,data_ldlc.fer,'-*','Color',[1,0,0],'linewidth',lw)

% LDA
semilogy(data_lda.VNR_dB,data_lda.fer,'-x','Color',[1,0.5,0],'linewidth',lw);

% Poltyrev limit
fer1 = get(gca,'ylim');
semilogy([0 0],fer1,'k--','LineWidth',2);

% Graph specifications
xlabel('VNR (dB)','Interpreter','Latex');
ylabel('Word Error Rate (WER)','Interpreter','Latex');
legend({'LDPC (Generalized D$''$)','Level 0 ($R_0 = 0.4094$)','Level 1 ($R_1 = 0.973$)','Level 2 (uncoded)',...
        'LDLC [29]',...
        'LDA [28]',...
        'Poltyrev limit'},...
        'units','pixels',...
        'Location','southwest',...
        'Interpreter','Latex');
% set(gca,'XTick',0:0.25:2);
axis([0 1.4 -inf inf]);
ax = gca;
ax.GridColor = 'k';
ax.GridAlpha = 0.18;
grid on;
set(gca,'XMinorGrid','off')
hold off;
ylim([1e-5-eps,1e0])
figure(gcf)

% export
set(gcf,'Color','w')
% uses export_fig library to export to pdf
%addpath('../Library/export_fig')
%export_fig('fer2-n10000.pdf');
