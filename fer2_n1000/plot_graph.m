% Authors: P. R. Branco da Silva and D. Silva
% Paper: Multilevel LDPC lattices with efficient encoding and decoding and a generalization of construction D'
% Journal: IEEE Transactions on Information Theory, vol. 65, no. 5, pp. 3246–3260, May 2019
% DOI: 10.1109/TIT.2018.2883119

% This script plots Fig. 3 of the paper (see caption for details)
% Design for P_e <= 10^-2 and n = 1000


% line width
lw = 1.0;

%% Load data

% Generalized Construction D' LDPC (2-level, n=1000) with multistage decoding
results = load('results_ldpc_cs.mat');

% QC-LDPC n=1190
% Kindly provided by H. Khodaiemehr and M.-R. Sadeghi
data_qcldpc = load('data_qcldpc.mat');

% LDLC N=1000, alpha=0.8571, M=2, d=7, TTD Decoder - Quantized
% Obtained from simulation files kindly provided by Brian Kurkoski
data_ldlc = load('data_ldlc.mat');

% LDA n=1000
% Kindly provided by N. di Pietro and J. J. Boutros
data_lda = load('data_lda.mat');

% GLD N=1000, L=250, n=4, k=3, p=11, seedperm=17
% Kindly provided by N. di Pietro and J. J. Boutros
data_gld = load('data_gld.mat');

%% Plotting
hold off;

% QC-LDPC
semilogy(data_qcldpc.VNR_dB,data_qcldpc.fer,'-+','Color','m','linewidth',lw);
hold on;

% LDPC Gen. Const. D'
semilogy(results.VNR_dB,results.fer,'-','Marker','o','Color',[0,0,1],'linewidth',lw);
semilogy(results.VNR_dB,results.fer0,'--','Color',[0 0.5 1],'LineWidth',lw);
semilogy(results.VNR_dB,results.fer1,'-.','Color',[0,0.5,0.5],'LineWidth',lw);
semilogy(results.VNR_dB,results.fer2,'-','Color',[0.5 0.5 0.5],'LineWidth',lw);

% LDLC
semilogy(data_ldlc.VNR_dB, data_ldlc.fer, '-','Marker','*','Color',[1 0 0],'LineWidth',lw);

% LDA
semilogy(data_lda.VNR_dB, data_lda.fer,'-','Marker','x','Color',[1,0.5,0],'linewidth',lw);

% GLD
semilogy(data_gld.VNR_dB, data_gld.fer,'-','Marker','v','Color',[0.5 0 0],'LineWidth',lw);

% Poltyrev limit
y1 = get(gca,'ylim');
semilogy([0 0],y1,'k--','LineWidth',2);

% Graph specifications
xlabel('VNR (dB)','Interpreter','Latex');
ylabel('Word Error Rate (WER)','Interpreter','Latex');
legend({'QC-LDPC [37]'....
        'LDPC (Generalized D$''$)','Level 0 ($R_0 = 0.5$)','Level 1 ($R_1 = 0.978$)','Level 2 (uncoded)',...
        'LDLC [29]',...
        'LDA [28]',...
        'GLD [4]',...
        'Poltyrev limit'},...
        'units','pixels',...
        'Location','NorthEast',...
        'Location','SouthWest',...
        'Interpreter','Latex');
set(gca,'XTick',0:0.2:3);
axis([0 2. -inf inf]);
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
%export_fig('fer2-n1000.pdf');
