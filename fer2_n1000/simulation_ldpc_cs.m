% Authors: P. R. Branco da Silva and D. Silva
% Paper: Multilevel LDPC lattices with efficient encoding and decoding and a generalization of construction D'
% Journal: IEEE Transactions on Information Theory, vol. 65, no. 5, pp. 3246–3260, May 2019
% DOI: 10.1109/TIT.2018.2883119

% Simulation of a Generalized Construction D' 2-level LDPC lattice
% Constructed with check splitting and decoded with multistage decoding
% n=1000, R0=0.5, R1=0.978 and a triangular construction with g=20

% Computes level 0, level 1, level 2 and the joint error probability

addpath('../Library')

n = 1000;
m0 = round(n*(1-0.5)); % R0 = 0.5
m1 = round(n*(1-0.978)); % R1 = 0.978

k0 = n-m0;
k1 = n-m1;

R1 = k1/n;
R0 = k0/n;
R = R0+R1;


VNR_dB         =  [0    0.25  0.5  0.75   1    1.25   1.3   1.356 1.383   1.4   1.6   1.7   1.85  1.9  ];
MaxFrameErrors =  [500  500   500  500    500  500    500   400   300     200   200   100   100   100];
VNR = 10.^(VNR_dB/10);
sigma = sqrt((2^(4-2*R))./(2*pi*exp(1).*VNR));

% triangular PEG and triangular check splitting with g=20
g = 20;
H1 = PEG_tri(3*ones(1,n),m1,g);
H0 = PEG_CS_tri(H1,m0);
H{1} = H0;
H{2} = H1;
H{3} = [];

result_vector = [];

for i = 1:length(sigma)
    rng(0); 
    [result] = ferldpc(H,'sigma',sigma(i),...
                               'ModuloFrontEnd',true,'Scaled',false,...
                               'Verbose',true,'MaxFrameErrors',MaxFrameErrors(i),...
                               'ErrorPropagation',false); 
    result.R0 = R0;
    result.R1 = R1;
    result.R = R;
    result.n = n;
    result.m0 = m0;
    result.m1 = m1;
    result.sigma = sigma(i);
    result.VNR_dB = VNR_dB(i);
    result_vector = [result_vector result];
    save('full_results_ldpc_cs.mat','result_vector');
end

% simplify .mat
load('full_results_ldpc_cs.mat');
n = length(result_vector);
fer = zeros(1,n); fer0 = zeros(1,n); fer1 = zeros(1,n); fer2 = zeros(1,n);
ser = zeros(1,n); ser0 = zeros(1,n); ser1 = zeros(1,n); ser2 = zeros(1,n);
VNR_dB = zeros(1,n);
for i=1:n
    fer(i) = result_vector(i).fer_joint;
    fer0(i) = result_vector(i).fer(1); fer1(i) = result_vector(i).fer(2); fer2(i) = result_vector(i).fer(3);
    ser(i) = result_vector(i).ser_joint;
    ser0(i) = result_vector(i).ser(1); ser1(i) = result_vector(i).ser(2); ser2(i) = result_vector(i).ser(3);
    VNR_dB(i) = result_vector(i).VNR_dB;
end
save('results_ldpc_cs.mat','fer','fer0','fer1','fer2',...
    'ser','ser1','ser2','VNR_dB');

