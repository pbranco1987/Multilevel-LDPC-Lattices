% Authors: P. R. Branco da Silva and D. Silva
% Paper: Multilevel LDPC lattices with efficient encoding and decoding and a generalization of construction D'
% Journal: IEEE Transactions on Information Theory, vol. 65, no. 5, pp. 3246–3260, May 2019
% DOI: 10.1109/TIT.2018.2883119

% Simulation of a Generalized Construction D' 2-level LDPC lattice
% Constructed with check splitting and decoded with multistage decoding
% n=10000, R0=0.4094, R1=0.9730 and a triangular construction with g=22

% Computes level 1 error probability

addpath('../Library')

n = 10000;
m0 = round(n*(1-0.4094)); % R0 = 0.4094
m1 = round(n*(1-0.9730)); % R1 = 0.9730
k0 = n-m0;
k1 = n-m1;

R1 = k1/n;
R0 = k0/n;
R = R0+R1;


VNR_dB         =  [0   0.2   0.4  0.5  0.6  0.7   0.8  0.884  0.9   1    1.1  1.2  1.3  1.4   1.5];
MaxFrameErrors =  [500 500   500  500  500  500   500  300    200   100  100  100  100  200   100];
VNR = 10.^(VNR_dB/10);
sigma = sqrt((2^(4-2*R))./(2*pi*exp(1).*VNR));

% Loads the H matrix constructed on the level 0 script
load('H.mat');
result_vector = [];

for i = 1:length(sigma)
    display(VNR_dB(i));
    rng(0); 
    [result] = ferldpc_C1_GPU(H,'sigma',sigma(i),...
                               'ModuloFrontEnd',true,'Scaled',false,...
                               'Verbose',true,'MaxFrameErrors',MaxFrameErrors(i),...
                               'ErrorPropagation',false,...
                               'Rate',R,'CodeLength',n);
    result.R0 = R0;
    result.R1 = R1;
    result.R = R;
    result.n = n;
    result.m0 = m0;
    result.m1 = m1;
    result.sigma = sigma(i);
    result.VNR_dB = VNR_dB(i);
    result_vector(i) = result;
    save('full_results_ldpc_cs_C1.mat','result_vector');
end

% processing to a simplified .mat
load('full_results_ldpc_cs_C1.mat');
n = length(result_vector);
VNR_dB = zeros(1,n);
fer = zeros(1,n);
ser = zeros(1,n);
for i=1:n
    VNR_dB(i) = result_vector(i).VNR_dB;
    fer(i) = result_vector(i).fer;
    ser(i) = result_vector(i).ser;
end
save('results_ldpc_cs_C1.mat','VNR_dB','ser','fer');
