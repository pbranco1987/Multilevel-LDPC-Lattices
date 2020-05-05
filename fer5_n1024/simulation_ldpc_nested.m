% Authors: P. R. Branco da Silva and D. Silva
% Paper: Multilevel LDPC lattices with efficient encoding and decoding and a generalization of construction D'
% Journal: IEEE Transactions on Information Theory, vol. 65, no. 5, pp. 3246–3260, May 2019
% DOI: 10.1109/TIT.2018.2883119

% Simulation of an Original Construction D' 2-level LDPC lattice
% Constructed with EPEG [9] and decoded with multistage decoding
% n=1024, R0=0.0967, R1=0.9043

addpath('../Library')
n = 1024;
m1 = 98;
k1 = n-m1;
R1 = k1/n;
m0 = 925;
k0 = n-m0;
R0 = k0/n;
R = R0+R1;

dv1 = 3*ones(1,n);
dv0 = 6*ones(1,n);
Dv = [dv0; dv1];
M = [m0 m1];
A = EPEG(Dv,M);
H0 = A(1:m0,:);
H1 = A(1:m1,:);
H{1} = H0;
H{2} = H1;
H{3} = [];

VNR_dB =           [0    0.5   1    1.5   2.0   2.25  2.5  2.75   2.90  3.1133];
MaxFrameErrors =   [300  200   200  100   100   100   100  100    100   60];
VNR = 10.^(VNR_dB/10);
sigma = sqrt((2^(4-2*R))./(2*pi*exp(1).*VNR));
result_vector = [];

for i = 1:length(sigma)
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
    save('full_results_ldpc_nested.mat','result_vector');

end

% simplify .mat
load('full_results_ldpc_nested.mat');
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
save('results_ldpc_nested.mat','fer','fer0','fer1','fer2',...
    'ser','ser1','ser2','VNR_dB');
