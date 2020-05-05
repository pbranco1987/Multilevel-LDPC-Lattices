% Authors: P. R. Branco da Silva and D. Silva
% Paper: Multilevel LDPC lattices with efficient encoding and decoding and a generalization of construction D'
% Journal: IEEE Transactions on Information Theory, vol. 65, no. 5, pp. 3246–3260, May 2019
% DOI: 10.1109/TIT.2018.2883119

% Simulation of a Generalized Construction D' 2-level LDPC lattice
% Constructed with check splitting and decoded with multistage decoding
% n=1024, R0=0.23 and R1=0.90 and a triangular construction with g=22

% Computes level 0 error probability

addpath('../Library')

n = 1024;

m1 = 103; % k1 = 921 / R1 = 0.8994
k1 = n-m1; 
R1 = k1/n;
m0 = 788; % k0 = 236 / R0 = 0.2305
k0 = n-m0;
R0 = k0/n;
R = R0+R1;
fprintf('\r\n Total rate = %d',R);


g = 22;
fprintf('\r\n Gap = %d \r\n',g);

H1 = PEG_tri(3*ones(1,n),m1,g);
[hasgirth4_1,~] = girth4(H1);
fprintf('\r\n Size of H1: ');
disp(size(H1));
if hasgirth4_1
    fprintf('\r\n Has girth 4\r\n');
    return
else
    fprintf('\r\n Does not have girth 4\r\n');
end
H0 = PEG_CS_tri(H1,m0);
[hasgirth4_0,~] = girth4(H0);
fprintf('\r\n Size of H0: ');
disp(size(H0));
if hasgirth4_0
    fprintf('\r\n Has girth 4\r\n');
    return
else
    fprintf('\r\n Does not have girth 4\r\n');
end
H{1} = H0;
H{2} = H1;
H{3} = [];

disp(H);

% VNR_dB 2.6113 is where the polar lattices attain FER = 1e-6
% (where the uncoded level reaches FER = 1e-6)
VNR_dB         = [0    0.5   1   1.5  1.7   1.9   2.1   2.3   2.5];
disp(VNR_dB);
MaxFrameErrors = 100;
disp(MaxFrameErrors);
VNR = 10.^(VNR_dB/10);
sigma = sqrt((2^(4-2*R))./(2*pi*exp(1).*VNR));

result_vector = [];
for i = 1:length(sigma)
    rng(0); 
    [result] = ferldpc_C0_GPU(H,'sigma',sigma(i),...
                               'ModuloFrontEnd',true,'Scaled',false,...
                               'Verbose',true,'MaxFrameErrors',MaxFrameErrors(i),...
                               'ErrorPropagation',false,...
                               'Rate',R,...
                               'CodeLength',n);
                           
    result.R0 = R0;
    result.R1 = R1;
    result.R = R;
    result.n = n;
    result.m0 = m0;
    result.m1 = m1;
    result.sigma = sigma(i);
    result.VNR_dB = VNR_dB(i);


    result_vector = [result_vector result];
    disp(result);
    save('full_results_ldpc_cs_C0.mat','result_vector');

end

% simplify .mat
load('full_results_ldpc_cs_C0.mat');
n = length(result_vector);
fer = zeros(1,n); 
ser = zeros(1,n); 
VNR_dB = zeros(1,n);
for i=1:n
    fer(i) = result_vector(i).fer;
    VNR_dB(i) = result_vector(i).VNR_dB;
end
save('results_ldpc_cs_C0.mat','fer','VNR_dB');
