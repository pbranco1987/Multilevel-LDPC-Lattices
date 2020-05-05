function [result] = ferldpc_C0_GPU(H,varargin)
% Examples:
%   ferldpc(H,'SNR_dB',10)
%   ferldpc({H0,H1},'SNR_dB',10,'ModuloFrontEnd',true,'Scaled','true')


p = inputParser;
%addParameter(p,'MyString','option1',@(x)any(validatestring(x,{'option1','option2'})))
addParameter(p,'sigma_ch',[],@isnumeric)
addParameter(p,'SNR_dB',[],@isnumeric)
addParameter(p,'EbN0_dB',[],@isnumeric)
addParameter(p,'INV_dB',[],@isnumeric)
ChannelParameterSet = p.Parameters;
addParameter(p,'ModuloFrontEnd',false,@islogical)
addParameter(p,'Scaled',false,@islogical)
addParameter(p,'MaxFrameErrors',30,@isnumeric)
addParameter(p,'MaximumIterationCount',50,@isnumeric)
addParameter(p,'Verbose',false,@islogical)
addParameter(p,'ErrorPropagation',false,@islogical)
addParameter(p,'MaxIterations',0,@isnumeric)
addParameter(p,'Rate',0,@isnumeric)
addParameter(p,'CodeLength',0,@isnumeric)

parse(p,varargin{:})
for s = p.Parameters
    eval([s{1} ' = p.Results.' s{1} ';'])
end
ChannelParameter = '';
for s = ChannelParameterSet
    if ~isempty(eval(s{1}))
        if isempty(ChannelParameter)
            ChannelParameter = s{1};
        else
            error('Only one type of channel parameter can be specified.')
        end
    end
end
if false
    % never executes, prevents compilation errors in parfor
    ModuloFrontEnd = false;
    Scaled = false;
    Verbose = false;
    ErrorPropagation = false;
end

if ModuloFrontEnd
    % Multi-level code, real channel, modulo
    Decoder = comm.gpu.LDPCDecoder(H{1},... % H{1} is the parity-check matrix of C0
                               'OutputValue','Whole codeword',...
                               'MaximumIterationCount',MaximumIterationCount,...
                               'IterationTerminationCondition','Maximum iteration count');

    L = 3;
    R = Rate;
    n = CodeLength;
else
    % Single-level code, real channel, no modulo operation
    Decoder = comm.LDPCDecoder(...
        'ParityCheckMatrix',H,...
        'OutputValue','Whole Codeword',...
        'IterationTerminationCondition','Parity check satisfied',...
        'MaximumIterationCount',MaximumIterationCount);
    n = size(H,2);
    R = 1-size(H,1)/n;
    L = 1;
end

q = 2^L;
Constellation = (0:q-1).';
d = -mean(Constellation);
P = mean(abs(Constellation + d).^2);
spectral_efficiency = R;  % information bits per symbol
switch ChannelParameter
    case 'sigma_ch'
        SNR = P/sigma_ch^2;
    case 'INV_dB'
        INV = 10^(INV_dB/10);
        sigma_ch = sqrt(1/INV);
        SNR = P/sigma_ch^2;
    case 'SNR_dB'
        SNR = 10^(SNR_dB/10);
        sigma_ch = sqrt(P/SNR);
    case 'EbN0_dB'
        EbN0 = 10^(EbN0_dB/10);
        SNR = 2*spectral_efficiency*EbN0;
        sigma_ch = sqrt(P/SNR);
end

if ~ModuloFrontEnd    
    L = 1;
end

%Standard: MaxIterations = 0
if MaxIterations == 0
    MaxIterFlag = 0;
else
    MaxIterFlag =1;
end

nfer = 0; 
nser = 0; 
frames = 0;
symbols = 0;
display(sigma_ch)

% while none of the following conditions take place:
    % frame count > max number of frames
    % error count > max number of errors

while ~( (MaxIterFlag == 1 && frames > MaxIterations) || (nfer > MaxFrameErrors))

    
    % variables internal to parfor
   
    infer = 0;
    inser = 0;
    
    iframes = 0;
    isymbols = 0;
    
    GPU_Len = 75;
    for par_loop = 1:50
        c = zeros(GPU_Len,n); % must have c == mod(c,q)
        c = gpuArray(c);
        if ModuloFrontEnd
            s = randi([0 q-1],GPU_Len,n);    
        else
            s = zeros(1,n);
        end
        x = mod(c + s, q) + d;
        z = sigma_ch*randn(GPU_Len,n);
        y = x + z;
        if Scaled
            alpha = 1/(1 + 1/SNR);
            noise_variance_factor = 1/(1 + 1/SNR);
        else
            alpha = 1;
            noise_variance_factor = 1;
        end
        r = alpha*y - s - d;
        sigmaeff = sigma_ch*sqrt(noise_variance_factor);

        if ModuloFrontEnd

            r0 = mod(r+1,2)-1;
            LLR = (1-2*abs(r0))./(2*sigmaeff.^2); % sigmaeff/2 --> C1  %%% faster 
            
            %LLR = LLR_mod2(ri,sigmaeff/2);  %%% slower
            c_hat = step(Decoder,LLR(:));
            err = ((c_hat)~=0); 
            err_matrix = reshape(err,[n,GPU_Len]); % each column is a received vector
            column_sum = sum(err_matrix); % sum of column elements
            
            inser = gather(inser + sum(err)); 
            infer = gather(infer + sum(column_sum~=0)); 
            iframes = gather(iframes + GPU_Len);
            isymbols = gather(isymbols + numel(err));
        else
            % Single-level code, real channel
            r0 = mod(r+1,2)-1;
    %             ri = r;
            LLR = (1 - 2*abs(r0))/(2*sigmaeff.^2);
            c_dec = step(Decoder,LLR(:))';
            err = (c_dec ~= 0);
            inser = inser + sum(sum(err));
            infer = infer + any(err>0);
            iframes = iframes + 1;
            isymbols = isymbols + numel(err);
        end
    end
    toc
    nfer = nfer + infer;
    nser = nser + inser;
    frames = frames + iframes;
    symbols = symbols + isymbols;
    
    fprintf('FER_0 = %d (%d/%d)\n', nfer/frames, nfer, frames); 
    fprintf('\n\n\n\n'); 

    kshow = 1e4; % division constant
    show_condition = frames/kshow;
    show_condition = mod(show_condition,1)==0;
    if ((infer > 0) || (show_condition))
        fprintf('FER_0 = %d (%d/%d)\n', nfer/frames, nfer, frames);
        fprintf('\n\n\n\n'); 
    end
    
    ksave = 1e4;
    save_condition = frames/ksave;
    save_condition = mod(save_condition,1)==0;
    if ((infer > 0) || (save_condition))
        save('provisional_result_c0.mat','sigma_ch','nfer','frames');
    end
end

result.fer = nfer./frames;
result.nfer = nfer;
result.nframes = frames;
result.ser = nser./symbols;
result.nser = nser;
result.nsymbols = symbols;
result.parameters = p.Results;
result.(ChannelParameter) = eval(ChannelParameter);
result.description = 'Error is only being measured for C0';

