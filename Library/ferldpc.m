function [result] = ferldpc(H,varargin)
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
    Decoder = LDPCMultiStageDecoder(H,MaximumIterationCount);
    n = Decoder.n;
    R = Decoder.Rsum;
    L = Decoder.L;
else
    % Single-level code, real channel, no modulo
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
Constellation = [0:q-1].';
d = -mean(Constellation);
P = mean(abs(Constellation + d).^2);
M = length(Constellation);

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

if ModuloFrontEnd    
    L = Decoder.L;
    if Decoder.uncoded_level
        coded_levels = L-1;
    else
        coded_levels = L;
    end
else
    L = 1;
    coded_levels = 1;
end

%Standard: MaxIterations = 0
if MaxIterations == 0
    MaxIterFlag = 0;
else
    MaxIterFlag =1;
end

%End variables
nfer_joint = 0;
nser_joint = 0;
nfer = zeros(1,L);
nser = zeros(1,L);
frames = 0;
symbols = 0;
display(sigma_ch)

% While none of the following conditions takes place:
    % Frame count > max number of frames
    % Error count > max number de errors
while ~( (MaxIterFlag == 1 && frames > MaxIterations) || (min(nfer(1:coded_levels) > MaxFrameErrors)))    
    
    error_propag = ErrorPropagation;
    
    % Internal errors to parfor
    infer = zeros(1,L);
    inser = zeros(1,L);
    infer_joint = 0;
    inser_joint = 0;
    iframes = 0;
    isymbols = 0;
    
    parfor par_loop = 1:1000
        c = zeros(1,n); % must have c == mod(c,q)
        if ModuloFrontEnd
            s = randi([0 q-1],1,n);
        else
            s = zeros(1,n);
        end
        x = mod(c + s, q) + d;
        z = sigma_ch*randn(1,n);
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
            % Multi-level code, real channel, modulo
            [c_dec, ci_dec] = step(Decoder, r, sigmaeff, error_propag);
            calc_ser_joint = zeros(1,n);  % aux variable for the joint ser
            inser_temp = zeros(1,L);
            infer_temp = zeros(1,L);
            for level=1:L
                err = (ci_dec(level,:) ~= c);
                calc_ser_joint = or(calc_ser_joint, err);  % if any level makes an error at the i-th position, then the i-th symbol is wrong 
                inser_temp(level) = inser_temp(level)  + sum(err);
                infer_temp(level) = infer_temp(level) + any(err);

            end

            inser = inser + inser_temp;
            infer = infer + infer_temp;

            infer_joint = infer_joint + any(c_dec~=0);
            inser_joint = inser_joint + sum(calc_ser_joint);
            iframes = iframes + 1;
            isymbols = isymbols + numel(err);
        else
            % Single-level code, real channel
            ri = mod(r+1,2)-1;
    %             ri = r;
            LLR = (1 - 2*abs(ri))/(2*sigmaeff.^2);
            c_dec = step(Decoder,LLR(:))';
            err = (c_dec ~= c);
            inser = inser + sum(sum(err));
            infer = infer + any(err>0);
            
            infer_joint = infer_joint + any(err>0);
            inser_joint = inser_joint + sum(sum(err));
            
            iframes = iframes + 1;
            isymbols = isymbols + numel(err);
        end
    end
    nfer = nfer + infer;
    nser = nser + inser;
    frames = frames + iframes;
    symbols = symbols + isymbols;
    nfer_joint = nfer_joint + infer_joint;
    nser_joint = nser_joint + inser_joint;

    fprintf('joint FER = %d (%d/%d)\n', nfer_joint/frames, nfer_joint, frames);
    for nivel = 1:L
        fprintf('Level %d: FER = %d  (%d/%d)\n',nivel-1, nfer(nivel)/frames, nfer(nivel), frames);
    end
    fprintf('\n\n\n\n');
    
end

result.fer_joint = nfer_joint./frames;
result.fer = nfer./frames;
result.nfer = nfer;
result.nfer_joint = nfer_joint;
result.nframes = frames;

result.ser_joint = nser_joint./symbols;
result.ser = nser./symbols;
result.nser = nser;
result.nser_joint = nser_joint;
result.nsymbols = symbols;

result.parameters = p.Results;
result.(ChannelParameter) = eval(ChannelParameter);
