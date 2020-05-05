function [result] = ferldpc_C1_GPU(H,varargin)
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
    
    
    % Decoder = LDPCMultiStageDecoder(H,MaximumIterationCount); %%%%
    %%%% orginal line of code
    
    
    %%% Alterations (probably will need to be changed later on)
    Decoder = comm.gpu.LDPCDecoder(H{2},... % alteration // H{1} is the parity-check matrix of C0
                               'OutputValue','Whole codeword',...
                               'MaximumIterationCount',MaximumIterationCount,...
                               'IterationTerminationCondition','Maximum iteration count');
    
    
                           
    
    %%%% Original lines of code:
%     n = Decoder.n;
%     R = Decoder.Rsum;
%     L = Decoder.L;




    %%%% Alteration
    L = 3; % added line of code
    R = Rate; % added line of code
    n = CodeLength; % added line of code
    


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
Constellation = (0:q-1).';
d = -mean(Constellation);
P = mean(abs(Constellation + d).^2);


%   M = length(Constellation);  %%%% original line of code



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
    
%%%% Alteration
%%% Original Lines of code:

%     L = Decoder.L;
%     if Decoder.uncoded_level
%         coded_levels = L-1;
%     else
%         coded_levels = L;
%     end

else
    L = 1;
    
    %%%% Alteration
    %%%coded_levels = 1; % original line of code
end

%Padrao MaxIterations = 0
if MaxIterations == 0
    MaxIterFlag = 0;
else
    MaxIterFlag =1;
end


%Variaveis finais


%nfer_conjunta = 0;  %%% original line of code
%nser_conjunta = 0;  %%% original line of code

%nfer = zeros(1,L); %%% original line of code
%nser = zeros(1,L); %%% original line of code





%%%% Alterations
nfer = 0; % alteration
nser = 0; % alteration






frames = 0;
symbols = 0;
display(sigma_ch)



%k = 5; %%% original line of code
%min_nfer = min(nfer); %%% original line of code






% Enquanto nenhuma das condicoes acontece:
    % Numero de frames > max numero de frames
    % Numero de erros > max numero de erros
% while ~( (MaxIterFlag == 1 && frames > MaxIterations) || (min(nfer(1:coded_levels) > MaxFrameErrors)))
% while ~( (MaxIterFlag == 1 && frames > MaxIterations) || (min_nfer > MaxFrameErrors))    %%% original line of code


while ~( (MaxIterFlag == 1 && frames > MaxIterations) || (nfer > MaxFrameErrors))
    
    %error_propag = ErrorPropagation;  %%% original line of code
    
    
    
    %Variaveis internas do parfor
    
    
    
    %infer = zeros(1,L);  %%% original line of code
    %inser = zeros(1,L);  %%% original line of code
    
    
    %%%% alterations
    infer = 0;
    inser = 0;
    
    
    
    %infer_conjunta = 0; %%% original line of code
    %inser_conjunta = 0; %%% original line of code
    
    
    
    iframes = 0;
    isymbols = 0;
    
    
    %%%% change back to parfor after tests!!!!!
  
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
        
        
        %%%% Alteration:
        % c1 = zeros(1,n); % added line of code
        
        
        
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

            %%%% Alterations:
            %%% must use LLR still !!!!
            r1 = mod((r/2)+1,2)-1;
            LLR = (1-2*abs(r1))./(2*(sigmaeff/2).^2); % sigmaeff/2 --> C1  %%% faster 
            
            c_hat = step(Decoder,LLR(:)); % added line of code
            err = ((c_hat)~=0); % added line of code
            err_matrix = reshape(err,[n,GPU_Len]); %Cada coluna é um vetor recebido
            column_sum = sum(err_matrix); % Soma os elementos da coluna, deveriam ser zero
            
            inser = gather(inser + sum(err)); 
            infer = gather(infer + sum(column_sum~=0)); 
            
           
            
            
            
            
            iframes = gather(iframes + GPU_Len);
            isymbols = gather(isymbols + numel(err));
        else
            % Single-level code, real channel
            r1 = mod(r+1,2)-1;
    %             ri = r;
            LLR = (1 - 2*abs(r1))/(2*sigmaeff.^2);
            c_dec = step(Decoder,LLR(:))';
            err = (c_dec ~= c);
            inser = inser + sum(sum(err));
            infer = infer + any(err>0);
            
            



            
            iframes = iframes + GPU_Len;
            isymbols = isymbols + numel(err);
        end
    end
    nfer = nfer + infer;
    nser = nser + inser;
    frames = frames + iframes;
    symbols = symbols + isymbols;
    



    %%%% Alterations:
    k = 1e4; % division constant
    show_condition = frames/k;
    show_condition = mod(show_condition,1)==0;
    if (1)
        fprintf('FER_1 = %d (%d/%d)\n', nfer/frames, nfer, frames); % added line of code
        fprintf('\n\n\n\n'); % added line of code
    end
    
    
    ksave = 1e4;
    save_condition = frames/ksave;
    save_condition = mod(save_condition,1)==0;
    if ((infer > 0) || (save_condition))
        save('provisional_result_c1.mat','sigma_ch','nfer','frames');
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




%%%% Alteration
result.description = 'Error is only being measured for C1'; % added line of code

