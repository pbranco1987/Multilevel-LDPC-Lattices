classdef LDPCMultiStageDecoder
    properties (SetAccess = immutable)
        H
        MaximumIterationCount
        L
        Rsum
        n
        uncoded_level
    end
    properties (GetAccess = private, SetAccess = immutable)
        Decoder
        k
        R
    end
    methods
        function obj = LDPCMultiStageDecoder(H,MaximumIterationCount)
% LDPCMultiStageDecoder Multistage decoder for multilevel binary LPDC code
%   DEC = LDPCMultiStageDecoder creates a decoder System object, DEC,
%   for a multilevel binary low-density parity-check (LDPC) code. This
%   object performs multistage LDPC decoding based on the specified
%   parity-check matrices, one for each level, which are assumed to have
%   the same number of columns.
% 
%   DEC = LDPCMultiStageDecoder(PARITY_ARRAY) sets the H
%   property to PARITY_ARRAY, which must be a cell array consisting of the 
%   parity-check matrices of each level. The levels are decoded in the same
%   order as PARITY_ARRAY.
%
%   DEC = LDPCMultiStageDecoder(PARITY_ARRAY, MAXITERS) sets the 
%   MaximumIterationCount property to MAXITERS, which in turn specifies 
%   the maximum number of decoding iterations for each individual decoder.
%
%   LDPCMultiStageDecoder methods:
% 
%   step                    - Decode input data (click for help)
%
%   LDPCMultiStageDecoder properties:
%
%   H                       - Parity-check matrices for each level
%   MaximumIterationCount   - Maximum number of decoding iterations

            if nargin >= 1
                if ~iscell(H)
                    H = {H};
                end
                for i = 1:length(H)
                    if isempty(H{i})
                        obj.uncoded_level = 1;
                    else
                        obj.uncoded_level = 0;
                    end
                    H{i} = sparse(H{i});
                end
                obj.H = H;
            else
                obj.H = {dvbs2ldpc(1/2)};
            end
            if nargin >= 2
                obj.MaximumIterationCount = MaximumIterationCount;
            else
                obj.MaximumIterationCount = 50;
            end
            obj.n = size(obj.H{1},2);
            obj.L = length(obj.H);

             
            for i = 1:obj.L-obj.uncoded_level
                if size(obj.H{i},2) ~= obj.n
                    error('Parity-check matrices must have the same number of columns.')
                end
                obj.Decoder{i} = comm.LDPCDecoder(...
                    'ParityCheckMatrix',obj.H{i},...
                    'OutputValue','Whole Codeword',...
                    'IterationTerminationCondition','Parity check satisfied',...
                    'MaximumIterationCount',obj.MaximumIterationCount);
                obj.k{i} = obj.n - size(obj.H{i},1);
                obj.R{i} = obj.k{i}/obj.n;
            end
            obj.Rsum = sum(cell2mat(obj.R));
        end
        function [c,ci] = step(obj, r, sigma, ErrorPropagation)
% step      Multistage decoding of multilevel binary LPDC code 
%   c = step(DEC, r, sigma) decodes a received word, r, assuming noise
%   variance sigma.^2, where DEC is an LDPCMultiStageDecoder system object.
%   The input r is interpreted as obtained from r = c + z, where c is a
%   codeword with entries in the set {0, 1, 2, ..., 2^L-1}, L is the number
%   of levels, and z is a zero-mean Gaussian noise vector. If sigma is a
%   scalar, then z is interpreted as i.i.d. with variance sigma^2. If sigma
%   is a vector with the same size as r, then z is interpreted as having
%   covariance matrix diag(sigma.^2).

            ci = zeros(obj.L,length(r)); 
            c = zeros(size(r));
            % All levels, except for the last, are coded in the same way
            for i = 1:obj.L-1 
                ri = mod(r+1,2)-1;
                LLR = (1-2*abs(ri))./(2*sigma.^2);
                ci(i,:) = step(obj.Decoder{i},LLR(:));
                c = c + 2^(i-1)*ci(i,:);
                if ErrorPropagation               
                    r = (r - ci(i,:))/2; %Com propagação de erro
                else
                    r = r/2;               %Sem propagação de erro e com palavra codigo zero
                end
                sigma = sigma/2;
            end
            %O ultimo nível é decodificado de maneira diferente, dependendo
            %se o último nível é codificado ou não
            if obj.uncoded_level
                ri = mod(r+1,2)-1;
                LLR = (1-2*abs(ri))./(2*sigma.^2);
                ci(i+1,:) = (LLR < 0);
                c = c + 2^(i)*ci(i+1,:);
            else
                ri = mod(r+1,2)-1; %Com ou sem mod 
                LLR = (1-2*abs(ri))./(2*sigma.^2);
                ci(i+1,:) = step(obj.Decoder{i+1},LLR(:));
                c = c + 2^(i)*ci(i+1,:);
            end
        end
    end
end
