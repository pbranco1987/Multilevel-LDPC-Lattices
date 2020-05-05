function [Arref,varargout] = gf2rref(A,varargin)
% A_ech is the matrix in reduced row echelon form
% r is the rank of matrix A

% A = unique(A,'rows','sorted');
% A = flipud(A);
% assert(isa(A, 'double'));
% assert(all(size(A) <= [1024 1024]));
% [m,n] = size(A);
% pivot_counter = 1;
% if m~=n
%     jdx = 1;
%     while pivot_counter <= m
%         indices = find(A(:,jdx)); % rows to be cleared
%         indx = indices(indices>=pivot_counter);
%         if ~isempty(indx)
%             ind = indices(indices~=indx(1));
%             A(ind,:) = xor(A(ind,:),repmat(A(indx(1),:),length(ind),1));
%             A([pivot_counter indx(1)],:) = A([indx(1) pivot_counter],:); %swap first row with pivot
%             pivot_counter = pivot_counter + 1;
%         end
%         jdx = jdx + 1;
%     end
% elseif m==n

Arref = logical(A);
pivot_counter = 1;
[m,n] = size(Arref);
for jdx = 1:n
    indices = find(Arref(:,jdx)); % rows to be cleared
    indx = indices(indices>=pivot_counter);
    if ~isempty(indx)
        ind = indices(indices~=indx(1));
        Arref(ind,:) = xor(Arref(ind,:),repmat(Arref(indx(1),:),length(ind),1));
        Arref([pivot_counter indx(1)],:) = Arref([indx(1) pivot_counter],:); %swap first row with pivot
        pivot_counter = pivot_counter + 1;
    end
    if pivot_counter>m
        break
    end
end

if nargin==2
    if strcmp(varargin{1},'MakeInvertible')
        [~,indep_cols] = sort(sum(Arref,1));
        Aind = A(:,indep_cols);
        %         lindep_columns = sum(Arref,1);
        %         lindep_columns = (1:n).*(lindep_columns~=1);
        %         dep_inds = (lindep_columns<=m)&(lindep_columns>0);
        %         dep_cols = lindep_columns(dep_inds);
        %         numdepcols = length(dep_cols);
        %         ind_cols = zeros(1,numdepcols);
        %         for ii = 1:numdepcols
        %             ind_cols(ii) = find(Arref(dep_cols(ii),:),1,'first');
        %         end
        %         A(:,dep_cols) = A(:,ind_cols);
        %         Aind = A;
    end
end

switch nargout
    case 2
        A_rank = sum(sum(Arref,2)~=0);
        if (nargin==2) && (exist('Aind','var')==1)
            varargout{1} = Aind;
        else
            varargout{1} = A_rank;
        end
    case 3
        A_rank = sum(sum(Arref,2)~=0);
        varargout{1} = A_rank;
        if (nargin==2) && (exist('Aind','var')==1)
            varargout{2} = Aind;
        else
            lindep_columns = sum(Arref,1);
            lindep_columns = (1:n).*(lindep_columns~=1);
            varargout{2} = lindep_columns(lindep_columns~=0);
        end
    case 4
        A_rank = sum(sum(Arref,2)~=0);
        varargout{1} = A_rank;
        lindep_columns = sum(Arref,1);
        lindep_columns = (1:n).*(lindep_columns~=1);
        varargout{2} = lindep_columns(lindep_columns~=0);
        if (nargin==2) && (exist('Aind','var')==1)
            varargout{3} = Aind;
        else
            error('Too many output variables');
        end
end