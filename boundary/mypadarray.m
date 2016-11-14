function b = mypadarray(varargin)


[a, method, padSize, padVal, direction] = ParseInputs(varargin{:});

if isempty(a),% treat empty matrix similar for any method

   if strcmp(direction,'both')
      sizeB = size(a) + 2*padSize;
   else
      sizeB = size(a) + padSize;
   end

   b = mkconstarray(class(a), padVal, sizeB);
   
else
  switch method
    case 'constant'
        b = ConstantPad(a, padSize, padVal, direction);
        
    case 'circular'
        b = CircularPad(a, padSize, direction);
  
    case 'symmetric'
        b = SymmetricPad(a, padSize, direction);
        
    case 'antisymmetric'
        b = AntisymmetricPad(a, padSize, direction);
        
    case 'replicate'
        b = ReplicatePad(a, padSize, direction);
  end      
end

if (islogical(a))
    b = logical(b);
end

%%%
%%% ConstantPad
%%%
function b = ConstantPad(a, padSize, padVal, direction)

numDims = prod(size(padSize));

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
sizeB = zeros(1,numDims);
for k = 1:numDims
    M = size(a,k);
    switch direction
        case 'pre'
            idx{k}   = (1:M) + padSize(k);
            sizeB(k) = M + padSize(k);
            
        case 'post'
            idx{k}   = 1:M;
            sizeB(k) = M + padSize(k);
            
        case 'both'
            idx{k}   = (1:M) + padSize(k);
            sizeB(k) = M + 2*padSize(k);
    end
end

% Initialize output array with the padding value.  Make sure the
% output array is the same type as the input.
b         = mkconstarray(class(a), padVal, sizeB);
b(idx{:}) = a;


%%%
%%% CircularPad
%%%
function b = CircularPad(a, padSize, direction)

numDims = prod(size(padSize));

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
for k = 1:numDims
  M = size(a,k);
  dimNums = [1:M];
  p = padSize(k);
    
  switch direction
    case 'pre'
       idx{k}   = dimNums(mod([-p:M-1], M) + 1);
    
    case 'post'
      idx{k}   = dimNums(mod([0:M+p-1], M) + 1);
    
    case 'both'
      idx{k}   = dimNums(mod([-p:M+p-1], M) + 1);
  
  end
end
b = a(idx{:});

%%%
%%% SymmetricPad
%%%
function b = SymmetricPad(a, padSize, direction)

numDims = prod(size(padSize));

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
for k = 1:numDims
  M = size(a,k);
  dimNums = [1:M M:-1:1];
  p = padSize(k);
    
  switch direction
    case 'pre'
      idx{k}   = dimNums(mod([-p:M-1], 2*M) + 1);
            
    case 'post'
      idx{k}   = dimNums(mod([0:M+p-1], 2*M) + 1);
            
    case 'both'
      idx{k}   = dimNums(mod([-p:M+p-1], 2*M) + 1);
  end
end
b = a(idx{:});

%%%
%%% AntisymmetricPad
%%%
function b = AntisymmetricPad(a, padSize, direction)

numDims = prod(size(padSize));

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
idy   = cell(1,numDims);
for k = 1:numDims
  M = size(a,k);
  dimNums = [1:M];
  p = padSize(k);

  if (M == 1) 
      idx{k} = 1;
      idy{k} = 1;
  else
    switch direction
        case 'pre'
            dimNums = [1:M M:-1:1];
            idx{k}  = dimNums(mod([-p:M-1], 2*M-1) + 1);
            idy{k}  = [ones(1,p) 1:M];
        case 'post'
            dimNums = [1:M M-1:-1:1];
            idx{k}  = dimNums(mod([0:M+p-1], 2*M-1) + 1);
            idy{k}  = [1:M ones(1,p)*M];
        case 'both'
            dimNums = [1:M M:-1:1];
            idx{k}  = dimNums(mod([-p:M-1], 2*M-1) + 1);
            dimNums = [1:M M-1:-1:1];
            idx{k}  = [idx{k}, dimNums(mod([M:M+p-1], 2*M-1) + 1)];
            idy{k}  = [ones(1,p) 1:M ones(1,p)*M];
        end
    end
end
b = 2*a(idy{:}) - a(idx{:});

%%%
%%% ReplicatePad
%%%
function b = ReplicatePad(a, padSize, direction)

numDims = prod(size(padSize));

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
for k = 1:numDims
  M = size(a,k);
  p = padSize(k);
  onesVector = ones(1,p);
    
  switch direction
    case 'pre'
      idx{k}   = [onesVector 1:M];
            
    case 'post'
      idx{k}   = [1:M M*onesVector];
            
    case 'both'
      idx{k}   = [onesVector 1:M M*onesVector];
  end
end
 b = a(idx{:});

%%%
%%% ParseInputs
%%%
function [a, method, padSize, padVal, direction] = ParseInputs(varargin)

% default values
a         = [];
method    = 'constant';
padSize   = [];
padVal    = 0;
direction = 'both';

%checknargin(2,4,nargin,mfilename);

a = varargin{1};

padSize = varargin{2};
%checkinput(padSize, {'double'}, {'real' 'vector' 'nonnan' 'nonnegative' ...
%                    'integer'}, mfilename, 'PADSIZE', 2);

% Preprocess the padding size
if (prod(size(padSize)) < ndims(a))
    padSize           = padSize(:);
    padSize(ndims(a)) = 0;
end

if nargin > 2

    firstStringToProcess = 3;
    
    if ~ischar(varargin{3})
        % Third input must be pad value.
        padVal = varargin{3};
%        checkinput(padVal, {'numeric' 'logical'}, {'scalar'}, ...
%                   mfilename, 'PADVAL', 3);
        
        firstStringToProcess = 4;
        
    end
    
    for k = firstStringToProcess:nargin
%        validStrings = {'circular' 'replicate' 'symmetric' 'antisymmetric' ...
%                        'pre' 'post' 'both'};
%        string = checkstrs(varargin{k}, validStrings, mfilename, ...
%                           'METHOD or DIRECTION', k);
        string = varargin{k};
        switch string
         case {'circular' 'replicate' 'symmetric' 'antisymmetric'}
          method = string;
          
         case {'pre' 'post' 'both'}
          direction = string;
          
         otherwise
          error('Images:padarray:unexpectedError', '%s', ...
                'Unexpected logic error.')
        end
    end
end
    
% Check the input array type
if strcmp(method,'constant') & ~(isnumeric(a) | islogical(a))
    id = sprintf('Images:%s:badTypeForConstantPadding', mfilename);
    msg1 = sprintf('Function %s expected A (argument 1)',mfilename);
    msg2 = 'to be numeric or logical for constant padding.';
    error(id,'%s\n%s',msg1,msg2);
end
