function out=blkdiagMat(A,N) 
Ar = repmat(A, 1, N);                                   % Repeat Matrix
Atemp = mat2cell(Ar, size(A,1), repmat(size(A,2),1,N));    % Create Cell Array Of Orignal Repeated Matrix
out = blkdiag(Atemp{:});  
