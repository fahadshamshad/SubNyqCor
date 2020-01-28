function M = AVMM(X,tol) % for generating orthogonal matrix

    if nargin==1
	  tol=1e-6;
    end
    n = size(X,1);
    M = zeros(n); % prealloc

    vi = randn(n,1);  


    M(:,1) = vi ./ norm(vi);
    
    for i=2:n
	  nrm = 0;
	  while nrm<tol
		vi = randn(n,1);
		vi = vi -  M(:,1:i-1)  * ( M(:,1:i-1).' * vi )  ;
		nrm = norm(vi);
	  end
	  M(:,i) = vi ./ nrm;

    end %i
    %C = M*X;
end

