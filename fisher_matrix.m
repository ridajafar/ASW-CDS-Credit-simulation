function I=fisher_matrix(lambda,theta)
%Construction of the fisher information matrix
I=zeros(2,2);
I(1,1)=(1-exp(-lambda(1)*theta))/lambda(1)^2;
I(2,2)=(exp(-lambda(1)*theta))/lambda(2)^2;

end

