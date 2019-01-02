
function [samples] = generateHaltonSamples(dim,n)
maxPrime = 1000;
primeNums = primes(maxPrime);
if dim > length(primeNums)
    disp('Change the maxPrime number, in generateHaltonSamples()');
end
skipPrime = 1;
base = primeNums(skipPrime:dim+skipPrime-1);

samples = zeros(n,dim);
for d = 1:dim
    samples(:,d) = HaltonSequence(n,base(d));
end
end

function hs = HaltonSequence(n,b)
% Preallocate the output
hs = zeros(n,1);
% Generate the numbers
for idx = 1:n
    hs(idx) = localHaltonSingleNumber(idx,b);
end
end

function hn = localHaltonSingleNumber(n,b)
% This function generates the n-th number in Halton's low
% discrepancy sequence.
n0 = n;
hn = 0;
f = 1/b;
while (n0>0)
    n1 = floor(n0/b);
    r = n0-n1*b;
    hn = hn + f*r;
    f = f/b;
    n0 = n1;
end
end