function loglike = DPloglike( alpha, K, m )

loglike = K*log(alpha) + sum(gammaln(m)) + gammaln(alpha) - gammaln(sum(m)+alpha);
