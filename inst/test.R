require(Rcpp)
require(predictimportance)
sourceCpp("discreteLV_C.cpp")

simtime = 200
out <- OneRun('ImportanceRevision/Original/Data/MPN25/MPN25-web-1-run-3-mat.txt',
              'ImportanceRevision/Original/Data/MPN25/MPN25-web-1-run-3-pop.txt',
              simtime = simtime, deltat = .01)
mycopy <- out$n0s
ptm <- proc.time()
print(out$n0s)
rout <- discreteLV(out$rmat, out$alphas, out$n0s, out$deltat, out$simtime)
print(proc.time() - ptm)
ptm <- proc.time()
print(out$n0s)
acout <- discreteLV_C_arma(out$rmat, out$alphas, mycopy, out$deltat, out$simtime)
print(proc.time() - ptm)
ptm <- proc.time()
print(out$n0s)
cout <- discreteLV_C(out$rmat, out$alphas, out$n0s, out$deltat, out$simtime)
print(proc.time() - ptm)

plot(1:simtime, rout[1,], col = 1, type = 'l', ylim = c(min(rout), max(rout)))
lines(1:simtime, acout[1,], col = 1, lty = 'dashed')
for(i in 2:50){
    lines(1:simtime, rout[i,], col = i)
    lines(1:simtime, acout[i,], col = i, lty = 'dashed')
}

## dev.new()
## hist(rout - acout)
