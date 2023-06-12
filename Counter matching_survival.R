library("survival")

pvalue.match = c()
pvalue.cmatch = c()
nsim = 500

for(z in 1:nsim){
  set.seed(z)
  
  # x1: proxy
  # x2: treatment
  # conf: confounder
  sim <- function(n=10000, a = 2, med0 = 140,
                  #p.x2 = 0.5, sen = 0.9, spe = 0.9,
                  p.x1 = 0.5, p.conf = 0.5,
                  hr.x2 = 1,
                  hr.conf = 1.5) {
    
    ## simulate x1, x2, x3
    #x2 <- rbinom(n, 1, p.x2)
    #x1 <- x2
    #x1[which(x2==1)] = rbinom(sum(x2), 1, sen)
    #x1[which(x2==0)] = rbinom(n - sum(x2), 1, 1 - spe)
    
    x1 = rbinom(n,1,p.x1)
    conf = rbinom(n,1,p.conf)
    x2 = rbinom(n,1, exp(-2 + 4*x1 + 1*conf)/(1+exp(-2 + 4*x1 + 1*conf)))
    
    ## Calculate coefficients
    b0=med0*log(2)^(-1/a)
    b.x2 = log (hr.x2)  
    b.conf = log(hr.conf)
    
    ## simulate non-censored event time
    t.event <- rweibull(n = n, shape = a,
                        scale = b0 * exp((-b.x2*x2 - b.conf*conf)/a))
    
    ## simulate non-informative censoring
    t.cens=rweibull (n = n, shape = a, scale = 50)
    ## get final simulated data
    tobs <- pmin(t.event, t.cens)
    out <- data.frame(x1, x2, conf,
                      tobs,
                      stt = as.numeric(tobs == t.event), 
                      id = 1:n)
    
    ## return output
    return(out)
  }
  
  data = sim()
  
  #-----------------------
  # No countermatching
  #-----------------------
  
  n = 10000
  source.noc = data[order(data$tobs),]
  source.noc$case = 0
  source.noc$ctrl = 0
  source.noc$sel = 0
  source.noc$pair = 0
  j = 1
  
  # We only analyze the first 200 cases
  ind = which(source.noc$stt ==1)
  ind = ind[1:200]
  
  for (i in ind){
    
    if (source.noc$sel[i] == 0){
      
      source.noc$case[i] = 1
      
      # All non-cases at the time t.obs of the case i under consideration
      all.ctrl = source.noc[(i+1):n,]
      
      # All non-cases eligible to be matched with the case i
      ctrl = all.ctrl[which(all.ctrl$sel == 0 & all.ctrl$conf == source.noc$conf[i]),]
      
      # Select one candidate from ctrl set
      ctrl.i.id = sample(ctrl$id, size= 1)
      
      # Record the selected ctrl
      source.noc$ctrl[i] = ctrl.i.id 
      
      # Record selection status
      source.noc$sel[i] = 1
      source.noc$sel[which(source.noc$id == ctrl.i.id)] = 1
      
      # Record pair order
      source.noc$pair[i] = j
      source.noc$pair[which(source.noc$id == ctrl.i.id)] = j
      j=j+1
    }
  }
  
  case.noc = source.noc[source.noc$sel == 1 & source.noc$case ==1,]
  ctrl.noc = source.noc[source.noc$sel == 1 & source.noc$case ==0,]
  ctrl.noc = ctrl.noc[order(ctrl.noc$pair),]
  
  # MH test on the first 200 matched sets
  pvalue.match[z] = 
    sens.analysis.mh(cases.exposed = case.noc$x2,
                     referents.exposed = ctrl.noc$x2,
                     no.referents = 1,
                     Gamma = 1)$lower.bound.pval
  
  
  #----------------------
  # With countermatching
  #----------------------
  
  source = data[order(data$tobs),]
  source$case = 0
  source$ctrl = 0
  source$sel = 0
  source$pair = 0
  j = 1
  
  ind = which(source$stt ==1)
  ind = ind[1:200]
  
  for (i in ind){
    
    if (source$sel[i] == 0){
      
      source$case[i] = 1
      
      # All non-cases at the time t.obs of the case i under consideration
      all.ctrl = source[(i+1):n,]
      
      # All non-cases eligible to be counter-matched with the case i
      ctrl = all.ctrl[which(all.ctrl$x1 == 1-source$x1[i] & 
                              all.ctrl$conf == source$conf[i] &
                              all.ctrl$sel == 0),]
      
      # Select one candidate from ctrl set
      ctrl.i.id = sample(ctrl$id, size= 1)
      
      # Record the selected ctrl
      source$ctrl[i] = ctrl.i.id 
      
      # Record selection status
      source$sel[i] = 1
      source$sel[which(source$id == ctrl.i.id)] = 1
      
      # Record pair order
      source$pair[i] = j
      source$pair[which(source$id == ctrl.i.id)] = j
      j=j+1
    }
  }
  
  # Stt = 0, case = 1: Never happens
  # Stt = 0, case = 0: control
  # Stt = 1, case = 0: control matched to a case at time t1 (that will also become a case at time t2 > t1)
  # Stt = 1, case = 1: case
  
  case = source[source$sel == 1 & source$case ==1,]
  ctrl = source[source$sel == 1 & source$case ==0,]
  ctrl = ctrl[order(ctrl$pair),]
  
  #all = rbind(case,ctrl)
  #all = all[order(all$pair),c("x1","x2","case","pair")]
  
  # MH test on the first 1000 cases
  pvalue.cmatch[z] = 
    sens.analysis.mh(cases.exposed = case$x2,
                     referents.exposed = ctrl$x2,
                     no.referents = 1,
                     Gamma = 1)$lower.bound.pval
}                                


sum(pvalue.cmatch<0.05)/500
sum(pvalue.match<0.05)/500

