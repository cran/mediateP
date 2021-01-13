mediate_contY_contM=function(data,
                             outcome="Y",
                             mediator="M",
                             exposure="X",
                             covariateY=c("X1","X2","X3","X4","X5","X6","X7","X8"),
                             covariateM=c("X1","X2","X3","X4","X5","X6","X7"),x0=0,x1=1) {


  #data=dat1;outcome=Outcome;mediator=Mediator;covariateY=covariateM=CovariateY;exposure=Exposure

  data = as.data.frame(data)
  ## (1) formula for outcome model and mediator model
  if (is.null(covariateY)) {
    formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,sep=""))
  } else {
    formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,"+",paste(covariateY,collapse="+"),sep=""))
  }
  if (is.null(covariateM)) {
    formula_M=as.formula(paste(mediator,"~",exposure,sep=""))
  } else {
    formula_M=as.formula(paste(mediator,"~",exposure,"+",paste(covariateM,collapse="+"),sep=""))
  }
  ## (2) estimate the outcome model and mediator model
  model_Y=summary(lm(formula_Y,data=data))
  model_M=summary(lm(formula_M,data=data))

  beta=model_Y$coef[,1];cov_beta=model_Y$cov.unscaled
  gamma=model_M$coef[,1];cov_gamma=model_M$cov.unscaled

  ## (3) covariance matrix of (beta,gamma)
  nbeta=dim(cov_beta)[1];ngamma=dim(cov_gamma)[1]
  S=matrix(0,ncol=nbeta+ngamma,nrow=nbeta+ngamma)
  S[1:ngamma,1:ngamma]=cov_gamma; S[(ngamma+1):(nbeta+ngamma),(ngamma+1):(nbeta+ngamma)]=cov_beta
  colnames(S)=rownames(S)=c(paste(names(gamma),"_gamma",sep=""),paste(names(beta),"_beta",sep=""))

  ## (4) approximate method
  ##### NIE TE and PM expressions
  if (is.null(covariateY)==0) {
    names(cY) = paste(names(beta),"_betafix",sep="")[-c(1:3)]
    beta_c=paste("beta_",covariateY,sep="")
  }
  if (is.null(covariateM)==0) {
    names(cM) = paste(names(gamma),"_gammafix",sep="")[-c(1:2)]
    gamma_c=paste("gamma_",covariateM,sep="")
  }

  NIEa_fun = function() {
    output = "beta2*gamma1*(x1-x0)"
    return(output)
  }
  variable=c("gamma0","gamma1",if(is.null(covariateM)==0) {gamma_c},"beta0","beta1","beta2",if(is.null(covariateY)==0) {beta_c})
  NIEa_D=deriv(parse(text=NIEa_fun()),variable)
  gamma0=gamma[1];gamma1=gamma[2];
  if(is.null(covariateM)==0) {
    for (i in (1:length(covariateM))) {assign(gamma_c[i],gamma[2+i])}
  }
  beta0=beta[1];beta1=beta[2];beta2=beta[3]
  if(is.null(covariateY)==0) {
    for (i in (1:length(covariateY))) {assign(beta_c[i],beta[3+i])}
  }

  TEa_fun = function() {
    output = "(beta2*gamma1+beta1)*(x1-x0)"
    return(output)
  }
  TEa_D=deriv(parse(text=TEa_fun()),variable)

  PMa_fun = function() {
    .UP = "beta2*gamma1"
    .BOT = "beta2*gamma1+beta1"
    output=paste("(",.UP,")/(",.BOT,")")
    return(output)
  }
  PMa_D=deriv(parse(text=PMa_fun()),variable)


  NIEa_D = eval(NIEa_D)
  NIEa_p = NIEa_D[1]
  lambda= t(attr(NIEa_D,"gradient"))
  V_NIEa = as.vector(t(lambda) %*% S %*% lambda)

  TEa_D = eval(TEa_D)
  TEa_p = TEa_D[1]
  lambda= t(attr(TEa_D,"gradient"))
  V_TEa = as.vector(t(lambda) %*% S %*% lambda)

  PMa_D = eval(PMa_D)
  PMa_p = PMa_D[1]
  lambda= t(attr(PMa_D,"gradient"))
  V_PMa = as.vector(t(lambda) %*% S %*% lambda)


  point_est = c(NIEa_p,TEa_p,PMa_p);
  names(point_est)=c("NIE","TE","PM")
  var_est = c(V_NIEa,V_TEa,V_PMa);
  names(var_est)=c("NIE","TE","PM")
  sd_est = sqrt(var_est)
  names(sd_est)=c("NIE","TE","PM")
  ci_est = rbind(point_est-1.96*sd_est,point_est+1.96*sd_est)
  rownames(ci_est) = c("Lower boundary","Upper boundary")

  return(list(point_est=point_est,var_est=var_est,sd_est=sd_est,ci_est=ci_est))
}
