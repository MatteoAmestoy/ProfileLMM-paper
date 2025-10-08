#include <RcppArmadillo.h>
#include <mvnorm.h>
#include <wishart.h>

// Package To Do:
//  - add shape checks
//  - add categorical U


using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

arma::sp_umat compute_Adj(arma::imat StoreZ,int No){
  int n = StoreZ.n_rows;
  arma::sp_umat adj(n,n);
  for (int i = 0; i < n; ++i) {
    for(int j = i+1; j < n; ++j) {
      int count = 0;
      for (int nobs = 0; nobs < No; ++nobs) {
        if (StoreZ(i,nobs) == StoreZ(j,nobs)) {
          count++;
        }
      }
      if(count>0){
        adj(i,j) = count;}
    }
  }
  return(adj);
}

arma::mat rdirichlet_cpp( arma::vec alpha_m) {
  int distribution_size = alpha_m.n_elem;
  // each row will be a draw from a Dirichlet
  arma::vec distribution = arma::zeros(distribution_size);

  double sum_term = 0;
  // loop through the distribution and draw Gamma variables
  for (int j = 0; j < distribution_size; ++j) {
    double cur = R::rgamma(alpha_m[j],1.0);
    distribution(j) = cur;
    sum_term += cur;
  }
  return(distribution/sum_term);
}


arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  bool success = false;
  int breakChol = 0;
  arma::mat R;
  while(success == false and breakChol<10)
  {
    success = arma::chol(R, sigma);

    if(success == false)
    {
      //Rcout<< "Csigma"<< sigma << "\n";
      Rcout << "Chol failed adding identity"<< "\n";
      breakChol += 1;
      sigma += arma::eye(sigma.n_rows,sigma.n_rows) * 1e-6;
    }
  }

  return arma::repmat(mu, 1, n).t() + Y * R;
}


// Create a data class
class DataObj {       // The class
  public:             // Access specifier
    arma::vec Y;        // outcome
    arma::mat XFE;        // covariates RE
    arma::mat XRE;        // covariates FE
    arma::mat XL;        // covariates latent
    arma::mat U;      // continuous clustered variables
    arma::vec ZRE;      // index of the RE membership
    arma::mat PX; //projection matrix on X
    arma::mat VmX; //inverse Variance of beta
    int n;
    int qFE;
    int qRE;
    int qL;
    int nRE;
    int qU;
    DataObj( arma::vec Y,// outcome
             arma::mat XFE,        // covariates RE
             arma::mat XRE,        // covariates FE
             arma::mat XL,        // covariates latent
             arma::mat U,      // continuous clustered variables
             arma::vec ZRE// index of the RE membership
    );
};

DataObj::DataObj( arma::vec Y_,// outcome
                  arma::mat XFE_,        // covariates RE
                  arma::mat XRE_,        // covariates FE
                  arma::mat XL_,        // covariates latent
                  arma::mat U_,      // continuous clustered variables
                  arma::vec ZRE_// index of the RE membership
){
  Y = Y_;
  XRE =XRE_;
  XFE= XFE_;
  U = U_;
  XL = XL_;
  ZRE = ZRE_;
  n  = Y.n_elem;
  qFE = XFE_.n_cols;
  qRE = XRE_.n_cols;
  qL = XL_.n_cols;
  qU = U_.n_cols;
  nRE = ZRE_.max() + 1;
};





class ParamLMM {        // The class
public:
  arma::vec beta; // FE parameters
  double sig2;
  arma::mat gamma;
  arma::mat W;

  double a; //prior invgamma
  double b; //prior invgamma
  double lambda; //prior beta
  arma::mat Phi; //prior Inverse wishart lat
  double eta; //prior Inverse wishart lat

  arma::vec YFE; // LMM (FE ) contribution to prediction
  arma::vec YLat; // LMM (Lat RE) contribution to prediction

  ParamLMM(arma::vec beta, // FE parameters
           double sig2,
           arma::mat gamma,
           arma::mat W,

           double a, //prior invgamma
           double b, //prior invgamma
           double lambda, //prior beta
           arma::mat Phi, //prior Inverse wishart lat
           double eta, //prior Inverse wishart lat
           arma::mat XFE,
           arma::vec YLat
  );
  void update(DataObj data, int nC, arma::vec Y, arma::ivec Z,  arma::vec cluster_count);
};
ParamLMM::ParamLMM(arma::vec beta_, // FE parameters
                   double sig2_,
                   arma::mat gamma_,
                   arma::mat W_,
                   double a_, //prior invgamma
                   double b_, //prior invgamma
                   double lambda_, //prior beta
                   arma::mat Phi_, //prior Inverse wishart lat
                   double eta_, //prior Inverse wishart lat
                   arma::mat XFE,
                   arma::vec YLat_

){
  beta = beta_;
  sig2 = sig2_;
  a = a_;
  b = b_;
  lambda = lambda_;
  gamma = gamma_;
  W = W_;
  Phi = Phi_;
  eta = eta_;
  YLat = YLat_;
  YFE = XFE * beta;// ADD LATENT EFFECT TO INITIALISATION?
};

void ParamLMM::update(DataObj data, int nC, arma::vec Y, arma::ivec Z,  arma::vec cluster_count){

  arma::mat W_m = arma::inv(W);// inverse of the current Lat RE covariance
  arma::mat omega = Phi; //initialisation of mean posterior of W
  arma::mat B = arma::eye(data.qFE,data.qFE);//initialisation of variance of FE posterior
  arma::vec bhat(data.qFE);// initialisation of mean of FE posterior
  B = B * lambda/sig2;

  arma::cube Plat(data.qL,data.qFE,nC);// storing projection matrix for FE correction

  for(int c = 0; c < nC; c++){
    // Rcout << "The value of ur : " << cluster_count(c) << "\n";
    if(cluster_count(c)>0){
      // pre-compute useful matrices
      arma::uvec idx =  arma::find(Z==c);
      arma::mat XL_ = data.XL.rows(idx);
      arma::mat X_ = data.XFE.rows(idx);
      arma::mat XXL = X_.t()*XL_;
      arma::mat XLXL = XL_.t()*XL_;

      arma::mat VmCore = arma::inv(W_m+XLXL/sig2)/pow(sig2,2);
      arma::mat XLY = XL_.t()*Y(idx);

      B = B + (X_.t()*X_/sig2-XXL*VmCore*XXL.t());
      bhat = bhat+(X_.t()*Y(idx)/sig2-XXL*VmCore*XLY);

      arma::mat Sig = W - W*XLXL*W/sig2+W*XLXL*VmCore*XLXL*W;



      // draw the random part of the latent RE.
      //Correction for the mean once beta is drawn
      gamma.col(c) = mvrnormArma(1,W*XLY/sig2-W*XLXL*VmCore*XLY,
                (Sig+Sig.t())/2.0).t();
      Plat.slice(c) = W*XXL.t()/sig2-W*XLXL*VmCore*XXL.t();
    }

  }

  // Rcout << "The value of ur : " << B << "\n";
  B = arma::inv(B);
  bhat = B*bhat;

  beta = mvrnormArma(1,bhat,(B+B.t())/2.0).t();
  YFE = data.XFE*beta;
  for(int c = 0; c < nC; c++){
    // Rcout << "The value of ur : " << cluster_count(c) << "\n";
    if(cluster_count(c)>0){

      gamma.col(c) = gamma.col(c) - Plat.slice(c)*beta;
      arma::uvec idx =  arma::find(Z==c);

      YLat(idx) =  data.XL.rows(idx) * gamma.col(c);
    }
    else{
      arma::vec mu(data.qL);
      gamma.col(c) = mvrnormArma(1,mu,(W+W.t())/2.0).t();
    }

    omega +=   gamma.col(c)*gamma.col(c).t();
  }

  //Rcout << "The value of ur : " << omega << "\n";
  W = riwish(nC + eta,(omega+omega.t())/2.0);
  double b_ = b+arma::as_scalar((Y-YFE-YLat).t()*(Y-YFE-YLat))/2.0;
  sig2 = 1.0/(Rcpp::rgamma( 1, a+data.n/2.0,1/b_)(0));

};



class ParamRE {       // The class
public:
  arma::mat Phi; //prior Inverse wishart RE
  double eta; //prior Inverse wishart RE
  arma::vec YRE; // random effect realisations
  arma::mat W; //random effect covariance matrix

  ParamRE(  arma::mat Phi, //prior Inverse wishart RE
            double eta, //prior Inverse wishart RE
            arma::vec YRE, // random effect realisations
            arma::mat W //random effect covariance matrix
  );

  void update(DataObj data, arma::vec Y, double sig2);

};
ParamRE::ParamRE(  arma::mat Phi_, //prior Inverse wishart RE
                   double eta_, //prior Inverse wishart RE
                   arma::vec YRE_, // random effect realisations
                   arma::mat W_ //random effect covariance matrix
){
  Phi = Phi_;
  eta = eta_;
  YRE = YRE_;
  W = W_;
};


void ParamRE::update(DataObj data, arma::vec Y, double sig2){

  arma::mat gamma(data.qRE,data.nRE);

  arma::mat omega = Phi;

  arma::mat W_m = inv(W);

  for(int ind = 0; ind < data.nRE; ind++){

    arma::uvec idx =  arma::find(data.ZRE==ind);
    arma::mat XRE_ = data.XRE.rows(idx);
    arma::mat XX = XRE_.t()*XRE_;

    arma::mat VmCore = W*XX*arma::inv(W_m+XX/sig2)/pow(sig2,2);
    arma::mat XY = XRE_.t()*Y(idx);


    arma::mat Sig = W - W*XX*W/sig2+VmCore*XX*W;
    gamma.col(ind) = mvrnormArma(1,W*XY/sig2-VmCore*XY,
              (Sig+Sig.t())/2).t();

    omega +=   gamma.col(ind)*gamma.col(ind).t();
    YRE(idx) = XRE_*gamma.col(ind);
  }
  W =  riwish(data.nRE+ eta,omega);
};

class ParamClus {       // The class
public:
  double lam0; // lambda prior NIW
  arma::vec mu0; // mu prior NIW
  double nu0; // nu prior NIW
  arma::mat  Psi0; // psi prior NIW
  arma::mat mu; //mean location of each latent class
  arma::cube Sigma; // covariance matrices of each latent cluster

  ParamClus( double lam0, // lambda prior NIW
             arma::vec mu0, // mu prior NIW
             double nu0, // nu prior NIW
             arma::mat  Psi0, // psi prior NIW
             arma::mat mu, //mean location of each latent class
             arma::cube Sigma // covariance matrices of each latent cluster
  );
  void update(DataObj data,int nC, arma::ivec Z);

};

ParamClus::ParamClus( double lam0_, // lambda prior NIW
                      arma::vec mu0_, // mu prior NIW
                      double nu0_, // nu prior NIW
                      arma::mat  Psi0_, // psi prior NIW
                      arma::mat mu_, //mean location of each latent class
                      arma::cube Sigma_ // covariance matrices of each latent cluster
){
  lam0 = lam0_;
  nu0 = nu0_;
  mu0 = mu0_;
  Psi0 = Psi0_;
  mu = mu_;
  Sigma = Sigma_;
};

void ParamClus::update(DataObj data,int nC, arma::ivec Z){
  arma::mat EUc(data.qU,nC); //average value of Uc per group
  arma::vec nCvec(nC); // nb of observation per group
  arma::cube vUc(data.qU,data.qU,nC); // E[Uc**2]
  EUc.fill(0);
  vUc.fill(0);

  for(int i = 0; i < data.n; i++){
    EUc.col(Z(i)) += data.U.row(i).t();
    nCvec(Z(i)) += 1.0 ;
    vUc.slice(Z(i)) += data.U.row(i).t()*data.U.row(i);

  };

  for(int c = 0; c < nC; c++){ // iterate over clusters
    arma::vec mun = mu0;
    double lambdan = lam0;
    double nun = nu0;
    arma::mat Psin = Psi0;

    if (nCvec(c) != 0 ){
      mun = (lam0*mu0+EUc.col(c))/(lam0+nCvec(c));
      lambdan +=  nCvec(c);
      nun += nCvec(c);
      Psin += ((EUc.col(c)/nCvec(c)-mu0)*(EUc.col(c)/nCvec(c)-mu0).t())*(lam0*nCvec(c))/(lam0+nCvec[c])+
        vUc.slice(c)-EUc.col(c)*EUc.col(c).t()/nCvec(c);
    }
    Sigma.slice(c) = riwish(nun, (Psin+Psin.t())/2.0);
    mu.col(c) = mvrnormArma(1,mun,(Sigma.slice(c)+Sigma.slice(c).t())/2.0/lambdan).t();
  }

};




class ParamAssign {       // The class
public:
  arma::ivec Z;
  arma::vec p0;
  double alpha;
  double scale;
  double shape;
  arma::vec cluster_count;
  int non_0_clust;
  arma::sp_umat adjMat;
  ParamAssign(
    arma::ivec Z,
    arma::vec p0,
    double scale,
    double shape
  );
  void update(DataObj data, int nC, arma::vec Y, double sig2,
              arma::cube SigmaGM, //current variance of gm
              arma::mat muGM, //current mean of gm
              arma::mat gammaL // current latent effect vector
  );

};
ParamAssign::ParamAssign(
  arma::ivec Z_,
  arma::vec p0_,
  double scale_,
  double shape_
){
  Z = Z_;
  p0 = p0_;
  int nC = p0.n_elem;
  int n = Z.n_elem;
  alpha = R::rgamma(shape_,scale_);
  scale = scale_;
  shape = shape_;
  arma::vec cluster_count_(nC);
  for(int i = 0; i < n; i++){
    cluster_count_(Z(i)) +=1;
  }
  cluster_count = cluster_count_;

  int non_0_clust_loc=0;
  for(int c = 0; c < nC; c++){
    if (cluster_count(c)>0){
      non_0_clust_loc += 1.0;
    }

  }
  non_0_clust = non_0_clust_loc;
  adjMat = arma::sp_umat(n,n);
};




void ParamAssign::update(DataObj data, int nC, arma::vec Y, double sig2 ,
                         arma::cube SigmaGM, //current variance of gm
                         arma::mat muGM, //current mean of gm
                         arma::mat gammaL // current latent effect vector
){

  // draw from a truncated dirichlet process Proof in Ishwaran p168
  double VarAcc; //auxilary accumulation variable

  arma::vec V(nC);
  VarAcc = 0;//sum of the counts for current class
  for(int c = 0; c < nC; c++){
    V(nC-1-c) = R::rbeta(1+cluster_count(nC-1-c),alpha+VarAcc);
    VarAcc += cluster_count(nC-1-c);
  }
  VarAcc = 1;//sum of the probability proportions still available
  for(int c = 0; c < nC; c++){
    p0(c) = V(c)*VarAcc;
    VarAcc = VarAcc*(1-V(c));
  }

  // arma::vec p0 = rdirichlet_cpp(cluster_count+alpha);
  // p0 = rdirichlet_cpp(cluster_count+alpha);
  arma::vec lp0(nC);

  arma::cube SigmaGMinvC(data.qU,data.qU,nC);
  arma::mat predC(data.n,nC);
  arma::vec cluster_count_loc(nC);

  // Rcout << "The value of ur : " << p0 << "\n";
  for(int c = 0; c < nC; c++){
    lp0(c) = log(p0(c))-arma::log_det_sympd(SigmaGM.slice(c))/2;
    SigmaGMinvC.slice(c) = arma::inv(SigmaGM.slice(c));
    predC.col(c) = Y -  data.XL*gammaL.col(c);
  }

  arma::vec logprob(nC);
  for(int o = 0; o < data.n; o++){

    logprob = lp0;
    double nconst = 0;

    for(int c = 0; c < nC; c++){
      logprob(c) += -arma::as_scalar((data.U.row(o)-muGM.col(c).t())*SigmaGMinvC.slice(c)*(data.U.row(o).t()-muGM.col(c))/2.0);
      logprob(c) += -pow(predC(o,c),2)/sig2/2.0;
      logprob(c) = exp(logprob(c));
      nconst += logprob(c);
    }

    double u = runif(1,0,1)(0)*nconst;
    for(int c = 0; c < nC; c++){
      if (u<logprob(c)){
        Z(o) = c;
        cluster_count_loc(c) += 1.0;
        u = nconst;
      }
      else{
        u = u-logprob(c);
      }
    }
  }
  int non_0_clust_loc;
  non_0_clust_loc = 0;
  for(int c = 0; c < nC; c++){
    if (cluster_count_loc(c)>0){
      non_0_clust_loc += 1.0;
    }

  }
  alpha = R::rgamma(shape+nC,1/(1/scale-log(VarAcc)));
  cluster_count = cluster_count_loc;
  non_0_clust = non_0_clust_loc;
}




// [[Rcpp::export]]
List GSLoopCPP(// parameters GS
    int nIt, // number of iterrations
    int nBurnIn, // number of burn in steps
    int nC, //number of clusters
    arma::vec Y,// outcome
    arma::mat XFE,        // covariates RE
    arma::mat XRE,        // covariates FE
    arma::mat XL,        // covariates latent
    arma::mat U,      // continuous clustered variables
    arma::vec ZRE,// index of the RE membership
    arma::vec beta, // FE parameters
    double sig2,
    arma::mat WRE,
    arma::mat muClus, //mean location of each latent class
    arma::cube SigmaClus, // covariance matrices of each latent cluster
    arma::mat WLat,
    arma::mat gammaLat,
    double a, //prior invgamma
    double b, //prior invgamma
    double lambdaFE, //prior beta FE
    arma::mat  PhiRE, //prior Inverse wishart RE
    double etaRE,//prior Inverse wishart RE
    double lam0, // lambda prior NIW
    arma::vec mu0, // mu prior NIW
    double nu0, // nu prior NIW
    arma::mat  Psi0, // psi prior NIW
    arma::mat PhiLat, //prior Inverse wishart lat
    double etaLat, //prior Inverse wishart lat
    double scale, // prior gamma on alpha (dirichlet prior)
    double shape // prior gamma on alpha (dirichlet prior)
)
{

  // Data object
  DataObj data(Y,// outcome
               XFE,        // covariates RE
               XRE,        // covariates FE
               XL,        // covariates latent
               U,      // continuous clustered variables
               ZRE// index of the RE membership
  );

  // Initialise theta

  arma::vec YLat(data.n); // random effect realisations
  ParamLMM thetaLMM(beta, // FE parameters
                    sig2,
                    gammaLat,
                    WLat,
                    a, //prior invgamma
                    b, //prior invgamma
                    lambdaFE, //projection matrix on X
                    PhiLat, //prior Inverse wishart lat
                    etaLat ,
                    data.XFE, //inverse Variance of beta
                    YLat
  );

  arma::vec YRE(data.n); // random effect realisations

  ParamRE thetaRE( PhiRE, //prior Inverse wishart RE
                   etaRE, //prior Inverse wishart RE
                   YRE, // random effect realisations
                   WRE //random effect covariance matrix
  );

  ParamClus thetaClus( lam0, // lambda prior NIW
                       mu0, // mu prior NIW
                       nu0, // nu prior NIW
                       Psi0, // psi prior NIW
                       muClus, //mean location of each latent class
                       SigmaClus // covariance matrices of each latent cluster
  );

  arma::vec p0(nC);
  arma::ivec Z(data.n);

  Z = arma::randi( data.n, arma::distr_param(0,nC-1) );

  ParamAssign thetaAssign(
      Z,
      p0,
      scale,
      shape
  );

  int nStore = nIt-nBurnIn-1;
  arma::cube StoreWL(data.qL,data.qL,nStore);
  arma::cube StoreWRE(data.qRE,data.qRE,nStore);
  arma::imat StoreZ(data.n,nStore);
  arma::cube StoremuClus(data.qU,nC,nStore);
  arma::mat StoreBeta(data.qFE,nStore);
  arma::vec Storesig2(nStore);
  arma::vec Storealpha(nStore);
  arma::cube StoreGamma(data.qL,nC,nStore);
  arma::field<arma::cube> StorePhiClus(nStore);


  // Rcout << "The value of count : " << thetaAssign.cluster_count << "\n";
  int storeIdx = 0;
  for(int it = 0; it < nIt; it++){
    thetaLMM.update(data, nC, data.Y-thetaRE.YRE, thetaAssign.Z, thetaAssign.cluster_count);
    // Rcout << "The value of Lat";
    thetaClus.update(data, nC, thetaAssign.Z);
    // Rcout << "The value of Clus";
    thetaRE.update(data, data.Y-thetaLMM.YFE-thetaLMM.YLat, thetaLMM.sig2);
    // Rcout << "The value of FE";
    thetaAssign.update(data, nC, data.Y-thetaRE.YRE-thetaLMM.YFE, thetaLMM.sig2, thetaClus.Sigma, //current variance of gm
                       thetaClus.mu, //current mean of gm
                       thetaLMM.gamma // current latent effect vector
    );
    // Rcout << "The value of ass";

    //
    if ((it%1000)==0){
      Rcout << "Iteration: " << it << "\n";
      Rcout << "Nb of non 0 clusters: " << thetaAssign.non_0_clust << "\n";
    }
    if(it>nBurnIn){
      StoreWL.slice(storeIdx) = thetaLMM.W;
      StoreWRE.slice(storeIdx) = thetaRE.W;
      StorePhiClus(storeIdx)  = thetaClus.Sigma;
      StoreZ.col(storeIdx) = thetaAssign.Z;
      StoremuClus.slice(storeIdx) = thetaClus.mu;
      StoreBeta.col(storeIdx) = thetaLMM.beta;
      StoreGamma.slice(storeIdx) = thetaLMM.gamma;
      Storesig2(storeIdx) = thetaLMM.sig2;
      Storealpha(storeIdx) = thetaAssign.alpha;
      storeIdx += 1;
    }

  }
  //Rcout << "The value of count : " << StorePhiClus(storeIdx-1) << "\n";
  List Store = List::create(Named("WL") = StoreWL ,
                            _["WRE"] = StoreWRE,
                            _["Z"] = StoreZ,
                            _["muClus"] = StoremuClus,
                            _["PhiClus"] = StorePhiClus,
                            _["beta"] = StoreBeta,
                            _["sig2"] = Storesig2,
                            _["alpha"] = Storealpha,
                            _["gamma"] = StoreGamma
                            // _["adjMat"] = compute_Adj(StoreZ,1500)
  );



  return Store;

}

