functions {
  real fiduciaJacobian_lp(int n,int d,real phi,real theta,vector dataX,vector dataY,vector dataZ){
    matrix[n*3,d+2] JacobianMatrix=rep_matrix(0,n*3,d+2);
    //The derivatives with respect to mu do not contribute
    for(k in 1:n){
      JacobianMatrix[k,1]=cos(theta)*cos(phi);
      JacobianMatrix[k,2]=-sin(theta)*sin(phi);
      JacobianMatrix[k,3]=dataX[k]-(cos(theta)*sin(phi));
    }
    for(k in 1:n){
      JacobianMatrix[n+k,1]=sin(theta)*cos(phi);
      JacobianMatrix[n+k,2]=cos(theta)*sin(phi);
      JacobianMatrix[n+k,4]=dataY[k]-(sin(theta)*sin(phi));
    }
    for(k in 1:n){
      JacobianMatrix[2*n+k,1]=-sin(phi);
      JacobianMatrix[2*n+k,5]=dataZ[k]-cos(phi)-1;
    }
    return (0.5)*log_determinant((JacobianMatrix')*JacobianMatrix);
  }
}

data {
  int<lower=1> n;
  int<lower=1> d;
  vector[3] y[n];
  matrix[d,d] sig;
  vector[n] dataX;
  vector[n] dataY;
  vector[n] dataZ;
}

parameters {
  real<lower=0,upper=pi()> phi;
  real<lower=0,upper=2*pi()> theta;
}

transformed parameters{
  vector[3] mu;
  real logJac;
  
  mu=rep_vector(0,3);
  mu[1]=cos(theta)*sin(phi);
  mu[2]=sin(theta)*sin(phi);
  mu[3]=cos(phi)+1;
  
  logJac=fiduciaJacobian_lp(n,d,phi,theta,dataX,dataY,dataZ);
}

model {
  // phi ~ uniform(0,pi());
  // theta ~ uniform(0,2*pi());
  y ~ multi_normal(mu, sig);
  // target+=log(sin(phi));
  target+=logJac; //adjusting for the rest of the Jacobian
}
