#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;
using namespace std;

//Some useful shortcuts
// [[Rcpp::export]]
double theta(vec p, int size, double type){
  double res;
  if (type==00)
  {
    res = pow(p(0),size);
  }
  else if (type==10)
  {
    res = pow(p(0)+p(1),size)-pow(p(0),size);
  }
  else if (type==01)
  {
    res = pow(p(0)+p(2),size)-pow(p(0),size);
  }
  else 
  {
    res = 1-pow(p(0)+p(1),size)-pow(p(0)+p(2),size)+pow(p(0),size);
  }
  return res;  
}

double f1(vec p, vec SE, vec SP,int s){
  double res = SE(0) + (theta(p,s,01)+theta(p,s,00))*(1-SE(0)-SP(0));
  return res;
}

double f2(vec p, vec SE, vec SP,int s){
  double res = SE(1) + (theta(p,s,10)+theta(p,s,00))*(1-SE(1)-SP(1));
  return res;
}

double f11(vec p, vec SE, vec SP,int s){
  double res = theta(p,s,00)*(1-SP(0))*(1-SP(1))+theta(p,s,10)*SE(0)*(1-SP(1))+theta(p,s,01)*SE(1)*(1-SP(0))+theta(p,s,11)*SE(0)*SE(1);
  return res;
}

double g1(vec p, int c , int s1, int s2){
  double res = R::choose(s2,c)*pow(theta(p,s1,10)+theta(p,s1,11),c)*pow(theta(p,s1,01)+theta(p,s1,00),s2-c);
  return res;
}

double g2(vec p, int c , int s1, int s2){
  double res = R::choose(s2,c)*pow(theta(p,s1,01)+theta(p,s1,11),c)*pow(theta(p,s1,10)+theta(p,s1,00),s2-c);
  return res;
}

double g12(int c1, int c2, int n1, double eta_00, double eta_10, double eta_01, double eta_11){
  double res = 0;
  for (int w=0; w<=min(c1,c2); w++){
    res += R::choose(n1,w)*R::choose(n1-w,c1-w)*R::choose(n1-c1,c2-w)*pow(eta_11,w)*pow(eta_10,c1-w)*pow(eta_01,c2-w)*pow(eta_00,n1-c1-c2+w);
  } 
  return res;
}


double f(vec p, int n1, int a, int n2, int b){
  double pre = pow(p(0)+p(1), (n2-b)*a);
  double prf = pow(p(0)+p(2), (n1-a)*b);
  double mid1 = 1- pow(1- pow(p(0)+p(2),b)*pow(p(0)/(p(0)+p(1)),n2-b),a);
  double mid2 = 1- pow(1- pow(p(0)+p(1),a)*pow(p(0)/(p(0)+p(2)),n1-a),b);
  double mid3 = 0;
  for (int l=1 ; l<=a;l++){
    for (int m=1; m<=b; m++){
      mid3 += pow(-1,m+l)*R::choose(a,l)*R::choose(b,m)*pow(p(0)/(p(0)+p(1)),l*(n2-b))*pow(p(0)/(p(0)+p(2)),m*(n1-a))*pow(p(0),l*m)*pow(p(0)+p(1),(a-l)*m)*pow(p(0)+p(2),(b-m)*l);
    }
  }
  double res = R::choose(n1,a)*R::choose(n2,b)*pow(p(0),(n1-a)*(n2-b))*pre*prf*(1-(mid1+mid2-mid3));
  return res;
}

//First probability
// [[Rcpp::export]]
double eff_nomaster_1(vec p, vec SE, vec SP, int n){
  double res_1 = (p(1)+p(3))*pow(SE(0),2)+ (p(0)+p(2))*pow(f1(p,SE,SP,n-1),2);
  double res_2 = 0; 
  for (int c1=0; c1<=n;c1++){
    res_2+= (1-SP(0))*pow(1-SE(0),c1)*pow(SP(0),n-c1)*(theta(p,n,00)+theta(p,n,01))*g1(p,c1,n-1,n)+ SE(0)*pow((1-SE(0)),c1)*pow(SP(0),n-c1)*(g1(p,c1,n,n)-(theta(p,n,00)+theta(p,n,01))*g1(p,c1,n-1,n));
  }
  double res = res_1+2*res_2;
  return res;
}

//Second probability
// [[Rcpp::export]]
double eff_nomaster_2(vec p, vec SE, vec SP, int n){
  double res_1 = (p(2)+p(3))*pow(SE(1),2)+ (p(0)+p(1))*pow(f2(p,SE,SP,n-1),2);
  double res_2 = 0; 
  for (int c2=0; c2<=n;c2++){
    res_2+= (1-SP(1))*pow(1-SE(1),c2)*pow(SP(1),n-c2)*(theta(p,n,00)+theta(p,n,10))*g2(p,c2,n-1,n)+ SE(1)*pow((1-SE(1)),c2)*pow(SP(1),n-c2)*(g2(p,c2,n,n)-(theta(p,n,00)+theta(p,n,10))*g2(p,c2,n-1,n));
  }
  double res = res_1+2*res_2;
  return res;
}

// Intersection probability: 1
double eff_nomaster_3_1(vec p, vec SE, vec SP, int n){
  double res = p(0)*pow(f11(p,SE,SP,n-1),2)+p(1)*pow(SE(0),2)*pow(f2(p,SE,SP,n-1),2)+p(2)*pow(SE(1),2)*pow(f1(p,SE,SP,n-1),2)+p(3)*pow(SE(0)*SE(1),2);
  return res;
}

// Intersection probability: 2
double diag_3_2(int r1, int r2, int c1, int c2, int n, vec SE, vec SP){
  double res = pow(SE(0),r1+c1)*pow(1-SP(0),2-r1-c1)*pow(SE(1),r2)*pow(1-SP(1),1-r2)*pow(1-SE(1),c2)*pow(SP(1),n-c2);
  return res;
}

double eff_nomaster_3_2(vec p, vec SE, vec SP, int n){
  double res = 0; 
  for (int c2=0; c2<=n; c2++){
    double h_111n_000c2 = theta(p,2*n-1,00)*g2(p,c2,n-1,n-1)+ theta(p,n,00)*theta(p,n-1,01)*g2(p,c2-1,n-1,n-1);
    double h_111n_000c2_d = h_111n_000c2*diag_3_2(0, 0, 0, c2,  n,  SE,  SP);
    
    double h_111n_001c2 = theta(p,n,00)*theta(p,n-1,10)*g2(p,c2,n-1,n-1)+ theta(p,n,00)*theta(p,n-1,11)*g2(p,c2-1,n-1,n-1);
    double h_111n_001c2_d = h_111n_001c2*diag_3_2(0, 0, 1, c2,  n,  SE,  SP);
    
    double h_111n_010c2 = theta(p,n,00)*R::choose(n-1,c2)*pow(p(2),c2)*pow(p(0),n-c2-1)*g2(p,0,n-1,n-c2-1)-h_111n_000c2;
    for (int w=0; w<=c2-1; w++)
    {
      h_111n_010c2+=theta(p,n,00)*R::choose(n-1,w)*pow(p(2),w)*pow(p(0),n-w-1)*g2(p,c2-w,n-1,n-w-1)+ theta(p,n,01)*R::choose(n-1,w)*pow(p(2),w)*pow(p(0),n-w-1)*g2(p,c2-w-1,n-1,n-w-1);
    }
    double h_111n_010c2_d = h_111n_010c2*diag_3_2(0, 1, 0, c2,  n,  SE,  SP);
    
    double h_111n_100c2 = theta(p,n,00)*theta(p,n-1,10)*g2(p,c2,n-1,n-1)+p(0)*theta(p,n-1,01)*theta(p,n-1,10)*g2(p,c2-1,n-1,n-1);
    double h_111n_100c2_d = h_111n_100c2*diag_3_2(1, 0, 0, c2,  n,  SE,  SP);
    
    double h_111n_011c2 =  -h_111n_000c2 - h_111n_001c2 - h_111n_010c2;
    for (int w=0; w<=c2; w++){
      h_111n_011c2 += R::choose(n,w)*pow(p(2),w)*pow(p(0),n-w)*g2(p,c2-w,n-1,n-w);
    }
    double h_111n_011c2_d = h_111n_011c2*diag_3_2(0, 1, 1, c2,  n,  SE,  SP);
    
    double h_111n_101c2 = (theta(p,n,10)+theta(p,n,00))*g2(p,c2,n-1,n)-h_111n_000c2-h_111n_001c2-h_111n_100c2;
    double h_111n_101c2_d = h_111n_101c2*diag_3_2(1,0, 1, c2,  n,  SE,  SP);
    
    double h_111n_110c2 = theta(p,n,00)*g2(p,c2,n,n-1)+theta(p,n,01)*g2(p,c2-1,n,n-1)-h_111n_000c2-h_111n_100c2-h_111n_010c2;
    double h_111n_110c2_d =  h_111n_110c2*diag_3_2(1,1, 0, c2,  n,  SE,  SP);
    
    double h_111n_111c2 = g2(p,c2,n,n) - h_111n_000c2 -h_111n_001c2-h_111n_010c2-h_111n_100c2-h_111n_011c2-h_111n_101c2-h_111n_110c2;
    double h_111n_111c2_d = h_111n_111c2*diag_3_2(1,1, 1, c2,  n,  SE,  SP);
    
    res+=h_111n_000c2_d+h_111n_001c2_d+h_111n_010c2_d+h_111n_100c2_d+h_111n_011c2_d+h_111n_101c2_d+h_111n_110c2_d+h_111n_111c2_d;
  }
  return res;
}


// Intersection probability: 3
double diag_3_3(int r1, int r2, int c1, int c2, int n, vec SE, vec SP){
  double res = pow(SE(1),r2+c2)*pow(1-SP(1),2-r2-c2)*pow(SE(0),r1)*pow(1-SP(0),1-r1)*pow(1-SE(0),c1)*pow(SP(0),n-c1);
  return res;
}

double eff_nomaster_3_3(vec p, vec SE, vec SP, int n){
  double res = 0; 
  for (int c1=0; c1<=n; c1++){
    double h_11n1_00c10 = theta(p,2*n-1,00)*g1(p,c1,n-1,n-1)+ theta(p,n,00)*theta(p,n-1,10)*g1(p,c1-1,n-1,n-1);
    double h_11n1_00c10_d = h_11n1_00c10*diag_3_3(0, 0, c1, 0,  n,  SE,  SP);
    
    double h_11n1_00c11 = theta(p,n,00)*theta(p,n-1,01)*g1(p,c1,n-1,n-1)+ theta(p,n,00)*theta(p,n-1,11)*g1(p,c1-1,n-1,n-1);
    double h_11n1_00c11_d = h_11n1_00c11*diag_3_3(0, 0, c1, 1,  n,  SE,  SP);
    
    double h_11n1_10c10 = theta(p,n,00)*R::choose(n-1,c1)*pow(p(1),c1)*pow(p(0),n-c1-1)*g1(p,0,n-1,n-c1-1)-h_11n1_00c10;
    for (int w=0; w<=c1-1; w++)
    {
      h_11n1_10c10+=theta(p,n,00)*R::choose(n-1,w)*pow(p(1),w)*pow(p(0),n-w-1)*g1(p,c1-w,n-1,n-w-1)+ theta(p,n,10)*R::choose(n-1,w)*pow(p(1),w)*pow(p(0),n-w-1)*g1(p,c1-w-1,n-1,n-w-1);
    }
    double h_11n1_10c10_d = h_11n1_10c10*diag_3_3(1, 0, c1, 0,  n,  SE,  SP);
    
    double h_11n1_01c10 = theta(p,n,00)*theta(p,n-1,01)*g1(p,c1,n-1,n-1)+p(0)*theta(p,n-1,01)*theta(p,n-1,10)*g1(p,c1-1,n-1,n-1);
    double h_11n1_01c10_d = h_11n1_01c10*diag_3_3(0, 1, c1, 0,  n,  SE,  SP);
    
    double h_11n1_10c11 =  -h_11n1_00c10 - h_11n1_10c10 - h_11n1_00c11;
    for (int w=0; w<=c1; w++){
      h_11n1_10c11 += R::choose(n,w)*pow(p(1),w)*pow(p(0),n-w)*g1(p,c1-w,n-1,n-w);
    }
    double h_11n1_10c11_d = h_11n1_10c11*diag_3_3(1, 0, c1, 1,  n,  SE,  SP);
    
    double h_11n1_01c11 = (theta(p,n,01)+theta(p,n,00))*g1(p,c1,n-1,n)-h_11n1_00c10-h_11n1_00c11-h_11n1_01c10;
    double h_11n1_01c11_d = h_11n1_01c11*diag_3_3(0, 1, c1, 1,  n,  SE,  SP);
    
    double h_11n1_11c10 = theta(p,n,00)*g1(p,c1,n,n-1)+theta(p,n,10)*g1(p,c1-1,n,n-1)-h_11n1_00c10-h_11n1_10c10-h_11n1_01c10;
    double h_11n1_11c10_d =  h_11n1_11c10*diag_3_3(1,1, c1, 0,  n,  SE,  SP);
    
    double h_11n1_11c11 = g1(p,c1,n,n) - h_11n1_00c10 -h_11n1_10c10-h_11n1_01c10-h_11n1_00c11-h_11n1_11c10-h_11n1_01c11-h_11n1_10c11;
    double h_11n1_11c11_d = h_11n1_11c11*diag_3_3(1,1, c1, 1,  n,  SE,  SP);
    
    res+=h_11n1_00c10_d+h_11n1_00c11_d+h_11n1_01c10_d+h_11n1_10c10_d+h_11n1_01c11_d+h_11n1_10c11_d+h_11n1_11c10_d+h_11n1_11c11_d;
  }
  return res;
}

// Intersection probability: 4
double diag_3_4(int r1, int r2, int c1, int c2, int n, vec SE, vec SP){
  double res = pow(SE(0),r1)*pow(1-SP(0),1-r1)*pow(SE(1),r2)*pow(1-SP(1),1-r2)*pow(1-SE(0),c1)*pow(SP(0),n-c1)*pow(1-SE(1),c2)*pow(SP(1),n-c2);
  return res;
}

double eff_nomaster_3_4(vec p, vec SE, vec SP, int n){
  double res = 0; 
  for (int c1=0 ; c1<=n; c1++){
    for (int c2=0; c2<=n; c2++){
      double h_11nn_00c1c2 = g12(c1,c2,n,theta(p,n,00), p(0)*theta(p,n-1,10),p(0)*theta(p,n-1,01), p(0)*theta(p,n-1,11));
      double h_11nn_00c1c2_d = h_11nn_00c1c2*diag_3_4(0,0,c1,c2,n,SE,SP);
      
      double h_11nn_10c1c2 = g12(c1,c2,n,theta(p,n,00), theta(p,n,10), p(0)*theta(p,n-1,01),p(0)*theta(p,n-1,11)+p(1)*(theta(p,n-1,01)+theta(p,n-1,11))) - h_11nn_00c1c2;
      double h_11nn_10c1c2_d = h_11nn_10c1c2*diag_3_4(1,0,c1,c2,n,SE,SP);
      
      double h_11nn_01c1c2 = g12(c1,c2,n,theta(p,n,00), p(0)*theta(p,n-1,10), theta(p,n,01), p(0)*theta(p,n-1,11)+p(2)*(theta(p,n-1,10)+theta(p,n-1,11))) - h_11nn_00c1c2;
      double h_11nn_01c1c2_d = h_11nn_01c1c2*diag_3_4(0,1,c1,c2,n,SE,SP);
      
      double h_11nn_11c1c2 = g12(c1,c2,n,theta(p,n,00),theta(p,n,10),theta(p,n,01),theta(p,n,11))-h_11nn_00c1c2-h_11nn_10c1c2-h_11nn_01c1c2;
      double h_11nn_11c1c2_d = h_11nn_11c1c2*diag_3_4(1,1,c1,c2,n,SE,SP);
      
      res+= h_11nn_00c1c2_d + h_11nn_10c1c2_d + h_11nn_01c1c2_d + h_11nn_11c1c2_d;
    }
    return res;
  }
  return res;
}

// Intersection probability: 5
double diag_3_5(int r1, int r2, int c1, int c2, int n, vec SE, vec SP){
  double res = pow(SE(0),r1)*pow(1-SP(0),1-r1)*pow(1-SE(1),r2)*pow(SP(1),n-r2)*pow(1-SE(0),c1)*pow(SP(0),n-c1)*pow(SE(1),c2)*pow(1-SP(1),1-c2);
  return res;
}

double eff_nomaster_3_5(vec p, vec SE, vec SP, int n){
  double res = 0; 
  for (int r2=0 ; r2<=n; r2++){
    for (int c1=0; c1<=n; c1++){
      double h_1nn1_0r2c10 = theta(p,2*n-1,00)*f(p,n-1,c1,n-1,r2) + theta(p,n,00)*theta(p,n-1,01)*f(p,n-1,c1,n-1,r2-1) + theta(p,n,00)*theta(p,n-1,10)*f(p,n-1,c1-1,n-1,r2) + p(0)*theta(p,n-1,10)*theta(p,n-1,01)*f(p,n-1,c1-1,n-1,r2-1);
      double h_1nn1_0r2c10_d = h_1nn1_0r2c10*diag_3_5(0,r2,c1,0,n,SE,SP);
      
      double h_1nn1_0r2c11 = theta(p,n,00)*f(p,n,c1,n-1,r2) + theta(p,n,01)*f(p,n,c1,n-1,r2-1) - h_1nn1_0r2c10;
      double h_1nn1_0r2c11_d = h_1nn1_0r2c11* diag_3_5(0,r2,c1,1,n,SE,SP);
      
      double h_1nn1_1r2c10 = theta(p,n,00)*f(p,n-1,c1,n,r2) + theta(p,n,10)*f(p,n-1,c1-1,n,r2) - h_1nn1_0r2c10;
      double h_1nn1_1r2c10_d = h_1nn1_1r2c10* diag_3_5(1,r2,c1,0,n,SE,SP);
      
      double h_1nn1_1r2c11 = f(p,n,c1,n,r2) - h_1nn1_0r2c10 - h_1nn1_1r2c10 - h_1nn1_0r2c11 ;
      double h_1nn1_1r2c11_d = h_1nn1_1r2c11* diag_3_5(1,r2,c1,1,n,SE,SP);
      
      res += h_1nn1_0r2c10_d+h_1nn1_0r2c11_d+h_1nn1_1r2c10_d+h_1nn1_1r2c11_d;
    }
    return res;
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////PSe1//////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

double pse1_nomaster_1(vec p, vec SE, vec SP, int n){
  double res = pow(SE(0),3)+2*pow(SE(0),2)*(1-SE(0))*pow(1-f1(p,SE,SP,n),n-1);
  return res;
}

double pse1_nomaster_2(vec p, vec SE, vec SP, int n){
  double res_1 = SE(0)*pow(p(1)+p(3),-1)*(p(1)*pow(f2(p,SE,SP,n-1),2)+p(3)*pow(SE(1),2));
  double res_2 = 0;
  for (int c2=0; c2<=n; c2++){
    double mid1 = p(1)*(theta(p,n-1,10)+theta(p,n-1,00))*g2(p,c2,n-1,n);
    double mid1_C = mid1 * (1-SP(1))*pow(1-SE(1),c2)*pow(SP(1),n-c2);
    
    double mid2 = p(1)*(theta(p,n-1,10)+theta(p,n-1,00))*g2(p,c2,n,n-1)+p(1)*(theta(p,n-1,01)+theta(p,n-1,11))*g2(p,c2-1,n,n-1)-mid1;
    double mid2_C = mid2 *SE(1)*pow(1-SE(1),c2)*pow(SP(1),n-c2);
    
    res_2 += mid1_C + mid2_C; 
  }
  res_2 = pow(p(1)+p(3),-1)*(p(3)*SE(1)*(1-SE(1))*pow(1-f2(p,SE,SP,n),n-1)+res_2);
  double res  = res_1 + 2* SE(0)*res_2; 
  return res;
}

double pse1_nomaster_3_1(vec p, vec SE, vec SP, int n){
  double res = pow(SE(0),3)*pow(p(1)+p(3),-1)*(p(3)*pow(SE(1),2)+p(1)*pow(f2(p,SE,SP,n-1),2));
  return res;
}

double pse1_nomaster_3_2(vec p, vec SE, vec SP, int n){
  double res = 0;
  for (int c2=0; c2<=n; c2++){
    double mid1 = p(1)*(theta(p,n-1,10)+theta(p,n-1,00))*g2(p,c2,n-1,n);
    double mid1_C = mid1 * (1-SP(1))*pow(1-SE(1),c2)*pow(SP(1),n-c2);
    
    double mid2 = p(1)*(theta(p,n-1,10)+theta(p,n-1,00))*g2(p,c2,n,n-1)+p(1)*(theta(p,n-1,01)+theta(p,n-1,11))*g2(p,c2-1,n,n-1)-mid1;
    double mid2_C = mid2 *SE(1)*pow(1-SE(1),c2)*pow(SP(1),n-c2);
    
    res += mid1_C + mid2_C; 
  }
  res = pow(SE(0),3)*pow(p(1)+p(3),-1)*(p(3)*SE(1)*(1-SE(1))*pow(1-f2(p,SE,SP,n),n-1)+res);
  return res;
}

double pse1_nomaster_3_3(vec p, vec SE, vec SP, int n){
  double res = 0;
  for (int c1=0; c1<=n-1;c1++){
    double mid1 = 0;
    for (int w=0; w<=c1;w++){
      mid1 += R::choose(n-1,w)*pow(p(1),w)*pow(p(0),n-1-w)*g1(p,c1-w,n-1,n-w-1);
    }
    mid1 = mid1*p(1);
    double mid1_c = mid1*(1-SP(1))*pow(SP(0),n-1-c1)*pow(1-SE(0),c1);
    
    double mid2 = p(1)*g1(p,c1,n,n-1)-mid1;
    double mid2_c = mid2*SE(1)*pow(SP(0),n-1-c1)*pow(1-SE(0),c1);
    
    res += mid1_c+mid2_c;
  }
  res = pow(p(1)+p(3),-1)*pow(SE(0),2)*(1-SE(0))*(p(3)*pow(SE(1),2)*pow(1-f1(p,SE,SP,n),n-1)+f2(p,SE,SP,n-1)*res);
  return res;
}

double pse1_nomaster_3_4(vec p, vec SE, vec SP, int n){
  double res = 0;
  for (int c1=0; c1<=n-1;c1++){
    for (int c2=0; c2<=n; c2++){
      double mid1 = p(1)*(theta(p,n-1,10)+theta(p,n-1,00))*g12(c1,c2,n-1,theta(p,n,00), theta(p,n,10), p(0)*theta(p,n-1,01),p(0)*theta(p,n-1,11)+p(1)*(theta(p,n-1,01)+theta(p,n-1,11))) + 
        p(1)*(theta(p,n-1,01)+theta(p,n-1,11))*g12(c1,c2-1,n-1,theta(p,n,00), theta(p,n,10), p(0)*theta(p,n-1,01),p(0)*theta(p,n-1,11)+p(1)*(theta(p,n-1,01)+theta(p,n-1,11)));
      double mid1_c = mid1 * pow(1-SE(0),c1)*pow(SP(0),n-1-c1)*(1-SP(1))*pow(1-SE(1),c2)*pow(SP(1),n-c2); 
      
      double mid2 = p(1)*(theta(p,n-1,10)+theta(p,n-1,00))*g12(c1,c2,n-1,theta(p,n,00), theta(p,n,10), theta(p,n,01), theta(p,n,11)) + 
        (p(3)+p(1)*(theta(p,n-1,01)+theta(p,n-1,11)))*g12(c1,c2-1,n-1,theta(p,n,00), theta(p,n,10), theta(p,n,01), theta(p,n,11)) - mid1;
      double mid2_c = mid2 * pow(1-SE(0),c1)*pow(SP(0),n-1-c1)*SE(1)*pow(1-SE(1),c2)*pow(SP(1),n-c2); 
      
      res += mid1_c+mid2_c;
    }
  }
  res = pow(SE(0),2)* (1-SE(0))*pow(p(1)+p(3),-1)*res;
  return res;
}

double pse1_nomaster_3_5(vec p, vec SE, vec SP, int n){
  double res1 = 0;
  for (int c1=0;c1<=n-1;c1++){
    for (int r2=0;r2<=n-1;r2++){
      double mid = 0;
      for (int w1=0;w1<=c1;w1++){
        for (int w2=0;w2<=r2;w2++){
          mid += R::choose(n-1-w1,c1-w1)*pow(p(1)+p(3),c1-w1)*pow(p(0)+p(2),n-c1-1)*R::choose(n-1-w2,r2-w2)*pow(p(2)+p(3),r2-w2)*pow(p(0)+p(1),n-r2-1)*f(p,n-1,w1,n-1,w2);
        }
      }
      res1 += mid *pow(1-SE(0),c1)*pow(SP(0),n-1-c1)*pow(1-SE(1),r2)*pow(SP(1),n-1-r2);
    }
  }
  res1 = res1 * p(3)*pow(SE(0),2)*(1-SE(0))*SE(1)*(1-SE(1));
  
  double res2 = 0; 
  for (int r2=0; r2<=n;r2++){
    for (int c1=0; c1<=n-1;c1++){
      double mid1 = p(1)*(theta(p,n-1,10)+theta(p,n-1,00))*f(p,n-1,c1,n,r2);
      double mid1_c = mid1* pow(1-SE(0),c1)*pow(SP(0),n-1-c1)*pow(1-SE(1),r2)*pow(SP(1),n-r2)*(1-SP(1));
      
      double mid2 = 0;
      for (int w1=0;w1<=c1;w1++){
        for (int w2=1; w2<=min(n-1,r2); w2++){
          mid2 += R::choose(n-1,w2)*pow(p(2)+p(3),w2)*pow(p(0)+p(1),n-w2-1)*f(p,n-1,w1,n-w2,r2-w2)*R::choose(n-w1-1,c1-w1)*pow(theta(p,w2,10)+theta(p,w2,11),c1-w1)*pow(theta(p,w2,00)+theta(p,w2,01),n-1-c1);
        }
      }
      mid2 = mid2 *p(1);
      double mid2_c = mid2 * pow(1-SE(0),c1)*pow(SP(0),n-1-c1)*pow(1-SE(1),r2)*pow(SP(1),n-r2)*SE(1);
      res2 += pow(SE(0),2)*(1-SE(0))*(mid1_c+mid2_c);
    }
  }
  double res = pow(p(1)+p(3),-1)*(res1+res2);
  return res;
}

double pse1_nomaster(vec p, vec SE, vec SP, int n){
  double res = pse1_nomaster_1(p,SE,SP,n) + pse1_nomaster_2(p,SE,SP,n) - (pse1_nomaster_3_1(p,SE,SP,n)+2*pse1_nomaster_3_2(p,SE,SP,n)+2*pse1_nomaster_3_3(p,SE,SP,n)+2*pse1_nomaster_3_4(p,SE,SP,n)+2*pse1_nomaster_3_5(p,SE,SP,n));
  return res;
}


/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////PSe2//////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


double pse2_nomaster_1(vec p, vec SE, vec SP, int n){
  double res = pow(SE(1),3)+2*pow(SE(1),2)*(1-SE(1))*pow(1-f2(p,SE,SP,n),n-1);
  return res;
}

double pse2_nomaster_2(vec p, vec SE, vec SP, int n){
  double res_1 = SE(1)*pow(p(2)+p(3),-1)*(p(2)*pow(f1(p,SE,SP,n-1),2)+p(3)*pow(SE(0),2));
  double res_2 = 0;
  for (int c1=0; c1<=n; c1++){
    double mid1 = p(2)*(theta(p,n-1,01)+theta(p,n-1,00))*g1(p,c1,n-1,n);
    double mid1_C = mid1 * (1-SP(0))*pow(1-SE(0),c1)*pow(SP(0),n-c1);
    
    double mid2 = p(2)*(theta(p,n-1,01)+theta(p,n-1,00))*g1(p,c1,n,n-1)+p(2)*(theta(p,n-1,10)+theta(p,n-1,11))*g1(p,c1-1,n,n-1)-mid1;
    double mid2_C = mid2 *SE(0)*pow(1-SE(0),c1)*pow(SP(0),n-c1);
    
    res_2 += mid1_C + mid2_C; 
  }
  res_2 = pow(p(2)+p(3),-1)*(p(3)*SE(0)*(1-SE(0))*pow(1-f1(p,SE,SP,n),n-1)+res_2);
  double res  = res_1 + 2* SE(1)*res_2; 
  return res;
}

double pse2_nomaster_3_1(vec p, vec SE, vec SP, int n){
  double res = pow(SE(1),3)*pow(p(2)+p(3),-1)*(p(3)*pow(SE(0),2)+p(2)*pow(f1(p,SE,SP,n-1),2));
  return res;
}


double pse2_nomaster_3_2(vec p, vec SE, vec SP, int n){
  double res = 0;
  for (int c1=0; c1<=n; c1++){
    double mid1 = p(2)*(theta(p,n-1,01)+theta(p,n-1,00))*g1(p,c1,n-1,n);
    double mid1_C = mid1 * (1-SP(0))*pow(1-SE(0),c1)*pow(SP(0),n-c1);
    
    double mid2 = p(2)*(theta(p,n-1,01)+theta(p,n-1,00))*g1(p,c1,n,n-1)+p(2)*(theta(p,n-1,10)+theta(p,n-1,11))*g1(p,c1-1,n,n-1)-mid1;
    double mid2_C = mid2 *SE(0)*pow(1-SE(0),c1)*pow(SP(0),n-c1);
    
    res += mid1_C + mid2_C; 
  }
  res = pow(SE(1),3)*pow(p(2)+p(3),-1)*(p(3)*SE(0)*(1-SE(0))*pow(1-f1(p,SE,SP,n),n-1)+res);
  return res;
}

double pse2_nomaster_3_3(vec p, vec SE, vec SP, int n){
  double res = 0;
  for (int c2=0; c2<=n-1;c2++){
    double mid1 = 0;
    for (int w=0; w<=c2;w++){
      mid1 += R::choose(n-1,w)*pow(p(2),w)*pow(p(0),n-1-w)*g2(p,c2-w,n-1,n-w-1);
    }
    mid1 = mid1*p(2);
    double mid1_c = mid1*(1-SP(0))*pow(SP(1),n-1-c2)*pow(1-SE(1),c2);
    
    double mid2 = p(2)*g2(p,c2,n,n-1)-mid1;
    double mid2_c = mid2*SE(0)*pow(SP(1),n-1-c2)*pow(1-SE(1),c2);
    
    res += mid1_c+mid2_c;
  }
  res = pow(p(2)+p(3),-1)*pow(SE(1),2)*(1-SE(1))*(p(3)*pow(SE(0),2)*pow(1-f2(p,SE,SP,n),n-1)+f1(p,SE,SP,n-1)*res);
  return res;
}


double pse2_nomaster_3_4(vec p, vec SE, vec SP, int n){
  double res = 0;
  for (int c1=0; c1<=n;c1++){
    for (int c2=0; c2<=n-1; c2++){
      double mid1 = p(2)*(theta(p,n-1,01)+theta(p,n-1,00))*g12(c1,c2,n-1,theta(p,n,00), p(0)*theta(p,n-1,10),theta(p,n,01), p(0)*theta(p,n-1,11)+p(2)*(theta(p,n-1,10)+theta(p,n-1,11))) + 
        p(2)*(theta(p,n-1,10)+theta(p,n-1,11))*g12(c1-1,c2,n-1,theta(p,n,00), p(0)*theta(p,n-1,10),theta(p,n,01),p(0)*theta(p,n-1,11)+p(2)*(theta(p,n-1,10)+theta(p,n-1,11)));
      double mid1_c = mid1 * pow(1-SE(1),c2)*pow(SP(1),n-1-c2)*(1-SP(0))*pow(1-SE(0),c1)*pow(SP(0),n-c1); 
      
      double mid2 = p(2)*(theta(p,n-1,01)+theta(p,n-1,00))*g12(c1,c2,n-1,theta(p,n,00), theta(p,n,10), theta(p,n,01), theta(p,n,11)) + 
        (p(3)+p(2)*(theta(p,n-1,10)+theta(p,n-1,11)))*g12(c1-1,c2,n-1,theta(p,n,00), theta(p,n,10), theta(p,n,01), theta(p,n,11)) - mid1;
      double mid2_c = mid2 * pow(1-SE(1),c2)*pow(SP(1),n-1-c2)*SE(0)*pow(1-SE(0),c1)*pow(SP(0),n-c1); 
      
      res += mid1_c+mid2_c;
    }
  }
  res = pow(SE(1),2)* (1-SE(1))*pow(p(2)+p(3),-1)*res;
  return res;
}


double pse2_nomaster_3_5(vec p, vec SE, vec SP, int n){
  double res1 = 0;
  for (int c1=0;c1<=n-1;c1++){
    for (int r2=0;r2<=n-1;r2++){
      double mid = 0;
      for (int w1=0;w1<=c1;w1++){
        for (int w2=0;w2<=r2;w2++){
          mid += R::choose(n-1-w1,c1-w1)*pow(p(1)+p(3),c1-w1)*pow(p(0)+p(2),n-c1-1)*R::choose(n-1-w2,r2-w2)*pow(p(2)+p(3),r2-w2)*pow(p(0)+p(1),n-r2-1)*f(p,n-1,w1,n-1,w2);
        }
      }
      res1 += mid *pow(1-SE(0),c1)*pow(SP(0),n-1-c1)*pow(1-SE(1),r2)*pow(SP(1),n-1-r2);
    }
  }
  res1 = res1 * p(3)*pow(SE(1),2)*(1-SE(1))*SE(0)*(1-SE(0));
  
  double res2 = 0; 
  for (int r2=0; r2<=n-1;r2++){
    for (int c1=0; c1<=n;c1++){
      double mid1 = p(2)*(theta(p,n-1,01)+theta(p,n-1,00))*f(p,n,c1,n-1,r2);
      double mid1_c = mid1* pow(1-SE(1),r2)*pow(SP(1),n-1-r2)*pow(1-SE(0),c1)*pow(SP(0),n-c1)*(1-SP(0));
      
      double mid2 = 0;
      for (int w2=0;w2<=r2;w2++){
        for (int w1=1; w1<=min(n-1,c1); w1++){
          mid2 += R::choose(n-1,w1)*pow(p(1)+p(3),w1)*pow(p(0)+p(2),n-w1-1)*f(p,n-w1,c1-w1,n-1,w2)*R::choose(n-w2-1,r2-w2)*pow(theta(p,w1,01)+theta(p,w1,11),r2-w2)*pow(theta(p,w1,10)+theta(p,w1,00),n-1-r2);
        }
      }
      mid2 = mid2 *p(2);
      double mid2_c = mid2 * pow(1-SE(1),r2)*pow(SP(1),n-1-r2)*pow(1-SE(0),c1)*pow(SP(0),n-c1)*SE(0);
      res2 += pow(SE(1),2)*(1-SE(1))*(mid1_c+mid2_c);
    }
  }
  double res = pow(p(2)+p(3),-1)*(res1+res2);
  return res;
}


////////////////////////////////////////////////////////////////////////////////////////
/////////////////// Efficiency with master pool ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

double g1n(vec p, int c , int s1, int s2){
  double res = R::choose(s2,c)*pow(theta(p,s1,10),c)*pow(theta(p,s1,00),s2-c);
  return res;
}

double g2n(vec p, int c , int s1, int s2){
  double res = R::choose(s2,c)*pow(theta(p,s1,01),c)*pow(theta(p,s1,00),s2-c);
  return res;
}

double alpha2(vec p, vec SE,vec SP, int n){
  double res = p(1)*pow(SE(0),2)*pow(p(0)+p(1),n*n-1) + p(0)*pow(theta(p,n-1,00)*(1-SP(0))+theta(p,n-1,10)*SE(0),2)*pow(p(0)+p(1),n*n-2*n+1);
  return res;
}

double alpha1(vec p, vec SE,vec SP, int n){
  double res = p(2)*pow(SE(1),2)*pow(p(0)+p(2),n*n-1) + p(0)*pow(theta(p,n-1,00)*(1-SP(1))+theta(p,n-1,01)*SE(1),2)*pow(p(0)+p(2),n*n-2*n+1);
  return res;
}

double beta2(vec p, vec SE, vec SP, int n){
  double res = 0;
  for (int c1=0;c1<=n;c1++){
    double mid1 = theta(p,n,00)*g1n(p,c1,n-1,n);
    double mid1_c = mid1*(1-SP(0))*pow(1-SE(0), c1)*pow(SP(0),n-c1);
    
    double mid2 = g1n(p,c1,n,n) - mid1;
    double mid2_c = mid2*SE(0)*pow(1-SE(0), c1)*pow(SP(0),n-c1);
    
    res += mid1_c + mid2_c;
  }
  return res;
}

double beta1(vec p, vec SE, vec SP, int n){
  double res = 0;
  for (int c2=0;c2<=n;c2++){
    double mid1 = theta(p,n,00)*g2n(p,c2,n-1,n);
    double mid1_c = mid1*(1-SP(1))*pow(1-SE(1), c2)*pow(SP(1),n-c2);
    
    double mid2 = g2n(p,c2,n,n) - mid1;
    double mid2_c = mid2*SE(1)*pow(1-SE(1), c2)*pow(SP(1),n-c2);
    
    res += mid1_c + mid2_c;
  }
  return res;
}


//First probability
double eff_master_1(vec p, vec SE, vec SP, int n, double et_nomaster_1){
  double mid1 = theta(p,n*n,00)*(pow(1-SP(0),2)+2*(1-SP(0))*pow(SP(0),n));
  double mid2 = alpha2(p,SE,SP,n) + 2*beta2(p,SE,SP,n) - mid1;
  double mid3 = theta(p,n*n,01)*(pow(1-SP(0),2)+2*(1-SP(0))*pow(SP(0),n));
  double mid4 = et_nomaster_1 - mid1 - mid2 - mid3; 
  double res = mid1*(1-SP(0)*SP(1)) + mid2*(1-SP(1)*(1-SE(0))) + mid3*(1-SP(0)*(1-SE(1)))+mid4*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//Second probability
double eff_master_2(vec p, vec SE, vec SP, int n, double et_nomaster_2){
  double mid1 = theta(p,n*n,00)*(pow(1-SP(1),2)+2*(1-SP(1))*pow(SP(1),n));
  double mid2 = alpha1(p,SE,SP,n) + 2*beta1(p,SE,SP,n) - mid1;
  double mid3 = theta(p,n*n,10)*(pow(1-SP(1),2)+2*(1-SP(1))*pow(SP(1),n));
  double mid4 = et_nomaster_2 - mid1 - mid2 - mid3; 
  double res = mid1*(1-SP(0)*SP(1)) + mid3*(1-SP(1)*(1-SE(0))) + mid2*(1-SP(0)*(1-SE(1)))+mid4*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//Intersection probability 3-1
double eff_master_3_1(vec p, vec SE, vec SP, int n, double et_nomaster_3_1){
  double mid1 = theta(p,n*n,00)*pow((1-SP(0))*(1-SP(1)),2);
  double mid2 = pow((1-SP(1)),2)*alpha2(p,SE,SP,n) - mid1;
  double mid3 = pow((1-SP(0)),2)*alpha1(p,SE,SP,n) - mid1;
  double mid4 = et_nomaster_3_1 - mid1 - mid2 - mid3; 
  double res = mid1*(1-SP(0)*SP(1)) + mid2*(1-SP(1)*(1-SE(0))) + mid3*(1-SP(0)*(1-SE(1)))+mid4*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//Intersection probability 3-2
double eff_master_3_2(vec p, vec SE, vec SP, int n, double et_nomaster_3_2){
  double mid1 = theta(p,n*n,00)*pow(1-SP(0),2)*(1-SP(1))*pow(SP(1),n);
  double mid2 = (1-SP(1))*pow(SP(1),n)*alpha2(p,SE,SP,n) - mid1;
  double mid3 = pow(1-SP(0),2)*beta1(p,SE,SP,n) - mid1;
  double mid4 = et_nomaster_3_2 - mid1 - mid2 - mid3; 
  double res = mid1*(1-SP(0)*SP(1)) + mid2*(1-SP(1)*(1-SE(0))) + mid3*(1-SP(0)*(1-SE(1)))+mid4*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//Intersection probability 3-3
double eff_master_3_3(vec p, vec SE, vec SP, int n, double et_nomaster_3_3){
  double mid1 = theta(p,n*n,00)*pow(1-SP(1),2)*(1-SP(0))*pow(SP(0),n);
  double mid2 = (1-SP(0))*pow(SP(0),n)*alpha1(p,SE,SP,n) - mid1;
  double mid3 = pow(1-SP(1),2)*beta2(p,SE,SP,n) - mid1;
  double mid4 = et_nomaster_3_3 - mid1 - mid2 - mid3; 
  double res = mid1*(1-SP(0)*SP(1)) + mid3*(1-SP(1)*(1-SE(0))) + mid2*(1-SP(0)*(1-SE(1)))+mid4*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//Intersection probability 3-4
double eff_master_3_4(vec p, vec SE, vec SP, int n, double et_nomaster_3_4){
  double mid1 = theta(p,n*n,00)*(1-SP(0))*pow(SP(0),n)*(1-SP(1))*pow(SP(1),n);
  double mid2 = (1-SP(1))*pow(SP(1),n)*beta2(p,SE,SP,n) - mid1;
  double mid3 = (1-SP(0))*pow(SP(0),n)*beta1(p,SE,SP,n) - mid1;
  double mid4 = et_nomaster_3_4 - mid1 - mid2 - mid3; 
  double res = mid1*(1-SP(0)*SP(1)) + mid2*(1-SP(1)*(1-SE(0))) + mid3*(1-SP(0)*(1-SE(1)))+mid4*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//Intersection probability 3-5
double eff_master_3_5(vec p, vec SE, vec SP, int n, double et_nomaster_3_5){
  double mid1 = theta(p,n*n,00)*(1-SP(0))*pow(SP(0),n)*(1-SP(1))*pow(SP(1),n);
  double mid2 = (1-SP(1))*pow(SP(1),n)*beta2(p,SE,SP,n) - mid1;
  double mid3 = (1-SP(0))*pow(SP(0),n)*beta1(p,SE,SP,n) - mid1;
  double mid4 = et_nomaster_3_5 - mid1 - mid2 - mid3; 
  double res = mid1*(1-SP(0)*SP(1)) + mid2*(1-SP(1)*(1-SE(0))) + mid3*(1-SP(0)*(1-SE(1)))+mid4*(1-(1-SE(0))*(1-SE(1)));
  return res;
}


////////////////////////////////////////////////////////////////////////////////////////
/////////////////// PSE1 with master pool ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

//First 
double pse1_master_1(vec p, vec SE, vec SP, int n, double PSE1_nomaster_1){
  double mid1 = pow(p(1)+p(3),-1)*p(1)*pow(p(0)+p(1),n*n-1)*pow(SE(0),3) + 2*pow(p(1)+p(3),-1)*p(1)*pow(p(0)+p(1),n-1)*pow(theta(p,n,10)*(1-SE(0))+theta(p,n,00)*SP(0),n-1)*pow(SE(0),2)*(1-SE(0));
  double mid2 = PSE1_nomaster_1-mid1; 
  double res = mid1*(1-(1-SE(0))*SP(1))+mid2*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//second 
double pse1_master_2(vec p, vec SE, vec SP, int n, double PSE1_nomaster_2){
  double mid1 = pow(p(1)+p(3),-1)*p(1)*pow(p(0)+p(1),n*n-1)*(SE(0)*pow(1-SP(1),2)+2*SE(0)*(1-SP(1))*pow(SP(1),n));
  double mid2 = PSE1_nomaster_2-mid1; 
  double res = mid1*(1-(1-SE(0))*SP(1))+mid2*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//Intersection 3-1 
double pse1_master_3_1(vec p, vec SE, vec SP, int n, double PSE1_nomaster_3_1){
  double mid1 = pow(p(1)+p(3),-1)*p(1)*pow(p(0)+p(1),n*n-1)*pow(SE(0),3)*pow(1-SP(1),2);
  double mid2 = PSE1_nomaster_3_1-mid1; 
  double res = mid1*(1-(1-SE(0))*SP(1))+mid2*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//Intersection 3-2 
double pse1_master_3_2(vec p, vec SE, vec SP, int n, double PSE1_nomaster_3_2){
  double mid1 = pow(p(1)+p(3),-1)*p(1)*pow(p(0)+p(1),n*n-1)*pow(SE(0),3)*pow(SP(1),n)*(1-SP(1));
  double mid2 = PSE1_nomaster_3_2-mid1; 
  double res = mid1*(1-(1-SE(0))*SP(1))+mid2*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//Intersection 3-3
double pse1_master_3_3(vec p, vec SE, vec SP, int n, double PSE1_nomaster_3_3){
  double mid1 = pow(p(1)+p(3),-1)*p(1)*pow(p(0)+p(1),n-1)*pow(theta(p,n,10)*(1-SE(0))+theta(p,n,00)*SP(0),n-1)*pow(SE(0),2)*(1-SE(0))*pow(1-SP(1),2);
  double mid2 = PSE1_nomaster_3_3-mid1; 
  double res = mid1*(1-(1-SE(0))*SP(1))+mid2*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//Intersection 3-4
double pse1_master_3_4(vec p, vec SE, vec SP, int n, double PSE1_nomaster_3_4){
  double mid1 = pow(p(1)+p(3),-1)*p(1)*pow(p(0)+p(1),n-1)*pow(theta(p,n,10)*(1-SE(0))+theta(p,n,00)*SP(0),n-1)*pow(SE(0),2)*(1-SE(0))*(1-SP(1))*pow(SP(1),n);
  double mid2 = PSE1_nomaster_3_4-mid1; 
  double res = mid1*(1-(1-SE(0))*SP(1))+mid2*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//Intersection 3-5
double pse1_master_3_5(vec p, vec SE, vec SP, int n, double PSE1_nomaster_3_5){
  double mid1 = pow(p(1)+p(3),-1)*p(1)*pow(p(0)+p(1),n-1)*pow(theta(p,n,10)*(1-SE(0))+theta(p,n,00)*SP(0),n-1)*pow(SE(0),2)*(1-SE(0))*(1-SP(1))*pow(SP(1),n);
  double mid2 = PSE1_nomaster_3_5-mid1; 
  double res = mid1*(1-(1-SE(0))*SP(1))+mid2*(1-(1-SE(0))*(1-SE(1)));
  return res;
}


////////////////////////////////////////////////////////////////////////////////////////
/////////////////// PSE2 with master pool ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

//First 
double pse2_master_1(vec p, vec SE, vec SP, int n, double PSE2_nomaster_1){
  double mid1 = pow(p(2)+p(3),-1)*p(2)*pow(p(0)+p(2),n*n-1)*pow(SE(1),3) + 2*pow(p(2)+p(3),-1)*p(2)*pow(p(0)+p(2),n-1)*pow(theta(p,n,01)*(1-SE(1))+theta(p,n,00)*SP(1),n-1)*pow(SE(1),2)*(1-SE(1));
  double mid2 = PSE2_nomaster_1-mid1; 
  double res = mid1*(1-(1-SE(1))*SP(0))+mid2*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//second 
double pse2_master_2(vec p, vec SE, vec SP, int n, double PSE2_nomaster_2){
  double mid1 = pow(p(2)+p(3),-1)*p(2)*pow(p(0)+p(2),n*n-1)*(SE(1)*pow(1-SP(0),2)+2*SE(1)*(1-SP(0))*pow(SP(0),n));
  double mid2 = PSE2_nomaster_2-mid1; 
  double res = mid1*(1-(1-SE(1))*SP(0))+mid2*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//Intersection 3-1 
double pse2_master_3_1(vec p, vec SE, vec SP, int n, double PSE2_nomaster_3_1){
  double mid1 = pow(p(2)+p(3),-1)*p(2)*pow(p(0)+p(2),n*n-1)*pow(SE(1),3)*pow(1-SP(0),2);
  double mid2 = PSE2_nomaster_3_1-mid1; 
  double res = mid1*(1-(1-SE(1))*SP(0))+mid2*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//Intersection 3-2 
double pse2_master_3_2(vec p, vec SE, vec SP, int n, double PSE2_nomaster_3_2){
  double mid1 = pow(p(2)+p(3),-1)*p(2)*pow(p(0)+p(2),n*n-1)*pow(SE(1),3)*pow(SP(0),n)*(1-SP(0));
  double mid2 = PSE2_nomaster_3_2-mid1; 
  double res = mid1*(1-(1-SE(1))*SP(0))+mid2*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//Intersection 3-3
double pse2_master_3_3(vec p, vec SE, vec SP, int n, double PSE2_nomaster_3_3){
  double mid1 = pow(p(2)+p(3),-1)*p(2)*pow(p(0)+p(2),n-1)*pow(theta(p,n,01)*(1-SE(1))+theta(p,n,00)*SP(1),n-1)*pow(SE(1),2)*(1-SE(1))*pow(1-SP(0),2);
  double mid2 = PSE2_nomaster_3_3-mid1; 
  double res = mid1*(1-(1-SE(1))*SP(0))+mid2*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//Intersection 3-4
double pse2_master_3_4(vec p, vec SE, vec SP, int n, double PSE2_nomaster_3_4){
  double mid1 = pow(p(2)+p(3),-1)*p(2)*pow(p(0)+p(2),n-1)*pow(theta(p,n,01)*(1-SE(1))+theta(p,n,00)*SP(1),n-1)*pow(SE(1),2)*(1-SE(1))*(1-SP(0))*pow(SP(0),n);
  double mid2 = PSE2_nomaster_3_4-mid1; 
  double res = mid1*(1-(1-SE(1))*SP(0))+mid2*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

//Intersection 3-5
double pse2_master_3_5(vec p, vec SE, vec SP, int n, double PSE2_nomaster_3_5){
  double mid1 = pow(p(2)+p(3),-1)*p(2)*pow(p(0)+p(2),n-1)*pow(theta(p,n,01)*(1-SE(1))+theta(p,n,00)*SP(1),n-1)*pow(SE(1),2)*(1-SE(1))*(1-SP(0))*pow(SP(0),n);
  double mid2 = PSE2_nomaster_3_5-mid1; 
  double res = mid1*(1-(1-SE(1))*SP(0))+mid2*(1-(1-SE(0))*(1-SE(1)));
  return res;
}

////////////////////////////////////////////////////////////////////////////////////////
/////////////////// Output /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List ARRAY(vec p, vec SE, vec SP, double n){
  double et_nomaster_1 = eff_nomaster_1(p,SE,SP,n);
  double et_nomaster_2 = eff_nomaster_2(p,SE,SP,n);
  double et_nomaster_3_1 = eff_nomaster_3_1(p,SE,SP,n);
  double et_nomaster_3_2 = eff_nomaster_3_2(p,SE,SP,n);
  double et_nomaster_3_3 = eff_nomaster_3_3(p,SE,SP,n);
  double et_nomaster_3_4 = eff_nomaster_3_4(p,SE,SP,n);
  double et_nomaster_3_5 = eff_nomaster_3_5(p,SE,SP,n);
  
  double PSE1_nomaster_1 = pse1_nomaster_1(p,SE,SP,n);     double PSE2_nomaster_1 = pse2_nomaster_1(p,SE,SP,n); 
  double PSE1_nomaster_2 = pse1_nomaster_2(p,SE,SP,n);     double PSE2_nomaster_2 = pse2_nomaster_2(p,SE,SP,n);
  double PSE1_nomaster_3_1 = pse1_nomaster_3_1(p,SE,SP,n); double PSE2_nomaster_3_1 = pse2_nomaster_3_1(p,SE,SP,n);
  double PSE1_nomaster_3_2 = pse1_nomaster_3_2(p,SE,SP,n); double PSE2_nomaster_3_2 = pse2_nomaster_3_2(p,SE,SP,n);
  double PSE1_nomaster_3_3 = pse1_nomaster_3_3(p,SE,SP,n); double PSE2_nomaster_3_3 = pse2_nomaster_3_3(p,SE,SP,n);
  double PSE1_nomaster_3_4 = pse1_nomaster_3_4(p,SE,SP,n); double PSE2_nomaster_3_4 = pse2_nomaster_3_4(p,SE,SP,n);
  double PSE1_nomaster_3_5 = pse1_nomaster_3_5(p,SE,SP,n); double PSE2_nomaster_3_5 = pse2_nomaster_3_5(p,SE,SP,n);
  
  double et_master_1 = eff_master_1(p,SE,SP,n,et_nomaster_1);
  double et_master_2 = eff_master_2(p,SE,SP,n,et_nomaster_2);
  double et_master_3_1 = eff_master_3_1(p,SE,SP,n,et_nomaster_3_1);
  double et_master_3_2 = eff_master_3_2(p,SE,SP,n,et_nomaster_3_2);
  double et_master_3_3 = eff_master_3_3(p,SE,SP,n,et_nomaster_3_3);
  double et_master_3_4 = eff_master_3_4(p,SE,SP,n,et_nomaster_3_4);
  double et_master_3_5 = eff_master_3_5(p,SE,SP,n,et_nomaster_3_5);
  
  double PSE1_master_1 = pse1_master_1(p,SE,SP,n,PSE1_nomaster_1);     double PSE2_master_1 = pse2_master_1(p,SE,SP,n,PSE2_nomaster_1); 
  double PSE1_master_2 = pse1_master_2(p,SE,SP,n,PSE1_nomaster_2);     double PSE2_master_2 = pse2_master_2(p,SE,SP,n,PSE2_nomaster_2);
  double PSE1_master_3_1 = pse1_master_3_1(p,SE,SP,n,PSE1_nomaster_3_1); double PSE2_master_3_1 = pse2_master_3_1(p,SE,SP,n,PSE2_nomaster_3_1);
  double PSE1_master_3_2 = pse1_master_3_2(p,SE,SP,n,PSE1_nomaster_3_2); double PSE2_master_3_2 = pse2_master_3_2(p,SE,SP,n,PSE2_nomaster_3_2);
  double PSE1_master_3_3 = pse1_master_3_3(p,SE,SP,n,PSE1_nomaster_3_3); double PSE2_master_3_3 = pse2_master_3_3(p,SE,SP,n,PSE2_nomaster_3_3);
  double PSE1_master_3_4 = pse1_master_3_4(p,SE,SP,n,PSE1_nomaster_3_4); double PSE2_master_3_4 = pse2_master_3_4(p,SE,SP,n,PSE2_nomaster_3_4);
  double PSE1_master_3_5 = pse1_master_3_5(p,SE,SP,n,PSE1_nomaster_3_5); double PSE2_master_3_5 = pse2_master_3_5(p,SE,SP,n,PSE2_nomaster_3_5);
  
  double efficiency_nomaster = 2/n+et_nomaster_1+et_nomaster_2 - (et_nomaster_3_1 + 2*et_nomaster_3_2+2*et_nomaster_3_3+2*et_nomaster_3_4 + 2*et_nomaster_3_5);
  double PSE1_nomaster = PSE1_nomaster_1 + PSE1_nomaster_2 - (PSE1_nomaster_3_1+2*PSE1_nomaster_3_2+2*PSE1_nomaster_3_3+2*PSE1_nomaster_3_4+2*PSE1_nomaster_3_5);
  double PSE2_nomaster = PSE2_nomaster_1 + PSE2_nomaster_2 - (PSE2_nomaster_3_1+2*PSE2_nomaster_3_2+2*PSE2_nomaster_3_3+2*PSE2_nomaster_3_4+2*PSE2_nomaster_3_5);
  double PSP1_nomaster = 1-(1-SP(0))*pow(p(0)+p(2),-1)*(efficiency_nomaster-2/n-(p(1)+p(3))*PSE1_nomaster/SE(0));
  double PSP2_nomaster = 1-(1-SP(1))*pow(p(0)+p(1),-1)*(efficiency_nomaster-2/n-(p(2)+p(3))*PSE2_nomaster/SE(1));
  double PPV1_nomaster = (p(1)+p(3))*PSE1_nomaster/((p(1)+p(3))*PSE1_nomaster+(p(0)+p(2))*(1-PSP1_nomaster));
  double PPV2_nomaster = (p(2)+p(3))*PSE2_nomaster/((p(2)+p(3))*PSE2_nomaster+(p(0)+p(1))*(1-PSP2_nomaster));
  double NPV1_nomaster = (p(0)+p(2))*PSP1_nomaster/((p(0)+p(2))*PSP1_nomaster+(p(1)+p(3))*(1-PSE1_nomaster));
  double NPV2_nomaster = (p(0)+p(1))*PSP2_nomaster/((p(0)+p(1))*PSP2_nomaster+(p(2)+p(3))*(1-PSE2_nomaster));
  
  double efficiency_master = 1/n/n + 2/n*(1-theta(p,n*n,00)*SP(0)*SP(1)-theta(p,n*n,10)*(1-SE(0))*SP(1)-theta(p,n*n,01)*SP(0)*(1-SE(1))-theta(p,n*n,11)*(1-SE(0))*(1-SE(1))) + et_master_1+et_master_2 - (et_master_3_1 + 2*et_master_3_2+2*et_master_3_3+2*et_master_3_4 + 2*et_master_3_5);
  double PSE1_master = PSE1_master_1 + PSE1_master_2 - (PSE1_master_3_1+2*PSE1_master_3_2+2*PSE1_master_3_3+2*PSE1_master_3_4+2*PSE1_master_3_5);
  double PSE2_master = PSE2_master_1 + PSE2_master_2 - (PSE2_master_3_1+2*PSE2_master_3_2+2*PSE2_master_3_3+2*PSE2_master_3_4+2*PSE2_master_3_5);
  double PSP1_master = 1- (1-SP(0))*pow(p(0)+p(2),-1)*(efficiency_master-1/n/n - 2/n*(1-theta(p,n*n,00)*SP(0)*SP(1)-theta(p,n*n,10)*(1-SE(0))*SP(1)-theta(p,n*n,01)*SP(0)*(1-SE(1))-theta(p,n*n,11)*(1-SE(0))*(1-SE(1))) -(p(1)+p(3))*PSE1_master/SE(0));
  double PSP2_master = 1- (1-SP(1))*pow(p(0)+p(1),-1)*(efficiency_master-1/n/n - 2/n*(1-theta(p,n*n,00)*SP(0)*SP(1)-theta(p,n*n,10)*(1-SE(0))*SP(1)-theta(p,n*n,01)*SP(0)*(1-SE(1))-theta(p,n*n,11)*(1-SE(0))*(1-SE(1))) -(p(2)+p(3))*PSE1_master/SE(1));
  double PPV1_master = (p(1)+p(3))*PSE1_master/((p(1)+p(3))*PSE1_master+(p(0)+p(2))*(1-PSP1_master));
  double PPV2_master = (p(2)+p(3))*PSE2_master/((p(2)+p(3))*PSE2_master+(p(0)+p(1))*(1-PSP2_master));
  double NPV1_master = (p(0)+p(2))*PSP1_master/((p(0)+p(2))*PSP1_master+(p(1)+p(3))*(1-PSE1_master));
  double NPV2_master = (p(0)+p(1))*PSP2_master/((p(0)+p(1))*PSP2_master+(p(2)+p(3))*(1-PSE2_master));
  
  return Rcpp::List::create(Rcpp::Named("ET_AT")=efficiency_nomaster,
                            Rcpp::Named("PSE1_AT")=PSE1_nomaster,
                            Rcpp::Named("PSE2_AT")=PSE2_nomaster,
                            Rcpp::Named("PSP1_AT")=PSP1_nomaster,
                            Rcpp::Named("PSP2_AT")=PSP2_nomaster,
                            Rcpp::Named("PPV1_AT")=PPV1_nomaster,
                            Rcpp::Named("PPV2_AT")=PPV2_nomaster,
                            Rcpp::Named("NPV1_AT")=NPV1_nomaster,
                            Rcpp::Named("NPV2_AT")=NPV2_nomaster,
                            Rcpp::Named("ET_ATM")=efficiency_master,
                            Rcpp::Named("PSE1_ATM")=PSE1_master,
                            Rcpp::Named("PSE2_ATM")=PSE2_master,
                            Rcpp::Named("PSP1_ATM")=PSP1_master,
                            Rcpp::Named("PSP2_ATM")=PSP2_master,
                            Rcpp::Named("PPV1_ATM")=PPV1_master,
                            Rcpp::Named("PPV2_ATM")=PPV2_master,
                            Rcpp::Named("NPV1_ATM")=NPV1_master,
                            Rcpp::Named("NPV2_ATM")=NPV2_master);
  
}


//Simulation AT
// [[Rcpp::export]]
List AT_SIM(int size, vec p, vec SE, vec SP,int row, int col, int iter){
  
  int num_pools = floor(size/row/col);
  int remainder = size - num_pools*row*col;
  vec T(iter); T.fill(num_pools*(row+col)+remainder);
  vec pse1(iter,fill::zeros); vec ppv1(iter,fill::zeros);
  vec pse2(iter,fill::zeros); vec ppv2(iter,fill::zeros);
  vec psp1(iter,fill::zeros); vec npv1(iter,fill::zeros);
  vec psp2(iter,fill::zeros); vec npv2(iter,fill::zeros);
  vec x(4);  x(0)=1; x(1)=2; x(2)=3; x(3)=4;
  mat tmp1(row,col,fill::zeros);
  mat tmp2(row,col,fill::zeros);
  mat tmp1_t(row,col,fill::zeros);
  mat tmp2_t(row,col,fill::zeros);
  vec res(size,fill::zeros);
  mat status(size,2,fill::zeros);
  mat test(size,2,fill::zeros);
  bool replace=true;
  
  for (int i=0; i<iter; i++){
    vec p_sim = p;
    res = Rcpp::RcppArmadillo::sample(x, size, replace, p_sim);
    status.fill(0);
    test.fill(0);
    
    for (int j=0; j<size;j++){
      if (res(j)==1){
        status(j,0)=0;
        status(j,1)=0;
      }
      if (res(j)==2){
        status(j,0)=1;
        status(j,1)=0;
      }
      if (res(j)==3){
        status(j,0)=0;
        status(j,1)=1;
      }
      if (res(j)==4){
        status(j,0)=1;
        status(j,1)=1;
      }
    } //end of data generation
    
    double c_pos = sum(status(span(0,size-1),0));
    double c_neg = size- c_pos;
    double g_pos = sum(status(span(0,size-1),1));
    double g_neg = size- g_pos;
    
    for (int j=1; j<=num_pools; j++){
      tmp1_t.fill(0);
      tmp2_t.fill(0);
      mat tmp(row*col, 2,fill::zeros);
      tmp = status(span((j-1)*row*col, j*col*row-1),span(0,1));
      
      for (int j1=0;j1<row; j1++){
        for (int j2=0;j2<col;j2++){
          tmp1(j1,j2)=tmp(j1*col+j2,0);
          tmp2(j1,j2)=tmp(j1*col+j2,1);
        }
      }//end of generating one matrix
      
      vec row_test1(row,fill::zeros);
      vec row_test2(row,fill::zeros);
      for (int k=0;k<row;k++){
        row_test1(k)=sum(tmp1(k,span(0,col-1)))>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
        row_test2(k)=sum(tmp2(k,span(0,col-1)))>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
      }
      
      vec col_test1(col,fill::zeros);
      vec col_test2(col,fill::zeros);
      for (int k=0;k<col;k++){
        col_test1(k)=sum(tmp1(span(0,row-1),k))>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
        col_test2(k)=sum(tmp2(span(0,row-1),k))>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
      }
      
      for (int j1=0;j1<row; j1++){
        for (int j2=0;j2<col;j2++){
          if (((row_test1(j1)==1) && (col_test1(j2)==1)) ||  ((row_test1(j1)==1) && (sum(col_test1)==0)) || ((col_test1(j2)==1) && (sum(row_test1)==0)) || ((row_test2(j1)==1) && (col_test2(j2)==1)) ||  ((row_test2(j1)==1) && (sum(col_test2)==0)) || ((col_test2(j2)==1) && (sum(row_test2)==0))){
            T(i) = T(i)+1;
            tmp1_t(j1,j2) = tmp1(j1,j2)>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
            tmp2_t(j1,j2) = tmp2(j1,j2)>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
          }
          test((j-1)*row*col + j1*col + j2, 0) = tmp1_t(j1,j2);
          test((j-1)*row*col + j1*col + j2, 1) = tmp2_t(j1,j2);
        }
      }
      
      if (remainder>0)
      {
        for (int k=size-remainder; k<size;k++)
        {
          test(k,0)=status(k,0)>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
          test(k,1)=status(k,1)>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
        }
      }
    } 
    
    double pse_1c=0; double pse_2c=0; double psp_1c=0; double psp_2c=0; double testp1=0; double testp2=0; double testn1=0; double testn2=0; 
    
    for (int m=0;m<size;m++)
    {
      if ((status(m,0)==1) && (test(m,0)==1)) pse_1c=pse_1c+1;
      if ((status(m,1)==1) && (test(m,1)==1)) pse_2c=pse_2c+1;
      if ((status(m,0)==0) && (test(m,0)==0)) psp_1c=psp_1c+1;
      if ((status(m,1)==0) && (test(m,1)==0)) psp_2c=psp_2c+1;
      if (test(m,0)==1) testp1++;
      if (test(m,1)==1) testp2++;
      if (test(m,0)==0) testn1++;
      if (test(m,1)==0) testn2++;
    }
    pse1(i)=pse_1c/c_pos;
    pse2(i)=pse_2c/g_pos;
    psp1(i)=psp_1c/c_neg;
    psp2(i)=psp_2c/g_neg;
    ppv1(i)=pse_1c/testp1;
    ppv2(i)=pse_2c/testp2;
    npv1(i)=psp_1c/testn1;
    npv2(i)=psp_2c/testn2;
    
  }
  
  return Rcpp::List::create(Rcpp::Named("Number of tests")=T,
                            Rcpp::Named("PSE1")=pse1,
                            Rcpp::Named("PSE2")=pse2,
                            Rcpp::Named("PSP1")=psp1,
                            Rcpp::Named("PSP2")=psp2,
                            Rcpp::Named("PPV1")=ppv1,
                            Rcpp::Named("PPV2")=ppv2,
                            Rcpp::Named("NPV1")=npv1,
                            Rcpp::Named("NPV2")=npv2,
                            Rcpp::Named("status")=status);
}


//Simulation ATM
// [[Rcpp::export]]
List ATM_SIM(int size, vec p, vec SE, vec SP, int row, int col, int iter){
  
  int num_pools = floor(size/row/col);
  int remainder = size - num_pools*row*col;
  vec T(iter); T.fill(num_pools+remainder);
  vec pse1(iter,fill::zeros); vec ppv1(iter,fill::zeros);
  vec pse2(iter,fill::zeros); vec ppv2(iter,fill::zeros);
  vec psp1(iter,fill::zeros); vec npv1(iter,fill::zeros);
  vec psp2(iter,fill::zeros); vec npv2(iter,fill::zeros);
  vec x(4);  x(0)=1; x(1)=2; x(2)=3; x(3)=4;
  mat tmp1(row,col);
  mat tmp2(row,col);
  mat tmp1_t(row,col);
  mat tmp2_t(row,col);
  vec res(size);
  mat status(size,2,fill::zeros);
  mat test(size,2,fill::zeros);
  bool replace=true;
  
  for (int i=0; i<iter; i++){
    vec p_sim = p;
    res = Rcpp::RcppArmadillo::sample(x, size, replace, p_sim);
    status.fill(0);
    test.fill(0);
    
    for (int j=0; j<size;j++){
      if (res(j)==1){
        status(j,0)=0;
        status(j,1)=0;
      }
      if (res(j)==2){
        status(j,0)=1;
        status(j,1)=0;
      }
      if (res(j)==3){
        status(j,0)=0;
        status(j,1)=1;
      }
      if (res(j)==4){
        status(j,0)=1;
        status(j,1)=1;
      }
    } //end of data generation
    
    double c_pos = sum(status(span(0,size-1),0));
    double c_neg = size- c_pos;
    double g_pos = sum(status(span(0,size-1),1));
    double g_neg = size- g_pos;
    
    for (int j=1; j<=num_pools; j++){
      tmp1_t.fill(0);
      tmp2_t.fill(0);
      
      mat tmp(row*col, 2,fill::zeros);
      tmp = status(span((j-1)*row*col, j*col*row-1),span(0,1));
      
      for (int j1=0;j1<row; j1++){
        for (int j2=0;j2<col;j2++){
          tmp1(j1,j2)=tmp(j1*col+j2,0);
          tmp2(j1,j2)=tmp(j1*col+j2,1);
        }
      }//end of generating one matrix
      
      int M1 = sum(tmp(span(0,col*row-1),0))>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
      int M2 = sum(tmp(span(0,col*row-1),1))>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
      
      if (M1+M2==0) {
        
        for (int j1=0;j1<row; j1++){
          for (int j2=0;j2<col;j2++){
            test((j-1)*row*col + j1*col + j2, 0) = 0;
            test((j-1)*row*col + j1*col + j2, 1) = 0;
          }
        }
      }
      
      else {
        T(i) = T(i)+col+row;
        
        vec row_test1(row,fill::zeros);
        vec row_test2(row,fill::zeros);
        for (int k=0;k<row;k++){
          row_test1(k)=sum(tmp1(k,span(0,col-1)))>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
          row_test2(k)=sum(tmp2(k,span(0,col-1)))>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
        }
        
        vec col_test1(col);
        vec col_test2(col);
        for (int k=0;k<col;k++){
          col_test1(k)=sum(tmp1(span(0,row-1),k))>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
          col_test2(k)=sum(tmp2(span(0,row-1),k))>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
        }
        
        for (int j1=0;j1<row; j1++){
          for (int j2=0;j2<col;j2++){
            if (((row_test1(j1)==1) && (col_test1(j2)==1)) ||  ((row_test1(j1)==1) && (sum(col_test1)==0)) || ((col_test1(j2)==1) && (sum(row_test1)==0)) || ((row_test2(j1)==1) && (col_test2(j2)==1)) ||  ((row_test2(j1)==1) && (sum(col_test2)==0)) || ((col_test2(j2)==1) && (sum(row_test2)==0))){
              T(i) = T(i)+1;
              tmp1_t(j1,j2) = tmp1(j1,j2)>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
              tmp2_t(j1,j2) = tmp2(j1,j2)>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
            }
            test((j-1)*row*col + j1*col + j2, 0) = tmp1_t(j1,j2);
            test((j-1)*row*col + j1*col + j2, 1) = tmp2_t(j1,j2);
          }
        }
      }
      
      if (remainder>0)
      {
        for (int k=size-remainder; k<size;k++)
        {
          test(k,0)=status(k,0)>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
          test(k,1)=status(k,1)>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
        }
      }
      
    }
    double pse_1c=0; double pse_2c=0; double psp_1c=0; double psp_2c=0; double testp1=0; double testp2=0; double testn1=0; double testn2=0; 
    
    for (int m=0;m<size;m++)
    {
      if ((status(m,0)==1) && (test(m,0)==1)) pse_1c=pse_1c+1;
      if ((status(m,1)==1) && (test(m,1)==1)) pse_2c=pse_2c+1;
      if ((status(m,0)==0) && (test(m,0)==0)) psp_1c=psp_1c+1;
      if ((status(m,1)==0) && (test(m,1)==0)) psp_2c=psp_2c+1;
      if (test(m,0)==1) testp1++;
      if (test(m,1)==1) testp2++;
      if (test(m,0)==0) testn1++;
      if (test(m,1)==0) testn2++;
    }
    pse1(i)=pse_1c/c_pos;
    pse2(i)=pse_2c/g_pos;
    psp1(i)=psp_1c/c_neg;
    psp2(i)=psp_2c/g_neg;
    ppv1(i)=pse_1c/testp1;
    ppv2(i)=pse_2c/testp2;
    npv1(i)=psp_1c/testn1;
    npv2(i)=psp_2c/testn2;
    
  }
  
  return Rcpp::List::create(Rcpp::Named("Number of tests")=T,
                            Rcpp::Named("PSE1")=pse1,
                            Rcpp::Named("PSE2")=pse2,
                            Rcpp::Named("PSP1")=psp1,
                            Rcpp::Named("PSP2")=psp2,
                            Rcpp::Named("PPV1")=ppv1,
                            Rcpp::Named("PPV2")=ppv2,
                            Rcpp::Named("NPV1")=npv1,
                            Rcpp::Named("NPV2")=npv2,
                            Rcpp::Named("STATUS")=status,
                            Rcpp::Named("test")=test);
}



//Simulation AT for three infections
// [[Rcpp::export]]
List AT_3_SIM(int size, vec p, vec SE, vec SP,int row, int col, int iter){
  
  int num_pools = floor(size/row/col);
  int remainder = size - num_pools*row*col;
  vec T(iter); T.fill(num_pools*(row+col)+remainder);
  vec pse1(iter,fill::zeros); vec ppv1(iter,fill::zeros);
  vec pse2(iter,fill::zeros); vec ppv2(iter,fill::zeros);
  vec pse3(iter,fill::zeros); vec ppv3(iter,fill::zeros);
  vec psp1(iter,fill::zeros); vec npv1(iter,fill::zeros);
  vec psp2(iter,fill::zeros); vec npv2(iter,fill::zeros);
  vec psp3(iter,fill::zeros); vec npv3(iter,fill::zeros);
  
  vec x(8);  x(0)=1; x(1)=2; x(2)=3; x(3)=4; x(4)=5; x(5)=6; x(6)=7; x(7)=8;
  mat tmp1(row,col,fill::zeros);
  mat tmp2(row,col,fill::zeros);
  mat tmp3(row,col,fill::zeros);
  mat tmp1_t(row,col,fill::zeros);
  mat tmp2_t(row,col,fill::zeros);
  mat tmp3_t(row,col,fill::zeros);
  vec res(size,fill::zeros);
  mat status(size,3,fill::zeros);
  mat test(size,3,fill::zeros);
  bool replace=true;
  
  for (int i=0; i<iter; i++){
    vec p_sim = p;
    res = Rcpp::RcppArmadillo::sample(x, size, replace, p_sim);
    status.fill(0);
    test.fill(0);
    
    for (int j=0; j<size;j++){
      if (res(j)==1){
        status(j,0)=0;
        status(j,1)=0;
        status(j,2)=0;
      }
      if (res(j)==2){
        status(j,0)=1;
        status(j,1)=0;
        status(j,2)=0;
      }
      if (res(j)==3){
        status(j,0)=0;
        status(j,1)=1;
        status(j,2)=0;
      }
      if (res(j)==4){
        status(j,0)=0;
        status(j,1)=0;
        status(j,2)=1;
      }
      if (res(j)==5){
        status(j,0)=1;
        status(j,1)=1;
        status(j,2)=0;
      }
      if (res(j)==6){
        status(j,0)=1;
        status(j,1)=0;
        status(j,2)=1;
      }
      if (res(j)==7){
        status(j,0)=0;
        status(j,1)=1;
        status(j,2)=1;
      }
      if (res(j)==8){
        status(j,0)=1;
        status(j,1)=1;
        status(j,2)=1;
      }
    } //end of data generation
    
    double c_pos = sum(status(span(0,size-1),0));
    double c_neg = size- c_pos;
    double g_pos = sum(status(span(0,size-1),1));
    double g_neg = size- g_pos;
    double h_pos = sum(status(span(0,size-1),2));
    double h_neg = size- h_pos;
    
    for (int j=1; j<=num_pools; j++){
      tmp1_t.fill(0);
      tmp2_t.fill(0);
      tmp3_t.fill(0);
      mat tmp(row*col,3,fill::zeros);
      tmp = status(span((j-1)*row*col, j*col*row-1),span(0,2));
      
      for (int j1=0;j1<row; j1++){
        for (int j2=0;j2<col;j2++){
          tmp1(j1,j2)=tmp(j1*col+j2,0);
          tmp2(j1,j2)=tmp(j1*col+j2,1);
          tmp3(j1,j2)=tmp(j1*col+j2,2);
        }
      }//end of generating one matrix
      
      vec row_test1(row,fill::zeros);
      vec row_test2(row,fill::zeros);
      vec row_test3(row,fill::zeros);
      for (int k=0;k<row;k++){
        row_test1(k)=sum(tmp1(k,span(0,col-1)))>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
        row_test2(k)=sum(tmp2(k,span(0,col-1)))>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
        row_test3(k)=sum(tmp3(k,span(0,col-1)))>0? (R::rbinom(1,SE(2))):(R::rbinom(1,1-SP(2)));
      }
      
      vec col_test1(col,fill::zeros);
      vec col_test2(col,fill::zeros);
      vec col_test3(col,fill::zeros);
      for (int k=0;k<col;k++){
        col_test1(k)=sum(tmp1(span(0,row-1),k))>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
        col_test2(k)=sum(tmp2(span(0,row-1),k))>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
        col_test3(k)=sum(tmp3(span(0,row-1),k))>0? (R::rbinom(1,SE(2))):(R::rbinom(1,1-SP(2)));
        
      }
      
      for (int j1=0;j1<row; j1++){
        for (int j2=0;j2<col;j2++){
          if (((row_test1(j1)==1) && (col_test1(j2)==1)) ||  ((row_test1(j1)==1) && (sum(col_test1)==0)) || ((col_test1(j2)==1) && (sum(row_test1)==0)) || ((row_test2(j1)==1) && (col_test2(j2)==1)) ||  ((row_test2(j1)==1) && (sum(col_test2)==0)) || ((col_test2(j2)==1) && (sum(row_test2)==0)) || ((row_test3(j1)==1) && (col_test3(j2)==1)) ||  ((row_test3(j1)==1) && (sum(col_test3)==0)) || ((col_test3(j2)==1) && (sum(row_test3)==0))){
            T(i) = T(i)+1;
            tmp1_t(j1,j2) = tmp1(j1,j2)>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
            tmp2_t(j1,j2) = tmp2(j1,j2)>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
            tmp3_t(j1,j2) = tmp3(j1,j2)>0? (R::rbinom(1,SE(2))):(R::rbinom(1,1-SP(2)));
          }
          test((j-1)*row*col + j1*col + j2, 0) = tmp1_t(j1,j2);
          test((j-1)*row*col + j1*col + j2, 1) = tmp2_t(j1,j2);
          test((j-1)*row*col + j1*col + j2, 2) = tmp3_t(j1,j2);
        }
      }
      
      if (remainder>0)
      {
        for (int k=size-remainder; k<size;k++)
        {
          test(k,0)=status(k,0)>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
          test(k,1)=status(k,1)>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
          test(k,2)=status(k,2)>0? (R::rbinom(1,SE(2))):(R::rbinom(1,1-SP(2)));
        }
      }
    } 
    
    double pse_1c=0; double pse_2c=0; double pse_3c=0; double psp_1c=0; double psp_2c=0; double psp_3c=0; 
    double testp1=0; double testp2=0; double testp3=0; double testn1=0; double testn2=0; double testn3=0;
    
    for (int m=0;m<size;m++)
    {
      if ((status(m,0)==1) && (test(m,0)==1)) pse_1c=pse_1c+1;
      if ((status(m,1)==1) && (test(m,1)==1)) pse_2c=pse_2c+1;
      if ((status(m,2)==1) && (test(m,2)==1)) pse_3c=pse_3c+1;
      if ((status(m,0)==0) && (test(m,0)==0)) psp_1c=psp_1c+1;
      if ((status(m,1)==0) && (test(m,1)==0)) psp_2c=psp_2c+1;
      if ((status(m,2)==0) && (test(m,2)==0)) psp_3c=psp_3c+1;
      if (test(m,0)==1) testp1++;
      if (test(m,1)==1) testp2++;
      if (test(m,2)==1) testp3++;
      if (test(m,0)==0) testn1++;
      if (test(m,1)==0) testn2++;
      if (test(m,2)==0) testn3++;
    }
    pse1(i)=pse_1c/c_pos;
    pse2(i)=pse_2c/g_pos;
    pse3(i)=pse_3c/h_pos;
    psp1(i)=psp_1c/c_neg;
    psp2(i)=psp_2c/g_neg;
    psp3(i)=psp_3c/h_neg;
    ppv1(i)=pse_1c/testp1;
    ppv2(i)=pse_2c/testp2;
    ppv3(i)=pse_3c/testp3;
    npv1(i)=psp_1c/testn1;
    npv2(i)=psp_2c/testn2;
    npv3(i)=psp_3c/testn3;
  }
  
  return Rcpp::List::create(Rcpp::Named("Number of tests")=T,
                            Rcpp::Named("PSE1")=pse1,
                            Rcpp::Named("PSE2")=pse2,
                            Rcpp::Named("PSE3")=pse3,
                            Rcpp::Named("PSP1")=psp1,
                            Rcpp::Named("PSP2")=psp2,
                            Rcpp::Named("PSP3")=psp3,
                            Rcpp::Named("PPV1")=ppv1,
                            Rcpp::Named("PPV2")=ppv2,
                            Rcpp::Named("PPV3")=ppv3,
                            Rcpp::Named("NPV1")=npv1,
                            Rcpp::Named("NPV2")=npv2,
                            Rcpp::Named("NPV3")=npv3,
                            Rcpp::Named("status")=status);
}





//Simulation ATM for three infections
// [[Rcpp::export]]
List ATM_3_SIM(int size, vec p, vec SE, vec SP,int row, int col, int iter){
  
  int num_pools = floor(size/row/col);
  int remainder = size - num_pools*row*col;
  vec T(iter); T.fill(num_pools+remainder);
  vec pse1(iter,fill::zeros); vec ppv1(iter,fill::zeros);
  vec pse2(iter,fill::zeros); vec ppv2(iter,fill::zeros);
  vec pse3(iter,fill::zeros); vec ppv3(iter,fill::zeros);
  vec psp1(iter,fill::zeros); vec npv1(iter,fill::zeros);
  vec psp2(iter,fill::zeros); vec npv2(iter,fill::zeros);
  vec psp3(iter,fill::zeros); vec npv3(iter,fill::zeros);
  vec x(8);  x(0)=1; x(1)=2; x(2)=3; x(3)=4; x(4)=5; x(5)=6; x(6)=7; x(7)=8;
  mat tmp1(row,col);
  mat tmp2(row,col);
  mat tmp3(row,col);
  mat tmp1_t(row,col);
  mat tmp2_t(row,col);
  mat tmp3_t(row,col);
  vec res(size);
  mat status(size,3,fill::zeros);
  mat test(size,3,fill::zeros);
  bool replace=true;
  
  for (int i=0; i<iter; i++){
    vec p_sim = p;
    res = Rcpp::RcppArmadillo::sample(x, size, replace, p_sim);
    status.fill(0);
    test.fill(0);
    
    for (int j=0; j<size;j++){
      if (res(j)==1){
        status(j,0)=0;
        status(j,1)=0;
        status(j,2)=0;
      }
      if (res(j)==2){
        status(j,0)=1;
        status(j,1)=0;
        status(j,2)=0;
      }
      if (res(j)==3){
        status(j,0)=0;
        status(j,1)=1;
        status(j,2)=0;
      }
      if (res(j)==4){
        status(j,0)=0;
        status(j,1)=0;
        status(j,2)=1;
      }
      if (res(j)==5){
        status(j,0)=1;
        status(j,1)=1;
        status(j,2)=0;
      }
      if (res(j)==6){
        status(j,0)=1;
        status(j,1)=0;
        status(j,2)=1;
      }
      if (res(j)==7){
        status(j,0)=0;
        status(j,1)=1;
        status(j,2)=1;
      }
      if (res(j)==8){
        status(j,0)=1;
        status(j,1)=1;
        status(j,2)=1;
      }
    } //end of data generation
    
    double c_pos = sum(status(span(0,size-1),0));
    double c_neg = size- c_pos;
    double g_pos = sum(status(span(0,size-1),1));
    double g_neg = size- g_pos;
    double h_pos = sum(status(span(0,size-1),2));
    double h_neg = size- h_pos;
    
    for (int j=1; j<=num_pools; j++){
      tmp1_t.fill(0);
      tmp2_t.fill(0);
      tmp3_t.fill(0);
      
      mat tmp(row*col,3,fill::zeros);
      tmp = status(span((j-1)*row*col, j*col*row-1),span(0,2));
      
      for (int j1=0;j1<row; j1++){
        for (int j2=0;j2<col;j2++){
          tmp1(j1,j2)=tmp(j1*col+j2,0);
          tmp2(j1,j2)=tmp(j1*col+j2,1);
          tmp3(j1,j2)=tmp(j1*col+j2,2);
        }
      }//end of generating one matrix
      
      int M1 = sum(tmp(span(0,col*row-1),0))>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
      int M2 = sum(tmp(span(0,col*row-1),1))>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
      int M3 = sum(tmp(span(0,col*row-1),2))>0? (R::rbinom(1,SE(2))):(R::rbinom(1,1-SP(2)));
      
      if (M1+M2+M3==0) {
        
        for (int j1=0;j1<row; j1++){
          for (int j2=0;j2<col;j2++){
            test((j-1)*row*col + j1*col + j2, 0) = 0;
            test((j-1)*row*col + j1*col + j2, 1) = 0;
            test((j-1)*row*col + j1*col + j2, 2) = 0;
          }
        }
      }
      
      else {
        T(i) = T(i)+col+row;
        
        vec row_test1(row,fill::zeros);
        vec row_test2(row,fill::zeros);
        vec row_test3(row,fill::zeros);
        
        for (int k=0;k<row;k++){
          row_test1(k)=sum(tmp1(k,span(0,col-1)))>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
          row_test2(k)=sum(tmp2(k,span(0,col-1)))>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
          row_test3(k)=sum(tmp3(k,span(0,col-1)))>0? (R::rbinom(1,SE(2))):(R::rbinom(1,1-SP(2)));
        }
        
        vec col_test1(col);
        vec col_test2(col);
        vec col_test3(col);
        
        for (int k=0;k<col;k++){
          col_test1(k)=sum(tmp1(span(0,row-1),k))>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
          col_test2(k)=sum(tmp2(span(0,row-1),k))>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
          col_test3(k)=sum(tmp3(span(0,row-1),k))>0? (R::rbinom(1,SE(2))):(R::rbinom(1,1-SP(2)));
        }
        
        for (int j1=0;j1<row; j1++){
          for (int j2=0;j2<col;j2++){
            if (((row_test1(j1)==1) && (col_test1(j2)==1)) ||  ((row_test1(j1)==1) && (sum(col_test1)==0)) || ((col_test1(j2)==1) && (sum(row_test1)==0)) || ((row_test2(j1)==1) && (col_test2(j2)==1)) ||  ((row_test2(j1)==1) && (sum(col_test2)==0)) || ((col_test2(j2)==1) && (sum(row_test2)==0)) || ((row_test3(j1)==1) && (col_test3(j2)==1)) ||  ((row_test3(j1)==1) && (sum(col_test3)==0)) || ((col_test3(j2)==1) && (sum(row_test3)==0))){
              T(i) = T(i)+1;
              tmp1_t(j1,j2) = tmp1(j1,j2)>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
              tmp2_t(j1,j2) = tmp2(j1,j2)>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
              tmp3_t(j1,j2) = tmp3(j1,j2)>0? (R::rbinom(1,SE(2))):(R::rbinom(1,1-SP(2)));
            }
            test((j-1)*row*col + j1*col + j2, 0) = tmp1_t(j1,j2);
            test((j-1)*row*col + j1*col + j2, 1) = tmp2_t(j1,j2);
            test((j-1)*row*col + j1*col + j2, 2) = tmp3_t(j1,j2);
          }
        }
      }
      
      if (remainder>0)
      {
        for (int k=size-remainder; k<size;k++)
        {
          test(k,0)=status(k,0)>0? (R::rbinom(1,SE(0))):(R::rbinom(1,1-SP(0)));
          test(k,1)=status(k,1)>0? (R::rbinom(1,SE(1))):(R::rbinom(1,1-SP(1)));
          test(k,2)=status(k,2)>0? (R::rbinom(1,SE(2))):(R::rbinom(1,1-SP(2)));
        }
      }
      
    }
    double pse_1c=0; double pse_2c=0; double pse_3c=0; double psp_1c=0; double psp_2c=0; double psp_3c=0; 
    double testp1=0; double testp2=0; double testp3=0; double testn1=0; double testn2=0; double testn3=0; 
    
    for (int m=0;m<size;m++)
    {
      if ((status(m,0)==1) && (test(m,0)==1)) pse_1c=pse_1c+1;
      if ((status(m,1)==1) && (test(m,1)==1)) pse_2c=pse_2c+1;
      if ((status(m,2)==1) && (test(m,2)==1)) pse_3c=pse_3c+1;
      if ((status(m,0)==0) && (test(m,0)==0)) psp_1c=psp_1c+1;
      if ((status(m,1)==0) && (test(m,1)==0)) psp_2c=psp_2c+1;
      if ((status(m,2)==0) && (test(m,2)==0)) psp_3c=psp_3c+1;
      if (test(m,0)==1) testp1++;
      if (test(m,1)==1) testp2++;
      if (test(m,2)==1) testp3++;
      if (test(m,0)==0) testn1++;
      if (test(m,1)==0) testn2++;
      if (test(m,2)==0) testn3++;
    }
    pse1(i)=pse_1c/c_pos;
    pse2(i)=pse_2c/g_pos;
    pse3(i)=pse_3c/h_pos;
    psp1(i)=psp_1c/c_neg;
    psp2(i)=psp_2c/g_neg;
    psp3(i)=psp_3c/h_neg;
    ppv1(i)=pse_1c/testp1;
    ppv2(i)=pse_2c/testp2;
    ppv3(i)=pse_3c/testp3;
    npv1(i)=psp_1c/testn1;
    npv2(i)=psp_2c/testn2;
    npv3(i)=psp_3c/testn3;
  }
  
  return Rcpp::List::create(Rcpp::Named("Number of tests")=T,
                            Rcpp::Named("PSE1")=pse1,
                            Rcpp::Named("PSE2")=pse2,
                            Rcpp::Named("PSE3")=pse3,
                            Rcpp::Named("PSP1")=psp1,
                            Rcpp::Named("PSP2")=psp2,
                            Rcpp::Named("PSP3")=psp3,
                            Rcpp::Named("PPV1")=ppv1,
                            Rcpp::Named("PPV2")=ppv2,
                            Rcpp::Named("PPV3")=ppv3,
                            Rcpp::Named("NPV1")=npv1,
                            Rcpp::Named("NPV2")=npv2,
                            Rcpp::Named("NPV3")=npv3,
                            Rcpp::Named("STATUS")=status);
}
