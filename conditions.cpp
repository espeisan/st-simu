#include <Fepic/CustomEigen>
#include <cmath>
#include <iostream>
using namespace Eigen;
using namespace std;

// space
typedef Matrix<double, Dynamic,1,0,6,1>              Vector;
typedef Matrix<double, Dynamic,Dynamic,RowMajor,3,3> Tensor;
typedef Matrix<double, 3, 3> Tensor3;
typedef Matrix<double, Dynamic,Dynamic,RowMajor,6,6> TensorZ;
const double pi  = 3.141592653589793;
const double pi2 = pi*pi;

double pho(Vector const& X, int tag);
double gama(Vector const& X, double t, int tag);
double muu(int tag);
Vector force(Vector const& X, double t, int tag);   //density*gravity (force/vol)
Vector u_exact(Vector const& X, double t, int tag);
double p_exact(Vector const& X, double t, int tag);
Vector z_exact(Vector const& X, double t, int tag);
Tensor grad_u_exact(Vector const& X, double t, int tag);
Vector grad_p_exact(Vector const& X, double t, int tag);
Vector traction(Vector const& X, Vector const& normal, double t, int tag);
Vector u_initial(Vector const& X, int tag);
double p_initial(Vector const& X, int tag);
Vector z_initial(Vector const& X, int tag);
Vector solid_normal(Vector const& X, double t, int tag);
Tensor feature_proj(Vector const& X, double t, int tag);
Vector gravity(Vector const& X, int dim);
Vector force_pp(Vector const& Xi, Vector const& Xj, double Ri, double Rj,
                 double ep1, double ep2, double zeta);
Vector force_pw(Vector const& Xi, Vector const& Xj, double Ri,
                 double ew1, double ew2, double zeta);
Vector force_ppl(Vector const& Xi, Vector const& Xj, double ep, double zeta);
Vector force_rga(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 Vector const& Gr, double const masj, double ep, double zeta);
Vector force_rgb(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 Vector const& Gr, double const rhoj, double const rhof, double ep, double zeta);
Vector force_rgc(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 double ep, double zeta);
TensorZ MI_tensor(double M, double R, int dim, Tensor3 TI);
Tensor RotM(double theta, int dim);
Vector SlipVel(Vector const& X, Vector const& XG, Vector const& normal, int dim, int tag, double theta);

double Dif_coeff(int tag);
double nuB_coeff(int tag);
double sig_coeff(int tag);
double bbG_coeff(Vector const& X, int tag);

Vector traction_maxwell(Vector const& E, Vector const& normal, double eps, int tag);
double per_Elect(int tag);

Vector exact_tangent_ellipse(Vector const& X, Vector const& Xc, double theta, double R1, double R2, int dim);
Vector exact_normal_ellipse(Vector const& X, Vector const& Xc, double theta, double R1, double R2, int dim);

Vector cubic_ellipse(double yb, Vector const& X0, Vector const& X2, Vector const& T0, Vector const& T2, int dim);
Vector Dcubic_ellipse(double yb, Vector const& X0, Vector const& X2, Vector const& T0, Vector const& T2, int dim);
Vector curved_Phi(double yb, Vector const& X0, Vector const& X2, Vector const& T0, Vector const& T2, int dim);
Vector Dcurved_Phi(double yb, Vector const& X0, Vector const& X2, Vector const& T0, Vector const& T2, int dim);

// gota estática 2d/////////////////////////////////////////////////////////////
#if (false)

double pho(Vector const& X, int tag)
{
  if (tag == 102)
  {
    return 6.4e1; //1.0e3;
  }
  else
  {
    return 1.0; //8.0e2;
  }
}

double cos_theta0()
{
  return 0.0;
}

double zeta(double u_norm, double angle)
{
  return 0.0;
}

double beta_diss()
{
  return 0.0;
}

double gama(Vector const& X, double t, int tag)
{
  return 1.0;
}

double muu(int tag)
{
  if (tag == 102)
  {
    return 1.0e3;
  }
  else
  {
    return 1.0;
  }
}

Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(X.size()));
  //f(1) = 10;
  return f;
}

Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double un = 998*9.8*2.85e-4*2.85e-4/(3*1e-3);
  double w1 = 0.0, w2 = 1.0, R1 = .25, R2 = 1;
  //Vector v(Vector::Ones(X.size()));  v << 1 , 2;
  Vector v(Vector::Zero(X.size()));
  double r  = sqrt(x*x+y*y);
  double Fr = (R2*w2-R1*w1)*(r-R1)/(R2-R1) + R1*w1;
  double Lr = r*w2;
  double nu = muu(tag)/pho(X,tag);
  double C  = 78.0;
  double uc = exp(-C*nu*t)*(Fr-Lr)+Lr;
  v(0) = -y*Lr/r;
  v(1) =  x*Lr/r;
  if ( t == 0 ){
    if (tag == 1){
      v(0) = -w2*y;
      v(1) =  w2*x;
    }
    else if (tag == 1 && t >= 3 && false){
      v(0) =  w2*y;
      v(1) = -w2*x;
    }
    else{
      v(0) = 0;
      v(1) = 0;
    }
  }
  return v;
}

Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double w2 = 2.0;
  Tensor dxU(Tensor::Zero(X.size(), X.size()));
  dxU(0,0) = 0; dxU(0,1) = -w2;
  dxU(1,0) = w2; dxU(1,1) = 0;

  return dxU;
}

double p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  return 0.0;
}

Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector dxP(X.size());

  return dxP;
}

Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  Vector T(Vector::Zero(X.size()));
  //T(0) = -p_exact(X,t,tag);
  //T(1) = muu(tag)*(cos(w_*t) + sin(w_*t));
  Tensor dxU(grad_u_exact(X,t,tag));
  Tensor I(Tensor::Identity(2,2));
  T = (- p_exact(X,t,tag)*I +  muu(tag)*(dxU + dxU.transpose()))*normal;
  return T;
}

Vector u_initial(Vector const& X, int tag)
{
  return u_exact(X,0,tag);
}

double p_initial(Vector const& X, int tag)
{
  return p_exact(X,0,tag);
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  return N;
}

Vector v_exact(Vector const& X, double t, int tag) //(X,t,tag)
{
  double const tet = 10. * (pi/180.);
  double const x = X(0);
  double const y = X(1);
  Vector v(Vector::Zero(X.size()));

  return v;
}

// posição do contorno
Vector x_exact(Vector const& X, double t, int tag)
{
  Vector r(Vector::Zero(X.size()));
  return r;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  //if (tag == 2) {N(1) = 1;};
  return N;
}

Tensor feature_proj(Vector const& X, double t, int tag)
{
  Tensor f(Tensor::Zero(X.size(), X.size()));
  return f;
}

Vector z_exact(Vector const& X, double t, int tag)
{
  double w2 = 2.0;
  int dim = X.size();
  int LZ = 3*(dim-1);
  Vector v(Vector::Zero(LZ)); //v << 0, 0, 10;
  if (t > 0){
    v(2) = w2;
  }
  //Vector v(Vector::Ones(LZ));
  return v;
}

Vector gravity(Vector const& X, int dim){
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(3*(dim-1)));
  return f;
}

Vector z_initial(Vector const& X, int tag)
{
  return z_exact(X,0,tag);
}

#endif


// solid asc 2d/////////////////////////////////////////////////////////////
#if (false)

double pho(Vector const& X, int tag)
{
//  if (tag == 15)
//  {
    return 1.0;//e3;///1e4;
//  }
//  else
//  {
//    return 0.0;
//  }
}

double cos_theta0()
{
  return 0.0;
}

double zeta(double u_norm, double angle)
{
  return 0.0;
}

double beta_diss()
{
  return 0.0;
}

double gama(Vector const& X, double t, int tag)
{
  return 0.5;
}

double muu(int tag)
{
//  if (tag == 15)
//  {
    return 1e-1;//1.0*0.1;
//  }
//  else
//  {
//    return 0.0;
//  }
}

Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(X.size()));
//  if (tag == 15)
//  {
    f(1) = -980.0*1.0;//pho(X,tag);//*1e4;//*1e3;
//
//  else
//  {
//    f(1) = 0.0;  //-8e-4*1e4;
//  }
  return f;
}

Vector gravity(Vector const& X, int dim){
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(3*(dim-1)));
  if (dim == 2){
    f(1) = -980.0;//-8e-4;  //*1e3;
  }
  else if (dim == 3){
    f(2) = -980.0;  //-8e-4*1e4;
  }
  return f;
}


Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double un = 998*9.8*2.85e-4*2.85e-4/(3*1e-3);
  //Vector v(Vector::Ones(X.size()));  v << 1 , 2;
  Vector v(Vector::Zero(X.size()));

  return v;
}

Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Tensor dxU(Tensor::Zero(X.size(), X.size()));

  return dxU;
}

Vector z_exact(Vector const& X, double t, int tag)
{
  int dim = X.size();
  int LZ = 3*(dim-1);
  Vector v(Vector::Zero(LZ)); //v << .1, .2, .3;
  //Vector v(Vector::Ones(LZ));
  return v;
}

double p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  
  return 0.0;
}

Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector dxP(X.size());

  return dxP;
}

Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  Vector T(Vector::Zero(X.size()));
  //T(0) = -p_exact(X,t,tag);
  //T(1) = muu(tag)*(cos(w_*t) + sin(w_*t));
  Tensor dxU(grad_u_exact(X,t,tag));
  Tensor I(Tensor::Identity(2,2));
  T = (- p_exact(X,t,tag)*I +  muu(tag)*(dxU + dxU.transpose()))*normal;
  return T;
}

Vector u_initial(Vector const& X, int tag)
{
  return u_exact(X,0,tag);
}

Vector z_initial(Vector const& X, int tag)
{
  return z_exact(X,0,tag);
}

double p_initial(Vector const& X, int tag)
{
  return p_exact(X,0,tag);
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  return N;
}

Vector v_exact(Vector const& X, double t, int tag) //(X,t,tag)
{
  double const x = X(0);
  double const y = X(1);
  Vector v(Vector::Zero(X.size()));

  return v;
}

// posição do contorno
Vector x_exact(Vector const& X, double t, int tag)
{
  Vector r(Vector::Zero(X.size()));
  return r;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  //if (tag == 2) {N(1) = 1;};
  return N;
}

Tensor feature_proj(Vector const& X, double t, int tag)
{
  Tensor f(Tensor::Zero(X.size(), X.size()));
  return f;
}

Vector force_pp(Vector const& Xi, Vector const& Xj, double Ri, double Rj,
                 double ep1, double ep2, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();
//  if (dij > Ri+Rj+zeta){
//    return f;
//  }
  if (dij <= Ri+Rj){
    f = (1/ep1)*(Ri+Rj-dij)*(Xi - Xj);
  }
  else if((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    f = (1/ep2)*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj);
  }
  return f;
}

Vector force_pw(Vector const& Xi, Vector const& Xj, double Ri,
                 double ew1, double ew2, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double di = (Xi - Xj).norm();
//  if (dij > Ri+Rj+zeta){
//    return f;
//  }
  if (di <= 2*Ri){
    f = (1/ew1)*(2*Ri-di)*(Xi - Xj);
  }
  else if((2*Ri <= di) && (di <= 2*Ri+zeta)){
    f = (1/ew2)*(2*Ri+zeta-di)*(2*Ri+zeta-di)*(Xi - Xj);
  }
  return f;
}

Vector force_ppl(Vector const& Xi, Vector const& Xj, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  f = (zeta/ep)*(Xi - Xj)/(Xi - Xj).norm();
  return f;
}

Vector force_rga(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 Vector const& Gr, double const masj, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();;
  double g = 0.0;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    g   = Gr.norm();
    f   = masj*g*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/(ep*zeta*zeta*dij);
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}

Vector force_rgb(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 Vector const& Gr, double const rhoj, double const rhof, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = dij = (Xi - Xj).norm();
  double g = 0.0;
  double R = 0.0;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    g   = Gr.norm();
    R   = std::max(Ri,Rj);
    f   = (rhoj-rhof)*pi*g*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/(ep*zeta*zeta*dij);
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}

Vector force_rgc(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    f   = (Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/ep;
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}
#endif

// canal /////////////////////////////////////////////////////////////
#if (false)

double pho(Vector const& X, int tag)
{
  return 1.0;
}

double cos_theta0()
{
  return 0.0;
}

double zeta(double u_norm, double angle)
{
  return 0*5.e-1;
}

double beta_diss()
{
  return 0*1.e-4;
}

double gama(Vector const& X, double t, int tag)
{
  return 0;
}

double muu(int tag)
{
  return  1.0;
}

Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(X.size()));
  return f;
}

Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double umax = 1;
  Vector v(Vector::Zero(X.size()));

  if (tag == 3 || tag == 2)
  //if (tag == 2)
  {
    v(0) = umax*y*(2-y);
    v(1) = 0;
  }

  return v;
}

Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Tensor dxU(Tensor::Zero(X.size(), X.size()));

  return dxU;
}

double p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  return 0;
}

Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector dxP(X.size());

  return dxP;
}

Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  Vector T(Vector::Zero(X.size()));
  //T(0) = -p_exact(X,t,tag);
  //T(1) = muu(tag)*(cos(w_*t) + sin(w_*t));
  Tensor dxU(grad_u_exact(X,t,tag));
  Tensor I(Tensor::Identity(2,2));
  T = (- p_exact(X,t,tag)*I +  muu(tag)*(dxU + dxU.transpose()))*normal;
  return T;
}

Vector u_initial(Vector const& X, int tag)
{
  return u_exact(X,0,tag);
}

double p_initial(Vector const& X, int tag)
{
  return p_exact(X,0,tag);
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  return N;
}

Vector v_exact(Vector const& X, double t, int tag) //(X,t,tag)
{
  double const tet = 10. * (pi/180.);
  double const x = X(0);
  double const y = X(1);
  Vector v(Vector::Zero(X.size()));

  return v;
}

// posição do contorno
Vector x_exact(Vector const& X, double t, int tag)
{
  Vector r(Vector::Zero(X.size()));
  return r;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  //N(1) = 1;
  return N;
}

Tensor feature_proj(Vector const& X, double t, int tag)
{
  Tensor f(Tensor::Zero(X.size(), X.size()));
  if (tag == 3)
  {
    f(0,0) = 1;
  }
//  else if (tag == 1)
//  {
//    f(0,0) = 1;
//  }
  return f;
}

#endif

// rot solid 2d/////////////////////////////////////////////////////////////
#if (false)

double pho(Vector const& X, int tag)
{
//  if (tag == 15)
//  {
    return 1.0;//e3;///1e4;
//  }
//  else
//  {
//    return 0.0;
//  }
}

double cos_theta0()
{
  return 0.0;
}

double zeta(double u_norm, double angle)
{
  return 0.0;
}

double beta_diss()
{
  return 0.0;
}

double gama(Vector const& X, double t, int tag)
{
  return 0.5;
}

double muu(int tag)
{
//  if (tag == 15)
//  {
    return 1e-2;//1.0*0.1;
//  }
//  else
//  {
//    return 0.0;
//  }
}

Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(X.size()));
//  if (tag == 15)
//  {
//    f(1) = -980.0*1.0;//pho(X,tag);//*1e4;//*1e3;
//
//  else
//  {
//    f(1) = 0.0;  //-8e-4*1e4;
//  }
  return f;
}

Vector gravity(Vector const& X, int dim){
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(3*(dim-1)));
  return f;
}

Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double un = 998*9.8*2.85e-4*2.85e-4/(3*1e-3);
  double w1 = 0.0, w2 = 2.0, R1 = .125, R2 = .5;
  //Vector v(Vector::Ones(X.size()));  v << 1 , 2;
  Vector v(Vector::Zero(X.size()));
  double r  = sqrt(x*x+y*y);
  double Fr = (R2*w2-R1*w1)*(r-R1)/(R2-R1) + R1*w1;
  double Lr = r*w2;
  double nu = muu(tag)/pho(X,tag);
  double C  = 78.0;
  double uc = exp(-C*nu*t)*(Fr-Lr)+Lr;
  v(0) = -y*Lr/r;
  v(1) =  x*Lr/r;
  if ( t == 0 ){
    if (tag == 1){
      v(0) = -w2*y;
      v(1) =  w2*x;
    }
    else if (tag == 1 && t >= 3 && false){
      v(0) =  w2*y;
      v(1) = -w2*x;
    }
    else{
      v(0) = 0;
      v(1) = 0;
    }
  }
  return v;
}

Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double w2 = 2.0;
  Tensor dxU(Tensor::Zero(X.size(), X.size()));
  dxU(0,0) = 0; dxU(0,1) = -w2;
  dxU(1,0) = w2; dxU(1,1) = 0;

  return dxU;
}

Vector z_exact(Vector const& X, double t, int tag)
{
  double w2 = 2.0;
  int dim = X.size();
  int LZ = 3*(dim-1);
  Vector v(Vector::Zero(LZ)); //v << 0, 0, 10;
  if (t > 0){
    v(2) = w2;
  }
  //Vector v(Vector::Ones(LZ));
  return v;
}

double p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  return 0.0;
}

Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector dxP(X.size());

  return dxP;
}

Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  Vector T(Vector::Zero(X.size()));
  //T(0) = -p_exact(X,t,tag);
  //T(1) = muu(tag)*(cos(w_*t) + sin(w_*t));
  Tensor dxU(grad_u_exact(X,t,tag));
  Tensor I(Tensor::Identity(2,2));
  T = (- p_exact(X,t,tag)*I +  muu(tag)*(dxU + dxU.transpose()))*normal;
  return T;
}

Vector u_initial(Vector const& X, int tag)
{
  return u_exact(X,0,tag);
}

Vector z_initial(Vector const& X, int tag)
{
  return z_exact(X,0,tag);
}

double p_initial(Vector const& X, int tag)
{
  return p_exact(X,0,tag);
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  return N;
}

Vector v_exact(Vector const& X, double t, int tag) //(X,t,tag)
{
  double const x = X(0);
  double const y = X(1);
  Vector v(Vector::Zero(X.size()));

  return v;
}

// posição do contorno
Vector x_exact(Vector const& X, double t, int tag)
{
  Vector r(Vector::Zero(X.size()));
  return r;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  //if (tag == 2) {N(1) = 1;};
  return N;
}

Tensor feature_proj(Vector const& X, double t, int tag)
{
  Tensor f(Tensor::Zero(X.size(), X.size()));
  return f;
}

Vector force_pp(Vector const& Xi, Vector const& Xj, double Ri, double Rj,
                 double ep1, double ep2, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();
//  if (dij > Ri+Rj+zeta){
//    return f;
//  }
  if (dij <= Ri+Rj){
    f = (1/ep1)*(Ri+Rj-dij)*(Xi - Xj);
  }
  else if((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    f = (1/ep2)*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj);
  }
  return f;
}

Vector force_pw(Vector const& Xi, Vector const& Xj, double Ri,
                 double ew1, double ew2, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double di = (Xi - Xj).norm();
//  if (dij > Ri+Rj+zeta){
//    return f;
//  }
  if (di <= 2*Ri){
    f = (1/ew1)*(2*Ri-di)*(Xi - Xj);
  }
  else if((2*Ri <= di) && (di <= 2*Ri+zeta)){
    f = (1/ew2)*(2*Ri+zeta-di)*(2*Ri+zeta-di)*(Xi - Xj);
  }
  return f;
}

Vector force_ppl(Vector const& Xi, Vector const& Xj, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  f = (zeta/ep)*(Xi - Xj)/(Xi - Xj).norm();
  return f;
}

Vector force_rga(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 Vector const& Gr, double const masj, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();;
  double g = 0.0;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    g   = Gr.norm();
    f   = masj*g*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/(ep*zeta*zeta*dij);
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}

Vector force_rgb(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 Vector const& Gr, double const rhoj, double const rhof, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = dij = (Xi - Xj).norm();
  double g = 0.0;
  double R = 0.0;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    g   = Gr.norm();
    R   = std::max(Ri,Rj);
    f   = (rhoj-rhof)*pi*g*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/(ep*zeta*zeta*dij);
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}

Vector force_rgc(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    f   = (Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/ep;
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}
/*
TensorZ MI_tensor(double M, double R, int dim)
{
  TensorZ MI(TensorZ::Zero(3*(dim-1),3*(dim-1)));
  if (dim == 2){
    MI(0,0) = M; MI(1,1) = M; MI(2,2) = 0.5*M*R*R;
  }
  else if(dim == 3){
    MI(0,0) = M; MI(1,1) = M; MI(2,2) = M;
    MI(3,3) = 0.4*M*R*R; MI(4,4) = 0.4*M*R*R; MI(5,5) = 0.4*M*R*R;
  }
  return MI;
}

Vector SlipVel(Vector const& X, Vector const& XG, int dim, int tag)
{
  Vector V(Vector::Zero(dim));
  Vector X3(Vector::Zero(3));
  Vector Xp(Vector::Zero(dim));
  double alp = 0.1;
  double bet = 0.1;

  if (dim == 2)
  {
    X3(0) = X(0); X3(1) = X(1);
    X3 = X3 - XG;
    X3.normalize();
    Xp(0) = X3(1); Xp(1) = -X3(0);
    if (tag == 104){
      V = alp*Xp;
    }
    else if (tag == 105){
      V = -bet*Xp;
    }
  }

  return V;
}
*/
#endif

// rot solid 2d slip vel/////////////////////////////////////////////////////////////
#if (false)

double pho(Vector const& X, int tag)
{
//  if (tag == 15)
//  {
    return 1.0;//e3;///1e4;
//  }
//  else
//  {
//    return 0.0;
//  }
}

double cos_theta0()
{
  return 0.0;
}

double zeta(double u_norm, double angle)
{
  return 0.0;
}

double beta_diss()
{
  return 0.0;
}

double gama(Vector const& X, double t, int tag)
{
  return 0.5;
}

double muu(int tag)
{
//  if (tag == 15)
//  {
    return 1;//1.0*0.1;
//  }
//  else
//  {
//    return 0.0;
//  }
}

Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(X.size()));
//  if (tag == 15)
//  {
//    f(1) = -980.0*1.0;//pho(X,tag);//*1e4;//*1e3;
//
//  else
//  {
//    f(1) = 0.0;  //-8e-4*1e4;
//  }
  return f;
}

Vector gravity(Vector const& X, int dim){
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(3*(dim-1)));
  return f;
}

Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double un = 998*9.8*2.85e-4*2.85e-4/(3*1e-3);
  double w1 = 0.0, w2 = 1.0, R1 = .25, R2 = 1;
  //Vector v(Vector::Ones(X.size()));  v << 1 , 2;
  Vector v(Vector::Zero(X.size()));
/*  double r  = sqrt(x*x+y*y);
  double Fr = (R2*w2-R1*w1)*(r-R1)/(R2-R1) + R1*w1;
  double Lr = r*w2;
  double nu = muu(tag)/pho(X,tag);
  double C  = 78.0;
  double uc = exp(-C*nu*t)*(Fr-Lr)+Lr;
  v(0) = -y*Lr/r;
  v(1) =  x*Lr/r;
  if ( t == 0 ){
    if (tag == 1){
      v(0) = -w2*y;
      v(1) =  w2*x;
    }
    else if (tag == 1 && t >= 3 && false){
      v(0) =  w2*y;
      v(1) = -w2*x;
    }
    else{
      v(0) = 0;
      v(1) = 0;
    }
  }*/
  return v;
}

Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double w2 = 2.0;
  Tensor dxU(Tensor::Zero(X.size(), X.size()));
  dxU(0,0) = 0; dxU(0,1) = -w2;
  dxU(1,0) = w2; dxU(1,1) = 0;

  return dxU;
}

Vector z_exact(Vector const& X, double t, int tag)
{
  double w2 = 2.0;
  int dim = X.size();
  int LZ = 3*(dim-1);
  Vector v(Vector::Zero(LZ)); //v << 0, 0, 10;
  if (t > 0){
    v(2) = w2;
  }
  //Vector v(Vector::Ones(LZ));
  return v;
}

double p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  return 0.0;
}

Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector dxP(X.size());

  return dxP;
}

Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  Vector T(Vector::Zero(X.size()));
  //T(0) = -p_exact(X,t,tag);
  //T(1) = muu(tag)*(cos(w_*t) + sin(w_*t));
  Tensor dxU(grad_u_exact(X,t,tag));
  Tensor I(Tensor::Identity(2,2));
  T = (- p_exact(X,t,tag)*I +  muu(tag)*(dxU + dxU.transpose()))*normal;
  return T;
}

Vector u_initial(Vector const& X, int tag)
{
  return u_exact(X,0,tag);
}

Vector z_initial(Vector const& X, int tag)
{
  return z_exact(X,0,tag);
}

double p_initial(Vector const& X, int tag)
{
  return p_exact(X,0,tag);
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  return N;
}

Vector v_exact(Vector const& X, double t, int tag) //(X,t,tag)
{
  double const x = X(0);
  double const y = X(1);
  Vector v(Vector::Zero(X.size()));

  return v;
}

// posição do contorno
Vector x_exact(Vector const& X, double t, int tag)
{
  Vector r(Vector::Zero(X.size()));
  return r;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  //if (tag == 2) {N(1) = 1;};
  return N;
}

Tensor feature_proj(Vector const& X, double t, int tag)
{
  Tensor f(Tensor::Zero(X.size(), X.size()));
  return f;
}

Vector force_pp(Vector const& Xi, Vector const& Xj, double Ri, double Rj,
                 double ep1, double ep2, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();
//  if (dij > Ri+Rj+zeta){
//    return f;
//  }
  if (dij <= Ri+Rj){
    f = (1/ep1)*(Ri+Rj-dij)*(Xi - Xj);
  }
  else if((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    f = (1/ep2)*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj);
  }
  return f;
}

Vector force_pw(Vector const& Xi, Vector const& Xj, double Ri,
                 double ew1, double ew2, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double di = (Xi - Xj).norm();
//  if (dij > Ri+Rj+zeta){
//    return f;
//  }
  if (di <= 2*Ri){
    f = (1/ew1)*(2*Ri-di)*(Xi - Xj);
  }
  else if((2*Ri <= di) && (di <= 2*Ri+zeta)){
    f = (1/ew2)*(2*Ri+zeta-di)*(2*Ri+zeta-di)*(Xi - Xj);
  }
  return f;
}

Vector force_ppl(Vector const& Xi, Vector const& Xj, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  f = (zeta/ep)*(Xi - Xj)/(Xi - Xj).norm();
  return f;
}

Vector force_rga(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 Vector const& Gr, double const masj, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();;
  double g = 0.0;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    g   = Gr.norm();
    f   = masj*g*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/(ep*zeta*zeta*dij);
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}

Vector force_rgb(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 Vector const& Gr, double const rhoj, double const rhof, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = dij = (Xi - Xj).norm();
  double g = 0.0;
  double R = 0.0;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    g   = Gr.norm();
    R   = std::max(Ri,Rj);
    f   = (rhoj-rhof)*pi*g*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/(ep*zeta*zeta*dij);
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}

Vector force_rgc(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    f   = (Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/ep;
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}

#endif

// rot solid 3d slip vel/////////////////////////////////////////////////////////////
#if (false)

double pho(Vector const& X, int tag)
{
//  if (tag == 15)
//  {
    return 1.0;//e3;///1e4;
//  }
//  else
//  {
//    return 0.0;
//  }
}

double cos_theta0()
{
  return 0.0;
}

double zeta(double u_norm, double angle)
{
  return 0.0;
}

double beta_diss()
{
  return 0.0;
}

double gama(Vector const& X, double t, int tag)
{
  return 0.5;
}

double muu(int tag)
{
//  if (tag == 15)
//  {
    return 0.1;//1.0*0.1;
//  }
//  else
//  {
//    return 0.0;
//  }
}

Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(X.size()));
//  if (tag == 15)
//  {
//    f(1) = -980.0*1.0;//pho(X,tag);//*1e4;//*1e3;
//
//  else
//  {
//    f(1) = 0.0;  //-8e-4*1e4;
//  }
  return f;
}

Vector gravity(Vector const& X, int dim){
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(3*(dim-1)));
  return f;
}

Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double un = 998*9.8*2.85e-4*2.85e-4/(3*1e-3);
  double w1 = 0.0, w2 = 1.0, R1 = .25, R2 = 1;
  //Vector v(Vector::Ones(X.size()));  v << 1 , 2;
  Vector v(Vector::Zero(X.size()));
  if (tag == 101 || tag == 102 || tag == 103){
    v(0) = 0.333328;
  }
/*  double r  = sqrt(x*x+y*y);
  double Fr = (R2*w2-R1*w1)*(r-R1)/(R2-R1) + R1*w1;
  double Lr = r*w2;
  double nu = muu(tag)/pho(X,tag);
  double C  = 78.0;
  double uc = exp(-C*nu*t)*(Fr-Lr)+Lr;
  v(0) = -y*Lr/r;
  v(1) =  x*Lr/r;
  if ( t == 0 ){
    if (tag == 1){
      v(0) = -w2*y;
      v(1) =  w2*x;
    }
    else if (tag == 1 && t >= 3 && false){
      v(0) =  w2*y;
      v(1) = -w2*x;
    }
    else{
      v(0) = 0;
      v(1) = 0;
    }
  }*/
  return v;
}

Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double w2 = 2.0;
  Tensor dxU(Tensor::Zero(X.size(), X.size()));
  dxU(0,0) = 0; dxU(0,1) = -w2;
  dxU(1,0) = w2; dxU(1,1) = 0;

  return dxU;
}

Vector z_exact(Vector const& X, double t, int tag)
{
  double w2 = 2.0;
  int dim = X.size();
  int LZ = 3*(dim-1);
  Vector v(Vector::Zero(LZ)); v << 0.333328, 0, 0, 0, 0, 0;
  //Vector v(Vector::Ones(LZ));
  return v;
}

double p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  return 0.0;
}

Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector dxP(X.size());

  return dxP;
}

Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  Vector T(Vector::Zero(X.size()));
  //T(0) = -p_exact(X,t,tag);
  //T(1) = muu(tag)*(cos(w_*t) + sin(w_*t));
  Tensor dxU(grad_u_exact(X,t,tag));
  Tensor I(Tensor::Identity(2,2));
  T = (- p_exact(X,t,tag)*I +  muu(tag)*(dxU + dxU.transpose()))*normal;
  return T;
}

Vector u_initial(Vector const& X, int tag)
{
  return u_exact(X,0,tag);
}

Vector z_initial(Vector const& X, int tag)
{
  return z_exact(X,0,tag);
}

double p_initial(Vector const& X, int tag)
{
  return p_exact(X,0,tag);
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  return N;
}

Vector v_exact(Vector const& X, double t, int tag) //(X,t,tag)
{
  double const x = X(0);
  double const y = X(1);
  Vector v(Vector::Zero(X.size()));

  return v;
}

// posição do contorno
Vector x_exact(Vector const& X, double t, int tag)
{
  Vector r(Vector::Zero(X.size()));
  return r;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  //if (tag == 2) {N(1) = 1;};
  return N;
}

Tensor feature_proj(Vector const& X, double t, int tag)
{
  Tensor f(Tensor::Zero(X.size(), X.size()));
  return f;
}

Vector force_pp(Vector const& Xi, Vector const& Xj, double Ri, double Rj,
                 double ep1, double ep2, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();
//  if (dij > Ri+Rj+zeta){
//    return f;
//  }
  if (dij <= Ri+Rj){
    f = (1/ep1)*(Ri+Rj-dij)*(Xi - Xj);
  }
  else if((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    f = (1/ep2)*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj);
  }
  return f;
}

Vector force_pw(Vector const& Xi, Vector const& Xj, double Ri,
                 double ew1, double ew2, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double di = (Xi - Xj).norm();
//  if (dij > Ri+Rj+zeta){
//    return f;
//  }
  if (di <= 2*Ri){
    f = (1/ew1)*(2*Ri-di)*(Xi - Xj);
  }
  else if((2*Ri <= di) && (di <= 2*Ri+zeta)){
    f = (1/ew2)*(2*Ri+zeta-di)*(2*Ri+zeta-di)*(Xi - Xj);
  }
  return f;
}

Vector force_ppl(Vector const& Xi, Vector const& Xj, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  f = (zeta/ep)*(Xi - Xj)/(Xi - Xj).norm();
  return f;
}

Vector force_rga(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 Vector const& Gr, double const masj, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();;
  double g = 0.0;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    g   = Gr.norm();
    f   = masj*g*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/(ep*zeta*zeta*dij);
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}

Vector force_rgb(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 Vector const& Gr, double const rhoj, double const rhof, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = dij = (Xi - Xj).norm();
  double g = 0.0;
  double R = 0.0;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    g   = Gr.norm();
    R   = std::max(Ri,Rj);
    f   = (rhoj-rhof)*pi*g*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/(ep*zeta*zeta*dij);
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}

Vector force_rgc(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    f   = (Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/ep;
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}

#endif

// rot solid 2d trap and channel/////////////////////////////////////////////////////////////
#if (false)

double pho(Vector const& X, int tag)
{
//  if (tag == 15)
//  {
    return 1000;//e3;///1e4;
//  }
//  else
//  {
//    return 0.0;
//  }
}

double cos_theta0()
{
  return 0.0;
}

double zeta(double u_norm, double angle)
{
  return 0.0;
}

double beta_diss()
{
  return 0.0;
}

double gama(Vector const& X, double t, int tag)
{
  return 0.5;
}

double muu(int tag)
{
//  if (tag == 15)
//  {
    return 1e-3;//1.0*0.1;
//  }
//  else
//  {
//    return 0.0;
//  }
}

Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(X.size()));
//  if (tag == 15)
//  {
//    f(1) = -980.0*1.0;//pho(X,tag);//*1e4;//*1e3;
//
//  else
//  {
//    f(1) = 0.0;  //-8e-4*1e4;
//  }
  return f;
}

Vector gravity(Vector const& X, int dim){
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(3*(dim-1)));
  return f;
}

Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(Vector::Zero(X.size()));
  double Um = 10e-6, H = 200e-6;
  if (tag == 4){
    //Um = Um*(1-exp(-t*t*t*t*1e10));//Um*(1-exp(t*t*t*t/1e-14));//Um*(-1-exp(t*t*t*t/1e-14)*0);
    v(0) = Um*4*(H-y)*y/(H*H);
  }
  return v;
}

Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double w2 = 2.0, Um = 0.1, H = 2e-6;
  Tensor dxU(Tensor::Zero(X.size(), X.size()));
  dxU(0,0) = 0; dxU(0,1) = Um/H - Um*2.0*y/H;
  dxU(1,0) = 0; dxU(1,1) = 0;

  return dxU;
}

Vector z_exact(Vector const& X, double t, int tag)
{
  double w2 = 0.0;
  int dim = X.size();
  int LZ = 3*(dim-1);
  Vector v(Vector::Zero(LZ)); //v << 0, 0, 10;
  if (t > 0){
    v(0) = 0.1;
  }
  //Vector v(Vector::Ones(LZ));
  return v;
}

double p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  return 0.0;
}

Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector dxP(X.size());

  return dxP;
}

Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  Vector T(Vector::Zero(X.size()));
  //T(0) = -p_exact(X,t,tag);
  //T(1) = muu(tag)*(cos(w_*t) + sin(w_*t));
  Tensor dxU(grad_u_exact(X,t,tag));
  Tensor I(Tensor::Identity(2,2));
  T = (- p_exact(X,t,tag)*I +  muu(tag)*(dxU + dxU.transpose()))*normal;
  return T;
}

Vector u_initial(Vector const& X, int tag)
{
  return u_exact(X,0,tag);
}

Vector z_initial(Vector const& X, int tag)
{
  return z_exact(X,0,tag);
}

double p_initial(Vector const& X, int tag)
{
  return p_exact(X,0,tag);
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  return N;
}

Vector v_exact(Vector const& X, double t, int tag) //(X,t,tag)
{
  double const x = X(0);
  double const y = X(1);
  Vector v(Vector::Zero(X.size()));

  return v;
}

// posição do contorno
Vector x_exact(Vector const& X, double t, int tag)
{
  Vector r(Vector::Zero(X.size()));
  return r;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  //if (tag == 2) {N(1) = 1;};
  return N;
}

Tensor feature_proj(Vector const& X, double t, int tag)
{
  Tensor f(Tensor::Zero(X.size(), X.size()));
  if (tag == 3){
    f(0,0) = 1;
  }
  return f;
}

Vector force_pp(Vector const& Xi, Vector const& Xj, double Ri, double Rj,
                 double ep1, double ep2, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();
//  if (dij > Ri+Rj+zeta){
//    return f;
//  }
  if (dij <= Ri+Rj){
    f = (1/ep1)*(Ri+Rj-dij)*(Xi - Xj);
  }
  else if((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    f = (1/ep2)*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj);
  }
  return f;
}

Vector force_pw(Vector const& Xi, Vector const& Xj, double Ri,
                 double ew1, double ew2, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double di = (Xi - Xj).norm();
//  if (dij > Ri+Rj+zeta){
//    return f;
//  }
  if (di <= 2*Ri){
    f = (1/ew1)*(2*Ri-di)*(Xi - Xj);
  }
  else if((2*Ri <= di) && (di <= 2*Ri+zeta)){
    f = (1/ew2)*(2*Ri+zeta-di)*(2*Ri+zeta-di)*(Xi - Xj);
  }
  return f;
}

Vector force_ppl(Vector const& Xi, Vector const& Xj, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  f = (zeta/ep)*(Xi - Xj)/(Xi - Xj).norm();
  return f;
}

Vector force_rga(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 Vector const& Gr, double const masj, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();;
  double g = 0.0;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    g   = Gr.norm();
    f   = masj*g*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/(ep*zeta*zeta*dij);
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}

Vector force_rgb(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 Vector const& Gr, double const rhoj, double const rhof, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = dij = (Xi - Xj).norm();
  double g = 0.0;
  double R = 0.0;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    g   = Gr.norm();
    R   = std::max(Ri,Rj);
    f   = (rhoj-rhof)*pi*g*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/(ep*zeta*zeta*dij);
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}

Vector force_rgc(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    f   = (Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/ep;
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}

#endif

// rot solid 2d trap and channel/////////////////////////////////////////////////////////////
#if (true)

double pho(Vector const& X, int tag)
{
//  if (tag == 15)
//  {
    return 1;//e3;///1e4;
//  }
//  else
//  {
//    return 0.0;
//  }
}

double cos_theta0()
{
  return 0.0;
}

double zeta(double u_norm, double angle)
{
  return 0.0;
}

double beta_diss()
{
  return 0.0;
}

double gama(Vector const& X, double t, int tag)
{
  return 0.5;
}

double muu(int tag)
{
//  if (tag == 15)
//  {
    return 1e-2;//1.0*0.1;
//  }
//  else
//  {
//    return 0.0;
//  }
}

Vector force(Vector const& X, double t, int tag)//gravity*pho
{
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(X.size()));
//  if (tag == 15)
//  {
    f(1) = 0; //-1;//-980.0*1.0;//pho(X,tag);//*1e4;//*1e3;
//
//  else
//  {
//    f(1) = 0.0;  //-8e-4*1e4;
//  }
  return f;
}

Vector gravity(Vector const& X, int dim){
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(3*(dim-1)));
  if (dim == 2){
    f(1) = 0; //-1;//-980.0;//-8e-4;  //*1e3;
  }
  else if (dim == 3){
    f(2) = -980.0;  //-8e-4*1e4;
  }
  return f;
}

Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(Vector::Zero(X.size()));
  double Um = 10e-6, H = 200e-6;
  if (false && tag == 4){
    //Um = Um*(1-exp(-t*t*t*t*1e10));//Um*(1-exp(t*t*t*t/1e-14));//Um*(-1-exp(t*t*t*t/1e-14)*0);
    v(0) = Um*4*(H-y)*y/(H*H);
  }
  return v;
}

Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double w2 = 2.0, Um = 0.1, H = 2e-6;
  Tensor dxU(Tensor::Zero(X.size(), X.size()));
  //dxU(0,0) = 0; dxU(0,1) = Um/H - Um*2.0*y/H;
  //dxU(1,0) = 0; dxU(1,1) = 0;

  return dxU;
}

Vector z_exact(Vector const& X, double t, int tag)
{
  double w2 = 0.0;
  int dim = X.size();
  int LZ = 3*(dim-1);
  Vector v(Vector::Zero(LZ)); //v << 0, 0, 10;
  //Vector v(Vector::Ones(LZ));
  return v;
}

double p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  return 0.0;
}

Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector dxP(X.size());

  return dxP;
}

Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  Vector T(Vector::Zero(X.size()));
  //T(0) = -p_exact(X,t,tag);
  //T(1) = muu(tag)*(cos(w_*t) + sin(w_*t));
  Tensor dxU(grad_u_exact(X,t,tag));
  Tensor I(Tensor::Identity(2,2));
  T = (- p_exact(X,t,tag)*I +  muu(tag)*(dxU + dxU.transpose()))*normal;
  return T;
}

Vector u_initial(Vector const& X, int tag)
{
  return u_exact(X,0,tag);
}

Vector z_initial(Vector const& X, int tag)
{
  return z_exact(X,0,tag);
}

double p_initial(Vector const& X, int tag)
{
  return p_exact(X,0,tag);
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  return N;
}

Vector v_exact(Vector const& X, double t, int tag) //(X,t,tag)
{
  double const x = X(0);
  double const y = X(1);
  Vector v(Vector::Zero(X.size()));

  return v;
}

// posição do contorno
Vector x_exact(Vector const& X, double t, int tag)
{
  Vector r(Vector::Zero(X.size()));
  return r;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  //if (tag == 2) {N(1) = 1;};
  return N;
}

Tensor feature_proj(Vector const& X, double t, int tag)
{
  Tensor f(Tensor::Zero(X.size(), X.size()));
  if (false && tag == 3){
    f(0,0) = 1;
  }
  return f;
}

Vector force_pp(Vector const& Xi, Vector const& Xj, double Ri, double Rj,
                 double ep1, double ep2, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();
//  if (dij > Ri+Rj+zeta){
//    return f;
//  }
  if (dij <= Ri+Rj){
    f = (1/ep1)*(Ri+Rj-dij)*(Xi - Xj);
  }
  else if((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    f = (1/ep2)*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj);
  }
  return f;
}

Vector force_pw(Vector const& Xi, Vector const& Xj, double Ri,
                 double ew1, double ew2, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double di = (Xi - Xj).norm();
//  if (dij > Ri+Rj+zeta){
//    return f;
//  }
  if (di <= 2*Ri){
    f = (1/ew1)*(2*Ri-di)*(Xi - Xj);
  }
  else if((2*Ri <= di) && (di <= 2*Ri+zeta)){
    f = (1/ew2)*(2*Ri+zeta-di)*(2*Ri+zeta-di)*(Xi - Xj);
  }
  return f;
}

Vector force_ppl(Vector const& Xi, Vector const& Xj, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  f = (zeta/ep)*(Xi - Xj)/(Xi - Xj).norm();
  return f;
}

Vector force_rga(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 Vector const& Gr, double const masj, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();;
  double g = 0.0;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    g   = Gr.norm();
    f   = masj*g*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/(ep*zeta*zeta*dij);
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}

Vector force_rgb(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 Vector const& Gr, double const rhoj, double const rhof, double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = dij = (Xi - Xj).norm();
  double g = 0.0;
  double R = 0.0;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    g   = Gr.norm();
    R   = std::max(Ri,Rj);
    f   = (rhoj-rhof)*pi*g*(Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/(ep*zeta*zeta*dij);
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}

Vector force_rgc(Vector const& Xi, Vector const& Xj, double const Ri, double const Rj,
                 double ep, double zeta)
{
  Vector f(Vector::Zero(Xi.size()));
  double dij = (Xi - Xj).norm();;
  if ((Ri+Rj <= dij) && (dij <= Ri+Rj+zeta)){
    f   = (Ri+Rj+zeta-dij)*(Ri+Rj+zeta-dij)*(Xi - Xj)/ep;
  }
  //else if (dij < Ri+Rj){cout << "ERROR: penetration!!!!!!!!!!!!!!!!" << endl;}
  return f;
}

#endif




// General functions/////////////////////////////////////////////////////////////

TensorZ MI_tensor(double M, double R, int dim, Tensor3 TI)
{
  TensorZ MI(TensorZ::Zero(3*(dim-1),3*(dim-1)));
  if (dim == 2){
    MI(0,0) = M; MI(1,1) = M; MI(2,2) = TI(2,2); //MI(2,2) = 0.5*M*R*R;
  }
  else if (dim == 3){
    MI(0,0) = M; MI(1,1) = M; MI(2,2) = M;
    MI.block(3,3,5,5) = TI;
    //MI(3,3) = 0.4*M*R*R; MI(4,4) = 0.4*M*R*R; MI(5,5) = 0.4*M*R*R;
  }
  return MI;
}

Tensor RotM(double theta, int dim)
{
  Tensor M(Tensor::Zero(3,3));
  if (dim == 2){
    M(0,0) = cos(theta); M(0,1) = -sin(theta);
    M(1,0) = sin(theta); M(1,1) =  cos(theta);
  }
  if (dim == 3){
    cout << "Not yet" << endl;
  }
  return M;
}

Vector SlipVel(Vector const& X, Vector const& XG, Vector const& normal, int dim, int tag, double theta)
{
  Vector V(Vector::Zero(dim));
  Vector X3(Vector::Zero(3));
  Vector Xp(Vector::Zero(dim));

  double alp = 50.0;
  double bet = 1.0;

  Tensor I(dim,dim);
  I.setIdentity();
  Tensor Pr = I - normal*normal.transpose();

  if (false && dim == 3)
  {
    if (false && tag == 100){
      V(0) = -0.01; V(1) = -0.01;
    }
    else if (tag == 103 || tag == 203){
      V(0) = -1.0;
    }
    V = Pr*V;
  }
  //V.normalize();

  if (true && dim == 2) //this gives structure node swim
  {
    if (tag == 101 || tag == 102){//(tag == 103 || tag == 104)
      theta = pi/4; //5 grados appr
      V(0) = cos(theta); V(1) = sin(theta);
      V = alp*Pr*V;
    }
  }
  double the = 0;
  if (false && dim == 2) //this gives total random swim node
  {
    for (int I = 0; I < 20; I++)
    {
      if (tag == 121 + I)
      {
        the = pi;//rand() % (6);         // v1 in the range 0 to 99
        V(0) = cos(the); V(1) = sin(the);
        V = alp*Pr*V;
      }
    }
  }
  return V;
}

double Dif_coeff(int tag)
{
  return 1;
}
double nuB_coeff(int tag)
{/*
  for (int i = 121; i < 141; i++){
    if (tag == i){
      return 1;
    }
  }
  return 0;*/
  /*
  if (tag == 102)
    return 1;
  else
    return 0;*/
  if (tag == 103 || tag == 104)
    return 1;
  else
    return 0;
}
double sig_coeff(int tag)
{
  return 10;
}
double bbG_coeff(Vector const& X, int tag)
{
  return 1;
}
double per_Elect(int tag)
{
  return 1;
}

Vector traction_maxwell(Vector const& E, Vector const& normal, double eps, int tag)
{
  int d = E.size();
  Vector T(Vector::Zero(3));
  Tensor I(Tensor::Identity(3,3));
  Vector Ec(Vector::Zero(3)), normalc(Vector::Zero(3));
  if (d == 2){
    Ec(0) = E(0); Ec(1) = E(1); normalc(0) = normal(0); normalc(1) = normal(1);
  }
  else{
    Ec = E; normalc = normal;
  }
  T = eps*(Ec*Ec.transpose() - Ec.norm()*Ec.norm()*I)*normal;
  return T;
}

Vector exact_tangent_ellipse(Vector const& X, Vector const& Xc, double theta, double R1, double R2, int dim)
{
  Vector Xcan(Vector::Zero(3));
  Vector Tan(Vector::Zero(dim));
  Vector X3(Vector::Zero(3));
  X3(0) = X(0); X3(1) = X(1);
  Xcan = RotM(-theta, dim)*(X3-Xc);
  Tan(0) = -R1*Xcan(1)/R2;
  Tan(1) =  R2*Xcan(0)/R1;
  Tan.normalize();
  return Tan;
}

Vector exact_normal_ellipse(Vector const& X, Vector const& Xc, double theta, double R1, double R2, int dim)
{
  Vector Xcan(Vector::Zero(3));
  Vector Nrm(Vector::Zero(dim));
  Vector X3(Vector::Zero(3));
  X3(0) = X(0); X3(1) = X(1);
  Xcan = RotM(-theta, dim)*(X3-Xc);
  Nrm(0) = R2*Xcan(0)/R1;
  Nrm(1) = R1*Xcan(1)/R2;
  Nrm.normalize();
  return Nrm;
}

Vector cubic_ellipse(double yb, Vector const& X0, Vector const& X2, Vector const& T0, Vector const& T2, int dim)
{
  Vector Phib(Vector::Zero(dim));
  Phib = X0 + yb*(T0) + yb*yb*(3*X2-3*X0-2*T0-T2) + yb*yb*yb*(2*X0-2*X2+T0+T2);
  return Phib;
}

Vector Dcubic_ellipse(double yb, Vector const& X0, Vector const& X2, Vector const& T0, Vector const& T2, int dim)
{
  Vector DPhib(Vector::Zero(dim));
  DPhib = T0 + 2*yb*(3*X2-3*X0-2*T0-T2) + 3*yb*yb*(2*X0-2*X2+T0+T2);
  return DPhib;
}

Vector curved_Phi(double yb, Vector const& X0, Vector const& X2, Vector const& T0, Vector const& T2, int dim)
{
  Vector Phi(Vector::Zero(dim));
  if (yb != 1)
    Phi = (1.0/(1.0-yb))*(cubic_ellipse(yb,X0,X2,T0,T2,dim)+yb*(X0-X2)-X0);
  else
    Phi = -Dcubic_ellipse(1.0,X0,X2,T0,T2,dim)-X0+X2;

  return Phi;
}

Vector Dcurved_Phi(double yb, Vector const& X0, Vector const& X2, Vector const& T0, Vector const& T2, int dim)
{
  Vector DPhi(Vector::Zero(dim));
  if (yb != 1){
    DPhi = (1.0/((1.0-yb)*(1.0-yb)))*((Dcubic_ellipse(yb,X0,X2,T0,T2,dim)+X0-X2)*(1-yb)
           +cubic_ellipse(yb,X0,X2,T0,T2,dim)+yb*(X0-X2)-X0);}
  else
    {DPhi = -Dcubic_ellipse(1.0,X0,X2,T0,T2,dim)-X0+X2;}

  return DPhi;
}
