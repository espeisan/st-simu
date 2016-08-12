#include <Fepic/CustomEigen>
#include <cmath>
#include <iostream>
using namespace Eigen;
using namespace std;

// space
typedef Matrix<double, Dynamic,1,0,6,1>              Vector;
typedef Matrix<double, Dynamic,Dynamic,RowMajor,3,3> Tensor;
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
TensorZ MI_tensor(double M, double R, int dim);
Tensor RotM(double theta, int dim);
Vector SlipVel(Vector const& X, Vector const& XG, Vector const& normal, int dim, int tag);


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

TensorZ MI_tensor(double M, double R, int dim)
{
  TensorZ MI(TensorZ::Zero(3*(dim-1),3*(dim-1)));
  if (dim == 2){
    MI(0,0) = M; MI(1,1) = M; MI(2,2) = 0.4*M*R*R;
  }
  else if(dim == 3){
    MI(0,0) = M; MI(1,1) = M; MI(2,2) = M;
    MI(3,3) = 0.4*M*R*R; MI(4,4) = 0.4*M*R*R; MI(5,5) = 0.4*M*R*R;
  }
  return MI;
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
    return 1.0;//1.0*0.1;
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

#endif

// rot solid 2d slip vel/////////////////////////////////////////////////////////////
#if (true)

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
    return 1.0;//1.0*0.1;
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
  //Vector v(Vector::Ones(X.size()));  v << 1 , 2;
  Vector v(Vector::Zero(X.size()));
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
  Vector v(Vector::Zero(LZ)); v << 0.0, 0.0, 10.0;
  //if (t > 0){
  //  v(2) = w2;
  //}
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

Vector SlipVel(Vector const& X, Vector const& XG, Vector const& normal, int dim, int tag)
{
  Vector V(Vector::Zero(dim));
  Vector X3(Vector::Zero(3));
  Vector Xp(Vector::Zero(dim));

  double alp = 1;
  double bet = 1;

  Tensor I(dim,dim);
  I.setIdentity();
  Tensor Pr = I - normal*normal.transpose();

  if (dim == 2)
  {
    if (false && tag == 100){
      V(0) = -0.01; V(1) = -0.01;
    }
    else if (tag == 103){
      V(0) = -0.0; V(1) = 0.0;
    }
    V = Pr*V;
  }

  return V;
}

#endif
