//static char help[] = "Navier-Stokes.\n\n";

//
// ALE
//
// pho*( (Ut + (U-Umsh) · nabla)U ) + grad p = muu* div grad U + force
//
// LIMITAÇÕES:
// * triangulo e tetrahedro apenas.
// * no maximo 1 grau de liberdade associado a uma aresta
// * condição de dirichlet na normal só pode valer 0 .. procure por UNORMAL_AQUI
// * elementos isoparamétricos
//
// OBS:
// * a normal no ponto de contato é igual ao limite da normal pela superfície livre
//
//
//


#include "common.hpp"
#include <iomanip>

PetscErrorCode CheckSnesConvergence(SNES snes, PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx);
PetscErrorCode FormJacobian_mesh(SNES snes,Vec Vec_up_1,Mat *Mat_Jac, Mat *prejac, MatStructure *flag, void *ptr);
PetscErrorCode FormFunction_mesh(SNES snes, Vec Vec_up_1, Vec Vec_fun, void *ptr);
PetscErrorCode FormJacobian_fs(SNES snes,Vec Vec_up_1,Mat *Mat_Jac, Mat *prejac, MatStructure *flag, void *ptr);
PetscErrorCode FormFunction_fs(SNES snes, Vec Vec_up_1, Vec Vec_fun, void *ptr);

class AppCtx;
class Statistics;

class GetDataVelocity : public DefaultGetDataVtk
{
public:
  GetDataVelocity(double *q_array_, AppCtx const& user_) : user(user_), q_array(q_array_){}
  //double get_data_r(int nodeid) const;
  void get_vec(int id, Real * vec_out) const;
  int vec_ncomps() const { return user.mesh->spaceDim(); }
  AppCtx const& user;
  double *q_array;
  virtual ~GetDataVelocity() {}
};

class GetDataPressure : public DefaultGetDataVtk
{
public:
  GetDataPressure(double *q_array_, AppCtx const& user_) : user(user_),q_array(q_array_){}
  double get_data_r(int nodeid) const;
  AppCtx const& user;
  double *q_array;
  virtual ~GetDataPressure() {}
};

class GetDataPressCellVersion : public DefaultGetDataVtk
{
public:
  GetDataPressCellVersion(double *q_array_, AppCtx const& user_) : user(user_),q_array(q_array_){}
  double get_data_r(int cellid) const;
  AppCtx const& user;
  double *q_array;
  virtual ~GetDataPressCellVersion() {}
};

class GetDataNormal : public DefaultGetDataVtk
{
public:
  GetDataNormal(double *q_array_, AppCtx const& user_) : user(user_), q_array(q_array_){}
  void get_vec(int id, Real * vec_out) const;
  int vec_ncomps() const { return user.mesh->spaceDim(); }
  AppCtx const& user;
  double *q_array;
  virtual ~GetDataNormal() {}
};

class GetDataMeshVel : public DefaultGetDataVtk
{
public:
  GetDataMeshVel(double *q_array_, AppCtx const& user_) : user(user_), q_array(q_array_){}
  void get_vec(int id, Real * vec_out) const;
  int vec_ncomps() const { return user.mesh->spaceDim(); }
  AppCtx const& user;
  double *q_array;
  virtual ~GetDataMeshVel() {}
};

class GetDataCellTag : public DefaultGetDataVtk
{
public:
  GetDataCellTag(AppCtx const& user_) : user(user_){}
  int get_data_i(int cellid) const;
  AppCtx const& user;
  virtual ~GetDataCellTag() {}
};


AppCtx::AppCtx(int argc, char **argv, bool &help_return, bool &erro)
{
  setUpDefaultOptions();
  if (getCommandLineOptions(argc, argv) )
    help_return = true;
  else
    help_return = false;

  // create some other things
  erro = createFunctionsSpace(); if (erro) return;
  createQuadrature();
}

bool AppCtx::err_checks()
{
  const bool mesh_has_edge_nodes  = mesh->numNodesPerCell() > mesh->numVerticesPerCell();
  bool const u_has_edge_assoc_dof = shape_phi_c->numDofsAssociatedToFacet() +
                                    shape_phi_c->numDofsAssociatedToCorner() > 0;
  if (!mesh_has_edge_nodes && u_has_edge_assoc_dof)
  {
    printf("ERROR: cant be superparametric element\n");
    return true;
  }
  return false;
}

void AppCtx::loadMesh()
{
  mesh.reset( Mesh::create(ECellType(mesh_cell_type),dim) );
  msh_reader.readFileMsh(filename.c_str(), mesh.get());
  vtk_printer.attachMesh(mesh.get());
  vtk_printer.isFamily(family_files);
  vtk_printer.setOutputFileName(filename_out.c_str());
  vtk_printer.setBinaryOutput(true);  //Initial true

  meshAliasesUpdate();

}

void AppCtx::loadDofs()
{
  dofsCreate();
  timer.restart();
  dofsUpdate();
  timer.elapsed("CuthillMcKeeRenumber");

  n_dofs_u_per_cell   = dof_handler[DH_UNKM].getVariable(VAR_U).numDofsPerCell();
  n_dofs_u_per_facet  = dof_handler[DH_UNKM].getVariable(VAR_U).numDofsPerFacet();
  n_dofs_u_per_corner = dof_handler[DH_UNKM].getVariable(VAR_U).numDofsPerCorner();
  n_dofs_p_per_cell   = dof_handler[DH_UNKM].getVariable(VAR_P).numDofsPerCell();
  n_dofs_p_per_facet  = dof_handler[DH_UNKM].getVariable(VAR_P).numDofsPerFacet();
  n_dofs_p_per_corner = dof_handler[DH_UNKM].getVariable(VAR_P).numDofsPerCorner();
  n_dofs_v_per_cell   = dof_handler[DH_UNKM].getVariable(VAR_M).numDofsPerCell();
  n_dofs_v_per_facet  = dof_handler[DH_UNKM].getVariable(VAR_M).numDofsPerFacet();
  n_dofs_v_per_corner = dof_handler[DH_UNKM].getVariable(VAR_M).numDofsPerCorner();
  if (N_Solids){
    n_dofs_z_per_cell   = nodes_per_cell*LZ;   //dof_handler[DH_UNKM].getVariable(VAR_Z).numDofsPerCell();
    n_dofs_z_per_facet  = nodes_per_facet*LZ;  //dof_handler[DH_UNKM].getVariable(VAR_Z).numDofsPerFacet();
    n_dofs_z_per_corner = nodes_per_corner*LZ; //dof_handler[DH_UNKM].getVariable(VAR_Z).numDofsPerCorner();
  }
}

void AppCtx::setUpDefaultOptions()
{
/* global settings */
  dim                    = 2;
  mesh_cell_type         = TRIANGLE3;
  function_space         = 1; // P1P1
  behaviors              = BH_GLS;
  //Re                   = 0.0;
  dt                     = 0.1;
  unsteady               = PETSC_TRUE;
  boundary_smoothing     = PETSC_TRUE;
  steady_tol             = 1.e-6;
  utheta                 = 1;      // time step, theta method (momentum)
  vtheta                 = 1;      // time step, theta method (mesh velocity)
  maxts                  = 10;   // max num of time steps
  finaltime              = -1;
  force_pressure         = PETSC_FALSE;  // elim null space (auto)
  print_to_matlab        = PETSC_FALSE;  // imprime o jacobiano e a função no formato do matlab
  force_dirichlet        = PETSC_TRUE;   // impõe cc ? (para debug)
  full_diriclet          = PETSC_TRUE;
  force_mesh_velocity    = PETSC_FALSE;
  solve_the_sys          = true;   // for debug
  plot_exact_sol         = PETSC_FALSE;
  quadr_degree_cell      = (dim==2) ? 3 : 3;   // ordem de quadratura
  quadr_degree_facet     = 3;
  quadr_degree_corner    = 3;
  quadr_degree_err       = 8; // to compute error
  grow_factor            = 0.05;
  print_step             = 1;
  family_files           = PETSC_TRUE;
  has_convec             = PETSC_TRUE;
  renumber_dofs          = PETSC_FALSE;
  fprint_ca              = PETSC_FALSE;
  nonlinear_elasticity   = PETSC_FALSE;
  mesh_adapt             = PETSC_TRUE;
  fprint_hgv             = PETSC_FALSE;
  filename               = (dim==2 ? "malha/cavity2d-1o.msh" : "malha/cavity3d-1o.msh");

  //Luzia's part
  h_star = 0.02;
  Q_star = 0.8;
  beta_l = 4.00;
  L_min  = 0.02; //0.02;//h_star*sqrt(3)/3.0;
  L_max = 0.2;
  L_range = 0.02;
  L_low = 0.1;//1.0*L_min/L_max;
  L_sup = 1.3; //1.0/L_low;
  TOLad = 0.7;

  n_unknowns_z = 0; n_unknowns_f = 0; n_unknowns_s = 0;
  n_nodes_fsi  = 0; n_nodes_fo   = 0; n_nodes_so   = 0;
  n_dofs_z_per_cell   = 0; n_dofs_z_per_facet  = 0; n_dofs_z_per_corner = 0;
}

bool AppCtx::getCommandLineOptions(int argc, char **/*argv*/)
{
  PetscBool          flg_fin, flg_fout, flg_min, flg_hout;
  char               finaux[PETSC_MAX_PATH_LEN], minaux[PETSC_MAX_PATH_LEN];
  char               foutaux[PETSC_MAX_PATH_LEN], houtaux[PETSC_MAX_PATH_LEN];
  PetscBool          ask_help;

  if (argc == 1)
  {
    cout << "\nusage:\n";
    cout << "\n\t./main `cat args`\n\n";
    cout << "where args is a file with the command line parameters, or\n\n";
    cout << "\t./main -help\n\n";
    cout << "to show options.\n\n";
    return true;
  }

  /* opções do usuário */ //*Esto se imprime cuando se da -help
  PetscOptionsBegin(PETSC_COMM_WORLD, "", "Options for the Navier-Stokes", "ahn?");   //Mandatory to PetscOptions*
  PetscOptionsInt("-dim", "space dimension", "main.cpp", dim, &dim, PETSC_NULL);
  PetscOptionsInt("-mesh_type", "mesh type", "main.cpp", mesh_cell_type, &mesh_cell_type, PETSC_NULL);
  PetscOptionsInt("-function_space", "function_space", "main.cpp", 1, &function_space, PETSC_NULL);
  PetscOptionsInt("-maxts", "maximum number of time steps", "main.cpp", maxts, &maxts, PETSC_NULL);
  //PetscOptionsScalar("-Re", "Reynolds number", "main.cpp", Re, &Re, PETSC_NULL);
  PetscOptionsScalar("-dt", "time step", "main.cpp", dt, &dt, PETSC_NULL);
  PetscOptionsScalar("-utheta", "utheta value", "main.cpp", utheta, &utheta, PETSC_NULL);
  PetscOptionsScalar("-vtheta", "vtheta value", "main.cpp", vtheta, &vtheta, PETSC_NULL);
  PetscOptionsScalar("-sst", "steady state tolerance", "main.cpp", steady_tol, &steady_tol, PETSC_NULL);
  PetscOptionsBool("-print_to_matlab", "print jacobian to matlab", "main.cpp", print_to_matlab, &print_to_matlab, PETSC_NULL);
  PetscOptionsBool("-force_dirichlet", "force dirichlet bound cond", "main.cpp", force_dirichlet, &force_dirichlet, PETSC_NULL);
  PetscOptionsBool("-force_pressure", "force_pressure", "main.cpp", force_pressure, &force_pressure, PETSC_NULL);
  PetscOptionsBool("-plot_es", "plot exact solution", "main.cpp", plot_exact_sol, &plot_exact_sol, PETSC_NULL);
  PetscOptionsBool("-family_files", "plot family output", "main.cpp", family_files, &family_files, PETSC_NULL);
  PetscOptionsBool("-has_convec", "convective term", "main.cpp", has_convec, &has_convec, PETSC_NULL);
  PetscOptionsBool("-unsteady", "unsteady problem", "main.cpp", unsteady, &unsteady, PETSC_NULL);
  PetscOptionsBool("-boundary_smoothing", "boundary_smoothing", "main.cpp", boundary_smoothing, &boundary_smoothing, PETSC_NULL);
  PetscOptionsBool("-force_mesh_velocity", "force_mesh_velocity", "main.cpp", force_mesh_velocity, &force_mesh_velocity, PETSC_NULL);
  PetscOptionsBool("-renumber_dofs", "renumber dofs", "main.cpp", renumber_dofs, &renumber_dofs, PETSC_NULL);
  PetscOptionsBool("-fprint_ca", "print contact angle", "main.cpp", fprint_ca, &fprint_ca, PETSC_NULL);
  PetscOptionsBool("-nonlinear_elasticity", "put a non-linear term in the elasticity problem", "main.cpp", nonlinear_elasticity, &nonlinear_elasticity, PETSC_NULL);
  PetscOptionsBool("-mesh_adapt", "adapt the mesh during simulation", "main.cpp", mesh_adapt, &mesh_adapt, PETSC_NULL);
  PetscOptionsInt("-quadr_e", "quadrature degree (for calculating the error)", "main.cpp", quadr_degree_err, &quadr_degree_err, PETSC_NULL);
  PetscOptionsInt("-quadr_c", "quadrature degree", "main.cpp", quadr_degree_cell, &quadr_degree_cell, PETSC_NULL);
  PetscOptionsInt("-quadr_f", "quadrature degree (facet)", "main.cpp", quadr_degree_facet, &quadr_degree_facet, PETSC_NULL);
  PetscOptionsInt("-quadr_r", "quadrature degree (corner)", "main.cpp", quadr_degree_corner, &quadr_degree_corner, PETSC_NULL);
  PetscOptionsInt("-print_step", "print_step", "main.cpp", print_step, &print_step, PETSC_NULL);
  PetscOptionsScalar("-beta1", "par vel do fluido", "main.cpp", beta1, &beta1, PETSC_NULL);
  PetscOptionsScalar("-beta2", "par vel elastica", "main.cpp", beta2, &beta2, PETSC_NULL);
  PetscOptionsScalar("-finaltime", "the simulation ends at this time.", "main.cpp", finaltime, &finaltime, PETSC_NULL);
  PetscOptionsBool("-ale", "mesh movement", "main.cpp", ale, &ale, PETSC_NULL);
  PetscOptionsGetString(PETSC_NULL,"-fin",finaux,PETSC_MAX_PATH_LEN-1,&flg_fin);
  PetscOptionsGetString(PETSC_NULL,"-fout",foutaux,PETSC_MAX_PATH_LEN-1,&flg_fout);
  PetscOptionsGetString(PETSC_NULL,"-min",minaux,PETSC_MAX_PATH_LEN-1,&flg_min);
  PetscOptionsGetString(PETSC_NULL,"-hout",houtaux,PETSC_MAX_PATH_LEN-1,&flg_hout);
  PetscOptionsHasName(PETSC_NULL,"-help",&ask_help);

  is_mr_ab           = PETSC_TRUE;
  is_bdf3            = PETSC_FALSE;
  is_bdf2            = PETSC_FALSE;
  is_bdf2_bdfe       = PETSC_FALSE;
  is_bdf2_ab         = PETSC_FALSE;
  is_bdf_cte_vel     = PETSC_FALSE;
  is_bdf_euler_start = PETSC_FALSE;
  is_bdf_extrap_cte  = PETSC_FALSE;
  is_basic           = PETSC_FALSE;

  if ((is_bdf2 && utheta!=1) || (is_bdf3 && utheta!=1))
  {
    cout << "ERROR: BDF2/3 with utheta!=1" << endl;
    throw;
  }
  if (is_bdf_extrap_cte && !is_bdf_cte_vel)
  {
    cout << "ERROR: is_bdf_extrap_cte && !is_bdf_cte_vel" << endl;
    throw;
  }
  if (is_bdf2_bdfe && !is_bdf2)
  {
    cout << "ERROR: !is_bdf_bdf_extrap && !is_bdf2" << endl;
    throw;
  }
  if (is_bdf2_ab && !is_bdf2)
  {
    cout << "ERROR: !is_bdf_ab && !is_bdf2" << endl;
    throw;
  }
  if (is_bdf3 && is_bdf2)
  {
    cout << "ERROR: is_bdf3 && is_bdf2" << endl;
    throw;
  }


  if (finaltime < 0)
    finaltime = maxts*dt;
  else
    maxts = 1 + static_cast<int> (round ( finaltime/dt ));

  dirichlet_tags.resize(16);
  PetscBool flg_tags;
  int nmax = dirichlet_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-dir_tags", dirichlet_tags.data(), &nmax, &flg_tags); //looks for -dir_tags and find nmax and gives bool value to flg_tags
  if (flg_tags)
    dirichlet_tags.resize(nmax);
  else
    dirichlet_tags.clear();

  neumann_tags.resize(16);
  nmax = neumann_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-neum_tags", neumann_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
    neumann_tags.resize(nmax);
  else
    neumann_tags.clear();

  interface_tags.resize(16);
  nmax = interface_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-interf_tags", interface_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
    interface_tags.resize(nmax);
  else
    interface_tags.clear();

  solid_tags.resize(16);
  nmax = solid_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-solid_tags", solid_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
    solid_tags.resize(nmax);
  else
    solid_tags.clear();

  triple_tags.resize(16);
  nmax = triple_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-triple_tags", triple_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
    triple_tags.resize(nmax);
  else
    triple_tags.clear();

  periodic_tags.resize(16);
  nmax = periodic_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-periodic_tags", periodic_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
    periodic_tags.resize(nmax);
  else
    periodic_tags.clear();
  if (periodic_tags.size() % 2 != 0)
  {
    std::cout << "Invalid periodic tags\n";
    throw;
  }

  feature_tags.resize(16);
  nmax = feature_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-feature_tags", feature_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
    feature_tags.resize(nmax);
  else
    feature_tags.clear();

  flusoli_tags.resize(1e8);  //cout << flusol_tags.max_size() << endl;
  nmax = flusoli_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-fsi_tags", flusoli_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
	flusoli_tags.resize(nmax);
  else
	flusoli_tags.clear();

  fluidonly_tags.resize(16);  //cout << flusol_tags.max_size() << endl;
  nmax = fluidonly_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-fonly_tags", fluidonly_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
	fluidonly_tags.resize(nmax);
  else
	fluidonly_tags.clear();

  solidonly_tags.resize(16);
  nmax = solidonly_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-sonly_tags", solidonly_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
    solidonly_tags.resize(nmax);
  else
    solidonly_tags.clear();
  N_Solids = solidonly_tags.size();
  LZ = dim*(dim+1)/2; //3*(dim-1);

  slipvel_tags.resize(16);
  nmax = slipvel_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-slipv_tags", slipvel_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
    {slipvel_tags.resize(nmax); is_slipv = PETSC_TRUE;}
  else
    {slipvel_tags.clear(); is_slipv = PETSC_FALSE;}

  PetscOptionsEnd();   //Finish PetscOptions*

  switch (function_space)
  {
    case P1P1:
    {
      if(dim==2) mesh_cell_type = TRIANGLE3;
      else       mesh_cell_type = TETRAHEDRON4;
      break;
    }
    case P1bP1_c:
    {
      if(dim==2) mesh_cell_type = TRIANGLE3;
      else       mesh_cell_type = TETRAHEDRON4;
      break;
    }
    case P2bPm1_c:
    {
      if(dim==2) mesh_cell_type = TRIANGLE6;
      else       mesh_cell_type = TETRAHEDRON10;
      break;
    }
    case P2P1:
    {
      if(dim==2) mesh_cell_type = TRIANGLE6;
      else       mesh_cell_type = TETRAHEDRON10;
      break;
    }
    case P1bP1:
    {
      if(dim==2) mesh_cell_type = TRIANGLE3;
      else       mesh_cell_type = TETRAHEDRON4;
      break;
    }
    case P2P0:
    {
      if(dim==2) mesh_cell_type = TRIANGLE6;
      else       mesh_cell_type = TETRAHEDRON10;
      break;
    }
    case P2bPm1:
    {
      if(dim==2) mesh_cell_type = TRIANGLE6;
      else       mesh_cell_type = TETRAHEDRON10;
      break;
    }
    case P1P1unstable:
    {
      if(dim==2) mesh_cell_type = TRIANGLE3;
      else       mesh_cell_type = TETRAHEDRON4;
      break;
    }
    case P2bP1:
    {
      if(dim==2) mesh_cell_type = TRIANGLE6;
      else       mesh_cell_type = TETRAHEDRON10;
      break;
    }
    default:
    {
      cout << "invalid function space " << function_space << endl;
      throw;
    }
  };


  if(ask_help)
    return true;

  if (flg_fin)
    filename.assign(finaux);

  if (flg_fout)
    filename_out.assign(foutaux);

  if (flg_min){
    filemass.assign(minaux);
    int N = N_Solids;
    double mass = 0.0, rad = 0.0, xg = 0.0, yg = 0.0, zg = 0.0, vol = 0.0;
    Vector3d Xg;
    ifstream is;
    is.open(filemass.c_str(),ifstream::in);
    if (!is.good()) {cout << "mass file not found" << endl; throw;}
    for (; N > 0; N--){
      is >> mass; is >> rad; is >> vol;
      is >> xg; is >> yg; if(dim == 3){is >> zg;}
      Xg << xg, yg, zg;
      MV.push_back(mass); RV.push_back(rad); VV.push_back(vol);
      XG_1.push_back(Xg); XG_0.push_back(Xg); XG_aux.push_back(Xg);
      theta_1.push_back(0.0); theta_0.push_back(0.0); theta_aux.push_back(0.0);
    }
    is.close();
    MV.resize(N_Solids); RV.resize(N_Solids); VV.resize(N_Solids);
    XG_1.resize(N_Solids); XG_0.resize(N_Solids); XG_aux.resize(N_Solids);
    theta_1.resize(N_Solids); theta_0.resize(N_Solids); theta_aux.resize(N_Solids);
  }

  if (flg_hout){
    filehist_out.assign(houtaux);
    fprint_hgv = PETSC_TRUE;
  }

  //if (neumann_tags.size() + interface_tags.size() == 0 || force_pressure)
  if (force_pressure)
  {
    full_diriclet = PETSC_TRUE;
    force_pressure = PETSC_TRUE;
  }
  else
    force_pressure = PETSC_FALSE;

  if (!unsteady)
  {
    if (ale) {
      printf("ale with steady problem?\n");
      //throw;
    }

    if (utheta != 1)
    {
      cout << "WARNING!!!\n";
      cout << ">>> steady problem .. setting utheta to 1" << endl;
      utheta = 1;
    }

  }

  return false;
}

bool AppCtx::createFunctionsSpace()
{
  EShapeType velo_shape, pres_shape;
  bool is_simplex = ctype2cfamily(ECellType(mesh_cell_type)) == SIMPLEX;
  ECellType cell_type = ECellType(mesh_cell_type);

  switch (function_space)
  {
    case 1: // P1P1 (or Q1Q1) GLS stabilization
      {
        behaviors = BH_GLS;
        velo_shape = is_simplex ? P1 : Q1;
        pres_shape = is_simplex ? P1 : Q1;
      }
      break;
    case 2: // P1+P1 with bubble condensation
      {
        behaviors = BH_bble_condens_PnPn;
        velo_shape = P1;
        pres_shape = P1;
      }
      break;
    case 3: // P2+Pm1 with bubble condensation and pressure gradient elimination
      {
        behaviors = BH_bble_condens_CR;
        velo_shape = P2;
        pres_shape = P0;
      }
      break;
    case 4: // P2P1 (or Q2Q1)
      {
        behaviors = 0;
        velo_shape = is_simplex ? P2 : Q2;
        pres_shape = is_simplex ? P1 : Q1;
      }
      break;
    case 5: // P1+P1 traditional
      {
        behaviors = 0;
        velo_shape = P1ph ;
        pres_shape = P1;
      }
      break;
    case 6: // P2P0
      {
        behaviors = 0;
        velo_shape = P2;
        pres_shape = P0;
      }
      break;
    case 7: // P2+Pm1 full
      {
        behaviors = 0;
        velo_shape = P2ph;
        pres_shape = Pm1;
      }
      break;
    case 8: // P1P1 unstable
      {
        behaviors = 0;
        velo_shape = P1;
        pres_shape = P1;
      }
      break;
    case 9: // P2bP1 bubble condensation
      {
        behaviors = BH_bble_condens_PnPn;
        velo_shape = is_simplex ? P2 : Q2;
        pres_shape = is_simplex ? P1 : Q1;
      }
      break;
    default:
      {
        behaviors = BH_GLS;
        velo_shape = is_simplex ? P1 : Q1;
        pres_shape = is_simplex ? P1 : Q1;
      }

  }

  shape_phi_c.reset(ShapeFunction::create(cell_type, velo_shape));
  shape_psi_c.reset(ShapeFunction::create(cell_type, pres_shape));
  shape_qsi_c.reset(ShapeFunction::create(cell_type));
  shape_phi_f.reset(ShapeFunction::create(facetof(cell_type), facetof(velo_shape)));
  shape_psi_f.reset(ShapeFunction::create(facetof(cell_type), facetof(pres_shape)));
  shape_qsi_f.reset(ShapeFunction::create(facetof(cell_type)));
  shape_phi_r.reset(ShapeFunction::create(facetof(facetof(cell_type)), facetof(facetof(velo_shape))));
  shape_psi_r.reset(ShapeFunction::create(facetof(facetof(cell_type)), facetof(facetof(pres_shape))));
  shape_qsi_r.reset(ShapeFunction::create(facetof(facetof(cell_type))));

  shape_bble.reset(ShapeFunction::create(cell_type, BUBBLE));

  pres_pres_block = false;
  if (behaviors & (BH_bble_condens_PnPn | BH_GLS))
    pres_pres_block = true;

  return false;
}

void AppCtx::createQuadrature()
{
  ECellType ct = ECellType(mesh_cell_type);

  quadr_cell.reset( Quadrature::create(ECellType(ct)) );
  quadr_cell->setOrder(quadr_degree_cell);
  n_qpts_cell = quadr_cell->numPoints();

  quadr_facet.reset( Quadrature::create(facetof(ECellType(ct))) );
  quadr_facet->setOrder(quadr_degree_facet);
  n_qpts_facet = quadr_facet->numPoints();

  quadr_corner.reset( Quadrature::create(facetof(facetof(ECellType(ct)))) );
  quadr_corner->setOrder(quadr_degree_corner);
  n_qpts_corner = quadr_corner->numPoints();

  quadr_err.reset( Quadrature::create(ECellType(ct)) );
  quadr_err->setOrder(quadr_degree_err);
  n_qpts_err = quadr_err->numPoints();

}

void AppCtx::meshAliasesUpdate()
{
  n_nodes  = mesh->numNodes();
  n_cells  = mesh->numCells();
  n_facets = mesh->numFacets();
  n_corners= mesh->numCorners();
  n_nodes_total  = mesh->numNodesTotal();   // includes disableds
  n_cells_total  = mesh->numCellsTotal();   // includes disableds
  n_facets_total = mesh->numFacetsTotal();  // includes disableds
  n_corners_total= mesh->numCornersTotal(); // includes disableds
  nodes_per_cell  = mesh->numNodesPerCell();
  nodes_per_facet = mesh->numNodesPerFacet();
  nodes_per_corner= mesh->numNodesPerCorner();
}

void AppCtx::dofsCreate()
{
  // mesh velocity
  dof_handler[DH_MESH].setMesh(mesh.get());
  dof_handler[DH_MESH].addVariable("mesh_velo",  shape_qsi_c.get(), dim);
  //dof_handler[DH_MESH].setVariablesRelationship(blocks.data());

//  int const * fo_tag = &fluidonly_tags[0];
  dof_handler[DH_UNKM].setMesh(mesh.get());
//    dof_handler[DH_UNKM].addVariable("velo_fluid", shape_phi_c.get(), dim, fluidonly_tags.size(), fo_tag);
  dof_handler[DH_UNKM].addVariable("velo_fluid", shape_phi_c.get(), dim);
  dof_handler[DH_UNKM].addVariable("pres_fluid", shape_psi_c.get(), 1);
  dof_handler[DH_UNKM].getVariable(VAR_P).setType(SPLITTED_BY_REGION_CELL,0,0);
  if (N_Solids && false){
    int const * fsi_tag = &flusoli_tags[0];
    dof_handler[DH_UNKM].addVariable("velo_solid", shape_phi_c.get(), 3*(dim-1), flusoli_tags.size(), fsi_tag);

//    int const * fo_tag = &fluidonly_tags[0];
//    int const * so_tag = &solidonly_tags[0];
//    dof_handler[DH_FOON].setMesh(mesh.get());
//    dof_handler[DH_FOON].addVariable("vel_ofluid", shape_phi_c.get(), dim, fluidonly_tags.size(), fo_tag);
//    dof_handler[DH_FOON].addVariable("vel_osolid", shape_phi_c.get(), dim, fluidonly_tags.size(), so_tag);
  }

}

bool isPeriodic(Point const* p1, Point const* p2, int dim)
{
  double const TOL = 1.e-9;

  if (dim == 2)
  {
    if (fabs(p1->getCoord(0) - p2->getCoord(0)) < TOL)
      return true;
    if (fabs(p1->getCoord(1) - p2->getCoord(1)) < TOL)
      return true;
  }
  else // dim == 3
  {
    int c = 0;
    if (fabs(p1->getCoord(0) - p2->getCoord(0)) < TOL)
      c++;
    if (fabs(p1->getCoord(1) - p2->getCoord(1)) < TOL)
      c++;
    if (fabs(p1->getCoord(2) - p2->getCoord(2)) < TOL)
      c++;

    if (c==2)
      return true;
  }

  return false;
}

void AppCtx::dofsUpdate()
{
  dof_handler[DH_UNKM].SetUp();
  //if (renumber_dofs)
  //  dof_handler[DH_UNKS].CuthillMcKeeRenumber();
  n_nodes_fsi = 0; n_nodes_fo = 0; n_nodes_so = 0;
  std::vector<int> dofs1;
  std::vector<int> dofs2;
  NN_Solids.assign(N_Solids,0);
  int tag, nod_id;

//  int dofs_a[10];
//  int dofs_b[10];

  // apply periodic boundary conditions here
  {
    if ((mesh_cell_type != TRIANGLE3) && (mesh_cell_type != TETRAHEDRON4) && !periodic_tags.empty() )
    {
      std::cout << "Periodic boundary conditions is not allowed with high order mesh\n";
      throw;
    }

      point_iterator point1 = mesh->pointBegin();
      point_iterator point1_end = mesh->pointEnd();
      point_iterator point2, point2_end;

      for (; point1 != point1_end; ++point1)
      {
    	tag = point1->getTag();
        nod_id = is_in_id(tag,flusoli_tags); if (nod_id) {NN_Solids[nod_id-1]++; n_nodes_fsi++;}
        if (is_in(tag,interface_tags)) {n_nodes_fsi++;}
    	if (is_in(tag,fluidonly_tags)) {n_nodes_fo++;}
    	if (is_in(tag,solidonly_tags)) {n_nodes_so++;}

        for (int i = 0; i < (int)periodic_tags.size(); i+=2)
        {

          int const tag1 = periodic_tags[i];
          int const tag2 = periodic_tags[i+1];

          if (tag1 != point1->getTag())
            continue;

          point2 = mesh->pointBegin();
          point2_end = mesh->pointEnd();
          for (; point2 != point2_end; ++point2)
          {
            if (tag2 != point2->getTag())
              continue;

            if (isPeriodic(&*point1, &*point2, dim))
            {
              for (int var = 0; var < 2; ++var)
              {
//                int ndpv = dof_handler[DH_UNKS].getVariable(var).numDofsPerVertex();

//                dof_handler[DH_UNKS].getVariable(var).getVertexDofs(dofs_a, &*point1);
//                dof_handler[DH_UNKS].getVariable(var).getVertexDofs(dofs_b, &*point2);
//                for (int k = 0; k < ndpv; ++k)
//                {
//                  dofs1.push_back(dofs_a[k]);
//                  dofs2.push_back(dofs_b[k]);
//                }
              }

            }
          }
        }
      }
  }
//  dof_handler[DH_UNKS].linkDofs(dofs1.size(), dofs1.data(), dofs2.data());
//  n_unknowns = dof_handler[DH_UNKS].numDofs();

  dof_handler[DH_MESH].SetUp();
  n_dofs_v_mesh = dof_handler[DH_MESH].numDofs();

  dof_handler[DH_UNKM].linkDofs(dofs1.size(), dofs1.data(), dofs2.data());
  n_unknowns_u = dof_handler[DH_UNKM].getVariable(VAR_U).numPositiveDofs();
  n_unknowns_p = dof_handler[DH_UNKM].getVariable(VAR_P).numPositiveDofs();
  n_unknowns_z = N_Solids*3*(dim-1);
  //cout << n_unknowns_z << "   " <<  dof_handler[DH_UNKM].getVariable(VAR_Z).numPositiveDofs();

//    dof_handler[DH_FOON].SetUp();
//    n_unknowns_f = dof_handler[DH_FOON].getVariable(VAR_F).numPositiveDofs();
//    n_unknowns_s = dof_handler[DH_FOON].getVariable(VAR_S).numPositiveDofs();

  n_unknowns_fs = n_unknowns_u + n_unknowns_p + N_Solids*3*(dim-1); //- dim*n_nodes_fsi;

  //for (unsigned int l = 0; l < NN_Solids.size(); l++) cout << NN_Solids[l] << " ";
}

PetscErrorCode AppCtx::allocPetscObjs()
{
  printf( "Allocating PETSC objects... ");

  PetscErrorCode      ierr;
  PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INDEX);

  //Vec Vec_res;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_res_fs);                     CHKERRQ(ierr);  //creating Vec_res empty PETSc vector
  ierr = VecSetSizes(Vec_res_fs, PETSC_DECIDE, n_unknowns_fs);         CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_res_fs);                                CHKERRQ(ierr);

  //Vec Vec_uzp_0;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_uzp_0);                      CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_uzp_0, PETSC_DECIDE, n_unknowns_fs);          CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_uzp_0);                                 CHKERRQ(ierr);
  ierr = VecSetOption(Vec_uzp_0, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);

  //Vec Vec_uzp_1;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_uzp_1);                       CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_uzp_1, PETSC_DECIDE, n_unknowns_fs);           CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_uzp_1);                                  CHKERRQ(ierr);
  ierr = VecSetOption(Vec_uzp_1, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);

  //Vec Vec_duzp;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_duzp);                       CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_duzp, PETSC_DECIDE, n_unknowns_fs);           CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_duzp);                                  CHKERRQ(ierr);
  ierr = VecSetOption(Vec_duzp, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);

  //Vec Vec_uzp_m1;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_uzp_m1);                       CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_uzp_m1, PETSC_DECIDE, n_unknowns_fs);           CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_uzp_m1);                                  CHKERRQ(ierr);
  ierr = VecSetOption(Vec_uzp_m1, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);

  if (is_bdf3)
  {
    //Vec Vec_duzp_0;
    ierr = VecCreate(PETSC_COMM_WORLD, &Vec_duzp_0);                          CHKERRQ(ierr);
    ierr = VecSetSizes(Vec_duzp_0, PETSC_DECIDE, n_unknowns_fs);              CHKERRQ(ierr);
    ierr = VecSetFromOptions(Vec_duzp_0);                                     CHKERRQ(ierr);
    ierr = VecSetOption(Vec_duzp_0, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);

    //Vec Vec_uzp_m2;
    ierr = VecCreate(PETSC_COMM_WORLD, &Vec_uzp_m2);                          CHKERRQ(ierr);
    ierr = VecSetSizes(Vec_uzp_m2, PETSC_DECIDE, n_unknowns_fs);              CHKERRQ(ierr);
    ierr = VecSetFromOptions(Vec_uzp_m2);                                     CHKERRQ(ierr);
    ierr = VecSetOption(Vec_uzp_m2, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);
  }

  //Vec Vec_v_mid
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_v_mid);                CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_v_mid, PETSC_DECIDE, n_dofs_v_mesh);    CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_v_mid);                           CHKERRQ(ierr);

  //Vec Vec_v_1
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_v_1);                  CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_v_1, PETSC_DECIDE, n_dofs_v_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_v_1);                             CHKERRQ(ierr);

  //Vec Vec_x_0;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_x_0);                  CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_x_0, PETSC_DECIDE, n_dofs_v_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_x_0);                             CHKERRQ(ierr);

  //Vec Vec_x_1;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_x_1);                  CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_x_1, PETSC_DECIDE, n_dofs_v_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_x_1);                             CHKERRQ(ierr);

  if (is_bdf3)
  {
    //Vec Vec_x_aux;
    ierr = VecCreate(PETSC_COMM_WORLD, &Vec_x_aux);              CHKERRQ(ierr);
    ierr = VecSetSizes(Vec_x_aux, PETSC_DECIDE, n_dofs_v_mesh);  CHKERRQ(ierr);
    ierr = VecSetFromOptions(Vec_x_aux);                         CHKERRQ(ierr);
  }

  //Vec Vec_normal;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_normal);               CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_normal, PETSC_DECIDE, n_dofs_v_mesh);   CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_normal);                          CHKERRQ(ierr);

  //Vec Vec_res_m;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_res_m);                CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_res_m, PETSC_DECIDE, n_dofs_v_mesh);    CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_res_m);                           CHKERRQ(ierr);

  if (is_slipv)
  {
    //Vec Vec_slipvel_0;
    ierr = VecCreate(PETSC_COMM_WORLD, &Vec_slipv_0);               CHKERRQ(ierr);
    ierr = VecSetSizes(Vec_slipv_0, PETSC_DECIDE, n_dofs_v_mesh);   CHKERRQ(ierr);
    ierr = VecSetFromOptions(Vec_slipv_0);                          CHKERRQ(ierr);
    ierr = VecSetOption(Vec_slipv_0, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);

    //Vec Vec_slipvel_1;
    ierr = VecCreate(PETSC_COMM_WORLD, &Vec_slipv_1);               CHKERRQ(ierr);
    ierr = VecSetSizes(Vec_slipv_1, PETSC_DECIDE, n_dofs_v_mesh);   CHKERRQ(ierr);
    ierr = VecSetFromOptions(Vec_slipv_1);                          CHKERRQ(ierr);
    ierr = VecSetOption(Vec_slipv_1, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);

    //Vec Vec_slipvel_m1;
    ierr = VecCreate(PETSC_COMM_WORLD, &Vec_slipv_m1);               CHKERRQ(ierr);
    ierr = VecSetSizes(Vec_slipv_m1, PETSC_DECIDE, n_dofs_v_mesh);   CHKERRQ(ierr);
    ierr = VecSetFromOptions(Vec_slipv_m1);                          CHKERRQ(ierr);
    ierr = VecSetOption(Vec_slipv_m1, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);

    if (is_bdf3)
    {
      //Vec Vec_slipvel_m2;
      ierr = VecCreate(PETSC_COMM_WORLD, &Vec_slipv_m2);               CHKERRQ(ierr);
      ierr = VecSetSizes(Vec_slipv_m2, PETSC_DECIDE, n_dofs_v_mesh);   CHKERRQ(ierr);
      ierr = VecSetFromOptions(Vec_slipv_m2);                          CHKERRQ(ierr);
      ierr = VecSetOption(Vec_slipv_m2, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);
    }

    //Vec Vec_uzp_0_ns;
    ierr = VecCreate(PETSC_COMM_WORLD, &Vec_uzp_0_ns);                      CHKERRQ(ierr);
    ierr = VecSetSizes(Vec_uzp_0_ns, PETSC_DECIDE, n_unknowns_fs);          CHKERRQ(ierr);
    ierr = VecSetFromOptions(Vec_uzp_0_ns);                                 CHKERRQ(ierr);
    ierr = VecSetOption(Vec_uzp_0_ns, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);

    //Vec Vec_uzp_1;
    ierr = VecCreate(PETSC_COMM_WORLD, &Vec_uzp_1_ns);                       CHKERRQ(ierr);
    ierr = VecSetSizes(Vec_uzp_1_ns, PETSC_DECIDE, n_unknowns_fs);           CHKERRQ(ierr);
    ierr = VecSetFromOptions(Vec_uzp_1_ns);                                  CHKERRQ(ierr);
    ierr = VecSetOption(Vec_uzp_1_ns, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);
  }

  std::vector<int> nnz;

  //  ------------------------------------------------------------------
  //  ----------------------- Mesh -------------------------------------
  //  ------------------------------------------------------------------
  nnz.clear();
  int n_mesh_dofs = dof_handler[DH_MESH].numDofs();
  {
    std::vector<SetVector<int> > table;
    dof_handler[DH_MESH].getSparsityTable(table); // TODO: melhorar desempenho, função mt lenta

    nnz.resize(n_mesh_dofs, 0);

    //FEP_PRAGMA_OMP(parallel for)
    for (int i = 0; i < n_mesh_dofs; ++i)
      {nnz[i] = table[i].size();}  //cout << nnz[i] << " ";} //cout << endl;
  }

  //Mat Mat_Jac_m;
  ierr = MatCreate(PETSC_COMM_WORLD, &Mat_Jac_m);                                        CHKERRQ(ierr);
  ierr = MatSetType(Mat_Jac_m,MATSEQAIJ);                                                CHKERRQ(ierr);
  ierr = MatSetSizes(Mat_Jac_m, PETSC_DECIDE, PETSC_DECIDE, n_mesh_dofs, n_mesh_dofs);   CHKERRQ(ierr);
  //ierr = MatSetFromOptions(Mat_Jac_m);                                                 CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(Mat_Jac_m,  0, nnz.data());                           CHKERRQ(ierr);
  //ierr = MatSetOption(Mat_Jac_m,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);                  CHKERRQ(ierr);
  ierr = MatSetOption(Mat_Jac_m,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);             CHKERRQ(ierr);
  ierr = MatSetOption(Mat_Jac_m,MAT_SYMMETRIC,PETSC_TRUE);                               CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD, &snes_m);                                          CHKERRQ(ierr);
  ierr = SNESSetFunction(snes_m, Vec_res_m, FormFunction_mesh, this);                    CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes_m, Mat_Jac_m, Mat_Jac_m, FormJacobian_mesh, this);         CHKERRQ(ierr);

  SNESGetLineSearch(snes_m,&linesearch);
  SNESLineSearchSetType(linesearch,SNESLINESEARCHBASIC);

  ierr = SNESGetKSP(snes_m,&ksp_m);                                                  CHKERRQ(ierr);
  ierr = KSPGetPC(ksp_m,&pc_m);                                                      CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp_m,Mat_Jac_m,Mat_Jac_m,SAME_NONZERO_PATTERN);            CHKERRQ(ierr);
  //ierr = KSPSetType(ksp_m,KSPGMRES);                                               CHKERRQ(ierr);
  ierr = KSPSetType(ksp_m,KSPPREONLY);                                               CHKERRQ(ierr);
  ierr = PCSetType(pc_m,PCLU);                                                       CHKERRQ(ierr);
  //ierr = PCFactorSetMatOrderingType(pc_m, MATORDERINGNATURAL);                     CHKERRQ(ierr);
  //ierr = KSPSetTolerances(ksp_m,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);  CHKERRQ(ierr);
  //ierr = SNESSetApplicationContext(snes_m,this);

  if(!nonlinear_elasticity)  //for linear systems
  {
    ierr = SNESSetType(snes_m, SNESKSPONLY); CHKERRQ(ierr);
  }
  ierr = SNESSetFromOptions(snes_m); CHKERRQ(ierr);  //prints Newton iterations information SNES Function norm

//~ #ifdef PETSC_HAVE_MUMPS
  //~ PCFactorSetMatSolverPackage(pc_m,MATSOLVERMUMPS);
//~ #endif

  //ierr = SNESMonitorSet(snes_m, SNESMonitorDefault, 0, 0); CHKERRQ(ierr);
  //ierr = SNESMonitorSet(snes_m,Monitor,0,0);CHKERRQ(ierr);
  //ierr = SNESSetTolerances(snes_m,1.e-11,1.e-11,1.e-11,10,PETSC_DEFAULT);
  //ierr = SNESSetFromOptions(snes_m); CHKERRQ(ierr);
  //ierr = SNESLineSearchSet(snes_m, SNESLineSearchNo, &user); CHKERRQ(ierr);


  //  ------------------------------------------------------------------
  //  ----------------------- Fluid ------------------------------------
  //  ------------------------------------------------------------------

  nnz.clear();
  {
    nnz.resize(dof_handler[DH_UNKM].numDofs()+N_Solids*LZ,0);
    std::vector<SetVector<int> > tabla;
    dof_handler[DH_UNKM].getSparsityTable(tabla); // TODO: melhorar desempenho, função mt lenta
    //FEP_PRAGMA_OMP(parallel for)
      for (int i = 0; i < n_unknowns_fs - N_Solids*LZ; ++i)
        nnz[i] = tabla[i].size();
  }

    ierr = MatCreate(PETSC_COMM_WORLD, &Mat_Jac_fs);                                      CHKERRQ(ierr);
    ierr = MatSetSizes(Mat_Jac_fs, PETSC_DECIDE, PETSC_DECIDE, n_unknowns_fs, n_unknowns_fs);   CHKERRQ(ierr);
    ierr = MatSetFromOptions(Mat_Jac_fs);                                                 CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(Mat_Jac_fs, 0, nnz.data());                          CHKERRQ(ierr);
    ierr = MatSetOption(Mat_Jac_fs,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);           CHKERRQ(ierr);
    ierr = SNESCreate(PETSC_COMM_WORLD, &snes_fs);                                        CHKERRQ(ierr);
    ierr = SNESSetFunction(snes_fs, Vec_res_fs, FormFunction_fs, this);                   CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes_fs, Mat_Jac_fs, Mat_Jac_fs, FormJacobian_fs, this);       CHKERRQ(ierr);
    ierr = SNESSetConvergenceTest(snes_fs,CheckSnesConvergence,this,PETSC_NULL);          CHKERRQ(ierr);
    ierr = SNESGetKSP(snes_fs,&ksp_fs);                                                   CHKERRQ(ierr);
    ierr = KSPGetPC(ksp_fs,&pc_fs);                                                       CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp_fs,Mat_Jac_fs,Mat_Jac_fs,SAME_NONZERO_PATTERN);            CHKERRQ(ierr);

    ierr = SNESSetFromOptions(snes_fs); CHKERRQ(ierr);

  printf(" done.\n");
  PetscFunctionReturn(0);
}

void AppCtx::matrixColoring()
{
  printf("Matrix coloring... ");  //clock_t t; t = clock();

  VectorXi                     mapU_c(n_dofs_u_per_cell);
  VectorXi                     mapZ_c(nodes_per_cell*LZ);
  VectorXi                     mapP_c(n_dofs_p_per_cell);

  //Z Dofs re-organization
  int tag;
  int nod_id, nod_is;

  MatrixXd Afsloc = MatrixXd::Zero(n_dofs_u_per_cell, n_dofs_u_per_cell);
  MatrixXd Dfsloc = MatrixXd::Zero(n_dofs_p_per_cell, n_dofs_u_per_cell);
  MatrixXd Gfsloc = MatrixXd::Zero(n_dofs_u_per_cell, n_dofs_p_per_cell);
  MatrixXd Efsloc = MatrixXd::Zero(n_dofs_p_per_cell, n_dofs_p_per_cell);

  MatrixXd Z1fsloc = MatrixXd::Zero(n_dofs_u_per_cell,nodes_per_cell*LZ);
  MatrixXd Z2fsloc = MatrixXd::Zero(nodes_per_cell*LZ,n_dofs_u_per_cell);
  MatrixXd Z3fsloc = MatrixXd::Zero(nodes_per_cell*LZ,nodes_per_cell*LZ);
  MatrixXd Z4fsloc = MatrixXd::Zero(nodes_per_cell*LZ,n_dofs_p_per_cell);
  MatrixXd Z5fsloc = MatrixXd::Zero(n_dofs_p_per_cell,nodes_per_cell*LZ);

  std::vector<bool>   SV(N_Solids,false);  //solid visited history
  bool                SFI = false;         //solid-fluid interception

  cell_iterator cell = mesh->cellBegin();
  cell_iterator cell_end = mesh->cellEnd();

  for (; cell != cell_end; ++cell){
    // mapeamento do local para o global:
    dof_handler[DH_UNKM].getVariable(VAR_U).getCellDofs(mapU_c.data(), &*cell);  // U global ids for the current cell
    dof_handler[DH_UNKM].getVariable(VAR_P).getCellDofs(mapP_c.data(), &*cell);  // P global ids for the current cell
    // mapping correction for valocity unknowns case FS
    if (N_Solids){
      mapZ_c = -VectorXi::Ones(nodes_per_cell*LZ);
      SFI = false;
      for (int j = 0; j < nodes_per_cell; ++j){
        tag = mesh->getNodePtr(cell->getNodeId(j))->getTag();
        nod_id = is_in_id(tag,flusoli_tags);
        nod_is = is_in_id(tag,solidonly_tags);
        if (nod_id || nod_is){
          for (int l = 0; l < LZ; l++){
            mapZ_c(j*LZ + l) = n_unknowns_u + n_unknowns_p + LZ*(nod_id+nod_is) - LZ + l;
          }
          for (int l = 0; l < dim; l++){
            mapU_c(j*dim + l) = -1;
          }
          SFI = true;
        }
      }
    }
    //cout << endl << mapU_c.transpose() << "   "<< mapZ_c.transpose() << endl;
    MatSetValues(Mat_Jac_fs, mapU_c.size(), mapU_c.data(), mapU_c.size(), mapU_c.data(), Afsloc.data(), ADD_VALUES);
    MatSetValues(Mat_Jac_fs, mapU_c.size(), mapU_c.data(), mapP_c.size(), mapP_c.data(), Gfsloc.data(), ADD_VALUES);
    MatSetValues(Mat_Jac_fs, mapP_c.size(), mapP_c.data(), mapU_c.size(), mapU_c.data(), Dfsloc.data(), ADD_VALUES);
    MatSetValues(Mat_Jac_fs, mapP_c.size(), mapP_c.data(), mapP_c.size(), mapP_c.data(), Efsloc.data(), ADD_VALUES);
    if (SFI){
      MatSetValues(Mat_Jac_fs, mapU_c.size(), mapU_c.data(), mapZ_c.size(), mapZ_c.data(), Z1fsloc.data(), ADD_VALUES);
      MatSetValues(Mat_Jac_fs, mapZ_c.size(), mapZ_c.data(), mapU_c.size(), mapU_c.data(), Z2fsloc.data(), ADD_VALUES);
      MatSetValues(Mat_Jac_fs, mapZ_c.size(), mapZ_c.data(), mapZ_c.size(), mapZ_c.data(), Z3fsloc.data(), ADD_VALUES);
      MatSetValues(Mat_Jac_fs, mapZ_c.size(), mapZ_c.data(), mapP_c.size(), mapP_c.data(), Z4fsloc.data(), ADD_VALUES);
      MatSetValues(Mat_Jac_fs, mapP_c.size(), mapP_c.data(), mapZ_c.size(), mapZ_c.data(), Z5fsloc.data(), ADD_VALUES);
    }
  }//end for cell

  Assembly(Mat_Jac_fs); //View(Mat_Jac_fs,"MatColored.m","JJ");

  printf(" done. \n");  //t = clock() - t; printf ("It took me %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
}

void AppCtx::printMatlabLoader()
{
  FILE *fp = fopen("loadmat.m", "w");
  //fprintf(fp, "clear;\n"                   );
  fprintf(fp, "jacob;\n"                   );
  fprintf(fp, "clear zzz;\n"               );
  fprintf(fp, "B=Jac;\n"                   );
  fprintf(fp, "B(B!=0)=1;\n"               );
  fprintf(fp, "nU = %d;\n",dof_handler[DH_UNKM].getVariable(VAR_U).numPositiveDofs() );
  fprintf(fp, "nP = %d;\n",dof_handler[DH_UNKM].getVariable(VAR_P).numPositiveDofs() );
  fprintf(fp, "nT = nU + nP;\n"            );
  fprintf(fp, "K=Jac(1:nU,1:nU);\n"        );
  fprintf(fp, "G=Jac(1:nU,nU+1:nT);\n"     );
  fprintf(fp, "D=Jac(nU+1:nT,1:nU);\n"     );
  fprintf(fp, "E=Jac(nU+1:nT,nU+1:nT);\n"  );
  fprintf(fp, "\n"                        );
  fprintf(fp, "rhs;\n"                     );
  fprintf(fp, "f=res(1:nU);\n"             );
  fprintf(fp, "g=res(nU+1:nT);\n"          );
  fclose(fp);
}

// must be called after loadDofs
void AppCtx::evaluateQuadraturePts()
{
  // avaliando phi_c e as matrizes de transf nos pontos de quadratura.
  phi_c.resize(n_qpts_cell);
  psi_c.resize(n_qpts_cell);
  qsi_c.resize(n_qpts_cell);
//std::cout << std::endl << n_qpts_cell << " jeje " << phi_c.size() << phi_c.rows() << phi_c.cols();
  phi_f.resize(n_qpts_facet);
  psi_f.resize(n_qpts_facet);
  qsi_f.resize(n_qpts_facet);

  phi_r.resize(n_qpts_corner);
  psi_r.resize(n_qpts_corner);
  qsi_r.resize(n_qpts_corner);

  dLphi_c.resize(n_qpts_cell);
  dLpsi_c.resize(n_qpts_cell);
  dLqsi_c.resize(n_qpts_cell);

  dLphi_f.resize(n_qpts_facet);
  dLpsi_f.resize(n_qpts_facet);
  dLqsi_f.resize(n_qpts_facet);

  dLphi_r.resize(n_qpts_corner);
  dLpsi_r.resize(n_qpts_corner);
  dLqsi_r.resize(n_qpts_corner);

  bble.resize(n_qpts_cell);
  dLbble.resize(n_qpts_cell);


  for (int qp = 0; qp < n_qpts_cell; ++qp)
  {
    phi_c[qp].resize(n_dofs_u_per_cell/dim);
    psi_c[qp].resize(n_dofs_p_per_cell);
    qsi_c[qp].resize(nodes_per_cell);

    dLphi_c[qp].resize(n_dofs_u_per_cell/dim, dim);
    dLpsi_c[qp].resize(n_dofs_p_per_cell, dim);
    dLqsi_c[qp].resize(nodes_per_cell, dim);

    dLbble[qp].resize(dim);
    bble[qp] = shape_bble->eval(quadr_cell->point(qp), 0);

    for (int n = 0; n < n_dofs_u_per_cell/dim; ++n)
    {
      phi_c[qp][n] = shape_phi_c->eval(quadr_cell->point(qp), n); //std::cout << phi_c[qp][n] << "\n";
      for (int d = 0; d < dim; ++d)
      {
        /* dLphi_c nao depende de qp no caso de funcoes lineares */
        dLphi_c[qp](n, d) = shape_phi_c->gradL(quadr_cell->point(qp), n, d);
      }
    }

    for (int n = 0; n < n_dofs_p_per_cell; ++n)
    {
      psi_c[qp][n] = shape_psi_c->eval(quadr_cell->point(qp), n);
      for (int d = 0; d < dim; ++d)
      {
        /* dLpsi_c nao depende de qp no caso de funcoes lineares */
        dLpsi_c[qp](n, d) = shape_psi_c->gradL(quadr_cell->point(qp), n, d);
      }
    }

    for (int n = 0; n < nodes_per_cell; ++n)
    {
      qsi_c[qp][n] = shape_qsi_c->eval(quadr_cell->point(qp), n);
      for (int d = 0; d < dim; ++d)
      {
        /* dLqsi_c nao depende de qp no caso de funcoes lineares */
        dLqsi_c[qp](n, d) = shape_qsi_c->gradL(quadr_cell->point(qp), n, d);
      }
    }

    for (int d = 0; d < dim; ++d)
    {
      dLbble[qp](d) = shape_bble->gradL(quadr_cell->point(qp), 0, d);
    }

  } // end qp

  // pts de quadratura no contorno
  for (int qp = 0; qp < n_qpts_facet; ++qp)
  {
    phi_f[qp].resize(n_dofs_u_per_facet/dim);
    psi_f[qp].resize(n_dofs_p_per_facet);
    qsi_f[qp].resize(n_dofs_v_per_facet/dim);

    dLphi_f[qp].resize(n_dofs_u_per_facet/dim, dim-1);
    dLpsi_f[qp].resize(n_dofs_p_per_facet, dim-1);
    dLqsi_f[qp].resize(n_dofs_v_per_facet/dim, dim-1);

    for (int n = 0; n < n_dofs_u_per_facet/dim; ++n)
    {
      phi_f[qp][n] = shape_phi_f->eval(quadr_facet->point(qp), n);
      for (int d = 0; d < dim-1; ++d)
      {
        /* dLphi_f nao depende de qp no caso de funcoes lineares */
        dLphi_f[qp](n, d) = shape_phi_f->gradL(quadr_facet->point(qp), n, d);
      }
    }

    for (int n = 0; n < n_dofs_p_per_facet; ++n)
    {
      psi_f[qp][n] = shape_psi_f->eval(quadr_facet->point(qp), n);
      for (int d = 0; d < dim-1; ++d)
      {
        /* dLpsi_f nao depende de qp no caso de funcoes lineares */
        dLpsi_f[qp](n, d) = shape_psi_f->gradL(quadr_facet->point(qp), n, d);
      }
    }

    for (int n = 0; n < nodes_per_facet; ++n)
    {
      qsi_f[qp][n] = shape_qsi_f->eval(quadr_facet->point(qp), n);
      for (int d = 0; d < dim-1; ++d)
      {
        /* dLqsi_f nao depende de qp no caso de funcoes lineares */
        dLqsi_f[qp](n, d) = shape_qsi_f->gradL(quadr_facet->point(qp), n, d);
      }
    }

  } // end qp

  // pts de quadratura no contorno do contorno
  for (int qp = 0; qp < n_qpts_corner; ++qp)
  {
    phi_r[qp].resize(n_dofs_u_per_corner/dim);
    psi_r[qp].resize(n_dofs_p_per_corner);
    qsi_r[qp].resize(nodes_per_corner);

    dLphi_r[qp].resize(n_dofs_u_per_corner/dim, 1 /* = dim-2*/);
    dLpsi_r[qp].resize(n_dofs_p_per_corner, 1 /* = dim-2*/);
    dLqsi_r[qp].resize(nodes_per_corner, 1 /* = dim-2*/);

    for (int n = 0; n < n_dofs_u_per_corner/dim; ++n)
    {
      phi_r[qp][n] = shape_phi_r->eval(quadr_corner->point(qp), n);
      for (int d = 0; d < dim-2; ++d)
      {
        /* dLphi_r nao depende de qp no caso de funcoes lineares */
        dLphi_r[qp](n, d) = shape_phi_r->gradL(quadr_corner->point(qp), n, d);
      }
    }

    for (int n = 0; n < n_dofs_p_per_corner; ++n)
    {
      psi_r[qp][n] = shape_psi_r->eval(quadr_corner->point(qp), n);
      for (int d = 0; d < dim-2; ++d)
      {
        /* dLpsi_r nao depende de qp no caso de funcoes lineares */
        dLpsi_r[qp](n, d) = shape_psi_r->gradL(quadr_corner->point(qp), n, d);
      }
    }

    for (int n = 0; n < nodes_per_corner; ++n)
    {
      qsi_r[qp][n] = shape_qsi_r->eval(quadr_corner->point(qp), n);
      for (int d = 0; d < dim-2; ++d)
      {
        /* dLqsi_r nao depende de qp no caso de funcoes lineares */
        dLqsi_r[qp](n, d) = shape_qsi_r->gradL(quadr_corner->point(qp), n, d);
      }
    }
  } // end qp

  //     Quadrature
  //     to compute the error
  //
  // velocity
  phi_err.resize(n_qpts_err);         // shape function evaluated at quadrature points
  dLphi_err.resize(n_qpts_err);       // matriz de gradiente no elemento unitário
  // pressure
  psi_err.resize(n_qpts_err);         // shape function evaluated at quadrature points
  dLpsi_err.resize(n_qpts_err);       // matriz de gradiente no elemento unitário
  // mesh
  qsi_err.resize(n_qpts_err);         // shape function evaluated at quadrature points
  dLqsi_err.resize(n_qpts_err);       // matriz de gradiente no elemento unitário

  for (int qp = 0; qp < n_qpts_err; ++qp)
  {
    phi_err[qp].resize(n_dofs_u_per_cell/dim);
    psi_err[qp].resize(n_dofs_p_per_cell);
    qsi_err[qp].resize(nodes_per_cell);

    dLphi_err[qp].resize(n_dofs_u_per_cell/dim, dim);
    dLpsi_err[qp].resize(n_dofs_p_per_cell, dim);
    dLqsi_err[qp].resize(nodes_per_cell, dim);

    for (int n = 0; n < n_dofs_u_per_cell/dim; ++n)
    {
      phi_err[qp][n] = shape_phi_c->eval(quadr_err->point(qp), n);
      for (int d = 0; d < dim; ++d)
      {
        /* dLphi_err nao depende de qp no caso de funcoes lineares */
        dLphi_err[qp](n, d) = shape_phi_c->gradL(quadr_err->point(qp), n, d);
      }
    }

    for (int n = 0; n < n_dofs_p_per_cell; ++n)
    {
      psi_err[qp][n] = shape_psi_c->eval(quadr_err->point(qp), n);
      for (int d = 0; d < dim; ++d)
      {
        /* dLpsi_err nao depende de qp no caso de funcoes lineares */
        dLpsi_err[qp](n, d) = shape_psi_c->gradL(quadr_err->point(qp), n, d);
      }
    }

    for (int n = 0; n < nodes_per_cell; ++n)
    {
      qsi_err[qp][n] = shape_qsi_c->eval(quadr_err->point(qp), n);
      for (int d = 0; d < dim; ++d)
      {
        /* dLqsi_err nao depende de qp no caso de funcoes lineares */
        dLqsi_err[qp](n, d) = shape_qsi_c->gradL(quadr_err->point(qp), n, d);
      }
    }

  } // end qp

  Vector Xc = Vector::Zero(dim);
  bool is_simplex = ctype2cfamily(ECellType(mesh_cell_type)) == SIMPLEX;
  if (is_simplex)
    for (int i = 0; i < dim; ++i)
      Xc(i) = 1./(dim+1);
  qsi_c_at_center.resize(nodes_per_cell);
  for (int i = 0; i < nodes_per_cell; ++i)
    qsi_c_at_center(i) = shape_qsi_c->eval(Xc.data(), i);

  // facets func derivatives on your nodes
  if (dim==2)
  {
    std::vector<double> parametric_pts;
    parametric_pts = genLineParametricPts(  ctypeDegree(ECellType(mesh_cell_type))  );

    dLphi_nf.resize(parametric_pts.size());

    for (int k = 0; k < (int)parametric_pts.size(); ++k)
    {
      dLphi_nf[k].resize(n_dofs_u_per_facet/dim, dim-1);
      for (int n = 0; n < n_dofs_u_per_facet/dim; ++n)
        for (int d = 0; d < dim-1; ++d)
          dLphi_nf[k](n, d) = shape_phi_f->gradL(&parametric_pts[k], n, d);
    }
  }
  else
  if(dim==3)
  {
    std::vector<Eigen::Vector2d> parametric_pts;
    parametric_pts = genTriParametricPts(  ctypeDegree(ECellType(mesh_cell_type))  );

    dLphi_nf.resize(parametric_pts.size());

    for (int k = 0; k < (int)parametric_pts.size(); ++k)
    {
      dLphi_nf[k].resize(n_dofs_u_per_facet/dim, dim-1);
      for (int n = 0; n < n_dofs_u_per_facet/dim; ++n)
        for (int d = 0; d < dim-1; ++d)
          dLphi_nf[k](n, d) = shape_phi_f->gradL(parametric_pts[k].data(), n, d);
    }
  }

}

void AppCtx::onUpdateMesh()
{
  allocPetscObjs();
  matrixColoring();
}

PetscErrorCode AppCtx::setInitialConditions()
{
  PetscErrorCode      ierr(0);

  Vector    Uf(dim), Zf(LZ), Vs(dim), Nr(dim);
  Vector    X(dim);
  Tensor    R(dim,dim);
  VectorXi  dofs(dim), dofs_mesh(dim);  //global components unknowns enumeration from point dofs
  VectorXi  dofs_fs(LZ);
  Vector3d  Xg, XG_temp, Us;
  int       nod_id, nod_is, tag, nod_vs, nodsum;

  VecZeroEntries(Vec_v_mid);  //this size(V) = size(X) = size(U)
  VecZeroEntries(Vec_v_1);
  VecZeroEntries(Vec_x_0);
  VecZeroEntries(Vec_x_1);
  VecZeroEntries(Vec_normal);
  VecZeroEntries(Vec_res_fs); //this size is the [U,P,Z] sol vec size
  VecZeroEntries(Vec_uzp_0);
  VecZeroEntries(Vec_uzp_1);
  VecZeroEntries(Vec_uzp_m1);
  VecZeroEntries(Vec_duzp); //borrar
  if (is_bdf3){
    VecZeroEntries(Vec_x_aux);
    VecZeroEntries(Vec_uzp_m2);
    VecZeroEntries(Vec_duzp_0); //borrar
  }
  if (is_slipv){
    VecZeroEntries(Vec_slipv_0);
    VecZeroEntries(Vec_slipv_1);
    VecZeroEntries(Vec_slipv_m1);
    if (is_bdf3){
      VecZeroEntries(Vec_slipv_m2);
    }
    VecZeroEntries(Vec_uzp_0_ns);
    VecZeroEntries(Vec_uzp_1_ns);
  }

  copyMesh2Vec(Vec_x_0);  //copy initial mesh coordinates to Vec_x_0 = (a1,a2,a3,b1,b2,b3,c1,c2,c3,...)
  copyMesh2Vec(Vec_x_1);  //copy initial mesh coordinates to Vec_x_1

  // normals at boundary points (n01,n02,n03, n11,n12,n13, n21,n22,n23, ..., 0) following node order .geo
  getVecNormals(&Vec_x_0, Vec_normal);  //std::cout << VecView(Vec_normal,PETSC_VIEWER_STDOUT_WORLD)
  //Assembly(Vec_normal);  View(Vec_normal,"matrizes/nrm.m","nrmm");//VecView(Vec_uzp_0,PETSC_VIEWER_STDOUT_WORLD);

  // compute mesh sizes
  {
    mesh_sizes.reserve(n_nodes_total + 1*n_nodes_total/10);  //modify the capacity of mesh_sizes
    mesh_sizes.assign(n_nodes_total, 0.0);  //put zeros in the mesh_size vector (n_nodes_total components)

    // stores how much edges were connected until now at each vertex
    std::vector<int> counter(n_nodes_total, 0);

    VectorXi edge_nodes(3);
    Vector Xa(dim), Xb(dim);
    Real h;

    const int n_edges_total = (dim==2) ? mesh->numFacetsTotal() : mesh->numCornersTotal();
    CellElement const* edge(NULL);

    //FEP_PRAGMA_OMP(for nowait)
    for (int a = 0; a < n_edges_total; ++a)  //total edges ordered by element ID in .geo
    {
      if (dim == 2)
        edge = mesh->getFacetPtr(a);
      else
        edge = mesh->getCornerPtr(a);

      if (edge==NULL || edge->isDisabled())
        continue;

      if (dim == 2)
        mesh->getFacetNodesId(a, edge_nodes.data());
      else
        mesh->getCornerNodesId(a, edge_nodes.data());

      mesh->getNodePtr(edge_nodes[0])->getCoord(Xa.data(),dim);  //coords node1 of current edge
      mesh->getNodePtr(edge_nodes[1])->getCoord(Xb.data(),dim);  //coords node2 of current edge
      h = (Xa-Xb).norm();  //cout << h << "     " << edge_nodes.transpose() << "     ";

      mesh_sizes[edge_nodes[0]] = (mesh_sizes[edge_nodes[0]]*counter[edge_nodes[0]] + h)/(counter[edge_nodes[0]] + 1.0);  //cout << mesh_sizes[edge_nodes[0]] << "  ";
      mesh_sizes[edge_nodes[1]] = (mesh_sizes[edge_nodes[1]]*counter[edge_nodes[1]] + h)/(counter[edge_nodes[1]] + 1.0);  //cout << mesh_sizes[edge_nodes[1]] << endl;

      ++counter[edge_nodes[0]];
      ++counter[edge_nodes[1]];
    } // end n_edges_total //cout << mesh_sizes[2] << " " << mesh_sizes[3] << endl;

  } // end compute mesh sizes

  std::vector<bool>   SV(N_Solids,false);  //solid visited history

  // velocidade inicial e pressao inicial, guardados em Vec_uzp_0
  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    tag = point->getTag();
    point->getCoord(X.data(),dim);

    // press
    if (mesh->isVertex(&*point))
    {
      getNodeDofs(&*point,DH_UNKM,VAR_P,dofs.data());
      double p_in = p_initial(X, tag);
      VecSetValues(Vec_uzp_0,1,dofs.data(),&p_in,INSERT_VALUES);
      VecSetValues(Vec_uzp_1,1,dofs.data(),&p_in,INSERT_VALUES);
    }

    // vel
    getNodeDofs(&*point, DH_UNKM, VAR_U, dofs.data());
    nod_id = is_in_id(tag,flusoli_tags);
    nod_is = is_in_id(tag,solidonly_tags);
    nod_vs = is_in_id(tag,slipvel_tags);
    nodsum = nod_id+nod_is+nod_vs;
    if (nodsum){
      Zf = z_initial(X, tag);
      Uf = SolidVel(X, XG_0[nodsum-1], Zf, dim);
      VecSetValues(Vec_uzp_1, dim, dofs.data(), Uf.data(), INSERT_VALUES);
      if (nod_vs){
        getNodeDofs(&*point, DH_MESH, VAR_M, dofs_mesh.data());
        VecGetValues(Vec_normal, dim, dofs_mesh.data(), Nr.data());
        Vs = SlipVel(X, XG_0[nod_vs-1], Nr, dim, tag);
        VecSetValues(Vec_slipv_0, dim, dofs_mesh.data(), Vs.data(), INSERT_VALUES);
        Uf = Uf + Vs;  //cout << tag << "  " << X(0)-3 << " " << X(1)-3<< "  " << Uf.transpose() << endl;
      }
      VecSetValues(Vec_uzp_0, dim, dofs.data(), Uf.data(), INSERT_VALUES);
      if (!SV[nodsum-1]){
        for (int l = 0; l < LZ; l++){
          dofs_fs(l) = n_unknowns_u + n_unknowns_p + LZ*(nodsum-1) + l;
        }
        VecSetValues(Vec_uzp_0, LZ, dofs_fs.data(), Zf.data(), INSERT_VALUES);  //cout << dofs_fs.transpose() << endl;
        VecSetValues(Vec_uzp_1, LZ, dofs_fs.data(), Zf.data(), INSERT_VALUES);
        Xg = XG_0[nodsum-1];
        Xg(0) = Xg(0)+dt*Zf(0); Xg(1) = Xg(1)+dt*Zf(1); if (dim == 3){Xg(2) = Xg(2)+dt*Zf(2);}
        XG_1[nodsum-1] = Xg;
        theta_1[nodsum-1] = theta_0[nodsum-1] + dt*Zf(2);
        SV[nodsum-1] = true;
      }
    }
    else{
      Uf = u_initial(X, tag);
      VecSetValues(Vec_uzp_0, dim, dofs.data(), Uf.data(), INSERT_VALUES);  //cout << dofs.transpose() << endl;
      VecSetValues(Vec_uzp_1, dim, dofs.data(), Uf.data(), INSERT_VALUES);
    }

  } // end point loop

  //Assembly(Vec_uzp_0);  View(Vec_uzp_0,"matrizes/vuzp0.m","vuzp0m");//VecView(Vec_uzp_0,PETSC_VIEWER_STDOUT_WORLD);
  //Assembly(Vec_uzp_1);  View(Vec_uzp_1,"matrizes/vuzp1.m","vuzp1m");//VecView(Vec_uzp_0,PETSC_VIEWER_STDOUT_WORLD);
  //Assembly(Vec_slipv_0);  View(Vec_slipv_0,"matrizes/slip0.m","slipm");//VecView(Vec_uzp_0,PETSC_VIEWER_STDOUT_WORLD);
//  VecCopy(Vec_uzp_0,Vec_uzp_1);  //PetscInt size1; //u_unk+z_unk+p_unk //VecGetSize(Vec_up_0,&size1);
  //VecCopy(Vec_slipv_0, Vec_slipv_1);
  //plotFiles();
  // remember: Vec_normals follows the Vec_x_1
  calcMeshVelocity(Vec_x_0, Vec_uzp_0, Vec_uzp_1, 1.0, Vec_v_mid, 0.0); //Vec_up_0 = Vec_up_1, vtheta = 1.0, mean b.c. = Vec_up_1 = Vec_v_mid init. guess SNESSolve
  // move the mesh
  VecWAXPY(Vec_x_1, dt, Vec_v_mid, Vec_x_0); // Vec_x_1 = dt*Vec_v_mid + Vec_x_0 // for zero Dir. cond. solution lin. elast. is Vec_v_mid = 0
  if (N_Solids) updateSolidMesh();
  //VecView(Vec_v_mid,PETSC_VIEWER_STDOUT_WORLD); //VecView(Vec_x_0,PETSC_VIEWER_STDOUT_WORLD); VecView(Vec_x_1,PETSC_VIEWER_STDOUT_WORLD);
  //double Ar; Vector Gr = getAreaMassCenterSolid(1,Ar); cout << Gr.transpose() <<"  "<< XG_0[0].transpose() <<"  "<< XG_1[0].transpose() <<"  "<< Ar << endl;

  if (false && is_mr_ab){
    //copyVec2Mesh(Vec_x_1);
    //plotFiles();
    PetscFunctionReturn(0);
  }

  if (is_basic){
    if (N_Solids){moveCenterMass(0.0);}
    PetscFunctionReturn(0);
  }

  if (ale)
  {
    printf("Initial conditions:\n");
    for (int m=0; m<1+is_bdf3; ++m)
    {
      setUPInitialGuess();  //applies b.c. to Vec_uzp_1, applied once per time_step

      for (int i = 0; i < 1; ++i)  // predictor-corrector to initialize
      {
        printf("\tIterations %d\n", i);
        // * SOLVE THE SYSTEM *
        if (solve_the_sys)
        { //VecView(Vec_uzp_1,PETSC_VIEWER_STDOUT_WORLD); //VecView(Vec_uzp_0,PETSC_VIEWER_STDOUT_WORLD);
          ierr = SNESSolve(snes_fs,PETSC_NULL,Vec_uzp_1);  CHKERRQ(ierr);
          if (N_Solids){updateSolidVel();}
        }  //View(Vec_uzp_0,"matrizes/vuzp_0.m","vuzpm_0"); View(Vec_uzp_1,"matrizes/vuzp_1.m","vuzpm_1");
        // * SOLVE THE SYSTEM *

        // update
        if (ale && m==0)
        {
          //calcMeshVelocity(Vec_x_0, Vec_up_0, Vec_up_1, 1.0, Vec_v_mid, 0.0); // Euler (tem que ser esse no começo)
          calcMeshVelocity(Vec_x_0, Vec_uzp_0, Vec_uzp_1, 0.5, Vec_v_mid, 0.0); // Adams-Bashforth
          // move the mesh
          VecWAXPY(Vec_x_1, dt, Vec_v_mid, Vec_x_0); // Vec_x_1 = Vec_v_mid*dt + Vec_x_0
          if (N_Solids){
            moveCenterMass(0.5);
            updateSolidMesh();
          }
          //Gr = getAreaMassCenterSolid(1,Ar); cout << Gr.transpose() <<"  "<< XG_0[0].transpose() <<"  "<< XG_1[0].transpose() <<"  "<< Ar << endl;
          //compute normal for the next time step, at n+1/2
          {
            Vec Vec_x_mid;
            int xsize;
            double *xarray;
            VecGetSize(Vec_x_0, &xsize);  //VecView(Vec_res_fs,PETSC_VIEWER_STDOUT_WORLD);
            VecGetArray(Vec_res_fs, &xarray);
            //prototipo no petsc-dev: VecCreateSeqWithArray(MPI_Comm comm,PetscInt bs,PetscInt n,const PetscScalar array[],Vec *V)
            //to save the vector Vec_x_mid in the array xarray, Vec_res_fs is just a reference
            VecCreateSeqWithArray(MPI_COMM_SELF, 1, xsize, xarray, &Vec_x_mid);
            Assembly(Vec_x_mid);

            VecCopy(Vec_x_0, Vec_x_mid);
            VecAXPY(Vec_x_mid,1.,Vec_x_1);
            VecScale(Vec_x_mid, 0.5);
            getVecNormals(&Vec_x_mid, Vec_normal);

            VecDestroy(&Vec_x_mid);
            VecRestoreArray(Vec_res_fs, &xarray);
          }

        }//end if ale
      }// end for 10 loops

      if (is_bdf3)
      {
        current_time += dt;
        time_step += 1;
        if (m==0)
        {
          // saving uzp n-1 and n-2
          //VecCopy(Vec_uzp_1,Vec_uzp_m1);
          ////-----
          VecCopy(Vec_uzp_0,Vec_uzp_m2);
          /*///-----
          VecCopy(Vec_uzp_1, Vec_duzp_0);
          VecAXPY(Vec_duzp_0,-1.0,Vec_uzp_0);
          VecScale(Vec_duzp_0, 1./dt);  //duzp_0=(uzp_1-uzp_0)/dt
          /*///-----
          VecCopy(Vec_x_0, Vec_x_aux);
          copyVec2Mesh(Vec_x_1);
          VecCopy(Vec_x_1, Vec_x_0);
          if (N_Solids) {XG_aux = XG_0; XG_0 = XG_1;}
          VecCopy(Vec_uzp_1, Vec_uzp_0);
          //View(Vec_uzp_0,"matrizes/uzp_0.m","uzp0");
          //View(Vec_duzp_0,"matrizes/duzp_0.m","duzp0");
          //View(Vec_uzp_m2,"matrizes/uzp_m2.m","uzpm2");
        }
        else
        {
          ////-----
          VecCopy(Vec_uzp_0,Vec_uzp_m1);
          /*///-----
          VecCopy(Vec_uzp_1, Vec_duzp);
          VecAXPY(Vec_duzp,-1.0,Vec_uzp_0);
          VecScale(Vec_duzp, 1./dt);  //duzp=(uzp_1-uzp_0)/dt
          /*///-----
          copyVec2Mesh(Vec_x_1);
          //View(Vec_uzp_1,"matrizes/uzp_1.m","uzp1");
          //View(Vec_duzp,"matrizes/duzp_1.m","duzp");
          //View(Vec_uzp_m1,"matrizes/uzp_m1.m","uzpm1");

          // solving geometry before the (u,p)
          // extrapolation of Vec_x_1, \bar{X^{n+1}} = 3X^{n}-3X^{n-1}+\bar{X^{n-2}}
          VecScale(Vec_x_1, 3.0);
          VecAXPY(Vec_x_1,-3.0,Vec_x_0);
          VecAXPY(Vec_x_1, 1.0,Vec_x_aux);
          // copyMesh2Vec(Vec_x_0); // we need Vec_x_0 later, don't touch it
          // estraga Vec_uzp_0, usado porque não se precisa mais
          ////-----
          VecScale(Vec_uzp_0,-3.0);
          VecAXPY(Vec_uzp_0,3.0,Vec_uzp_1);
          VecAXPY(Vec_uzp_0,1.0,Vec_uzp_m2);
          //View(Vec_uzp_0,"matrizes/uzp_0.m","uzp0");
          /*///-----
          VecCopy(Vec_duzp, Vec_uzp_0);
          VecScale(Vec_uzp_0, 2.0);
          VecAXPY(Vec_uzp_0,-1.0,Vec_duzp_0);
          VecScale(Vec_uzp_0, dt);
          VecAXPY(Vec_uzp_0, 1.0,Vec_uzp_1); // Vec_up_0 tem agora a extrapolacao do futuro Vec_up_1
          //View(Vec_uzp_0,"matrizes/uzp_0o.m","uzp0o");
          /*///-----

          // calc V^{n+1} and update with D_{3}X^{n+1} = dt*V^{n+1}
          calcMeshVelocity(Vec_x_0, Vec_uzp_0, Vec_uzp_1, 0.0, Vec_v_1, current_time);  //vtheta = 0 porque Vec_uzp_0 JÁ É o extrapolado
          VecCopy(Vec_v_1, Vec_x_1);
          VecScale(Vec_x_1, 6./11.*dt);
          VecAXPY(Vec_x_1, 2./11.,Vec_x_aux);
          VecCopy(Vec_x_0, Vec_x_aux);
          VecAXPY(Vec_x_1,-9./11.,Vec_x_0);
          copyMesh2Vec(Vec_x_0);
          VecAXPY(Vec_x_1,18./11.,Vec_x_0);
          if (N_Solids){moveCenterMass(0.0);}// update center of mass D_{3}XG^{n+1} = dt*V^{n+1}
          VecCopy(Vec_uzp_1, Vec_uzp_0);
          VecCopy(Vec_v_1,Vec_v_mid);
        }
      }//end if bdf3
    }//end for bdf3 for m
  }// end if(ale)
  // normals
  getVecNormals(&Vec_x_1, Vec_normal);

  PetscFunctionReturn(0);
}

PetscErrorCode AppCtx::checkSnesConvergence(SNES snes, PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason)
{
  PetscErrorCode ierr;

  ierr = SNESConvergedDefault(snes,it,xnorm,pnorm,fnorm,reason,NULL); CHKERRQ(ierr);

  // se não convergiu, não terminou ou não é um método ale, retorna
  if (*reason<=0 || !ale)
  {
    return ierr;
  }
  else
  {
    // se não é a primeira vez que converge
    if (converged_times)
    {
      return ierr;
    }
    else
    {

      //Vec *Vec_up_k = &Vec_up_1;
      ////SNESGetSolution(snes,Vec_up_k);
      //
      //copyMesh2Vec(Vec_x_1);
      //calcMeshVelocity(*Vec_up_k, Vec_x_1, Vec_v_mid, current_time+dt);
      //moveMesh(Vec_x_1, Vec_vmsh_0, Vec_v_mid, 0.5);
      //
      //// mean mesh velocity
      //VecAXPY(Vec_v_mid,1,Vec_vmsh_0);
      //VecScale(Vec_v_mid, 0.5);
      //
      //*reason = SNES_CONVERGED_ITERATING;
      //++converged_times;

      return ierr;
    }
  }


/*
  converged
  SNES_CONVERGED_FNORM_ABS         =  2,  ||F|| < atol
  SNES_CONVERGED_FNORM_RELATIVE    =  3,  ||F|| < rtol*||F_initial||
  SNES_CONVERGED_PNORM_RELATIVE    =  4,  Newton computed step size small; || delta x || < stol
  SNES_CONVERGED_ITS               =  5,  maximum iterations reached
  SNES_CONVERGED_TR_DELTA          =  7,
   diverged
  SNES_DIVERGED_FUNCTION_DOMAIN    = -1,  the new x location passed the function is not in the domain of F
  SNES_DIVERGED_FUNCTION_COUNT     = -2,
  SNES_DIVERGED_LINEAR_SOLVE       = -3,  the linear solve failed
  SNES_DIVERGED_FNORM_NAN          = -4,
  SNES_DIVERGED_MAX_IT             = -5,
  SNES_DIVERGED_LINE_SEARCH        = -6,  the line search failed
  SNES_DIVERGED_LOCAL_MIN          = -8,  || J^T b || is small, implies converged to local minimum of F()
  SNES_CONVERGED_ITERATING         =  0
*/
}

PetscErrorCode AppCtx::setUPInitialGuess()
{
  // set U^{n+1} b.c.

  VectorXi    u_dofs_fs(dim), p_dofs(dim);
  VectorXi    x_dofs(dim);
  int         tag;
  Vector      X1(dim);
  Vector      U1(dim);

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    tag = point->getTag();

    getNodeDofs(&*point, DH_MESH, VAR_M, x_dofs.data());
    getNodeDofs(&*point, DH_UNKM, VAR_U, u_dofs_fs.data());

    VecGetValues(Vec_x_1, dim, x_dofs.data(), X1.data());

    if (is_in(tag, dirichlet_tags))
    {
      U1 = u_exact(X1, current_time+dt, tag);  //current_time+dt 'cause we want U^{n+1}
      VecSetValues(Vec_uzp_1, dim, u_dofs_fs.data(), U1.data(), INSERT_VALUES);
    }
    else if (is_in(tag, solid_tags) || is_in(tag, feature_tags) || is_in(tag, triple_tags))
    {
      U1.setZero();
      VecSetValues(Vec_uzp_1, dim, u_dofs_fs.data(), U1.data(), INSERT_VALUES);
    }

  } // end for point

  Assembly(Vec_uzp_1);

  PetscFunctionReturn(0);
}

PetscErrorCode AppCtx::solveTimeProblem()
{
  PetscErrorCode      ierr(0);

  double   initial_volume = getMeshVolume(), final_volume;
  Vector   X(dim), Xe(dim), U0(Vector::Zero(3)), U1(Vector::Zero(3)), XG_temp(3);
  double   x_error=0;
  VectorXi dofs(dim);

  if (print_to_matlab)
    printMatlabLoader();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Solve nonlinear system
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  printf("initial volume: %.15lf \n", initial_volume);
  current_time = 0;
  time_step = 0;
  double Qmax = 0;
  double steady_error = 1;
  setInitialConditions();  //called only here
  int its;

  printf("\nNum. of time iterations (maxts): %d\n",maxts);
  printf("Starting time loop ... \n");

  if (is_bdf2)
  {
    current_time += dt;
    time_step += 1;

    VecCopy(Vec_uzp_0,Vec_uzp_m1);
    copyVec2Mesh(Vec_x_1);

    if (is_bdf2_bdfe)
    {
      // extrapolation \tilde{X}^{n+1}=2X^{n}-X^{n-1}
      VecScale(Vec_x_1, 2.0);
      VecAXPY(Vec_x_1,-1.0,Vec_x_0);
      // calc V^{n+1} and update with D_{2}X^{n+1} = dt*V^{n+1}
      calcMeshVelocity(Vec_x_0, Vec_uzp_0, Vec_uzp_1, 2.0, Vec_v_1, current_time);
      VecCopy(Vec_v_1, Vec_x_1);
      VecScale(Vec_x_1, 2./3.*dt);
      VecAXPY(Vec_x_1,-1./3.,Vec_x_0);
      copyMesh2Vec(Vec_x_0);
      VecAXPY(Vec_x_1,4./3.,Vec_x_0);
      if (N_Solids){moveCenterMass(2.0);}
      VecCopy(Vec_v_1,Vec_v_mid);
    }
    else if (is_bdf2_ab)
    {
      // extrapolated geometry
      //VecScale(Vec_x_1, 2.0);
      //VecAXPY(Vec_x_1,-1.0,Vec_x_0);  // \bar{X}^(n+1/2)=2.0*X^(n)-1.0X^(n-1)
      VecScale(Vec_x_1, 1.5);
      VecAXPY(Vec_x_1,-0.5,Vec_x_0);  // \bar{X}^(n+1/2)=1.5*X^(n)-0.5X^(n-1)
      copyMesh2Vec(Vec_x_0);          //copy current mesh to Vec_x_0
      calcMeshVelocity(Vec_x_0, Vec_uzp_0, Vec_uzp_1, 1.5, Vec_v_1, current_time); // Adams-Bashforth
      VecWAXPY(Vec_x_1, dt, Vec_v_1, Vec_x_0); // Vec_x_1 = Vec_v_1*dt + Vec_x_0
      if (N_Solids){moveCenterMass(1.5);}
      // velocity at integer step
      //VecScale(Vec_v_1, 1.5);
      //VecAXPY(Vec_v_1,-.5,Vec_v_mid);
      VecScale(Vec_v_1, 2.0);
      VecAXPY(Vec_v_1,-1.0,Vec_v_mid);
    }

    VecCopy(Vec_uzp_1, Vec_uzp_0);
  }

  if (is_mr_ab && false){
    current_time += dt;
    time_step += 1;

    // extrapolated geometry
    VecScale(Vec_x_1, 2.0);         // look inside residue why this extrap is used
    VecAXPY(Vec_x_1,-1.0,Vec_x_0);  // \bar{X}^(n+1/2)=2.0*X^(n)-1.0X^(n-1)
    //VecScale(Vec_x_1, 1.5);
    //VecAXPY(Vec_x_1,-0.5,Vec_x_0);  // \bar{X}^(n+1/2)=1.5*X^(n)-0.5X^(n-1)
    copyMesh2Vec(Vec_x_0);          //copy current mesh to Vec_x_0
    calcMeshVelocity(Vec_x_0, Vec_uzp_0, Vec_uzp_1, 1.5, Vec_v_mid, current_time); // Adams-Bashforth
    VecWAXPY(Vec_x_1, dt, Vec_v_mid, Vec_x_0); // Vec_x_1 = Vec_v_mid*dt + Vec_x_0
    if (N_Solids){
      moveCenterMass(1.5);
      VecCopy(Vec_slipv_1,Vec_slipv_0);
      updateSolidMesh();
    }
  }

  //print solid's center information
  ofstream filg, filv;
  Vector3d Xgg;
  VectorXd v_coeffs_s(LZ*N_Solids);
  VectorXi mapvs(LZ*N_Solids);
  //int TT = 0;
  char gravc[PETSC_MAX_PATH_LEN], velc[PETSC_MAX_PATH_LEN];
  sprintf(gravc,"%s/HistGrav.txt",filehist_out.c_str());
  sprintf(velc,"%s/HistVel.txt",filehist_out.c_str());
  if (fprint_hgv){
    filg.open(gravc);
    filv.open(velc);
    filg.close();
    filv.close();
  }

  for(;;)  // equivalent to forever or while(true), must be a break inside
  {

    cout << endl;
    cout << "current time: " << current_time << endl;
    cout << "time step: "    << time_step  << endl;
    cout << "steady error: " << steady_error << endl;

    if (maxts == 0)
    {
      computeError(Vec_x_0, Vec_uzp_0,current_time);
      break;
    }

    bool const full_implicit = false;
    //bool const try2 = false; // extrapolate geometry Vec_x_1 <- 2*Vec_x_1 - Vec_x_0

    // * SOLVE THE SYSTEM *
    if (solve_the_sys)
    {

      setUPInitialGuess();  //setup Vec_uzp_1 for SNESSolve

      if (fprint_hgv){
        if ((time_step%print_step)==0 || time_step == (maxts-1)){
          getSolidVolume();
          for (int S = 0; S < LZ*N_Solids; S++){
            mapvs[S] = n_unknowns_u + n_unknowns_p + S;
          }
          VecGetValues(Vec_uzp_0,mapvs.size(),mapvs.data(),v_coeffs_s.data());
          filg.open(gravc,iostream::app);
          filv.open(velc,iostream::app);
          filg << current_time << " ";
          filv << current_time << " ";
          for (int S = 0; S < (N_Solids-1); S++){
            Xgg = XG_0[S];
            filg << Xgg(0) << " " << Xgg(1) << " "; if (dim == 3){filg << Xgg(2) << " ";} filg << " " << VV[S] << " " << theta_0[S];
            filv << v_coeffs_s(LZ*S) << " " << v_coeffs_s(LZ*S+1) << " " << v_coeffs_s(LZ*S+2) << " ";
          }
          Xgg = XG_0[N_Solids-1];
          filg << Xgg(0) << " " << Xgg(1); if (dim == 3){filg << " "<< Xgg(2);} filg << " " << VV[N_Solids-1] << " " << theta_0[N_Solids-1] << endl;
          filv << v_coeffs_s(LZ*(N_Solids-1)) << " " << v_coeffs_s(LZ*(N_Solids-1)+1) << " "
               << v_coeffs_s(LZ*(N_Solids-1)+2);
          if (dim == 3){
            filv << v_coeffs_s(LZ*(N_Solids-1)+3) << " " << v_coeffs_s(LZ*(N_Solids-1)+4) << " "
                << v_coeffs_s(LZ*(N_Solids-1)+5);
          }
          filv << endl;
          filg.close();  //TT++;
          filv.close();
        }
      }
      //char buf1[18], buf2[10];
      //sprintf(buf1,"matrizes/sol%d.m",time_step); sprintf(buf2,"solm%d",time_step);
      //View(Vec_uzp_0, buf1, buf2);
      //sprintf(buf1,"matrizes/mal%d.m",time_step); sprintf(buf2,"malm%d",time_step);
      //View(Vec_x_0, buf1, buf2);  //TT++;

      for (int kk = 0 ; kk < 1+3*full_implicit; kk++)
      {
        printf("\tFixed Point Iteration %d\n", kk);

        ierr = SNESSolve(snes_fs,PETSC_NULL,Vec_uzp_1);        CHKERRQ(ierr);
        if (N_Solids){updateSolidVel();}
      }

    } // end if (solve_the_system)

    if (is_bdf2)
    {
      if (plot_exact_sol)
        computeError(Vec_x_0, Vec_uzp_0,current_time);

      VecCopy(Vec_uzp_0,Vec_uzp_m1);
      if (is_bdf2_ab)
        VecCopy(Vec_v_1,Vec_v_mid);
    }
    else if (is_bdf3)
    {
      if (plot_exact_sol)
        computeError(Vec_x_0, Vec_uzp_0,current_time);
      ////-----
      VecCopy(Vec_uzp_m1,Vec_uzp_m2);
      VecCopy(Vec_uzp_0,Vec_uzp_m1);
      /*///-----
      VecCopy(Vec_duzp, Vec_duzp_0);
      VecCopy(Vec_uzp_1, Vec_duzp);
      VecAXPY(Vec_duzp,-1.0,Vec_uzp_0); // Vec_dup -= Vec_up_0
      VecScale(Vec_duzp, 1./dt);
      /*///-----
    }
    else if (is_mr_ab) //for MR-AB, modifies P from Vec_uzp_0
    {
      if (time_step == 0)
      {
        pressureTimeCorrection(Vec_uzp_0, Vec_uzp_1, 0., 1); // press(n) = press(n+1/2) - press(n-1/2)
        if (plot_exact_sol && maxts <= 1)
          computeError(Vec_x_0, Vec_uzp_0,current_time);
      }
      else
      {
        pressureTimeCorrection(Vec_uzp_0, Vec_uzp_1, .5, .5); // press(n) = press(n+1/2) - press(n-1/2)
        if (plot_exact_sol)
          computeError(Vec_x_0, Vec_uzp_0,current_time);//computeError(Vec_x_0, Vec_up_0,current_time);
      }
    }

    bool must_print = false;
    if (is_bdf2)
    {
      if (((time_step-1)%print_step)==0 || time_step == (maxts-1))
        must_print = true;
    }
    else
    if (is_bdf3)
    {
      if (((time_step-2)%print_step)==0 || time_step == (maxts-1))
        must_print = true;
    }
    else
      if ((time_step%print_step)==0 || time_step == (maxts-1))
        must_print = true;

    if (must_print)  // print files
    {
      if (family_files)
      {
        plotFiles();
        ierr = SNESGetIterationNumber(snes_fs,&its);     CHKERRQ(ierr);
        cout << "num snes iterations: " << its << endl;
      }

    }

    //printContactAngle(fprint_ca);
    //////if (plot_exact_sol) eh melhor depois da atualização
    //////  computeError(Vec_x_1, Vec_up_1,current_time+dt);
    current_time += dt;
    time_step += 1;

    // update
    if (ale)
    {
      // mesh adaption (it has topological changes) (2d only) it destroys Vec_normal and Vec_v_mid
      if (mesh_adapt)
        meshAdapt_s(); //meshAdapt()_l;

      copyVec2Mesh(Vec_x_1);

      // MESH FLIPPING
      if (mesh_adapt)
      {
        if(time_step%1 == 0)
        {
          if (dim==2)
          {
            Facet *f;
            Cell* ca, *cb;
            // Delaunay
            for (int i = 0; i < n_facets_total; ++i)
            {
              f = mesh->getFacetPtr(i);
              if (f->isDisabled())
                continue;
              if (!MeshToolsTri::inCircle2d(f, &*mesh))
              {
                ca = mesh->getCellPtr(f->getIncidCell());
                cb = mesh->getCellPtr(ca->getIncidCell(f->getPosition()));

                if (cb)
                  if (ca->getTag() == cb->getTag())
                    MeshToolsTri::flipEdge(f, &*mesh, true);
              }
            }
          } // end dim==2
        } // end ts%1
      } // end mesh adapt

      if (is_bdf2)
      {
        if (is_bdf2_bdfe)
        {
          // extrapolation \tilde{X}^{n+1}=2X^{n}-X^{n-1}
          VecScale(Vec_x_1, 2.0);
          VecAXPY(Vec_x_1,-1.0,Vec_x_0);
          // calc V^{n+1} and update with D_{2}X^{n+1} = dt*V^{n+1}
          calcMeshVelocity(Vec_x_0, Vec_uzp_0, Vec_uzp_1, 2.0, Vec_v_1, current_time);
          VecCopy(Vec_v_1, Vec_x_1);
          VecScale(Vec_x_1, 2./3.*dt);
          VecAXPY(Vec_x_1,-1./3.,Vec_x_0);
          copyMesh2Vec(Vec_x_0);
          VecAXPY(Vec_x_1,4./3.,Vec_x_0);
          if (N_Solids){moveCenterMass(2.0);}
          VecCopy(Vec_v_1,Vec_v_mid);
        }
        else if (is_bdf2_ab)
        {
          // extrapolated geometry
          //VecScale(Vec_x_1, 2.0);
          //VecAXPY(Vec_x_1,-1.0,Vec_x_0);  // \bar{X}^(n+1/2)=2.0*X^(n)-1.0X^(n-1)
          VecScale(Vec_x_1, 1.5);
          VecAXPY(Vec_x_1,-0.5,Vec_x_0);  // \bar{X}^(n+1/2)=1.5*X^(n)-0.5X^(n-1)
          copyMesh2Vec(Vec_x_0);          //copy current mesh to Vec_x_0
          calcMeshVelocity(Vec_x_0, Vec_uzp_0, Vec_uzp_1, 1.5, Vec_v_1, current_time); // Adams-Bashforth
          VecWAXPY(Vec_x_1, dt, Vec_v_1, Vec_x_0); // Vec_x_1 = Vec_v_1*dt + Vec_x_0
          if (N_Solids){moveCenterMass(1.5);}
          // velocity at integer step
          //VecScale(Vec_v_1, 1.5);
          //VecAXPY(Vec_v_1,-.5,Vec_v_mid);
          VecScale(Vec_v_1, 2.0);
          VecAXPY(Vec_v_1,-1.0,Vec_v_mid);

        }
      }
      else if (is_bdf3)
      {
        // extrapolation of Vec_x_1, \bar{X^{n+1}} = 3X^{n}-3X^{n-1}+\bar{X^{n-2}}
        VecScale(Vec_x_1, 3.0);
        VecAXPY(Vec_x_1,-3.0,Vec_x_0);
        VecAXPY(Vec_x_1, 1.0,Vec_x_aux);
        // copyMesh2Vec(Vec_x_0); // we need Vec_x_0 later, don't touch it
        // estraga Vec_up_0, usado porque não se precisa mais
        ////-----
        VecScale(Vec_uzp_0,-3.0);
        VecAXPY(Vec_uzp_0,3.0,Vec_uzp_1);
        VecAXPY(Vec_uzp_0,1.0,Vec_uzp_m2);
        /*///-----
        VecCopy(Vec_duzp, Vec_uzp_0);
        VecScale(Vec_uzp_0, 2.0);
        VecAXPY(Vec_uzp_0,-1.0,Vec_duzp_0);
        VecScale(Vec_uzp_0, dt);
        VecAXPY(Vec_uzp_0, 1.0,Vec_uzp_1); // Vec_up_0 tem agora a extrapolacao do futuro Vec_up_1
        /*///-----

        // calc V^{n+1} and update with D_{3}X^{n+1} = dt*V^{n+1}
        calcMeshVelocity(Vec_x_0, Vec_uzp_0, Vec_uzp_1, 0.0, Vec_v_1, current_time);  //comentado?
        VecCopy(Vec_v_1, Vec_x_1);
        VecScale(Vec_x_1, 6./11.*dt);
        VecAXPY(Vec_x_1, 2./11.,Vec_x_aux);
        VecCopy(Vec_x_0, Vec_x_aux);
        VecAXPY(Vec_x_1,-9./11.,Vec_x_0);
        copyMesh2Vec(Vec_x_0);
        VecAXPY(Vec_x_1,18./11.,Vec_x_0);
        if (N_Solids){moveCenterMass(0.0);}// update center of mass D_{3}XG^{n+1} = dt*V^{n+1}
        // VecCopy(Vec_uzp_1, Vec_uzp_0) is done after the end if ale
        VecCopy(Vec_v_1,Vec_v_mid);
      }
      else if (is_mr_ab) //for MR-AB
      {
        // extrapolated geometry
        VecScale(Vec_x_1, 2.0);         // look inside residue why this extrap is used
        VecAXPY(Vec_x_1,-1.0,Vec_x_0);  // \bar{X}^(n+1/2)=2.0*X^(n)-1.0X^(n-1)
        //VecScale(Vec_x_1, 1.5);
        //VecAXPY(Vec_x_1,-0.5,Vec_x_0);  // \bar{X}^(n+1/2)=1.5*X^(n)-0.5X^(n-1)
        copyMesh2Vec(Vec_x_0);          //copy current mesh to Vec_x_0
        velNoSlip(Vec_uzp_0,Vec_slipv_0,Vec_uzp_0_ns);
        velNoSlip(Vec_uzp_1,Vec_slipv_1,Vec_uzp_1_ns);
        calcMeshVelocity(Vec_x_0, Vec_uzp_0, Vec_uzp_1, 1.5, Vec_v_mid, current_time); // Adams-Bashforth
        VecWAXPY(Vec_x_1, dt, Vec_v_mid, Vec_x_0); // Vec_x_1 = Vec_v_mid*dt + Vec_x_0
        if (N_Solids){
          moveCenterMass(1.5);
          if (is_slipv) VecCopy(Vec_slipv_1,Vec_slipv_0);
          updateSolidMesh();
        }
      }
      else //for basic
      {
        //VecCopy(Vec_x_0,Vec_x_1);
        copyMesh2Vec(Vec_x_0);          //copy current mesh to Vec_x_0
        calcMeshVelocity(Vec_x_0, Vec_uzp_0, Vec_uzp_1, 1.0, Vec_v_mid, current_time); // Adams-Bashforth
        VecWAXPY(Vec_x_1, dt, Vec_v_mid, Vec_x_0); // Vec_x_1 = Vec_v_mid*dt + Vec_x_0
        if (N_Solids){moveCenterMass(1.0);}
      }

/*
      if ( !full_implicit && !is_bdf3 )
      {
        {
          calcMeshVelocity(Vec_x_0, Vec_uzp_0, Vec_uzp_1, 1.5, Vec_v_mid, current_time); // Adams-Bashforth
        }
        // move the mesh
        if (!try2 && !is_bdf_bdf_extrap && !is_bdf_ab){
          VecWAXPY(Vec_x_1, dt, Vec_v_mid, Vec_x_0); // Vec_x_1 = Vec_v_mid*dt + Vec_x_0
          moveCenterMass(1.5);
        }
      }//end !full_implicit && !is_bdf3
*/
      //compute normal for the next time step, at n+1/2
      {
        Vec Vec_x_mid;
        int xsize;
        double *xarray;
        VecGetSize(Vec_x_0, &xsize);
        VecGetArray(Vec_res_fs, &xarray); //VecGetArray(Vec_res, &xarray);
        //prototipo no petsc-dev: VecCreateSeqWithArray(MPI_Comm comm,PetscInt bs,PetscInt n,const PetscScalar array[],Vec *V)
        //VecCreateSeqWithArray(MPI_COMM_SELF, xsize, xarray, &Vec_x_mid);
        VecCreateSeqWithArray(MPI_COMM_SELF, 1,xsize, xarray, &Vec_x_mid);
        Assembly(Vec_x_mid);

        VecCopy(Vec_x_0, Vec_x_mid);
        VecAXPY(Vec_x_mid,1.,Vec_x_1);
        VecScale(Vec_x_mid, 0.5);
        getVecNormals(&Vec_x_mid, Vec_normal);

        VecDestroy(&Vec_x_mid);
        VecRestoreArray(Vec_res_fs, &xarray); //VecRestoreArray(Vec_res, &xarray);
      }

      VecCopy(Vec_uzp_1, Vec_uzp_0); //VecCopy(Vec_up_1, Vec_up_0);

    }//end if ale
    else // if it's not ALE
    {
      VecCopy(Vec_uzp_1, Vec_uzp_0); //VecCopy(Vec_up_1, Vec_up_0);
    }


    // compute steady error
    VecNorm(Vec_uzp_1, NORM_1, &Qmax); //VecNorm(Vec_up_1, NORM_1, &Qmax);
    VecCopy(Vec_uzp_0, Vec_res_fs); //VecCopy(Vec_up_0, Vec_res);
    VecAXPY(Vec_res_fs,-1.0,Vec_uzp_1); //VecAXPY(Vec_res,-1.0,Vec_up_1);
    //steady_error = VecNorm(Vec_res, NORM_1)/(Qmax==0.?1.:Qmax);
    VecNorm(Vec_res_fs, NORM_1, &steady_error); //VecNorm(Vec_res, NORM_1, &steady_error);
    steady_error /= (Qmax==0.?1.:Qmax);


    if(time_step >= maxts) {
      cout << "\n==========================================\n";
      cout << "stop reason:\n";
      cout << "maximum number of iterations reached. \n";
      break;
    }
    if (steady_error <= steady_tol) {
      cout << "\n==========================================\n";
      cout << "stop reason:\n";
      cout << "steady state reached. \n";
      break;
    }
  }
  //
  // END TIME LOOP
  //

  // PRINT AGAIN
  if (false)
  {
    if (family_files)
    {
      double  *q_array;
      double  *nml_array;
      double  *v_array;
      VecGetArray(Vec_uzp_0, &q_array);
      VecGetArray(Vec_normal, &nml_array);
      VecGetArray(Vec_v_mid, &v_array);
      vtk_printer.writeVtk();

      /* ---- nodes data ---- */
      vtk_printer.addNodeVectorVtk("u", GetDataVelocity(q_array, *this));
      //vtk_printer.addNodeVectorVtk("normal",  GetDataNormal(nml_array, *this));
      //vtk_printer.addNodeVectorVtk("v",  GetDataMeshVel(v_array, *this));
      vtk_printer.printPointTagVtk();

      if (!shape_psi_c->discontinuous())
        vtk_printer.addNodeScalarVtk("pressure", GetDataPressure(q_array, *this));
      else
        vtk_printer.addCellScalarVtk("pressure", GetDataPressCellVersion(q_array, *this));

      vtk_printer.addCellIntVtk("cell_tag", GetDataCellTag(*this));
      //vtk_printer.printPointTagVtk("point_tag");

      VecRestoreArray(Vec_uzp_0, &q_array);  //VecRestoreArray(Vec_up_0, &q_array);
      VecRestoreArray(Vec_normal, &nml_array);
      VecRestoreArray(Vec_v_mid, &v_array);

      ierr = SNESGetIterationNumber(snes_fs,&its);     CHKERRQ(ierr);
      cout << "num snes iterations: " << its << endl;
    }

  }//end if false


  cout << endl;

  final_volume = getMeshVolume();
  printf("final volume: %.15lf \n", final_volume);
  printf("volume error 100*(f-i)/i: %.15lf per percent\n", 100*abs(final_volume-initial_volume)/initial_volume);
  printf("x error : %.15lf \n", x_error);


  SNESConvergedReason reason;
  ierr = SNESGetIterationNumber(snes_fs,&its);     CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(snes_fs,&reason);  CHKERRQ(ierr);
  cout << "num snes iterations: " << its << endl;
  cout << "reason: " << reason << endl;
  switch (reason)
  {
  case( SNES_CONVERGED_FNORM_ABS       ): printf("SNES_CONVERGED_FNORM_ABS      /* ||F|| < atol */\n"); break;
  case( SNES_CONVERGED_FNORM_RELATIVE  ): printf("SNES_CONVERGED_FNORM_RELATIVE /* ||F|| < rtol*||F_initial|| */\n"); break;
  case( SNES_CONVERGED_SNORM_RELATIVE  ): printf("SNES_CONVERGED_PNORM_RELATIVE /* Newton computed step size small); break; || delta x || < stol */\n"); break;
  case( SNES_CONVERGED_ITS             ): printf("SNES_CONVERGED_ITS            /* maximum iterations reached */\n"); break;
  case( SNES_CONVERGED_TR_DELTA        ): printf("SNES_CONVERGED_TR_DELTA       \n"); break;
        /* diverged */
  case( SNES_DIVERGED_FUNCTION_DOMAIN  ): printf("SNES_DIVERGED_FUNCTION_DOMAIN /* the new x location passed the function is not in the domain of F */\n"); break;
  case( SNES_DIVERGED_FUNCTION_COUNT   ): printf("SNES_DIVERGED_FUNCTION_COUNT   \n"); break;
  case( SNES_DIVERGED_LINEAR_SOLVE     ): printf("SNES_DIVERGED_LINEAR_SOLVE    /* the linear solve failed */\n"); break;
  case( SNES_DIVERGED_FNORM_NAN        ): printf("SNES_DIVERGED_FNORM_NAN       \n"); break;
  case( SNES_DIVERGED_MAX_IT           ): printf("SNES_DIVERGED_MAX_IT          \n"); break;
  case( SNES_DIVERGED_LINE_SEARCH      ): printf("SNES_DIVERGED_LINE_SEARCH     /* the line search failed */ \n"); break;
  case( SNES_DIVERGED_LOCAL_MIN        ): printf("SNES_DIVERGED_LOCAL_MIN       /* || J^T b || is small, implies converged to local minimum of F() */\n"); break;
  case( SNES_CONVERGED_ITERATING       ): printf("SNES_CONVERGED_ITERATING      \n"); break;
  case( SNES_DIVERGED_INNER            ): printf("SNES_DIVERGED_INNER           \n"); break;
  }

  if (solve_the_sys)
    MatrixInfo(Mat_Jac_fs); //MatrixInfo(Mat_Jac);

  int lits;
  SNESGetLinearSolveIterations(snes_fs,&lits); //SNESGetLinearSolveIterations(snes,&lits);

  cout << "Greatest error reached during the simulation:" << endl;
  printf("%-21s %-21s %-21s %-21s %-21s %-21s %-21s %-21s %s\n", "# hmean", "u_L2_norm", "p_L2_norm", "grad_u_L2_norm", "grad_p_L2_norm", "u_L2_facet_norm", "u_inf_facet_norm", "u_inf_norm", "p_inf_norm" );
  printf("%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n\n",Stats.hmean, Stats.u_L2_norm, Stats.p_L2_norm, Stats.grad_u_L2_norm, Stats.grad_p_L2_norm, Stats.u_L2_facet_norm,  Stats.u_inf_facet_norm, Stats.u_inf_norm, Stats.p_inf_norm);

  ////if (unsteady)
  //{
  //  cout << "\nmean errors: \n";
  //  cout << "# hmean            u_L2_err         p_L2_err          grad_u_L2_err     grad_p_L2_err" << endl;
  //  printf( "%.15lf %.15lf %.15lf %.15lf %.15lf\n", Stats.mean_hmean() , Stats.mean_u_L2_norm(), Stats.mean_p_L2_norm(), Stats.mean_grad_u_L2_norm(), Stats.mean_grad_p_L2_norm());
  //}


  PetscFunctionReturn(0);
}

void AppCtx::computeError(Vec const& Vec_x, Vec &Vec_up, double tt)
{

  MatrixXd            u_coefs_c(n_dofs_u_per_cell/dim, dim);
  MatrixXd            u_coefs_c_trans(dim,n_dofs_u_per_cell/dim);
  MatrixXd            u_coefs_f(n_dofs_u_per_facet/dim, dim);
  MatrixXd            u_coefs_f_trans(dim,n_dofs_u_per_facet/dim);
  VectorXd            p_coefs_c(n_dofs_p_per_cell);
  MatrixXd            x_coefs_c(nodes_per_cell, dim);
  MatrixXd            x_coefs_c_trans(dim, nodes_per_cell);
  MatrixXd            x_coefs_f(nodes_per_facet, dim);
  MatrixXd            x_coefs_f_trans(dim, nodes_per_facet);
  Tensor              F_c(dim,dim), invF_c(dim,dim), invFT_c(dim,dim);
  Tensor              F_f(dim,dim-1), invF_f(dim-1,dim), fff(dim-1,dim-1);
  MatrixXd            dxphi_err(n_dofs_u_per_cell/dim, dim);
  MatrixXd            dxpsi_err(n_dofs_p_per_cell, dim);
  MatrixXd            dxqsi_err(nodes_per_cell, dim);
  Tensor              dxU(dim,dim); // grad u
  Vector              dxP(dim);     // grad p
  Vector              Xqp(dim);
  Vector              Uqp(dim);

  double              Pqp;
  VectorXi            cell_nodes(nodes_per_cell);
  double              J, JxW;
  double              weight;
  int                 tag;
  double              volume=0;

  double              p_L2_norm = 0.;
  double              u_L2_norm = 0.;
  double              grad_u_L2_norm = 0.;
  double              grad_p_L2_norm = 0.;
  double              p_inf_norm = 0.;
  double              u_inf_norm = 0.;
  double              hmean = 0.;
  double              u_L2_facet_norm = 0.;
  double              u_inf_facet_norm = 0.;

  VectorXi            mapU_c(n_dofs_u_per_cell);
  VectorXi            mapU_f(n_dofs_u_per_facet);
  VectorXi            mapU_r(n_dofs_u_per_corner);
  VectorXi            mapP_c(n_dofs_p_per_cell);
  VectorXi            mapP_r(n_dofs_p_per_corner);
  VectorXi            mapM_c(dim*nodes_per_cell);
  VectorXi            mapM_f(dim*nodes_per_facet);
  VectorXi            mapM_r(dim*nodes_per_corner);

  VecSetOption(Vec_up, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); //?

  cell_iterator cell = mesh->cellBegin();
  cell_iterator cell_end = mesh->cellEnd();
  for (; cell != cell_end; ++cell)
  {
    tag = cell->getTag();

    // mapeamento do local para o global:
    dof_handler[DH_MESH].getVariable(VAR_M).getCellDofs(mapM_c.data(), &*cell);
    dof_handler[DH_UNKM].getVariable(VAR_U).getCellDofs(mapU_c.data(), &*cell);
    dof_handler[DH_UNKM].getVariable(VAR_P).getCellDofs(mapP_c.data(), &*cell);

    /*  Pega os valores das variáveis nos graus de liberdade */
    VecGetValues(Vec_up, mapU_c.size(), mapU_c.data(), u_coefs_c.data());
    VecGetValues(Vec_up, mapP_c.size(), mapP_c.data(), p_coefs_c.data());
    VecGetValues(Vec_x,  mapM_c.size(), mapM_c.data(), x_coefs_c.data());

    u_coefs_c_trans = u_coefs_c.transpose();
    x_coefs_c_trans = x_coefs_c.transpose();

    for (int qp = 0; qp < n_qpts_err; ++qp)
    {
      F_c    = x_coefs_c_trans * dLqsi_err[qp];
      J      = F_c.determinant();
      invF_c = F_c.inverse();
      invFT_c= invF_c.transpose();

      dxphi_err = dLphi_err[qp] * invF_c;
      dxpsi_err = dLpsi_err[qp] * invF_c;
      dxqsi_err = dLqsi_err[qp] * invF_c;

      dxP  = dxpsi_err.transpose() * p_coefs_c;
      dxU  = u_coefs_c_trans * dxphi_err;       // n+utheta

      Xqp  = x_coefs_c_trans * qsi_err[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
      Uqp  = u_coefs_c_trans * phi_err[qp];
      Pqp  = p_coefs_c.dot(psi_err[qp]);

      weight = quadr_err->weight(qp);

      JxW = J*weight;

      //  note: norm(u, H1)^2 = norm(u, L2)^2 + norm(gradLphi_c, L2)^2
      u_L2_norm        += (u_exact(Xqp, tt, tag) - Uqp).squaredNorm()*JxW;
      p_L2_norm        += sqr(p_exact(Xqp, tt, tag) - Pqp)*JxW;
      grad_u_L2_norm   += (grad_u_exact(Xqp, tt, tag) - dxU).squaredNorm()*JxW;
      grad_p_L2_norm   += (grad_p_exact(Xqp, tt, tag) - dxP).squaredNorm()*JxW;
      u_inf_norm       = max(u_inf_norm, (u_exact(Xqp, tt, tag) - Uqp).norm());
      p_inf_norm       = max(p_inf_norm, fabs(p_exact(Xqp, tt, tag) - Pqp));

      volume += JxW;
    } // fim quadratura

  } // end elementos
#if (true)
  facet_iterator facet = mesh->facetBegin();
  facet_iterator facet_end = mesh->facetEnd();
  for (; facet != facet_end; ++facet)
  {
    tag = facet->getTag();

    if (!(is_in(tag, dirichlet_tags) || is_in(tag, neumann_tags) || is_in(tag, interface_tags) ||
          is_in(tag, solid_tags) || is_in(tag, periodic_tags) || is_in(tag, feature_tags)))
      continue;

    // mapeamento do local para o global:
    //
    dof_handler[DH_UNKM].getVariable(VAR_U).getFacetDofs(mapU_f.data(), &*facet);
    dof_handler[DH_MESH].getVariable(VAR_M).getFacetDofs(mapM_f.data(), &*facet);

    /*  Pega os valores das variáveis nos graus de liberdade */
    VecGetValues(Vec_up, mapU_f.size(), mapU_f.data(), u_coefs_f.data());
    VecGetValues(Vec_x,  mapM_f.size(), mapM_f.data(), x_coefs_f.data());

    u_coefs_f_trans = u_coefs_f.transpose();
    x_coefs_f_trans = x_coefs_f.transpose();

    for (int qp = 0; qp < n_qpts_facet; ++qp)
    {
      F_f   = x_coefs_f_trans * dLqsi_f[qp];

      fff    = F_f.transpose()*F_f;
      J      = sqrt(fff.determinant());

      Xqp  = x_coefs_f_trans * qsi_f[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
      Uqp  = u_coefs_f_trans * phi_f[qp];

      weight = quadr_facet->weight(qp);

      JxW = J*weight;

      Vector U_exact = u_exact(Xqp, tt, tag);

      u_L2_facet_norm += (U_exact - Uqp).squaredNorm()*JxW;
      //u_L2_facet_norm += JxW;

      double const diff = (U_exact - Uqp).norm();

      if (diff > u_inf_facet_norm)
        u_inf_facet_norm = diff;


    } // fim quadratura

  }
#endif

  u_L2_norm      = sqrt(u_L2_norm     );
  p_L2_norm      = sqrt(p_L2_norm     );
  grad_u_L2_norm = sqrt(grad_u_L2_norm);
  grad_p_L2_norm = sqrt(grad_p_L2_norm);
  u_L2_facet_norm = sqrt(u_L2_facet_norm);

  // ASSUME QUE SÓ POSSA TER NO MÁXIMO 1 NÓ POR ARESTA


  VectorXi edge_nodes(3);
  Vector Xa(dim), Xb(dim);
  int n_edges=0;

  if (dim==2)
  //FEP_PRAGMA_OMP(parallel default(none) shared(hmean))
  {
    const int n_edges_total = mesh->numFacetsTotal();
    Facet const* edge(NULL);

    //FEP_PRAGMA_OMP(for nowait)
    for (int a = 0; a < n_edges_total; ++a)
    {
      edge = mesh->getFacetPtr(a);
      if (edge->isDisabled())
        continue;

      mesh->getFacetNodesId(&*edge, edge_nodes.data());

      mesh->getNodePtr(edge_nodes[0])->getCoord(Xa.data(),dim);
      mesh->getNodePtr(edge_nodes[1])->getCoord(Xb.data(),dim);
      hmean += (Xa-Xb).norm();
      ++n_edges;
    }
  }
  else
  if (dim==3)
  //FEP_PRAGMA_OMP(parallel default(none) shared(cout,hmean))
  {
    const int n_edges_total = mesh->numCornersTotal();
    Corner const* edge(NULL);
//cout << "aqui" << endl;
    //FEP_PRAGMA_OMP(for nowait)
    for (int a = 0; a < n_edges_total; ++a)
    {
      edge = mesh->getCornerPtr(a);
      if (edge->isDisabled())
        continue;

      mesh->getCornerNodesId(&*edge, edge_nodes.data());

      mesh->getNodePtr(edge_nodes[0])->getCoord(Xa.data(),dim);
      mesh->getNodePtr(edge_nodes[1])->getCoord(Xb.data(),dim);
      hmean += (Xa-Xb).norm();
      ++n_edges;
    }
  }
  hmean /= n_edges;
  //hme = hmean;
  //if (time_step==maxts)
  {
    cout << endl;
    //cout << "errors computed at last time step: "<< endl;
    //cout << "# hmean               u_L2_norm         p_L2_norm         grad_u_L2_norm    grad_p_L2_norm   u_L2_facet_norm    u_inf_facet_norm" << endl;
    printf("%-21s %-21s %-21s %-21s %-21s %-21s %s\n", "# hmean", "u_L2_norm", "p_L2_norm", "grad_u_L2_norm", "grad_p_L2_norm", "u_L2_facet_norm", "u_inf_facet_norm" );
    printf("%.15e %.15e %.15e %.15e %.15e %.15e %.15e\n\n", hmean, u_L2_norm, p_L2_norm, grad_u_L2_norm, grad_p_L2_norm, u_L2_facet_norm, u_inf_facet_norm);
  }

  Stats.add_p_L2_norm        (p_L2_norm       );
  Stats.add_u_L2_norm        (u_L2_norm       );
  Stats.add_grad_u_L2_norm   (grad_u_L2_norm  );
  Stats.add_grad_p_L2_norm   (grad_p_L2_norm  );
  Stats.add_p_inf_norm       (p_inf_norm      );
  Stats.add_u_inf_norm       (u_inf_norm      );
  Stats.add_hmean            (hmean           );
  Stats.add_u_L2_facet_norm  (u_L2_facet_norm );
  Stats.add_u_inf_facet_norm (u_inf_facet_norm);


}


void AppCtx::pressureTimeCorrection(Vec &Vec_up_0, Vec &Vec_up_1, double a, double b) // p(n+1) = a*p(n+.5) + b* p(n)
{
  Vector    Uf(dim);
  Vector    X(dim);
  Tensor    R(dim,dim);
  int       dof;

  double P0, P1, P2;

  if (behaviors & BH_Press_grad_elim)
  {
    cout << "FIX ME: NOT SUPPORTED YET!!!!\n";
    throw;
  }

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    // press
    if (mesh->isVertex(&*point))
    {
      getNodeDofs(&*point,DH_UNKM,VAR_P,&dof);
      VecGetValues(Vec_up_1, 1, &dof, &P1);
      VecGetValues(Vec_up_0, 1, &dof, &P0);
      P2 = a*P1 + b*P0;
      VecSetValues(Vec_up_0, 1, &dof, &P2, INSERT_VALUES);
    }

  } // end point loop
  Assembly(Vec_up_0);
}

double AppCtx::getMaxVelocity()
{
  Vector     Uf(dim); // Ue := elastic velocity
  VectorXi   vtx_dofs_fluid(dim); // indices de onde pegar a velocidade
  double     Umax=0;

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    if (!mesh->isVertex(&*point))
        continue;
    dof_handler[DH_UNKM].getVariable(VAR_U).getVertexDofs(vtx_dofs_fluid.data(), &*point);
    //VecGetValues(Vec_up_1, dim, vtx_dofs_fluid.data(), Uf.data());
    Umax = max(Umax,Uf.norm());
  }
  return Umax;
}

double AppCtx::getMeshVolume()
{
  MatrixXd            x_coefs_c(nodes_per_cell, dim);
  MatrixXd            x_coefs_c_trans(dim, nodes_per_cell);
  Tensor              F_c(dim,dim), invF_c(dim,dim), invFT_c(dim,dim);
  VectorXi            cell_nodes(nodes_per_cell);
  double              Jx;
  //int                 tag;
  double              volume=0;

  cell_iterator cell = mesh->cellBegin();
  cell_iterator cell_end = mesh->cellEnd();
  for (; cell != cell_end; ++cell)
  {
    //tag = cell->getTag();

    mesh->getCellNodesId(&*cell, cell_nodes.data());
    mesh->getNodesCoords(cell_nodes.begin(), cell_nodes.end(), x_coefs_c.data());
    x_coefs_c_trans = x_coefs_c.transpose();

    for (int qp = 0; qp < n_qpts_err; ++qp)
    {
      F_c    = x_coefs_c_trans * dLqsi_err[qp];
      Jx     = F_c.determinant();
      volume += Jx * quadr_err->weight(qp);

    } // fim quadratura

  } // end elementos
  return volume;
}

void AppCtx::getSolidVolume()
{
  MatrixXd            x_coefs_c(nodes_per_cell, dim);
  MatrixXd            x_coefs_c_trans(dim, nodes_per_cell);
  Tensor              F_c(dim,dim), invF_c(dim,dim), invFT_c(dim,dim);
  VectorXi            cell_nodes(nodes_per_cell);
  double              Jx;
  int                 tag, nod_id;
  VV.assign(N_Solids,0.0);

  cell_iterator cell = mesh->cellBegin();
  cell_iterator cell_end = mesh->cellEnd();
  for (; cell != cell_end; ++cell)
  {
    tag = cell->getTag();
    nod_id = is_in_id(tag, solidonly_tags);
    if (nod_id){
      mesh->getCellNodesId(&*cell, cell_nodes.data());
      mesh->getNodesCoords(cell_nodes.begin(), cell_nodes.end(), x_coefs_c.data());
      x_coefs_c_trans = x_coefs_c.transpose();

      for (int qp = 0; qp < n_qpts_err; ++qp)
      {
        F_c    = x_coefs_c_trans * dLqsi_err[qp];
        Jx     = F_c.determinant();
        VV[nod_id-1] += Jx * quadr_err->weight(qp);

      } // fim quadratura
    }
  } // end elementos
}

Vector AppCtx::getAreaMassCenterSolid(int sol_id, double &A){
  std::vector<Vector2d> XP;
  Vector Xp(dim);
  int tag, nod_id;
  int NS = NN_Solids[sol_id-1];
  VectorXi mapM_r(dim);
  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  //capture the boundary nodes
  for (; point != point_end; ++point){
	tag = point->getTag();
	nod_id = is_in_id(tag,flusoli_tags);
	if (nod_id == sol_id){
	  dof_handler[DH_UNKM].getVariable(VAR_U).getVertexDofs(mapM_r.data(),&*point);
      //VecGetValues(Vec_x_0, mapM_r.size(), mapM_r.data(), Xp.data());
      //VecGetValues(Vec_x_1, mapM_r.size(), mapM_r.data(), Xp.data());
	  point->getCoord(Xp.data(),dim);
      XP.push_back(Xp);
	}
  }
  XP.resize(NS);
  //ordering the point cloud XP
  std::vector<Vector2d> XPO;
  XPO = ConvexHull2d(XP);
  XPO.push_back(XPO[0]);
  //for (int j = 0; j < NN_Solids[sol_id-1]+1; j++)
//    cout << XPO[j].transpose() << "   ";
//  cout << endl;

  //area
  A = 0;//double A = 0;
  for (int i = 0; i < NS; i++){
    A = A + XPO[i](0)*XPO[i+1](1) - XPO[i+1](0)*XPO[i](1);
  }
  A = 0.5*A;  //cout << A << endl;
  //mass center
  Vector CM(dim);
  CM = Vector::Zero(2);
  for (int i = 0; i < NS; i++){
    CM(0) = CM(0) + (XPO[i](0) + XPO[i+1](0))*(XPO[i](0)*XPO[i+1](1) - XPO[i+1](0)*XPO[i](1));
	CM(1) = CM(1) + (XPO[i](1) + XPO[i+1](1))*(XPO[i](0)*XPO[i+1](1) - XPO[i+1](0)*XPO[i](1));
  }
  CM(0) = CM(0)/(6*A);
  CM(1) = CM(1)/(6*A);
  //cout << CM.transpose() << endl;
  return CM;
}

std::vector<Vector2d> AppCtx::ConvexHull2d(std::vector<Vector2d> & LI){
  int N = LI.size(), k = 0;
  std::vector<Vector2d> H(2*N);
  //lexicographically
  std::sort(LI.begin(),LI.end(),lessVector);
//  for (unsigned int j = 0; j < LI.size(); j++)
//	  cout << LI[j].transpose() << "   ";
//  cout << endl;
  // Build lower hull
  for (int i = 0; i < N; ++i) {
    while (k >= 2 && cross2d(H[k-2], H[k-1], LI[i]) <= 0) k--;
    H[k++] = LI[i];
  }
  // Build upper hull
  for (int i = N-2, t = k+1; i >= 0; i--) {
    while (k >= t && cross2d(H[k-2], H[k-1], LI[i]) <= 0) k--;
    H[k++] = LI[i];
  }

  H.resize(k);
  return H;
}

void AppCtx::calcHmean(double &hmean, double &hmin, double &hmax){
  VectorXi edge_nodes(3);
  Vector Xa(dim), Xb(dim);
  int n_edges=0;
  double hmn = 0, hmx = 0;
  hmean = 0;

  if (dim==2)
  //FEP_PRAGMA_OMP(parallel default(none) shared(hmean))
  {
    const int n_edges_total = mesh->numFacetsTotal();
    Facet const* edge(NULL);
    //double hlist[n_edges_total];

    //FEP_PRAGMA_OMP(for nowait)
    for (int a = 0; a < n_edges_total; ++a)
    {
      edge = mesh->getFacetPtr(a);
      if (edge->isDisabled())
        continue;

      mesh->getFacetNodesId(&*edge, edge_nodes.data());

      mesh->getNodePtr(edge_nodes[0])->getCoord(Xa.data(),dim);
      mesh->getNodePtr(edge_nodes[1])->getCoord(Xb.data(),dim);
      hmean += (Xa-Xb).norm(); //hlist[a] = (Xa-Xb).norm();
      ++n_edges;
    }
    //cout << "Here" << endl;
    //hmn = *min_element(hlist,hlist+n_edges_total);
    //hmx = *max_element(hlist,hlist+n_edges_total);
  }
  else
  if (dim==3)
  //FEP_PRAGMA_OMP(parallel default(none) shared(cout,hmean))
  {
    const int n_edges_total = mesh->numCornersTotal();
    Corner const* edge(NULL);
    //double hlist[n_edges_total];
    //FEP_PRAGMA_OMP(for nowait)
    for (int a = 0; a < n_edges_total; ++a)
    {
      edge = mesh->getCornerPtr(a);
      if (edge->isDisabled())
        continue;

      mesh->getCornerNodesId(&*edge, edge_nodes.data());

      mesh->getNodePtr(edge_nodes[0])->getCoord(Xa.data(),dim);
      mesh->getNodePtr(edge_nodes[1])->getCoord(Xb.data(),dim);
      hmean += (Xa-Xb).norm(); //hlist[a] = (Xa-Xb).norm();
      ++n_edges;
    }
    //hmn = *min_element(hlist,hlist+n_edges_total);
    //hmx = *max_element(hlist,hlist+n_edges_total);
  }
  hmean /= n_edges;
  hmin = hmn; hmax = hmx;
}

bool AppCtx::proxTest(MatrixXd &ContP, MatrixXd &ContW, double const INF){
  ContP = MatrixXd::Constant(N_Solids,N_Solids,INF);
  ContW = MatrixXd::Constant(N_Solids,5,INF);  //5 tags for dirichlet

//  std::vector<bool> Contact(N_Solids+1,false);
//  std::vector<double> Hmin(N_Solids+1,0.0);
  VectorXi edge_nodes(3); // 3 nodes at most
  Vector Xa(dim), Xb(dim);
  const int n_edges_total = mesh->numFacetsTotal();
  Facet const* edge(NULL);
  int tag_a, tag_b, tag_e;
  int K, L;
  bool RepF = false;
  double hv = 0;

  //FEP_PRAGMA_OMP(for nowait)
  for (int a = 0; a < n_edges_total; ++a)
  {
    edge = mesh->getFacetPtr(a);
    if (edge->isDisabled())
      continue;

    mesh->getFacetNodesId(&*edge, edge_nodes.data());

    tag_a = mesh->getNodePtr(edge_nodes[0])->getTag();
    if (mesh->isVertex(mesh->getNodePtr(edge_nodes[1])))
      tag_b = mesh->getNodePtr(edge_nodes[1])->getTag();
    else
      tag_b = mesh->getNodePtr(edge_nodes[2])->getTag();
    tag_e = edge->getTag();
    //cout << tag_a << "  " << tag_b << "  " << tag_e << endl;

    if ((tag_a==tag_b) && (tag_b==tag_e))  //only edges outside the same type of domain
      continue;
    if ( (is_in(tag_a, fluidonly_tags) && !is_in(tag_a, dirichlet_tags)) ||
         (is_in(tag_b, fluidonly_tags) && !is_in(tag_b, dirichlet_tags)) )
      continue;
    if (is_in(tag_e, dirichlet_tags))                                  //avoid dirichlet edges
      continue;
    if (is_in(tag_a, dirichlet_tags) && is_in(tag_b, dirichlet_tags))  //avoid dirichlet edges
      continue;

    RepF = true;

    mesh->getNodePtr(edge_nodes[0])->getCoord(Xa.data(),dim);
    if (mesh->isVertex(mesh->getNodePtr(edge_nodes[1])))
      mesh->getNodePtr(edge_nodes[1])->getCoord(Xb.data(),dim);
    else
      mesh->getNodePtr(edge_nodes[2])->getCoord(Xb.data(),dim);
    hv = (Xa-Xb).norm();  //cout << Xa.transpose() << "   " << Xb.transpose() << endl;

    K = is_in_id(tag_a,flusoli_tags);  //K or L = 0, means the point is wall
    L = is_in_id(tag_b,flusoli_tags);
    if ((K != 0) && (L != 0)){
      if (hv < ContP(K-1,L-1)){
        ContP(K-1,L-1) = hv;
        ContP(L-1,K-1) = hv;
      }  //saves the smallest distance
    }
    else if (K == 0){
      K = is_in_id(tag_a,dirichlet_tags);
      cout << K << "  " << mesh->getPointId(mesh->getNodePtr(edge_nodes[0]))
                << "  " << mesh->getPointId(mesh->getNodePtr(edge_nodes[1])) << endl;
      if (hv < ContW(L-1,K-1)){
        ContW(L-1,K-1) = hv;
      }
    }
    else if (L == 0){
      L = is_in_id(tag_b,dirichlet_tags);
      cout << L << "  " << mesh->getPointId(mesh->getNodePtr(edge_nodes[0]))
                << "  " << mesh->getPointId(mesh->getNodePtr(edge_nodes[1])) << endl;
      if (hv < ContW(K-1,L-1)){
        ContW(K-1,L-1) = hv;
      }
    }
  }//end for edges
  return RepF;
}

void AppCtx::forceDirichlet(){
  Vector    U1(dim);
  Vector    X1(dim);
  VectorXi  u_dofs(dim);
  int tag;
  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for ( ; point != point_end; ++point)
  {
    tag = point->getTag();
    if (is_in(tag,dirichlet_tags)){
      point->getCoord(X1.data(),dim);
      getNodeDofs(&*point, DH_UNKM, VAR_U, u_dofs.data());
      U1 = u_exact(X1, current_time+dt, tag);
      VecSetValues(Vec_uzp_1, dim, u_dofs.data(), U1.data(), INSERT_VALUES);
    }
  }
  Assembly(Vec_uzp_1);
}

PetscErrorCode AppCtx::updateSolidVel()
{
  VectorXi  dofs(dim), mapM_r(dim);
  VectorXi  dofs_fs(LZ);
  Vector    X(dim);
  Vector    Uf(dim), Zf(LZ), Vs(dim);
  Vector3d  Xg;
  int       nod_id, nod_is, nod_vs, nodsum;
  int       tag;
  std::vector<bool>   SV(N_Solids,false);  //solid visited history

  //XG_0 = XG;

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    tag = point->getTag();
    // vel
    nod_id = is_in_id(tag,flusoli_tags);
    nod_is = is_in_id(tag,solidonly_tags);
    nod_vs = is_in_id(tag,slipvel_tags);
    nodsum = nod_id+nod_is+nod_vs;
    if (nodsum){
      for (int l = 0; l < LZ; l++){
        dofs_fs(l) = n_unknowns_u + n_unknowns_p + LZ*(nodsum) - LZ + l;
      }
      dof_handler[DH_UNKM].getVariable(VAR_U).getVertexDofs(mapM_r.data(),&*point);
      //VecGetValues(Vec_x_0, mapM_r.size(), mapM_r.data(), X.data());
      VecGetValues(Vec_x_1, mapM_r.size(), mapM_r.data(), X.data());
      //point->getCoord(X.data(),dim);
      VecGetValues(Vec_uzp_1, LZ, dofs_fs.data(), Zf.data());
      //update center of mass position
/*      if(!SV[nod_id+nod_is-1]){
        Xg = XG_0[nod_id+nod_is-1];
        Xg(0) = Xg(0)+dt*Zf(0); Xg(1) = Xg(1)+dt*Zf(1); if (dim == 3){Xg(2) = Xg(2)+dt*Zf(2);}
        XG[nod_id+nod_is-1] = Xg;  SV[nod_id+nod_is-1] = true;
      }*/
      Uf = SolidVel(X, XG_1[nodsum-1], Zf, dim);
      if (nod_vs){
        getNodeDofs(&*point, DH_MESH, VAR_M, dofs.data());
        VecGetValues(Vec_slipv_1, dim, dofs.data(), Vs.data());
        Uf = Uf + Vs;  //cout << tag << "  " << X(0)-3 << " " << X(1)-3<< "  " << Uf.transpose() << endl;
      }
      //cout << dofs_fs.transpose() << ", " << X.transpose() << ", " << XG_0[nod_id+nod_is-1].transpose() << ", " << Zf.transpose() << ", "<< Uf.transpose() << endl;
      getNodeDofs(&*point,DH_UNKM,VAR_U,dofs.data());
      VecSetValues(Vec_uzp_1, dim, dofs.data(), Uf.data(), INSERT_VALUES);
    }
  } // end point loop
  Assembly(Vec_uzp_1);
  //char buf1[18], buf2[10];
  //sprintf(buf1,"matrizes/sol%d.m",time_step); sprintf(buf2,"solm%d",time_step);
  //View(Vec_uzp_1, buf1, buf2);
  PetscFunctionReturn(0);  //cout << endl;
}

PetscErrorCode AppCtx::moveCenterMass(double vtheta)
{
  VectorXi    dofs(dim), dof(1);
  Vector      U0(Vector::Zero(3)), U1(Vector::Zero(3)), XG_temp(Vector::Zero(3)), theta0(1), theta1(1);
  double      theta_temp;

  for (int s = 0; s < N_Solids; s++){

    for (int l = 0; l < dim; l++){
      dofs(l) = n_unknowns_u + n_unknowns_p + LZ*(s+1) - LZ + l;
    }
    VecGetValues(Vec_uzp_0, dim, dofs.data(), U0.data());
    VecGetValues(Vec_uzp_1, dim, dofs.data(), U1.data());
    dof(0) = n_unknowns_u + n_unknowns_p + LZ*s + dim;
    VecGetValues(Vec_uzp_0, 1, dof.data(), theta0.data());
    VecGetValues(Vec_uzp_1, 1, dof.data(), theta1.data());

    if (is_bdf3 && time_step > 1){
      XG_temp   = XG_1[s];
      XG_1[s]   = (18./11.)*XG_1[s] - (9./11.)*XG_0[s] + (2./11.)*XG_aux[s] + (6./11.)*dt*U0;
      XG_aux[s] = XG_0[s];
      XG_0[s]   = XG_temp;
      theta_temp   = theta_1[s];
      theta_1[s]   = (18./11.)*theta_1[s] - (9./11.)*theta_0[s] + (2./11.)*theta_aux[s] + (6./11.)*dt*theta0(0);
      theta_aux[s] = theta_0[s];
      theta_0[s]   = theta_temp;
    }
    else if (is_bdf2 && time_step > 0){
      if (is_bdf2_bdfe){
        XG_temp = XG_1[s];
        XG_1[s] = (4./3.)*XG_1[s] - (1./3.)*XG_0[s] + (2./3.)*dt*(vtheta*U1 + (1.-vtheta)*U0);
        XG_0[s] = XG_temp;
        theta_temp = theta_1[s];
        theta_1[s] = (4./3.)*theta_1[s] - (1./3.)*theta_0[s] + (2./3.)*dt*(vtheta*theta1(0) + (1.-vtheta)*theta0(0));
        theta_0[s] = theta_temp;
      }
      else if (is_bdf2_ab){
        XG_temp = XG_1[s];
        XG_1[s] = dt*(vtheta*U1 + (1.-vtheta)*U0) + XG_0[s];
        XG_0[s] = XG_temp;
        theta_temp = theta_1[s];
        theta_1[s] = dt*(vtheta*theta1(0) + (1.-vtheta)*theta0(0)) + theta_0[s];
        theta_0[s] = theta_temp;
      }
    }
    else{  //for MR-AB and basic, and for all at time = t0
      if (time_step == 0){
        XG_1[s] = dt*(vtheta*U1 + (1.-vtheta)*U0) + XG_0[s];
        theta_1[s] = dt*(vtheta*theta1(0) + (1.-vtheta)*theta0(0)) + theta_0[s];
      }
      else{
        XG_temp = XG_1[s];
        XG_1[s] = dt*(vtheta*U1 + (1.-vtheta)*U0) + XG_0[s];
        XG_0[s] = XG_temp;
        theta_temp = theta_1[s];
        theta_1[s] = dt*(vtheta*theta1(0) + (1.-vtheta)*theta0(0)) + theta_0[s];
        theta_0[s] = theta_temp;
      }
    }

  }
  PetscFunctionReturn(0);
}

PetscErrorCode AppCtx::updateSolidMesh()
{
  int tag, nod_id, nod_is, nod_vs, nodsum;
  VectorXi  dofs(dim);
  Vector    X0(Vector::Zero(3)), Vs(Vector::Zero(3));

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    tag = point->getTag();

    nod_id = is_in_id(tag,flusoli_tags);
    nod_is = is_in_id(tag,solidonly_tags);
    nod_vs = is_in_id(tag,slipvel_tags);
    nodsum = nod_id+nod_is+nod_vs;
    if (nodsum){
      getNodeDofs(&*point, DH_MESH, VAR_M, dofs.data());
      VecGetValues(Vec_x_0, dofs.size(), dofs.data(), X0.data());
      X0 = RotM(theta_1[nodsum-1]-theta_0[nodsum-1],dim)*(X0 - XG_0[nodsum-1]) + XG_1[nodsum-1];// - XG_0[nodsum-1];
      VecSetValues(Vec_x_1, dofs.size(), dofs.data(), X0.data(), INSERT_VALUES);
      nod_vs = is_in_id(tag,slipvel_tags);
      if (nod_vs){
        VecGetValues(Vec_slipv_0, dim, dofs.data(), Vs.data());
        Vs = RotM(theta_1[nod_vs-1]-theta_0[nod_vs-1],dim)*Vs;// + XG_1[nod_id+nod_is-1]; // - XG_0[nod_id+nod_is-1];
        VecSetValues(Vec_slipv_1, dim, dofs.data(), Vs.data(), INSERT_VALUES);
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode AppCtx::velNoSlip(Vec const& Vec_uzp, Vec const& Vec_sv, Vec &Vec_uzp_ns)
{
  int tag;
  VectorXi dofs(dim), dofs_mesh(dim);
  Vector   U(dim), Uns(dim);

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)  //to calculate Vec_v_mid at each point (initial guess)
  {
    tag = point->getTag();
    if (is_in(tag,slipvel_tags)){
      getNodeDofs(&*point, DH_UNKM, VAR_U, dofs.data());
      getNodeDofs(&*point, DH_MESH, VAR_M, dofs_mesh.data());
      VecGetValues(Vec_uzp,  dim, dofs.data(), U.data());
      VecGetValues(Vec_sv,  dim, dofs_mesh.data(), Uns.data());
      Uns = U - Uns;
      VecSetValues(Vec_uzp_ns, dim, dofs.data(), Uns.data(), INSERT_VALUES);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode AppCtx::plotFiles()
{
  double  *q_array;
  double  *nml_array;
  double  *v_array, *vs_array;;
  VecGetArray(Vec_uzp_0, &q_array);  //VecGetArray(Vec_up_0, &q_array);
  VecGetArray(Vec_normal, &nml_array);
  VecGetArray(Vec_v_mid, &v_array);
  if (is_slipv) VecGetArray(Vec_slipv_0, &vs_array);
  vtk_printer.writeVtk();

  /* ---- nodes data ---- */
  vtk_printer.addNodeVectorVtk("u", GetDataVelocity(q_array, *this));
  vtk_printer.addNodeVectorVtk("n", GetDataNormal(nml_array, *this));
  vtk_printer.addNodeVectorVtk("v", GetDataMeshVel(v_array, *this));
  if (is_slipv) vtk_printer.addNodeVectorVtk("vs", GetDataMeshVel(vs_array, *this));
  vtk_printer.printPointTagVtk();

  if (!shape_psi_c->discontinuous())
    vtk_printer.addNodeScalarVtk("p", GetDataPressure(q_array, *this));
  else
    vtk_printer.addCellScalarVtk("p", GetDataPressCellVersion(q_array, *this));

  vtk_printer.addCellIntVtk("cell_tag", GetDataCellTag(*this));

  //vtk_printer.printPointTagVtk("point_tag");
  VecRestoreArray(Vec_uzp_0, &q_array);  //VecRestoreArray(Vec_up_0, &q_array);
  VecRestoreArray(Vec_normal, &nml_array);
  VecRestoreArray(Vec_v_mid, &v_array);
  if (is_slipv) VecRestoreArray(Vec_slipv_0, &vs_array);

  PetscFunctionReturn(0);
}

void AppCtx::freePetscObjs()
{
  Destroy(Mat_Jac_fs); //Destroy(Mat_Jac);
  Destroy(Mat_Jac_m);
  Destroy(Vec_res_fs); //Destroy(Vec_res);
  Destroy(Vec_res_m);
  Destroy(Vec_uzp_0); //Destroy(Vec_up_0);
  Destroy(Vec_uzp_1); //Destroy(Vec_up_1);
  Destroy(Vec_duzp); //Destroy(Vec_dup);
  Destroy(Vec_v_mid);
  Destroy(Vec_v_1);
  Destroy(Vec_x_0);
  Destroy(Vec_x_1);
  Destroy(Vec_normal);
  Destroy(Vec_uzp_m1);
  //Destroy(ksp);
  //Destroy(snes);
  //SNESDestroy(&snes);
  SNESDestroy(&snes_m);
  SNESDestroy(&snes_fs);
  if (is_bdf3)
  {
    Destroy(Vec_x_aux);
    Destroy(Vec_duzp_0);
    Destroy(Vec_uzp_m2);
  }
  if (is_slipv){
    Destroy(Vec_slipv_0);
    Destroy(Vec_slipv_1);
    Destroy(Vec_slipv_m1);
    if (is_bdf3)
      Destroy(Vec_slipv_m2);
  }
}

void GetDataVelocity::get_vec(int id, Real * vec_out) const
{
  Point const* point = user.mesh->getNodePtr(id);
  std::vector<int> dofs(user.dim);

  user.getNodeDofs(&*point, DH_UNKM, VAR_U, dofs.data());
  for (int i = 0; i < user.dim; ++i)
    vec_out[i] = q_array[*(dofs.data()+i)];
 /*
  int tag = point->getTag();
  int nod_id = is_in_id(tag,user.flusoli_tags);
  int nod_is = is_in_id(tag,user.solidonly_tags);
  Vector Xp(user.dim), Z(user.LZ), U(user.dim);
  Vector3d Xg;

  if (nod_id || nod_is){
    Xg = user.XG_1[(nod_id+nod_is)-1];
    point->getCoord(Xp.data(),user.dim);
    for (int l = 0; l < user.LZ; l++){
      Z(l) = q_array[user.n_unknowns_u + user.n_unknowns_p + user.LZ*(nod_id+nod_is) - user.LZ + l];
    }
    U = SolidVel(Xp,Xg,Z,user.dim);
    for (int i = 0; i < user.dim; ++i)
      vec_out[i] = U(i);
  }
  else{
    user.getNodeDofs(&*point, DH_UNKM, VAR_U, dofs.data());
    for (int i = 0; i < user.dim; ++i)
      vec_out[i] = q_array[*(dofs.data()+i)];
  }
*/
}

double GetDataPressure::get_data_r(int nodeid) const
{
  Point const*const point = user.mesh->getNodePtr(nodeid);
  if (!user.mesh->isVertex(point))
  {
    int dofs[3];
    // position of the node at edge
    const int m = point->getPosition() - user.mesh->numVerticesPerCell();
    Cell const*const cell = user.mesh->getCellPtr(point->getIncidCell());
    if (user.dim==3)
    {
      const int edge_id = cell->getCornerId(m);
      user.dof_handler[DH_UNKM].getVariable(VAR_P).getCornerDofs(dofs, user.mesh->getCornerPtr(edge_id));
    }
    else
    {
      const int edge_id = cell->getFacetId(m);
      user.dof_handler[DH_UNKM].getVariable(VAR_P).getFacetDofs(dofs, user.mesh->getFacetPtr(edge_id));
    }
    return (q_array[dofs[0]] + q_array[dofs[1]])/2.;
  }
  int dof;
  user.dof_handler[DH_UNKM].getVariable(VAR_P).getVertexDofs(&dof, point);
  return q_array[dof];
}

double GetDataPressCellVersion::get_data_r(int cellid) const
{
  // assume que só há 1 grau de liberdade na célula
  int dof[user.dof_handler[DH_UNKM].getVariable(VAR_P).numDofsPerCell()];
  user.dof_handler[DH_UNKM].getVariable(VAR_P).getCellDofs(dof, user.mesh->getCellPtr(cellid));
  return q_array[dof[0]];
}

void GetDataNormal::get_vec(int id, Real * vec_out) const
{
  Point const* point = user.mesh->getNodePtr(id);
  vector<int> dofs(user.dim);

  user.getNodeDofs(&*point, DH_MESH, VAR_M, dofs.data());

  for (int i = 0; i < user.dim; ++i)
    vec_out[i] = q_array[*(dofs.data()+i)];
}

void GetDataMeshVel::get_vec(int id, Real * vec_out) const
{
  Point const* point = user.mesh->getNodePtr(id);
  vector<int> dofs(user.dim);

  user.getNodeDofs(&*point, DH_MESH, VAR_M, dofs.data());

  for (int i = 0; i < user.dim; ++i)
    vec_out[i] = q_array[*(dofs.data()+i)];
}

int GetDataCellTag::get_data_i(int cellid) const
{
  // assume que só há 1 grau de liberdade na célula
  return user.mesh->getCellPtr(cellid)->getTag();
}

/* petsc functions*/
//extern PetscErrorCode FormJacobian(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
//extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
//extern PetscErrorCode Monitor(SNES,PetscInt,PetscReal,void *);
//extern PetscErrorCode FormJacobian(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
//extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
extern PetscErrorCode FormJacobian_fs(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
extern PetscErrorCode FormFunction_fs(SNES,Vec,Vec,void*);


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{

  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

  bool help_return;
  bool erro;
  AppCtx user(argc, argv, help_return, erro);

  if (help_return)
    return 0;

  if (erro)
    return 1;

#ifdef FEP_HAS_OPENMP
  {
    int nthreads = 0;
    #pragma omp parallel
    {
      #pragma omp critical
      nthreads = omp_get_num_threads();
    }
    printf("OpenMP version.\n");
    printf("Num. threads: %d\n", nthreads);
  }
#else
  {
    printf("Serial version.\n");
  }
#endif


  user.loadMesh();
  user.loadDofs();  //counts dofs for U,P,V in cell, facet, corner
  user.evaluateQuadraturePts();

  erro = user.err_checks(); if (erro) return 1;

  // print info
  cout << endl; cout << "mesh: " << user.filename << endl;
  user.mesh->printInfo();
  cout << "\n# velocity unknowns: " << user.n_unknowns_u;
  cout << "\n# pressure unknowns: " << user.n_unknowns_p << " = "
                                    << user.n_unknowns_u/user.dim << " + " << user.n_nodes_fsi;
//  cout << "\n# total unknowns: " << user.n_unknowns << endl;
  if (user.N_Solids){
    cout << "\n# velocity fluid only unknowns: " << user.n_unknowns_u - user.dim*(user.n_nodes_so+user.n_nodes_fsi);
    cout << "\n# velocity solid only unknowns: " << user.n_nodes_so*user.dim;
    cout << "\n# velocity solid int unknowns: " << user.n_unknowns_z << " (" <<
    		                                       user.n_nodes_fsi << ") = " <<
    		                                       user.n_unknowns_z*user.n_nodes_fsi;
    cout << "\n# total unknowns: " << user.n_unknowns_fs << " = " << user.n_unknowns_fs - user.dim*user.n_nodes_fsi
                                                         << " + " << user.dim*user.n_nodes_fsi;
    cout << "\n# solids: " << user.N_Solids << endl;
    cout << "unknowns distribution: " << 0 << "-" << user.n_unknowns_u-1 <<
    		", " << user.n_unknowns_u << "-" <<
			user.n_unknowns_u + user.n_unknowns_p - 1 <<
			", " << user.n_unknowns_u + user.n_unknowns_p <<
			"-"  << user.n_unknowns_u + user.n_unknowns_p +
			        user.n_unknowns_z - 1;
  }
  cout << endl;
  user.mesh->printStatistics();
  user.mesh->timer.printTimes();

  cout << endl;
  if (user.N_Solids){
    cout << "flusoli_tags =  ";
    for (int i = 0; i < (int)user.flusoli_tags.size(); ++i)
    {
      cout << user.flusoli_tags[i] << " ";
    }
  cout << endl;
  }

  printf("on update mesh\n");
  user.onUpdateMesh();
  user.solveTimeProblem();

  cout << "\n";
  user.timer.printTimes();

  user.freePetscObjs();
  PetscFinalize();

  cout << "\a" << endl;
  return 0.;
}

/* ------------------------------------------------------------------- */
#if (false)
#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
PetscErrorCode FormJacobian(SNES snes,Vec Vec_up_1,Mat *Mat_Jac, Mat *prejac, MatStructure *flag, void *ptr)
{
  AppCtx *user    = static_cast<AppCtx*>(ptr);
  user->formJacobian(snes,Vec_up_1,Mat_Jac,prejac,flag);
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
PetscErrorCode FormFunction(SNES snes, Vec Vec_up_1, Vec Vec_fun, void *ptr)
{
  AppCtx *user    = static_cast<AppCtx*>(ptr);
  user->formFunction(snes,Vec_up_1,Vec_fun);
  PetscFunctionReturn(0);
}
#endif
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "CheckSnesConvergence"
PetscErrorCode CheckSnesConvergence(SNES snes, PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx)
{
  AppCtx *user    = static_cast<AppCtx*>(ctx);
  user->checkSnesConvergence(snes, it, xnorm, pnorm, fnorm, reason);
  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormJacobian_mesh"
PetscErrorCode FormJacobian_mesh(SNES snes,Vec Vec_up_1,Mat *Mat_Jac, Mat *prejac, MatStructure *flag, void *ptr)
{
  AppCtx *user    = static_cast<AppCtx*>(ptr);
  user->formJacobian_mesh(snes,Vec_up_1,Mat_Jac,prejac,flag);
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "FormFunction_mesh"
PetscErrorCode FormFunction_mesh(SNES snes, Vec Vec_up_1, Vec Vec_fun, void *ptr)
{
  AppCtx *user    = static_cast<AppCtx*>(ptr);
  user->formFunction_mesh(snes,Vec_up_1,Vec_fun);
  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormJacobian_fs"
PetscErrorCode FormJacobian_fs(SNES snes,Vec Vec_up_1,Mat *Mat_Jac, Mat *prejac, MatStructure *flag, void *ptr)
{
  AppCtx *user    = static_cast<AppCtx*>(ptr);
  user->formJacobian_fs(snes,Vec_up_1,Mat_Jac,prejac,flag);
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
PetscErrorCode FormFunction_fs(SNES snes, Vec Vec_up_1, Vec Vec_fun, void *ptr)
{
  AppCtx *user    = static_cast<AppCtx*>(ptr);
  user->formFunction_fs(snes,Vec_up_1,Vec_fun);
  PetscFunctionReturn(0);
}
