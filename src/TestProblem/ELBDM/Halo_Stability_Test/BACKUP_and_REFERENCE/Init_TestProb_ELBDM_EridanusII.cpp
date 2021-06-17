#include "GAMER.h"
#include "TestProb.h"



// soliton-specific global variables
// =======================================================================================
       bool     FixDM;                                   // true --> do not evolve psidm at all
       double   Soliton_CoreRadius;                      // soliton core radius
static int      Soliton_InputMode;                       // soliton input mode: 1/2/3 -> table/approximate analytical form/none
static double   Soliton_OuterSlope;                      // soliton outer slope (only used by Soliton_InputMode=2)
static char     Soliton_DensProf_Filename[MAX_STRING];   // filename of the reference soliton density profile

static int      Soliton_DensProf_NBin;                   // number of radial bins of the soliton density profile
static double  *Soliton_DensProf   = NULL;               // soliton density profile [radius/density]
static double   Soliton_ScaleL     = NULL;               // L/D: length/density scale factors of each soliton
                                                         //      (defined as the ratio between the core radii/peak
                                                         //      density of the target and reference soliton profiles)
static double   Soliton_ScaleD     = NULL;
static double   Soliton_CM_MaxR;                         // maximum radius for determining CM
static double   Soliton_CM_TolErrR;                      // maximum allowed errors for determining CM

       bool     Tidal_Enabled;                           // enable tidal field (as an external potential)
       bool     Tidal_RotatingFrame;                     // true/false --> rotating/inertial frame
       double   Tidal_Mass;                              // point mass
       double   Tidal_R;                                 // point mass distance
       double   Tidal_Angle0;                            // initial angle of the point mass
       bool     Tidal_FixedPos;                          // Fix the point mass position
       bool     Tidal_Centrifugal;                       // Add the centrifugal pseudo force
       double   Tidal_CutoffR;                           // Cut-off radius

       double   Tidal_Vrot;                              // rotational velocity due to the point mass
       double   Tidal_CM[3] = { __DBL_MAX__, __DBL_MAX__, __DBL_MAX__ }; // Center of mass of the satellite

static int      Sponge_Mode;                             // 1/2/3: sponge BC/truncation BC/none
static double   Sponge_Width;                            // sponge width
static double   Sponge_Amp;                              // sponge amplitude

       double   Sponge_dt = 0.0;                         // evolution time-step
// =======================================================================================


// particle-specific global variables
// =======================================================================================
       int    Star_RSeed;           // random seed for setting particle position and velocity
       double Star_Rho0;            // peak density
       double Star_R0;              // scale radius
       double Star_MaxR;            // maximum radius for particles
       double Star_Center[3];       // center coordinates
       int    Star_MassProfNBin;    // number of radial bins in the mass profile table
       int    Star_SigmaMode;       // different modes for assigning stellar velocity dispersion
                                    // 0: self-bound; 1: constant soliton peak density

       bool   Star_AddParForRestart;         // add particles after restart
       long   Star_AddParForRestart_NPar;    // number of particles for Star_AddParForRestart
       double Star_AddParForRestart_PeakRho; // Peak density for Star_AddParForRestart and Star_SigmaMode=1

static double Star_FreeT;           // free-fall time at Star_R0

       bool   ParFileCM_Enabled;    // shift the CM of particles loaded from PAR_IC
       double ParFileCM_Dis[3];     // CM displacement/velocity for ParFileCM_Enabled
       double ParFileCM_Vel[3];

       int    DensRecMode;          // different modes for recording the peak density (1/2/3: ELBDM/particle/total density)
// =======================================================================================

// problem-specific function prototypes
#ifdef PARTICLE
void Par_Init_ByFunction_EridanusII( const long NPar_ThisRank, const long NPar_AllRank,
                                     real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                     real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                     real *AllAttribute[PAR_NATT_TOTAL] );
void Init_User_EridanusII();
#endif

// external potential routines
void Init_ExtPot_EridanusII();
void SetExtPotAuxArray_EridanusII( double AuxArray_Flt[], int AuxArray_Int[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


// errors
#  if ( MODEL != ELBDM )
   Aux_Error( ERROR_INFO, "MODEL != ELBDM !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef GRAVITY
   if ( OPT__BC_POT != BC_POT_ISOLATED )
      Aux_Error( ERROR_INFO, "must adopt isolated BC for gravity --> reset OPT__BC_POT !!\n" );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM  &&  defined GRAVITY )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",           &VARIABLE,                   DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "FixDM",                     &FixDM,                      false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "Soliton_CoreRadius",        &Soliton_CoreRadius,        -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Soliton_InputMode",         &Soliton_InputMode,          1,             1,                3                 );
   ReadPara->Add( "Soliton_OuterSlope",        &Soliton_OuterSlope,        -8.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_DensProf_Filename",  Soliton_DensProf_Filename,  Useless_str,   Useless_str,      Useless_str       );
   ReadPara->Add( "Soliton_CM_MaxR",           &Soliton_CM_MaxR,           -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Soliton_CM_TolErrR",        &Soliton_CM_TolErrR,        -1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Star_RSeed",                &Star_RSeed,                 123,           0,                NoMax_int         );
   ReadPara->Add( "Star_SigmaMode",            &Star_SigmaMode,             0,             0,                1                 );
   ReadPara->Add( "Star_Rho0",                 &Star_Rho0,                 -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Star_R0",                   &Star_R0,                   -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Star_MaxR",                 &Star_MaxR,                 -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Star_CenterX",              &Star_Center[0],             amr->BoxCenter[0], NoMin_double, NoMax_double      );
   ReadPara->Add( "Star_CenterY",              &Star_Center[1],             amr->BoxCenter[1], NoMin_double, NoMax_double      );
   ReadPara->Add( "Star_CenterZ",              &Star_Center[2],             amr->BoxCenter[2], NoMin_double, NoMax_double      );
   ReadPara->Add( "Star_MassProfNBin",         &Star_MassProfNBin,          1000,          2,                NoMax_int         );
   ReadPara->Add( "Star_AddParForRestart",     &Star_AddParForRestart,      false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "Star_AddParForRestart_NPar",&Star_AddParForRestart_NPar,-1L,            NoMin_long,       NoMax_long        );
   ReadPara->Add( "Star_AddParForRestart_PeakRho", &Star_AddParForRestart_PeakRho, -1.0,   NoMin_double,     NoMax_double      );
   ReadPara->Add( "Tidal_Enabled",             &Tidal_Enabled,              false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "Tidal_RotatingFrame",       &Tidal_RotatingFrame,        true,          Useless_bool,     Useless_bool      );
   ReadPara->Add( "Tidal_Mass",                &Tidal_Mass,                -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Tidal_R",                   &Tidal_R,                   -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Tidal_Angle0",              &Tidal_Angle0,               0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Tidal_FixedPos",            &Tidal_FixedPos,             false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "Tidal_Centrifugal",         &Tidal_Centrifugal,          false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "Tidal_CutoffR",             &Tidal_CutoffR,              __DBL_MAX__,   0.0,              NoMax_double      );

   ReadPara->Add( "Sponge_Mode",               &Sponge_Mode,                3,             1,                3                 );
   ReadPara->Add( "Sponge_Width",              &Sponge_Width,               10.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Sponge_Amp",                &Sponge_Amp,                 1.0,           0.0,              NoMax_double      );

   ReadPara->Add( "ParFileCM_Enabled",         &ParFileCM_Enabled,          false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "ParFileCM_DisX",            &ParFileCM_Dis[0],           0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "ParFileCM_DisY",            &ParFileCM_Dis[1],           0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "ParFileCM_DisZ",            &ParFileCM_Dis[2],           0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "ParFileCM_VelX",            &ParFileCM_Vel[0],           0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "ParFileCM_VelY",            &ParFileCM_Vel[1],           0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "ParFileCM_VelZ",            &ParFileCM_Vel[2],           0.0,           NoMin_double,     NoMax_double      );

   ReadPara->Add( "DensRecMode",               &DensRecMode,                3,             1,                3                 );

   ReadPara->Read( FileName );

   delete ReadPara;

// convert to code units
   Tidal_Mass    *= Const_Msun / UNIT_M;
   Tidal_R       *= Const_kpc  / UNIT_L;
   Tidal_CutoffR *= Const_kpc  / UNIT_L;
   Sponge_Width  *= Const_kpc  / UNIT_L;
   Sponge_Amp    *= (1.0/Const_Gyr) / (1.0/UNIT_T);
   for (int d=0; d<3; d++)
   {
      ParFileCM_Dis[d] *= Const_kpc        / UNIT_L;
      ParFileCM_Vel[d] *= Const_km/Const_s / (UNIT_L/UNIT_T);
   }

// (1-2) set the default values
   if ( Soliton_CM_TolErrR < 0.0 )  Soliton_CM_TolErrR = 1.0*amr->dh[MAX_LEVEL];

// (1-3) check the runtime parameters
   if ( FixDM  &&  OPT__FIXUP_FLUX )   Aux_Error( ERROR_INFO, "must disable OPT__FIXUP_FLUX for FixDM !!\n" );

#  ifndef PARTICLE
   if ( DensRecMode == 2 )    Aux_Error( ERROR_INFO, "DensRecMode == 2 must work with PARTICLE !!\n" );
#  endif


// (2) set the problem-specific derived parameters
#  ifdef GRAVITY
   Star_FreeT = sqrt( (3.0*M_PI*pow(2.0,1.5)) / (32.0*NEWTON_G*Star_Rho0) );
#  endif

   Tidal_Vrot = sqrt( NEWTON_G * Tidal_Mass / Tidal_R );


// (3) load the reference soliton density profile and evaluate the scale factors
   if ( OPT__INIT != INIT_BY_RESTART  &&  Soliton_InputMode == 1 )
   {
//    load the reference profile
      const bool RowMajor_No  = false;    // load data into the column-major order
      const bool AllocMem_Yes = true;     // allocate memory for Soliton_DensProf
      const int  NCol         = 2;        // total number of columns to load
      const int  Col[NCol]    = {0, 1};   // target columns: (radius, density)

      Soliton_DensProf_NBin = Aux_LoadTable( Soliton_DensProf, Soliton_DensProf_Filename, NCol, Col, RowMajor_No, AllocMem_Yes );


//    get the core radius of the reference profile
      const double *RadiusRef = Soliton_DensProf + 0*Soliton_DensProf_NBin;
      const double *DensRef   = Soliton_DensProf + 1*Soliton_DensProf_NBin;
      const double  DensCore  = 0.5*DensRef[0];   // define core radius as the half-density radius

      double CoreRadiusRef = NULL_REAL;

      for (int b=1; b<Soliton_DensProf_NBin-1; b++)
      {
         if ( DensRef[b] >= DensCore  &&  DensRef[b+1] <= DensCore )
         {
            CoreRadiusRef = 0.5*( RadiusRef[b] + RadiusRef[b+1] );
            break;
         }
      }

      if ( CoreRadiusRef == NULL_REAL )
         Aux_Error( ERROR_INFO, "cannot determine the reference core radius !!\n" );


//    evaluate the scale factors of each soliton
      Soliton_ScaleL = Soliton_CoreRadius / CoreRadiusRef;
      Soliton_ScaleD = 1.0 / ( 4.0*M_PI*NEWTON_G*SQR(ELBDM_ETA)*POW4(Soliton_ScaleL) );
   } // if ( OPT__INIT != INIT_BY_RESTART )


// (4) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 12.0*Const_Gyr/UNIT_T;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }

   if ( Tidal_Enabled )
   {
      if ( OPT__EXT_POT != EXT_POT_FUNC )
         Aux_Error( ERROR_INFO, "must set OPT__EXT_POT = EXT_POT_FUNC for Tidal_Enabled !!\n" );

      if ( !OPT__RECORD_USER )
         Aux_Error( ERROR_INFO, "must enable OPT__RECORD_USER for Tidal_Enabled !!\n" );
   }

   if ( Sponge_Mode != 3 )
   {
      if ( !OPT__RESET_FLUID )
         Aux_Error( ERROR_INFO, "must enable OPT__RESET_FLUID for Sponge_Mode != 3!!\n" );
   }


// (5) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "======================================================================================\n" );
      Aux_Message( stdout, "  test problem ID                           = %d\n",     TESTPROB_ID                );
      Aux_Message( stdout, "  fix dark matter                           = %d\n",     FixDM                      );
      Aux_Message( stdout, "  soliton core radius                       = %13.6e\n", Soliton_CoreRadius         );
      Aux_Message( stdout, "  soliton input mode                        = %d\n",     Soliton_InputMode          );
      if      ( Soliton_InputMode == 2 )
      Aux_Message( stdout, "  soliton outer slope                       = %13.6e\n", Soliton_OuterSlope         );
      else if ( Soliton_InputMode == 1 ) {
      Aux_Message( stdout, "  density profile filename                  = %s\n",     Soliton_DensProf_Filename  );
      Aux_Message( stdout, "  number of bins of the density profile     = %d\n",     Soliton_DensProf_NBin      ); }
      Aux_Message( stdout, "  soliton CM max radius                     = %13.6e\n", Soliton_CM_MaxR            );
      Aux_Message( stdout, "  soliton CM tolerated error                = %13.6e\n", Soliton_CM_TolErrR         );
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "  star cluster properties:\n" );
      Aux_Message( stdout, "  random seed for setting particle position = %d\n",     Star_RSeed );
      Aux_Message( stdout, "  velocity dispersion mode                  = %d\n",     Star_SigmaMode );
      Aux_Message( stdout, "  peak density                              = %13.7e\n", Star_Rho0 );
      Aux_Message( stdout, "  scale radius                              = %13.7e\n", Star_R0 );
      Aux_Message( stdout, "  maximum radius of particles               = %13.7e\n", Star_MaxR );
      Aux_Message( stdout, "  central coordinate x                      = %13.7e\n", Star_Center[0] );
      Aux_Message( stdout, "  central coordinate y                      = %13.7e\n", Star_Center[1] );
      Aux_Message( stdout, "  central coordinate z                      = %13.7e\n", Star_Center[2] );
      Aux_Message( stdout, "  number of radial bins in the mass profile = %d\n",     Star_MassProfNBin );
      Aux_Message( stdout, "  free-fall time at the scale radius        = %13.7e\n", Star_FreeT );
      Aux_Message( stdout, "  add particles after restart               = %d\n",     Star_AddParForRestart );
      Aux_Message( stdout, "     number of particles to be added        = %ld\n",    Star_AddParForRestart_NPar );
      Aux_Message( stdout, "     peak DM density for estimating sigma   = %14.7e\n", Star_AddParForRestart_PeakRho );
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "  Tidal_Enabled       = %d\n",            Tidal_Enabled                  );
      if ( Tidal_Enabled ) {
      Aux_Message( stdout, "  Tidal_RotatingFrame = %d\n",            Tidal_RotatingFrame            );
      Aux_Message( stdout, "  Tidal_Mass          = %13.7e Msun\n",   Tidal_Mass*UNIT_M/Const_Msun   );
      Aux_Message( stdout, "  Tidal_R             = %13.7e kpc\n",    Tidal_R*UNIT_L/Const_kpc       );
      Aux_Message( stdout, "  Tidal_Angle0        = %13.7e\n",        Tidal_Angle0                   );
      Aux_Message( stdout, "  Tidal_FixedPos      = %d\n",            Tidal_FixedPos                 );
      Aux_Message( stdout, "  Tidal_Centrifugal   = %d\n",            Tidal_Centrifugal              );
      Aux_Message( stdout, "  Tidal_CutoffR       = %13.7e kpc\n",    Tidal_CutoffR*UNIT_L/Const_kpc );
      Aux_Message( stdout, "  Tidal_Vrot          = %13.7e km/s\n",   Tidal_Vrot*UNIT_V/(Const_km)   ); }
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "  Sponge_Mode         = %d\n",            Sponge_Mode                    );
      if ( Sponge_Mode != 3 ) {
      Aux_Message( stdout, "  Sponge_Width        = %13.7e kpc\n",    Sponge_Width*UNIT_L/Const_kpc  );
      Aux_Message( stdout, "  Sponge_Amp          = %13.7e Gyr^-1\n", Sponge_Amp*Const_Gyr/UNIT_T    ); }
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "  ParFileCM_Enabled   = %d\n",            ParFileCM_Enabled              );
      if ( ParFileCM_Enabled ) {
      Aux_Message( stdout, "  Particle CM displacement x = %14.7e kpc\n",  ParFileCM_Dis[0]*UNIT_L/Const_kpc );
      Aux_Message( stdout, "  Particle CM displacement y = %14.7e kpc\n",  ParFileCM_Dis[1]*UNIT_L/Const_kpc );
      Aux_Message( stdout, "  Particle CM displacement z = %14.7e kpc\n",  ParFileCM_Dis[2]*UNIT_L/Const_kpc );
      Aux_Message( stdout, "  Particle CM velocity x     = %14.7e km/s\n", ParFileCM_Vel[0]*(UNIT_L/UNIT_T)/(Const_km/Const_s) );
      Aux_Message( stdout, "  Particle CM velocity y     = %14.7e km/s\n", ParFileCM_Vel[1]*(UNIT_L/UNIT_T)/(Const_km/Const_s) );
      Aux_Message( stdout, "  Particle CM velocity z     = %14.7e km/s\n", ParFileCM_Vel[2]*(UNIT_L/UNIT_T)/(Const_km/Const_s) ); }
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "  Density recording mode = %d\n", DensRecMode );
      Aux_Message( stdout, "======================================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

   const double Soliton_Center[3] = { amr->BoxCenter[0],
                                      amr->BoxCenter[1],
                                      amr->BoxCenter[2] };
   const double r_tar             = sqrt( SQR(x-Soliton_Center[0]) +
                                          SQR(y-Soliton_Center[1]) +
                                          SQR(z-Soliton_Center[2]) );

   if ( Soliton_InputMode == 1 )
   {
      const double *Table_Radius  = Soliton_DensProf + 0*Soliton_DensProf_NBin;  // radius
      const double *Table_Density = Soliton_DensProf + 1*Soliton_DensProf_NBin;  // density

      double r_ref, dens_ref;

//    rescale radius (target radius --> reference radius)
      r_ref = r_tar / Soliton_ScaleL;

//    linear interpolation
      dens_ref = Mis_InterpolateFromTable( Soliton_DensProf_NBin, Table_Radius, Table_Density, r_ref );

      if ( dens_ref == NULL_REAL )
      {
         if      ( r_ref <  Table_Radius[0] )
            dens_ref = Table_Density[0];

         else if ( r_ref >= Table_Radius[Soliton_DensProf_NBin-1] )
            dens_ref = Table_Density[Soliton_DensProf_NBin-1];

         else
            Aux_Error( ERROR_INFO, "interpolation failed at radius %13.7e (min/max radius = %13.7e/%13.7e) !!\n",
                       r_ref, Table_Radius[0], Table_Radius[Soliton_DensProf_NBin-1] );
      }

//    rescale density (reference density --> target density) and add to the fluid array
      fluid[DENS] = dens_ref*Soliton_ScaleD;
   }

   else if ( Soliton_InputMode == 2 )
   {
      const double m22      = ELBDM_MASS*UNIT_M/(Const_eV/SQR(Const_c))/1.0e-22;
      const double rc_kpc   = Soliton_CoreRadius*UNIT_L/Const_kpc;
      const double peak_rho = 1.945e7/SQR( m22*rc_kpc*rc_kpc )*Const_Msun/CUBE(Const_kpc)/(UNIT_M/CUBE(UNIT_L));

      fluid[DENS] = peak_rho*pow( 1.0+9.06e-2*SQR(r_tar/rc_kpc), Soliton_OuterSlope );
   }

   else if ( Soliton_InputMode == 3 )
      fluid[DENS] = 0.0;

   else
      Aux_Error( ERROR_INFO, "Unsupported Soliton_InputMode (%d) !!\n", Soliton_InputMode );


// set the real and imaginary parts
   fluid[REAL] = sqrt( fluid[DENS] );
   fluid[IMAG] = 0.0;                  // imaginary part is always zero --> no initial velocity

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  End_EridanusII
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_EridanusII()
{

   delete [] Soliton_DensProf;

} // FUNCTION : End_EridanusII



//-------------------------------------------------------------------------------------------------------
// Function    :  BC
// Description :  Set the extenral boundary condition
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//
// Parameter   :  fluid    : Fluid field to be set
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void BC_EridanusII( real fluid[], const double x, const double y, const double z, const double Time,
                    const int lv, double AuxArray[] )
{

   fluid[REAL] = (real)0.0;
   fluid[IMAG] = (real)0.0;
   fluid[DENS] = (real)0.0;

} // FUNCTION : BC



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCenterOfMass
// Description :  Record the center of mass (CM)
//
// Note        :  1. Invoked by Record_EridanusII() recursively
//                2. Only include cells within CM_MaxR from CM_Old[] when updating CM
//
// Parameter   :  CM_Old[] : Previous CM
//                CM_New[] : New CM to be returned
//                CM_MaxR  : Maximum radius to compute CM
//
// Return      :  CM_New[]
//-------------------------------------------------------------------------------------------------------
void GetCenterOfMass( const double CM_Old[], double CM_New[], const double CM_MaxR )
{

   const double CM_MaxR2          = SQR( CM_MaxR );
   const double HalfBox[3]        = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const bool   Periodic          = ( OPT__BC_FLU[0] == BC_FLU_PERIODIC );
   const bool   IntPhase_No       = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const real   MinTemp_No        = -1.0;
   const bool   DE_Consistency_No = false;
#  ifdef PARTICLE
   const bool   TimingSendPar_No  = false;
   const bool   PredictParPos_No  = false;
   const bool   JustCountNPar_No  = false;
#  ifdef LOAD_BALANCE
   const bool   SibBufPatch       = true;
   const bool   FaSibBufPatch     = true;
#  else
   const bool   SibBufPatch       = NULL_BOOL;
   const bool   FaSibBufPatch     = NULL_BOOL;
#  endif
#  endif // #ifdef PARTICLE

   int   *PID0List = NULL;
   double M_ThisRank, MR_ThisRank[3], M_AllRank, MR_AllRank[3];
   real (*TotalDens)[PS1][PS1][PS1];

   M_ThisRank = 0.0;
   for (int d=0; d<3; d++)    MR_ThisRank[d] = 0.0;


   for (int lv=0; lv<NLEVEL; lv++)
   {
//    initialize the particle density array (rho_ext) and collect particles to the target level
#     ifdef PARTICLE
      Prepare_PatchData_InitParticleDensityArray( lv );

      Par_CollectParticle2OneLevel( lv, PredictParPos_No, NULL_REAL, SibBufPatch, FaSibBufPatch, JustCountNPar_No,
                                    TimingSendPar_No );
#     endif

//    get the total density on grids
      TotalDens = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
      PID0List  = new int  [ amr->NPatchComma[lv][1]/8 ];

      for (int PID0=0, t=0; PID0<amr->NPatchComma[lv][1]; PID0+=8, t++)    PID0List[t] = PID0;

      Prepare_PatchData( lv, Time[lv], TotalDens[0][0][0], NULL, 0, amr->NPatchComma[lv][1]/8, PID0List, _TOTAL_DENS, _NONE,
                         OPT__RHO_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                         MinDens_No, MinPres_No, MinTemp_No, DE_Consistency_No );

      delete [] PID0List;


//    free memory for collecting particles from other ranks and levels, and free density arrays with ghost zones (rho_ext)
#     ifdef PARTICLE
      Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch, FaSibBufPatch );

      Prepare_PatchData_FreeParticleDensityArray( lv );
#     endif


//    calculate the center of mass
      const double dh = amr->dh[lv];
      const double dv = CUBE( dh );

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       skip non-leaf patches
         if ( amr->patch[0][lv][PID]->son != -1 )  continue;

         const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

         double x, y, z, dx, dy, dz;

         for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;  dz = z - CM_Old[2];
                                       if ( Periodic ) {
                                          if      ( dz > +HalfBox[2] )  {  z -= amr->BoxSize[2];  dz -= amr->BoxSize[2];  }
                                          else if ( dz < -HalfBox[2] )  {  z += amr->BoxSize[2];  dz += amr->BoxSize[2];  }
                                       }
         for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;  dy = y - CM_Old[1];
                                       if ( Periodic ) {
                                          if      ( dy > +HalfBox[1] )  {  y -= amr->BoxSize[1];  dy -= amr->BoxSize[1];  }
                                          else if ( dy < -HalfBox[1] )  {  y += amr->BoxSize[1];  dy += amr->BoxSize[1];  }
                                       }
         for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;  dx = x - CM_Old[0];
                                       if ( Periodic ) {
                                          if      ( dx > +HalfBox[0] )  {  x -= amr->BoxSize[0];  dx -= amr->BoxSize[0];  }
                                          else if ( dx < -HalfBox[0] )  {  x += amr->BoxSize[0];  dx += amr->BoxSize[0];  }
                                       }

//          only include cells within CM_MaxR
            const double R2 = SQR(dx) + SQR(dy) + SQR(dz);
            if ( R2 < CM_MaxR2 )
            {
               const double dm = TotalDens[PID][k][j][i]*dv;

               M_ThisRank     += dm;
               MR_ThisRank[0] += dm*x;
               MR_ThisRank[1] += dm*y;
               MR_ThisRank[2] += dm*z;
            }
         }}}
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

      delete [] TotalDens;
   } // for (int lv=0; lv<NLEVEL; lv++)


// collect data from all ranks to calculate the CM
// --> note that all ranks will get CM_New[]
   MPI_Allreduce( &M_ThisRank, &M_AllRank, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce( MR_ThisRank, MR_AllRank, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)    CM_New[d] = MR_AllRank[d] / M_AllRank;

// map the new CM back to the simulation domain
   if ( Periodic )
   for (int d=0; d<3; d++)
   {
      if      ( CM_New[d] >= amr->BoxSize[d] )  CM_New[d] -= amr->BoxSize[d];
      else if ( CM_New[d] < 0.0              )  CM_New[d] += amr->BoxSize[d];

   }

   for (int d=0; d<3; d++)
      if ( CM_New[d] >= amr->BoxSize[d]  ||  CM_New[d] < 0.0 )
         Aux_Error( ERROR_INFO, "CM_New[%d] = %14.7e lies outside the domain !!\n", d, CM_New[d] );

} // FUNCTION : GetCenterOfMass



//-------------------------------------------------------------------------------------------------------
// Function    :  Record_EridanusII
// Description :  Record the maximum density and center coordinates
//
// Note        :  1. It will also record the real and imaginary parts associated with the maximum density
//                2. For the center coordinates, it will record the position of maximum density, minimum potential,
//                   and center-of-mass
//                3. Output filename is fixed to "Record__Center"
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Record_EridanusII()
{

   const char filename_center  [] = "Record__Center";
   const int  CountMPI            = 10;

   double dens, max_dens_loc=-__DBL_MAX__, max_dens_pos_loc[3], real_loc, imag_loc;
   double pote, min_pote_loc=+__DBL_MAX__, min_pote_pos_loc[3];
   double send[CountMPI], (*recv)[CountMPI]=new double [MPI_NRank][CountMPI];

   const long   DensMode          = ( DensRecMode == 1 ) ? _DENS :
#                                   ifdef PARTICLE
                                    ( DensRecMode == 2 ) ? _PAR_DENS :
#                                   endif
                                                           _TOTAL_DENS;
   const bool   IntPhase_No       = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const real   MinTemp_No        = -1.0;
   const bool   DE_Consistency_No = false;
#  ifdef PARTICLE
   const bool   TimingSendPar_No  = false;
   const bool   PredictParPos_No  = false;
   const bool   JustCountNPar_No  = false;
#  ifdef LOAD_BALANCE
   const bool   SibBufPatch       = true;
   const bool   FaSibBufPatch     = true;
#  else
   const bool   SibBufPatch       = NULL_BOOL;
   const bool   FaSibBufPatch     = NULL_BOOL;
#  endif
#  endif // #ifdef PARTICLE


// collect local data
   for (int lv=0; lv<NLEVEL; lv++)
   {
//    initialize the particle density array (rho_ext) and collect particles to the target level
#     ifdef PARTICLE
      Prepare_PatchData_InitParticleDensityArray( lv );

      Par_CollectParticle2OneLevel( lv, PredictParPos_No, NULL_REAL, SibBufPatch, FaSibBufPatch, JustCountNPar_No,
                                    TimingSendPar_No );
#     endif

//    get the total density on grids
      real (*TotalDens)[PS1][PS1][PS1] = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
      int   *PID0List                  = new int  [ amr->NPatchComma[lv][1]/8 ];

      for (int PID0=0, t=0; PID0<amr->NPatchComma[lv][1]; PID0+=8, t++)    PID0List[t] = PID0;

      Prepare_PatchData( lv, Time[lv], TotalDens[0][0][0], NULL, 0, amr->NPatchComma[lv][1]/8, PID0List, DensMode, _NONE,
                         OPT__RHO_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                         MinDens_No, MinPres_No, MinTemp_No, DE_Consistency_No );

      delete [] PID0List;


//    free memory for collecting particles from other ranks and levels, and free density arrays with ghost zones (rho_ext)
#     ifdef PARTICLE
      Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch, FaSibBufPatch );

      Prepare_PatchData_FreeParticleDensityArray( lv );
#     endif

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       skip non-leaf patches
         if ( amr->patch[0][lv][PID]->son != -1 )  continue;

         for (int k=0; k<PS1; k++)  {  const double z = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*amr->dh[lv];
         for (int j=0; j<PS1; j++)  {  const double y = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*amr->dh[lv];
         for (int i=0; i<PS1; i++)  {  const double x = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*amr->dh[lv];

//          dens = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
            dens = TotalDens[PID][k][j][i];
            pote = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot[k][j][i];

            if ( dens > max_dens_loc )
            {
               max_dens_loc        = dens;
               real_loc            = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[REAL][k][j][i];
               imag_loc            = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[IMAG][k][j][i];
               max_dens_pos_loc[0] = x;
               max_dens_pos_loc[1] = y;
               max_dens_pos_loc[2] = z;
            }

            if ( pote < min_pote_loc )
            {
               min_pote_loc        = pote;
               min_pote_pos_loc[0] = x;
               min_pote_pos_loc[1] = y;
               min_pote_pos_loc[2] = z;
            }
         }}}
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

      delete [] TotalDens;
   } // for (int lv=0; lv<NLEVEL; lv++)


// gather data to the root rank
   send[0] = max_dens_loc;
   send[1] = real_loc;
   send[2] = imag_loc;
   send[3] = max_dens_pos_loc[0];
   send[4] = max_dens_pos_loc[1];
   send[5] = max_dens_pos_loc[2];
   send[6] = min_pote_loc;
   send[7] = min_pote_pos_loc[0];
   send[8] = min_pote_pos_loc[1];
   send[9] = min_pote_pos_loc[2];

   MPI_Gather( send, CountMPI, MPI_DOUBLE, recv[0], CountMPI, MPI_DOUBLE, 0, MPI_COMM_WORLD );


// record the maximum density and center coordinates
   double max_dens      = -__DBL_MAX__;
   double min_pote      = +__DBL_MAX__;
   int    max_dens_rank = -1;
   int    min_pote_rank = -1;

   if ( MPI_Rank == 0 )
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         if ( recv[r][0] > max_dens )
         {
            max_dens      = recv[r][0];
            max_dens_rank = r;
         }

         if ( recv[r][6] < min_pote )
         {
            min_pote      = recv[r][6];
            min_pote_rank = r;
         }
      }

      if ( max_dens_rank < 0  ||  max_dens_rank >= MPI_NRank )
         Aux_Error( ERROR_INFO, "incorrect max_dens_rank (%d) !!\n", max_dens_rank );

      if ( min_pote_rank < 0  ||  min_pote_rank >= MPI_NRank )
         Aux_Error( ERROR_INFO, "incorrect min_pote_rank (%d) !!\n", min_pote_rank );

      static bool FirstTime = true;

      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(filename_center) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_center );
         else
         {
            FILE *file_center = fopen( filename_center, "w" );
            fprintf( file_center, "#%19s  %10s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %10s  %14s  %14s  %14s\n",
                     "Time", "Step", "Dens", "Real", "Imag", "Dens_x", "Dens_y", "Dens_z", "Pote", "Pote_x", "Pote_y", "Pote_z",
                     "NIter", "CM_x", "CM_y", "CM_z" );
            fclose( file_center );
         }

         FirstTime = false;
      }

      FILE *file_center = fopen( filename_center, "a" );
      fprintf( file_center, "%20.14e  %10ld  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e",
               Time[0], Step, recv[max_dens_rank][0], recv[max_dens_rank][1], recv[max_dens_rank][2], recv[max_dens_rank][3],
                              recv[max_dens_rank][4], recv[max_dens_rank][5], recv[min_pote_rank][6], recv[min_pote_rank][7],
                              recv[min_pote_rank][8], recv[min_pote_rank][9] );
      fclose( file_center );
   } // if ( MPI_Rank == 0 )


// compute the center of mass until convergence
   const double TolErrR2 = SQR( Soliton_CM_TolErrR );
   const int    NIterMax = 10;

   double dR2, CM_Old[3], CM_New[3];
   int NIter = 0;

// set an initial guess by the peak density position
   if ( MPI_Rank == 0 )
      for (int d=0; d<3; d++)    CM_Old[d] = recv[max_dens_rank][3+d];

   MPI_Bcast( CM_Old, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   while ( true )
   {
      GetCenterOfMass( CM_Old, CM_New, Soliton_CM_MaxR );

      dR2 = SQR( CM_Old[0] - CM_New[0] )
          + SQR( CM_Old[1] - CM_New[1] )
          + SQR( CM_Old[2] - CM_New[2] );
      NIter ++;

      if ( dR2 <= TolErrR2  ||  NIter >= NIterMax )
         break;
      else
         memcpy( CM_Old, CM_New, sizeof(double)*3 );
   }

   if ( MPI_Rank == 0 )
   {
      if ( dR2 > TolErrR2 )
         Aux_Message( stderr, "WARNING : dR (%13.7e) > Soliton_CM_TolErrR (%13.7e) !!\n", sqrt(dR2), Soliton_CM_TolErrR );

      FILE *file_center = fopen( filename_center, "a" );
      fprintf( file_center, "  %10d  %14.7e  %14.7e  %14.7e\n", NIter, CM_New[0], CM_New[1], CM_New[2] );
      fclose( file_center );
   }


   for (int d=0; d<3; d++)    Tidal_CM[d] = CM_New[d];


   delete [] recv;

} // FUNCTION : Record_EridanusII



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExtPotAuxArray_EridanusII
// Description :  Set the auxiliary array ExtPot_AuxArray[] used by ExtPot_EridanusII()
//
// Note        :  1. To adopt this routine, link to the function pointer "Init_ExtPotAuxArray_Ptr"
//                   in a test problem initializer as follows:
//
//                      void Init_ExtPotAuxArray_EridanusII( double AuxArray[] );
//
//                      ...
//
//                      Init_ExtPotAuxArray_Ptr = Init_ExtPotAuxArray_EridanusII;
//
//                   --> Then it will be invoked by Init_ExtAccPot()
//                2. AuxArray[] has the size of EXT_POT_NAUX_MAX defined in Macro.h (default = 10)
//
// Parameter   :  AuxArray : Array to be filled up
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
void Init_ExtPotAuxArray_EridanusII( double AuxArray[] )
{

// ExtPot_AuxArray has the size of EXT_POT_NAUX_MAX (default = 10)
   if ( Tidal_RotatingFrame )
   {
      AuxArray[0] = Tidal_CM[0];
      AuxArray[1] = Tidal_CM[1];
      AuxArray[2] = Tidal_CM[2];
   }

   else
   {
      AuxArray[0] = amr->BoxCenter[0];
      AuxArray[1] = amr->BoxCenter[1];
      AuxArray[2] = amr->BoxCenter[2];
   }

   AuxArray[3] = NEWTON_G*Tidal_Mass;
   AuxArray[4] = Tidal_R;
   AuxArray[5] = Tidal_Vrot;
   AuxArray[6] = ( Tidal_FixedPos ) ? +1.0 : -1.0;
   AuxArray[7] = ( Tidal_Centrifugal ) ? +1.0 : -1.0;
   AuxArray[8] = Tidal_Angle0;
   AuxArray[9] = ( Tidal_RotatingFrame ) ? +1.0 : -1.0;

} // FUNCTION : Init_ExtPotAuxArray_EridanusII



//-------------------------------------------------------------------------------------------------------
// Function    :  Reset
// Description :  Reset the wave function outside a specific sphere to be zero
//
// Note        :  1. Linked to the function pointer "Flu_ResetByUser_Func_Ptr"
//
// Parameter   :  fluid    : Fluid array storing both the input (origial) and reset values
//                           --> Including both active and passive variables
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  true  : This cell has been reset
//                false : This cell has not been reset
//-------------------------------------------------------------------------------------------------------
bool Reset( real fluid[], const double x, const double y, const double z, const double Time,
            const int lv, double AuxArray[] )
{

   const real dr[3] = { x-Tidal_CM[0], y-Tidal_CM[1], z-Tidal_CM[2] };
   const real r     = SQRT( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );

// sponge BC
   if      ( Sponge_Mode == 1 )
   {
      const double v    = 0.5*Sponge_Amp*(  1.0 + tanh( (r-Tidal_CutoffR)/Sponge_Width )  );
      const double damp = exp( -v*Sponge_dt );

      fluid[REAL] *= damp;
      fluid[IMAG] *= damp;
      fluid[DENS]  = SQR( fluid[REAL] ) + SQR( fluid[IMAG] );

      return true;
   }

// truncation BC
   else if ( Sponge_Mode == 2 )
   {
      const real WaveFloor = 1.0e-3;
      const real DensFloor = SQR(WaveFloor);

      if ( r > Tidal_CutoffR )
      {
         fluid[REAL] = WaveFloor;
         fluid[IMAG] = (real)0.0;
         fluid[DENS] = DensFloor;

         return true;
      }

      else
         return false;
   }

// nothing
   else
   {
      return false;
   }

} // FUNCTION : Reset



//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_UserWorkBeforePoisson_EridanusII
// Description :  User-specified work before invoking the Poisson solver
//
// Note        :  1. Invoked by Gra_AdvanceDt() using the function pointer "Poi_UserWorkBeforePoisson_Ptr"
//
// Parameter   :  Time : Target physical time
//                lv   : Target refinement level
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Poi_UserWorkBeforePoisson_EridanusII( const double Time, const int lv )
{

   if ( OPT__EXT_POT )
   {
      SetExtPotAuxArray_EridanusII( ExtPot_AuxArray_Flt, ExtPot_AuxArray_Int );

#     ifdef GPU
      CUAPI_SetConstMemory_ExtAccPot();
#     endif
   }

} // FUNCTION : Poi_UserWorkBeforePoisson_EridanusII
#endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_EridanusII
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_EridanusII()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM  &&  defined GRAVITY )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr        = SetGridIC;
#  ifdef PARTICLE
   Init_User_Ptr                 = Init_User_EridanusII;
#  else
   Init_User_Ptr                 = NULL;
#  endif
   Flag_User_Ptr                 = NULL;
   Mis_GetTimeStep_User_Ptr      = NULL;
   BC_User_Ptr                   = BC_EridanusII;
   Flu_ResetByUser_Func_Ptr      = Reset;
   Output_User_Ptr               = NULL;
   Aux_Record_User_Ptr           = Record_EridanusII;
   End_User_Ptr                  = End_EridanusII;
#  ifdef GRAVITY
   Init_ExtPot_Ptr               = Init_ExtPot_EridanusII;
   Poi_UserWorkBeforePoisson_Ptr = Poi_UserWorkBeforePoisson_EridanusII;
#  endif // #ifdef GRAVITY
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr       = Par_Init_ByFunction_EridanusII;
#  endif
#  endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_EridanusII
