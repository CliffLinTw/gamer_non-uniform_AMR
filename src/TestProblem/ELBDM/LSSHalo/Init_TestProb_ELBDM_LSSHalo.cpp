#include "GAMER.h"
#include "TestProb.h"

// problem-specific global variables
// =======================================================================================

static FieldIdx_t Su_Idx_Dens0 = Idx_Undefined;  // field index for storing the **initial** density

static double Disc_Cen[3];
static double Disc_BulkVel[3];
static double Disc_Decay_R;
static double Disc_Radius;
static double Disc_Mass;
static int    Disc_RSeed;
static double Disc_CM_MaxR;      // maximum radius for determining CM
static double Disc_CM_TolErrR;   // maximum allowed errors for determining CM
static double VelDisp;
static int Idx_ParLabel = Idx_Undefined;
bool FixDM;
bool OutputWaveFunction;
bool AddParWhenRestart;
int AddParNPar;
int AddParNRing;

static RandomNumber_t *RNG = NULL;
static double pi = 3.1415926;
static double G = 6.67408E-8;
static void   Vec2_FixRadius ( const double r, double RanVec[], double NormVec[], const double RanV);
static double Disc_Interpolation( const double RanM, const double R, const double DiscMConst);
static void   Maxwell_Interpolation( double VD_Vec[], const double v_rms );
static double myErfInv2(double x);

# ifdef PARTICLE

static void Par_Disc_Heating(const long NPar_ThisRank, const long NPar_AllRank,
                                 real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                 real *ParVelX, real *ParVelY, real *ParVelZ,
                                 real *ParTime, real *AllAttribute[PAR_NATT_TOTAL]);


static void Par_AfterAcceleration( const long NPar_ThisRank, const long NPar_AllRank,
                                 real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                 real *ParVelX, real *ParVelY, real *ParVelZ,
                                 real *ParAccX, real *ParAccY, real *ParAccZ,
                                 real *ParTime, real *AllAttribute[PAR_NATT_TOTAL] );

// this function pointer may be overwritten by various test problem initializers

void (*Par_Disc_Heating_Ptr)(const long NPar_ThisRank, const long NPar_AllRank,
                                 real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                 real *ParVelX, real *ParVelY, real *ParVelZ,
                                 real *ParTime, real *AllAttribute[PAR_NATT_TOTAL])= Par_Disc_Heating;

void (*Par_AfterAcceleration_Ptr)( const long NPar_ThisRank, const long NPar_AllRank,
                                 real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                 real *ParVelX, real *ParVelY, real *ParVelZ,
                                 real *ParAccX, real *ParAccY, real *ParAccZ,
                                 real *ParTime, real *AllAttribute[PAR_NATT_TOTAL] ) = Par_AfterAcceleration;
# endif //ifdef PARTICLE

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
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n");
#  endif

#  ifndef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n");
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
// ReadPara->Add( "KEY_IN_THE_FILE",      &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************


   ReadPara->Add( "Disc_Cen_X",           &Disc_Cen[0],           1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_Cen_Y",           &Disc_Cen[1],           1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_Cen_Z",           &Disc_Cen[2],           1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_BulkVel_X",       &Disc_BulkVel[0],       0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_BulkVel_Y",       &Disc_BulkVel[1],       0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_BulkVel_Z",       &Disc_BulkVel[2],       0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_Radius",          &Disc_Radius,           1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_Decay_Radius",    &Disc_Decay_R,          1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_Mass",            &Disc_Mass,             1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Disc_Random_Seed",     &Disc_RSeed,            123,             0,              NoMax_int         );
   ReadPara->Add( "Disc_CM_MaxR",         &Disc_CM_MaxR,         -1.0,    Eps_double,              NoMax_double      );
   ReadPara->Add( "Disc_CM_TolErrR",      &Disc_CM_TolErrR,      -1.0,  NoMin_double,              NoMax_double      );
   ReadPara->Add( "Velocity_Dispersion",  &VelDisp,               0.0,           0.0,              1.0               );
   ReadPara->Add( "Fix_DM",               &FixDM,                true,  Useless_bool,              Useless_bool      );
   ReadPara->Add( "Output_Wave_Function", &OutputWaveFunction,   true,  Useless_bool,              Useless_bool      );
   ReadPara->Add( "Add_Par_When_Restart", &AddParWhenRestart,    true,  Useless_bool,              Useless_bool      );
   ReadPara->Add( "Add_Par_NPar",         &AddParNPar,        2000000,             0,              NoMax_int         );
   ReadPara->Add( "Add_Par_NRing",        &AddParNRing,          1000,             1,              NoMax_int         );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   if ( Disc_CM_TolErrR < 0.0 )  Disc_CM_TolErrR = 1.0*amr->dh[MAX_LEVEL];

// (1-3) check the runtime parameters
   if ( OPT__INIT == INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "OPT__INIT=1 is not supported for this test problem !!\n" );


// (2) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__;   // ~7 Gyr

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n"  );
      Aux_Message( stdout, "  test problem ID             = %d\n"                        , TESTPROB_ID        );
      Aux_Message( stdout, "  disc center                 = %13.7e\n, %13.7e\n, %13.7e\n", Disc_Cen[0],
                                                                                           Disc_Cen[1],
                                                                                           Disc_Cen[2]        );
      Aux_Message( stdout, "  disc bulk velocity          = %13.7e\n, %13.7e\n, %13.7e\n", Disc_BulkVel[0],
                                                                                           Disc_BulkVel[1],
                                                                                           Disc_BulkVel[2]    );
      Aux_Message( stdout, "  disc decay radius           = %13.7e\n",                     Disc_Decay_R       );
      Aux_Message( stdout, "  disc radius                 = %13.7e\n",                     Disc_Radius        );
      Aux_Message( stdout, "  disc mass                   = %13.7e\n",                     Disc_Mass          );
      Aux_Message( stdout, "  disc random seed            = %d\n",                         Disc_RSeed         );
      Aux_Message( stdout, "  CM max radius               = %13.6e\n",                     Disc_CM_MaxR       );
      Aux_Message( stdout, "  CM tolerated error          = %13.6e\n",                     Disc_CM_TolErrR    );
      Aux_Message( stdout, "  velocity dispersion         = %13.7e\n",                     VelDisp            );
      Aux_Message( stdout, "  fix DM                      = %d\n",                         FixDM              );
      Aux_Message( stdout, "  output wavefunction         = %d\n",                         OutputWaveFunction );
      Aux_Message( stdout, "  add particles when restart  = %d\n",                         AddParWhenRestart  );
      Aux_Message( stdout, "  number of particles         = %d\n",                         AddParNPar         );
      Aux_Message( stdout, "  number of rings             = %d\n",                         AddParNRing        );

      Aux_Message( stdout, "=============================================================================\n"  );
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

   Aux_Error( ERROR_INFO, "OPT__INIT=1 is not supported for this test problem !!\n" );

} // FUNCTION : SetGridIC

void BC_Disc_Heating( real fluid[], const double x, const double y, const double z, const double Time,
         const int lv, double AuxArray[] )
{

   fluid[REAL] = (real)0.0;
   fluid[IMAG] = (real)0.0;
   fluid[DENS] = (real)0.0;

} // FUNCTION : BC_Disc_Heating

//-------------------------------------------------------------------------------------------------------
//// Function    :  Init_Disc
//// Description :  Set the initial condition ogf velicity of particles by their acceleration and position
////
//// Note        :  1. Linked to the function pointer "Init_User_Ptr"
////
//// Parameter   :  None
////
//// Return      :  None
////-------------------------------------------------------------------------------------------------------
void GetCenterOfMass_Disc_Heating( const double CM_Old[], double CM_New[], const double CM_MaxR )
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

      Par_CollectParticle2OneLevel( lv, _PAR_MASS|_PAR_POSX|_PAR_POSY|_PAR_POSZ ,PredictParPos_No, NULL_REAL, SibBufPatch, FaSibBufPatch, JustCountNPar_No,
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

} // FUNCTION : GetCenterOfMass_Disc_Heating

/*

//-------------------------------------------------------------------------------------------------------
// Function    :  GetCenterOfMass_Disc_Heating
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


void GetCenterOfMass_Disc_Heating( const double CM_Old[], double CM_New[], const double CM_MaxR )
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

//      Prepare_PatchData( lv, Time[lv], TotalDens[0][0][0], NULL, 0, amr->NPatchComma[lv][1]/8, PID0List, _TOTAL_DENS, _NONE,
//                         OPT__RHO_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
//                         MinDens_No, MinPres_No, MinTemp_No, DE_Consistency_No );

      Prepare_PatchData( lv, Time[lv], TotalDens[0][0][0], 0, amr->NPatchComma[lv][1]/8, PID0List, _TOTAL_DENS,
                         OPT__RHO_INT_SCHEME, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                         MinDens_No, MinPres_No, DE_Consistency_No );

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

} // FUNCTION : GetCenterOfMass_Disc_Heating


*/

//-------------------------------------------------------------------------------------------------------
// Function    :  Record_Disc_Heating
// Description :  Record the maximum density and center coordinates
//
// Note        :  1. It will also record the real and imaginary parts associated with the maximum density
//                2. For the center coordinates, it will record the position of maximum density, minimum potential,
//                   and center-of-mass
//                3. Output filenames are fixed to "Record__MaxDens" and "Record__Center"
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------

/*

void Record_Disc_Heating()
{

   const char filename_max_dens[] = "Record__MaxDens";
   const char filename_center  [] = "Record__Center";
   const int  CountMPI            = 10;

   double dens, max_dens_loc=-__DBL_MAX__, max_dens_pos_loc[3], real_loc, imag_loc;
   double pote, min_pote_loc=+__DBL_MAX__, min_pote_pos_loc[3];
   double send[CountMPI], (*recv)[CountMPI]=new double [MPI_NRank][CountMPI];


// collect local data
   for (int lv=0; lv<NLEVEL; lv++)
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
//    skip non-leaf patches
      if ( amr->patch[0][lv][PID]->son != -1 )  continue;

      for (int k=0; k<PS1; k++)  {  const double z = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*amr->dh[lv];
      for (int j=0; j<PS1; j++)  {  const double y = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*amr->dh[lv];
      for (int i=0; i<PS1; i++)  {  const double x = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*amr->dh[lv];

         dens = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
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
   }


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
         if ( Aux_CheckFileExist(filename_max_dens) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_max_dens );
         else
         {
            FILE *file_max_dens = fopen( filename_max_dens, "w" );
            fprintf( file_max_dens, "#%19s   %10s   %14s   %14s   %14s\n", "Time", "Step", "Dens", "Real", "Imag" );
            fclose( file_max_dens );
         }

         if ( Aux_CheckFileExist(filename_center) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_center );
         else
         {
            FILE *file_center = fopen( filename_center, "w" );
            fprintf( file_center, "#%19s  %10s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %10s  %14s  %14s  %14s\n",
                     "Time", "Step", "Dens", "Dens_x", "Dens_y", "Dens_z", "Pote", "Pote_x", "Pote_y", "Pote_z",
                     "NIter", "CM_x", "CM_y", "CM_z" );
            fclose( file_center );
         }

         FirstTime = false;
      }

      FILE *file_max_dens = fopen( filename_max_dens, "a" );
      fprintf( file_max_dens, "%20.14e   %10ld   %14.7e   %14.7e   %14.7e\n",
               Time[0], Step, recv[max_dens_rank][0], recv[max_dens_rank][1], recv[max_dens_rank][2] );
      fclose( file_max_dens );

      FILE *file_center = fopen( filename_center, "a" );
      fprintf( file_center, "%20.14e  %10ld  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e",
               Time[0], Step, recv[max_dens_rank][0], recv[max_dens_rank][3], recv[max_dens_rank][4], recv[max_dens_rank][5],
                              recv[min_pote_rank][6], recv[min_pote_rank][7], recv[min_pote_rank][8], recv[min_pote_rank][9] );
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
      GetCenterOfMass_Disc_Heating( CM_Old, CM_New, Soliton_CM_MaxR );

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

} // FUNCTION : Record_Disc_Heating

*/

void Record_Disc_Heating()
{

   const char filename_max_dens[] = "Record__MaxDens";
   const char filename_center  [] = "Record__Center";
   const int  CountMPI            = 10;

   double dens, max_dens_loc=-__DBL_MAX__, max_dens_pos_loc[3], real_loc, imag_loc;
   double pote, min_pote_loc=+__DBL_MAX__, min_pote_pos_loc[3];
   double send[CountMPI], (*recv)[CountMPI]=new double [MPI_NRank][CountMPI];


// collect local data
   for (int lv=0; lv<NLEVEL; lv++)
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
//    skip non-leaf patches
      if ( amr->patch[0][lv][PID]->son != -1 )  continue;

      for (int k=0; k<PS1; k++)  {  const double z = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*amr->dh[lv];
      for (int j=0; j<PS1; j++)  {  const double y = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*amr->dh[lv];
      for (int i=0; i<PS1; i++)  {  const double x = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*amr->dh[lv];

         dens = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
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
   }


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
         if ( Aux_CheckFileExist(filename_max_dens) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_max_dens );
         else
         {
            FILE *file_max_dens = fopen( filename_max_dens, "w" );
            fprintf( file_max_dens, "#%19s   %10s   %14s   %14s   %14s\n", "Time", "Step", "Dens", "Real", "Imag" );
            fclose( file_max_dens );
         }

         if ( Aux_CheckFileExist(filename_center) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_center );
         else
         {
            FILE *file_center = fopen( filename_center, "w" );
            fprintf( file_center, "#%19s  %10s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %10s  %14s  %14s  %14s\n",
                     "Time", "Step", "Dens", "Dens_x", "Dens_y", "Dens_z", "Pote", "Pote_x", "Pote_y", "Pote_z",
                     "NIter", "CM_x", "CM_y", "CM_z" );
            fclose( file_center );
         }

         FirstTime = false;
      }

      FILE *file_max_dens = fopen( filename_max_dens, "a" );
      fprintf( file_max_dens, "%20.14e   %10ld   %14.7e   %14.7e   %14.7e\n",
               Time[0], Step, recv[max_dens_rank][0], recv[max_dens_rank][1], recv[max_dens_rank][2] );
      fclose( file_max_dens );

      FILE *file_center = fopen( filename_center, "a" );
      fprintf( file_center, "%20.14e  %10ld  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e",
               Time[0], Step, recv[max_dens_rank][0], recv[max_dens_rank][3], recv[max_dens_rank][4], recv[max_dens_rank][5],
                              recv[min_pote_rank][6], recv[min_pote_rank][7], recv[min_pote_rank][8], recv[min_pote_rank][9] );
      fclose( file_center );
   } // if ( MPI_Rank == 0 )

   // compute the center of mass until convergence
      const double TolErrR2 = SQR( Disc_CM_TolErrR );
      const int    NIterMax = 10;

      double dR2, CM_Old[3], CM_New[3];
      int NIter = 0;

   // set an initial guess by the peak density position
      if ( MPI_Rank == 0 )
         for (int d=0; d<3; d++)    CM_Old[d] = recv[max_dens_rank][3+d];

      MPI_Bcast( CM_Old, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );

      while ( true )
      {
         GetCenterOfMass_Disc_Heating( CM_Old, CM_New, Disc_CM_MaxR );

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
            Aux_Message( stderr, "WARNING : dR (%13.7e) > Disc_CM_TolErrR (%13.7e) !!\n", sqrt(dR2), Disc_CM_TolErrR );

         FILE *file_center = fopen( filename_center, "a" );
         fprintf( file_center, "  %10d  %14.7e  %14.7e  %14.7e\n", NIter, CM_New[0], CM_New[1], CM_New[2] );
         fclose( file_center );
      }


   delete [] recv;

} // FUNCTION : Record_Disc_Heating



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Disc_Heating
// Description :  Particle IC
//
// Note        :  None
//
// Parameter   :  NPar_ThisRank, NPar_AllRank,*ParMass, *ParPosX, *ParPosY, *ParPosZ, *ParVelX, *ParVelY,
//                *ParVelZ, *ParTime, *AllAttribute[PAR_NATT_TOTAL]
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------

#  ifdef PARTICLE

void Par_Disc_Heating(const long NPar_ThisRank, const long NPar_AllRank, real *ParMass, real *ParPosX,
                                            real *ParPosY,      real *ParPosZ,           real *ParVelX, real *ParVelY,
                                            real *ParVelZ,      real *ParTime,           real *AllAttribute[PAR_NATT_TOTAL])
{

   if ( MPI_Rank == 0 ) Aux_Message( stdout, "%s ...\n", __FUNCTION__);

   real *Mass_AllRank = NULL;
   real *Pos_AllRank[3]  = {NULL, NULL, NULL};
   real *Vel_AllRank[3]  = {NULL, NULL, NULL};

   if ( MPI_Rank == 0 )
   {
      const double ParM = Disc_Mass / NPar_AllRank;
      double Ran, RanR, RanM, RanV, RanVec[3], NormVec[3];
      Aux_Message(stdout, " Particle Mass = %13.7e\n", ParM) ;
      Mass_AllRank = new real [NPar_AllRank];
      for (int d = 0; d < 3; d++)
      {
         Pos_AllRank[d] = new real [NPar_AllRank];
         Vel_AllRank[d] = new real [NPar_AllRank];
      } //for (int d = 0; d < 3; d++)

//   initialize the RNG
     RNG = new RandomNumber_t( 1 );
     RNG->SetSeed( 0, Disc_RSeed );

     const double DiscMConst = Disc_Mass / (Disc_Decay_R - (Disc_Decay_R + Disc_Radius) * exp( - Disc_Radius/Disc_Decay_R ) );

     for ( long p = 0; p < NPar_AllRank; p++)
     {
//      mass
        Mass_AllRank[p] = ParM;

//      position
        Ran  = RNG->GetValue( 0, 0.0, 1.0);
        RanM = Ran*Disc_Mass;
        RanR = Disc_Interpolation( RanM, Disc_Decay_R, DiscMConst );
        RanV = sqrt(G*RanM/RanR);
        Vec2_FixRadius( RanR, RanVec, NormVec, RanV );
        for (int d = 0; d < 3; d++) Pos_AllRank[d][p] = RanVec[d] + Disc_Cen[d];

//        if (RanVec[0] < 0)
//           Mass_AllRank[p] = ParM/100.0;

//      velocity
        for (int d = 0; d< 3; d++) Vel_AllRank[d][p] = 0;   //NormVec[d]+Disc_BulkVel[d];
     } // for ( long p = 0; p < NPar_AllRank; p++)
     Aux_Message( stdout, " Particle mass              = %13.7e\n", ParM );

   } //if ( MPI_Rank == 0 )
 // synchronize all particles to the physical time on the base level
   for (long p = 0; p<NPar_ThisRank; p++) ParTime[p] = Time[0];
    //get the number of particles in each rank and set the corresponding offsets
      if ( NPar_AllRank > (long)__INT_MAX__ )
         Aux_Error( ERROR_INFO, "NPar_Active_AllRank (%ld) exceeds the maximum integer (%ld) --> MPI will likely fail !!\n",
                    NPar_AllRank, (long)__INT_MAX__ );

      int NSend[MPI_NRank], SendDisp[MPI_NRank];
      int NPar_ThisRank_int = NPar_ThisRank;    // (i) convert to "int" and (ii) remove the "const" declaration
                                                // --> (ii) is necessary for OpenMPI version < 1.7

      MPI_Gather( &NPar_ThisRank_int, 1, MPI_INT, NSend, 1, MPI_INT, 0, MPI_COMM_WORLD );

      if ( MPI_Rank == 0 )
      {
         SendDisp[0] = 0;
         for (int r=1; r<MPI_NRank; r++)  SendDisp[r] = SendDisp[r-1] + NSend[r-1];
      }


 //  send particle attributes from the master rank to all ranks
      real *Mass   =   ParMass;
      real *Pos[3] = { ParPosX, ParPosY, ParPosZ };
      real *Vel[3] = { ParVelX, ParVelY, ParVelZ };

   #  ifdef FLOAT8
      MPI_Scatterv( Mass_AllRank, NSend, SendDisp, MPI_DOUBLE, Mass, NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );

      for (int d=0; d<3; d++)
      {
         MPI_Scatterv( Pos_AllRank[d], NSend, SendDisp, MPI_DOUBLE, Pos[d], NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );
         MPI_Scatterv( Vel_AllRank[d], NSend, SendDisp, MPI_DOUBLE, Vel[d], NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      }

   #  else
      MPI_Scatterv( Mass_AllRank, NSend, SendDisp, MPI_FLOAT,  Mass, NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );

      for (int d=0; d<3; d++)
      {
         MPI_Scatterv( Pos_AllRank[d], NSend, SendDisp, MPI_FLOAT,  Pos[d], NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );
         MPI_Scatterv( Vel_AllRank[d], NSend, SendDisp, MPI_FLOAT,  Vel[d], NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );
      }
   #  endif  //ifdef FLOAT8


      if ( MPI_Rank == 0 )
      {
         delete RNG;
         delete [] Mass_AllRank;

         for (int d=0; d<3; d++)
         {
            delete [] Pos_AllRank[d];
            delete [] Vel_AllRank[d];
         }
      }

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );
}

void Par_AfterAcceleration(const long NPar_ThisRank, const long NPar_AllRank, real *ParMass,
                                            real *ParPosX, real *ParPosY, real *ParPosZ,
                                            real *ParVelX, real *ParVelY, real *ParVelZ,
                                            real *ParAccX, real *ParAccY, real *ParAccZ,
                                            real *ParTime,  real *AllAttribute[PAR_NATT_TOTAL])
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// set other particle attributes
// ============================================================================================================
   real *ParPos[3] = { ParPosX, ParPosY, ParPosZ };
   real *ParVel[3] = { ParVelX, ParVelY, ParVelZ };
   real *ParAcc[3] = { ParAccX, ParAccY, ParAccZ };


   real ParRadius[2];
   real NormParRadius[2];
   double V_acc;
   double Ran[3];
   
   double *V_acc_total      = new double[2*AddParNRing];
   double *V_acc_total_recv = new double[2*AddParNRing];

  // double V_acc_total[2*AddParNRing], V_acc_total_recv[2*AddParNRing];
   //int RingID_list[NPar_ThisRank];
   int RingID;
   double Ran0, Ran1, Ran2, VD_Vec[3];


   int Count_Outliers = 0; 
//   const int TargetRank = (MPI_Rank+1)%2;
 //  const int Tag        = 123;  

  //   initialize the RNG
   RNG = new RandomNumber_t( 1 );
   RNG->SetSeed( 0, Disc_RSeed+100 );

   for (long p=0; p<2*AddParNRing; p++)
   {
      V_acc_total[p] = 0.0;
      V_acc_total_recv[p] = 0.0;
   }
   

   //double VD_Vec[3];

   for (long p=0; p<NPar_ThisRank; p++)
   {
      if ( ParMass[p] < 0.0 )  continue;
      ParRadius[0] =  ParPos[0][p]-Disc_Cen[0];
      ParRadius[1] =  ParPos[1][p]-Disc_Cen[1];

      RingID = int(pow(ParRadius[0]*ParRadius[0] + ParRadius[1]*ParRadius[1], 0.5)/(Disc_Radius/AddParNRing));
    //  RingID = 0;

    //  NormParRadius[0] =   ParRadius[0]/ sqrt( SQR(ParRadius[0]) + SQR(ParRadius[1]) );
    //  NormParRadius[1] =   ParRadius[1]/ sqrt( SQR(ParRadius[0]) + SQR(ParRadius[1]) );

      V_acc = sqrt(fabs(ParRadius[0]*ParAcc[0][p]+ParRadius[1]*ParAcc[1][p]));
      V_acc_total[RingID*2+0] += 1.0;
      V_acc_total[RingID*2+1] += V_acc;
   }
   



//   MPI_Allreduce( V_acc_total, V_acc_total_recv, 2*AddParNRing, MPI_DOUBLE, MPI_SUM , MPI_COMM_WORLD);

 //  if ( MPI_Rank == 0 )
 //  {
 //     MPI_Send( &V_acc_total,       Count, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD );
 //     MPI_Recv( &V_acc_total_recv,  Count, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE );

 //  }
 //  else
 //  {
 //     MPI_Recv( &V_acc_total_recv,  Count, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE );
 //     MPI_Send( &V_acc_total,       Count, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD );
 //  }
    

   //MPI_Finalize();
    MPI_Allreduce( V_acc_total, V_acc_total_recv, 2*AddParNRing, MPI_DOUBLE, MPI_SUM , MPI_COMM_WORLD);

    if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n"  );
      Aux_Message( stdout, "  Number of Particles = %13.7e\n, %13.7e\n, %13.7e\n,  %13.7e\n,  %13.7e\n,  %13.7e\n,  %13.7e\n,  %13.7e\n,  %13.7e\n,  %13.7e\n,  %13.7e\n,  %13.7e\n",
                                                                                           V_acc_total_recv[0],
                                                                                           V_acc_total_recv[2],
                                                                                           V_acc_total_recv[4],
                                                                                           V_acc_total_recv[6],
                                                                                           V_acc_total_recv[8],
                                                                                           V_acc_total_recv[10],
                                                                                           V_acc_total_recv[1988],
                                                                                           V_acc_total_recv[1990],
                                                                                           V_acc_total_recv[1992],
                                                                                           V_acc_total_recv[1994],
                                                                                           V_acc_total_recv[1996],
                                                                                           V_acc_total_recv[1998]     );
      Aux_Message( stdout, "  Ring sum            = %13.7e\n, %13.7e\n, %13.7e\n,  %13.7e\n,  %13.7e\n,  %13.7e\n,  %13.7e\n,  %13.7e\n,  %13.7e\n,  %13.7e\n,  %13.7e\n,  %13.7e\n", 
                                                                                           V_acc_total_recv[1],
                                                                                           V_acc_total_recv[3],
                                                                                           V_acc_total_recv[5],
                                                                                           V_acc_total_recv[7],
                                                                                           V_acc_total_recv[9],
                                                                                           V_acc_total_recv[11],
                                                                                           V_acc_total_recv[1989],
                                                                                           V_acc_total_recv[1991],
                                                                                           V_acc_total_recv[1993],
                                                                                           V_acc_total_recv[1995],
                                                                                           V_acc_total_recv[1997],
                                                                                           V_acc_total_recv[1999]     );

      Aux_Message( stdout, "=============================================================================\n"  );
   }

   MPI_Allreduce( V_acc_total, V_acc_total_recv, 2*AddParNRing, MPI_DOUBLE, MPI_SUM , MPI_COMM_WORLD);


   for (long p=0; p<NPar_ThisRank; p++)
   {
      if ( ParMass[p] < 0.0 )  continue;
      ParRadius[0] =  ParPos[0][p]-Disc_Cen[0];
      ParRadius[1] =  ParPos[1][p]-Disc_Cen[1];

      RingID = int(pow(ParRadius[0]*ParRadius[0] + ParRadius[1]*ParRadius[1], 0.5)/(Disc_Radius/AddParNRing));
    //  RingID = 0;
//      NormParRadius[0] =   ParRadius[0]/ sqrt( SQR(ParRadius[0]) + SQR(ParRadius[1]) );
//      NormParRadius[1] =   ParRadius[1]/ sqrt( SQR(ParRadius[0]) + SQR(ParRadius[1]) );

      NormParRadius[0] =   ParRadius[0]/ pow(ParRadius[0]*ParRadius[0] + ParRadius[1]*ParRadius[1], 0.5);
      NormParRadius[1] =   ParRadius[1]/ pow(ParRadius[0]*ParRadius[0] + ParRadius[1]*ParRadius[1], 0.5);

      //if (p%100000 ==0) Aux_Message( stdout, "PART_A\n", __FUNCTION__ );
       //Maxwell_Interpolation( VD_Vec, VelDisp*V_acc );
      //if (p%100000 ==0) Aux_Message( stdout, "PARTB\n", __FUNCTION__ );

//      ParVel[0][p] = - V_acc_total[RingID*2+1]/V_acc_total[RingID*2+0]*NormParRadius[1]+Disc_BulkVel[0];//+VD_Vec[0];
//      ParVel[1][p] =   V_acc_total[RingID*2+1]/V_acc_total[RingID*2+0]*NormParRadius[0]+Disc_BulkVel[1];//+VD_Vec[1];
//
//
//      Ran0 = 0.9998*(RNG->GetValue( 0, 0.0, 1.0)-0.5);
//      Ran1 = 0.9998*(RNG->GetValue( 0, 0.0, 1.0)-0.5);
//      Ran2 = 0.9998*(RNG->GetValue( 0, 0.0, 1.0)-0.5);
 

 //     VD_Vec[0] = fabs(Ran0)/Ran0 * pow(1.5, -0.5)*VelDisp*(V_acc_total[RingID*2+1]/V_acc_total[RingID*2+0])*myErfInv2(2*fabs(Ran0));
 //     VD_Vec[1] = fabs(Ran1)/Ran1 * pow(1.5, -0.5)*VelDisp*(V_acc_total[RingID*2+1]/V_acc_total[RingID*2+0])*myErfInv2(2*fabs(Ran1));
 //     VD_Vec[2] = fabs(Ran2)/Ran2 * pow(1.5, -0.5)*VelDisp*(V_acc_total[RingID*2+1]/V_acc_total[RingID*2+0])*myErfInv2(2*fabs(Ran2));


///      ParVel[0][p] = - V_acc_total[RingID*2+1]/V_acc_total[RingID*2+0]*NormParRadius[1]+Disc_BulkVel[0]; //+VD_Vec[0];
///      ParVel[1][p] =   V_acc_total[RingID*2+1]/V_acc_total[RingID*2+0]*NormParRadius[0]+Disc_BulkVel[1]; //+VD_Vec[1];
///      ParVel[2][p] =                                                                    Disc_BulkVel[2]; //+VD_Vec[2];



//     VD_Vec[0] = fabs(Ran0)/Ran0 * pow(1.5, -0.5)*VelDisp*(V_acc_total_recv[RingID*2+1]/V_acc_total_recv[RingID*2+0])*myErfInv2(2*fabs(Ran0));
//     VD_Vec[1] = fabs(Ran1)/Ran1 * pow(1.5, -0.5)*VelDisp*(V_acc_total_recv[RingID*2+1]/V_acc_total_recv[RingID*2+0])*myErfInv2(2*fabs(Ran1));
//     VD_Vec[2] = fabs(Ran2)/Ran2 * pow(1.5, -0.5)*VelDisp*(V_acc_total_recv[RingID*2+1]/V_acc_total_recv[RingID*2+0])*myErfInv2(2*fabs(Ran2));

     Ran0 = 0.9998*(RNG->GetValue( 0, 0.0, 1.0))+0.0001;
     Ran1 = 0.9998*(RNG->GetValue( 0, 0.0, 1.0))+0.0001;
     Ran2 = 0.9998*(RNG->GetValue( 0, 0.0, 1.0))+0.0001;



     VD_Vec[0] = pow(1.5, -0.5)*VelDisp*(V_acc_total_recv[RingID*2+1]/V_acc_total_recv[RingID*2+0])*myErfInv2(2*Ran0-1);
     VD_Vec[1] = pow(1.5, -0.5)*VelDisp*(V_acc_total_recv[RingID*2+1]/V_acc_total_recv[RingID*2+0])*myErfInv2(2*Ran1-1);
     VD_Vec[2] = pow(1.5, -0.5)*VelDisp*(V_acc_total_recv[RingID*2+1]/V_acc_total_recv[RingID*2+0])*myErfInv2(2*Ran2-1);

     
      /// pow(3/2, -0.5): see meeting 20210805

     if (isnan(VD_Vec[0]) == 0 && isnan(VD_Vec[1]) == 0 && isnan(VD_Vec[2]) == 0 ){   
        ParVel[0][p] = - V_acc_total_recv[RingID*2+1]/V_acc_total_recv[RingID*2+0]*NormParRadius[1]+Disc_BulkVel[0]+VD_Vec[0];
        ParVel[1][p] =   V_acc_total_recv[RingID*2+1]/V_acc_total_recv[RingID*2+0]*NormParRadius[0]+Disc_BulkVel[1]+VD_Vec[1];
        ParVel[2][p] =                                                                              Disc_BulkVel[2]+VD_Vec[2];
     }
     else{
        ParVel[0][p] = - V_acc_total_recv[RingID*2+1]/V_acc_total_recv[RingID*2+0]*NormParRadius[1]+Disc_BulkVel[0];
        ParVel[1][p] =   V_acc_total_recv[RingID*2+1]/V_acc_total_recv[RingID*2+0]*NormParRadius[0]+Disc_BulkVel[1];
        ParVel[2][p] =                                                                              Disc_BulkVel[2];
        Count_Outliers  += 1;
     }

   }

   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n"  );
      Aux_Message( stdout, "  number of outliers             = %d\n",                         Count_Outliers  );

      Aux_Message( stdout, "=============================================================================\n"  );
   }


   delete V_acc_total;
   delete V_acc_total_recv;





// ============================================================================================================


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction



double myErfInv2(double x){
   double tt1, tt2, lnx, sgn;
   sgn = (x < 0) ? -1.0f : 1.0f;

   x = (1 - x)*(1 + x);        // x = 1 - x*x;
   lnx = logf(x);

   tt1 = 2/(3.1415926*0.147) + 0.5f * lnx;
   tt2 = 1/(0.147) * lnx;

   return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
}


void Maxwell_Interpolation( double VD_Vec[], const double v_rms )
{
   double slope, TempCDF, c1, c2, Ran, RanV; //countwhile;

   Ran = RNG->GetValue( 0, 0.0, 1.0 );
   c1 = 4.0*pi*pow(1.5/pi, 1.5)*pow(v_rms, -3.0);
   c2 = 2.0*pow(v_rms, 2.0)/3.0;
   //countwhile = 0;

   RanV = v_rms*pow(2.0/3.0, 0.5);   //Turning point
   TempCDF =c1*c2/4.0* ( pow(pi*c2, 0.5)*erf(RanV/pow(c2, 0.5)) - 2*RanV*exp(-pow(RanV,2)/c2));

   while (fabs(Ran- TempCDF)> 0.00001)
   {
      slope = c1*pow(RanV,2)*exp(-pow(RanV,2)/c2);
      RanV    += (Ran- TempCDF) / slope;
      TempCDF =c1*c2/4.0* ( pow(pi*c2, 0.5)*erf(RanV/pow(c2, 0.5)) - 2*RanV*exp(-pow(RanV,2)/c2));
    //  countwhile += 1;

      //if ( MPI_Rank == 0 && countwhile<100)    Aux_Message( stdout, "  number of RanV = %13.7e\n",                  RanV );
   }

   /////////////////////

   double xx,yy,zz, radii;

   radii = 1.1;

   while (radii > 1.0)
      xx = RNG->GetValue( 0, 0.0, 1.0 );
      yy = RNG->GetValue( 0, 0.0, 1.0 );
      zz = RNG->GetValue( 0, 0.0, 1.0 );
      radii = pow(xx*xx+yy*yy+zz*zz, 0.5);

   VD_Vec[0] = RanV*xx/radii;
   VD_Vec[1] = RanV*yy/radii;
   VD_Vec[2] = RanV*zz/radii;

}



//-------------------------------------------------------------------------------------------------------
//// Function    :  Vec2_FixRadius
//// Description :  Generate a random vector on z = 0 with fixed Radius and generate the
//                  corresponging normal vector with magnitude RanV
////
//// Note        :  None
////
//// Parameter   :   r,  RanVec[], NormVec[], RanV
////
//// Return      :  None
////-------------------------------------------------------------------------------------------------------


void Vec2_FixRadius( const double r, double RanVec[], double NormVec[], const double RanV )
{
   double theta;
   theta = 2*pi*RNG->GetValue( 0, 0.0, 1.0 );
   RanVec[0] = r*cos(theta);
   RanVec[1] = r*sin(theta);
   RanVec[2] = 0;
   NormVec[0] = -RanV*sin(theta);
   NormVec[1] =  RanV*cos(theta);
   NormVec[2] = 0;
}


//-------------------------------------------------------------------------------------------------------
////// Function    :  Disc_Interpolation
////// Description :  Use RanM to find the corresponding RanR by the method of interpolation
//////
////// Note        :  None
//////
////// Parameter   :  RanM, R, DiscMConst
//////
////// Return      :  None
//////-------------------------------------------------------------------------------------------------------

double Disc_Interpolation( const double RanM, const double R, const double DiscMConst )
{
   double r, slope, TempM, ERRM;
   ERRM = 1;
   r = R;
   TempM = DiscMConst*(R - (2*R)*exp(-1));
   while (ERRM > 0.0001)
   {
      slope = DiscMConst * (r/R)*exp(-(r/R));
      r    += (RanM-TempM) / slope;
      TempM = DiscMConst * (R - (R+r)*exp(-r/R));
      ERRM = fabs((TempM-RanM)/RanM);
   }
   return r;
}

#  endif
# ifdef PARTICLE


//-------------------------------------------------------------------------------------------------------
//// Function    :  AddNewParticleAttribute_LSSHalo
//// Description :  Add the problem-specific particle attributes
////
//// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
////                2. Invoke AddParticleField() for each of the problem-specific particle attribute:
////                   --> Attribute label sent to AddParticleField() will be used as the output name of the attribute
////                   --> Attribute index returned by AddParticleField() can be used to access the particle attribute data
////                3. Pre-declared attribute indices are put in Field.h
////
//// Parameter   :  None
////
//// Return      :  None
////-------------------------------------------------------------------------------------------------------
void AddNewParticleAttribute_LSSHalo()
{

// "Idx_ParLabel" has been predefined in Field.h
   if ( Idx_ParLabel == Idx_Undefined )  Idx_ParLabel = AddParticleAttribute( "ParLabel" );

} // FUNCTION : AddNewParticleAttribute_LSSHalo

void Init_Disc()
{

   if ( amr->Par->Init != PAR_INIT_BY_RESTART  || !OPT__RESTART_RESET || !AddParWhenRestart )   return;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

   const long   NNewPar        = ( MPI_Rank == 0 ) ? AddParNPar : 0;
   const long   NPar_AllRank   = NNewPar;

   real *NewParAtt[PAR_NATT_TOTAL];

   for (int v=0; v<PAR_NATT_TOTAL; v++)   NewParAtt[v] = new real [NNewPar];


// set particle attributes
// ============================================================================================================
   real *Time_AllRank   = NewParAtt[PAR_TIME];
   real *Mass_AllRank   = NewParAtt[PAR_MASS];
   real *Pos_AllRank[3] = { NewParAtt[PAR_POSX], NewParAtt[PAR_POSY], NewParAtt[PAR_POSZ] };
   real *Vel_AllRank[3] = { NewParAtt[PAR_VELX], NewParAtt[PAR_VELY], NewParAtt[PAR_VELZ] };
   real *Label_AllRank   = NewParAtt[Idx_ParLabel];

   if ( MPI_Rank == 0 )
   {
      const double ParM = Disc_Mass / NPar_AllRank;
      double Ran, RanR, RanM, RanV, RanVec[3], NormVec[3];
      Aux_Message(stdout, " Particle Mass = %13.7e\n", ParM) ;

//    initialize the RNG
      RNG = new RandomNumber_t( 1 );
      RNG->SetSeed( 0, Disc_RSeed );

      const double DiscMConst = Disc_Mass / (Disc_Decay_R - (Disc_Decay_R + Disc_Radius) * exp( - Disc_Radius/Disc_Decay_R ) );

      for ( long p = 0; p < NPar_AllRank; p++)
      {
         Time_AllRank[p] = 0;
//       mass
         Mass_AllRank[p] = ParM;

//       label
         Label_AllRank[p] = 0;

//       position
         Ran  = RNG->GetValue( 0, 0.0, 1.0);
         RanM = Ran*Disc_Mass;
         RanR = Disc_Interpolation( RanM, Disc_Decay_R, DiscMConst );
         RanV = sqrt(G*RanM/RanR);
         Vec2_FixRadius( RanR, RanVec, NormVec, RanV );
         for (int d = 0; d < 3; d++) Pos_AllRank[d][p] = RanVec[d] + Disc_Cen[d];

         //if (RanVec[0] < 0)
         //   Mass_AllRank[p] = ParM/100.0;
         //   Label_AllRank[p] = 2;
//       velocity
         for (int d = 0; d< 3; d++) Vel_AllRank[d][p] = Disc_BulkVel[d];
      } // for ( long p = 0; p < NPar_AllRank; p++)
      Aux_Message( stdout, " Particle mass              = %13.7e\n", ParM );

   } //if ( MPI_Rank == 0 )

// add particles here
   Par_AddParticleAfterInit( NNewPar, NewParAtt );
// free memory
   for (int v=0; v<PAR_NATT_TOTAL; v++)   delete [] NewParAtt[v];



  #  ifdef GRAVITY
  if (  OPT__SELF_GRAVITY  ||  OPT__EXT_POT   )
  {
//    initialize the k-space Green's function for the isolated BC.
     if ( OPT__SELF_GRAVITY  &&  OPT__BC_POT == BC_POT_ISOLATED )  Init_GreenFuncK();


//    evaluate the initial average density if it is not set yet (may already be set in Init_ByRestart)
     if ( AveDensity_Init <= 0.0 )    Poi_GetAverageDensity();


//    evaluate the gravitational potential
     if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", "Calculating gravitational potential" );

     for (int lv=0; lv<NLEVEL; lv++)
     {
        if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Lv %2d ... ", lv );

        Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, NULL_INT, DATA_GENERAL, _DENS, _NONE, Rho_ParaBuf, USELB_YES );

        Gra_AdvanceDt( lv, Time[lv], NULL_REAL, NULL_REAL, NULL_INT, amr->PotSg[lv], true, false, false, false, true );

        if ( lv > 0 )

        Buf_GetBufferData( lv, NULL_INT, NULL_INT, amr->PotSg[lv], POT_FOR_POISSON, _POTE, _NONE, Pot_ParaBuf, USELB_YES );

        if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
     } // for (int lv=0; lv<NLEVEL; lv++)

     if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", "Calculating gravitational potential" );
  } // if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
#  endif // #ifdef GARVITY
#  endif // # ifdef PARTICLE

// initialize particle acceleration
#  if ( defined PARTICLE  &&  defined STORE_PAR_ACC )
  if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", "Calculating particle acceleration" );

  const bool StoreAcc_Yes    = true;
  const bool UseStoredAcc_No = false;

  for (int lv=0; lv<NLEVEL; lv++)
  Par_UpdateParticle( lv, amr->PotSgTime[lv][ amr->PotSg[lv] ], NULL_REAL, PAR_UPSTEP_ACC_ONLY, StoreAcc_Yes, UseStoredAcc_No );

  if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", "Calculating particle acceleration" );
#  endif
#  ifdef PARTICLE
  Par_AfterAcceleration_Ptr( amr->Par->NPar_AcPlusInac, amr->Par->NPar_Active_AllRank,
                                  amr->Par->Mass, amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ,
                                  amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ,
                                  amr->Par->AccX, amr->Par->AccY, amr->Par->AccZ,
                                  amr->Par->Time, amr->Par->Attribute );

} // FUNCTION : Init_Disc

#  endif //ifdef PARTICLE
#  endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_LSSHalo
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Su_AddNewField_ELBDM_Halo_Stability_Test(void)
{

#  if ( NCOMP_PASSIVE_USER > 0 )
   Su_Idx_Dens0 = AddField( "Dens0", NORMALIZE_NO, INTERP_FRAC_NO );
   if ( MPI_Rank == 0 )   printf("Su_Idx_Dens0 = %d \n", Su_Idx_Dens0);
   #  endif

} // FUNCTION : AddNewField_ELBDM_Halo_Stability_Test

void Init_TestProb_ELBDM_LSSHalo()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// validate the compilation flags and runtime parameters
   Validate();

#  if ( MODEL == ELBDM  &&  defined GRAVITY )
// set the problem-specific runtime parameters
   SetParameter();

//   Init_Function_User_Ptr      = SetGridIC;
   Init_Function_User_Ptr      = NULL;
   Init_Field_User_Ptr         = Su_AddNewField_ELBDM_Halo_Stability_Test;
   Flag_User_Ptr               = NULL;
   Mis_GetTimeStep_User_Ptr    = NULL;
   BC_User_Ptr                 = BC_Disc_Heating;
   Aux_Record_User_Ptr         = Record_Disc_Heating;
   End_User_Ptr                = NULL;
//#  ifdef GRAVITY
//   Init_ExternalAcc_Ptr        = NULL;
//   Init_ExternalPot_Ptr        = NULL;
//#  endif // #ifdef GRAVITY
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr     = Par_Disc_Heating;         // option: PAR_INIT=1;            example: Particle/Par_Init_ByFunction.cpp
   Par_Init_Attribute_User_Ptr = AddNewParticleAttribute_LSSHalo;    // set PAR_NATT_USER;             example: TestProblem/Hydro/AGORA_IsolatedGalaxy/Init_TestProb_Hydro_AGORA_IsolatedGalaxy.cpp --> AddNewParticleAttribute()
   Init_User_Ptr               = Init_Disc;
#  endif
#  endif // if ( MODEL == ELBDM  &&  defined GRAVITY )

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_LSSHalo
