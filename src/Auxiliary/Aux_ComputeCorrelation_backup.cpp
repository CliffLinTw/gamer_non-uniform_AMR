#include "GAMER.h"

extern void SetTempIntPara( const int lv, const int Sg_Current, const double PrepTime, const double Time0, const double Time1,
                            bool &IntTime, int &Sg, int &Sg_IntT, real &Weighting, real &Weighting_IntT );


//-------------------------------------------------------------------------------------------------------
// Function    :  Interpolate
// Description :  Interpolate the initial density for a given radius
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void InterpolateDensity(real *mean_inter, real *std_inter, const Profile_with_Sigma_t *prof_init[], const int NProf, const int bin_index, const double r)
{
   double delta_r, x;
   if (r>prof_init[0]->Radius[bin_index])
   {
      int bin_index_right = bin_index+1;
      delta_r             = prof_init[0]->Radius[bin_index_right] - prof_init[0]->Radius[bin_index];
      x                   = (r - prof_init[0]->Radius[bin_index]) / delta_r;
//    check x 
      if (x<(real)(0.0))
         Aux_Error( ERROR_INFO, "x (%14.7e) < 0.0 !!\n", x );
      else if (x>(real)(1.0))
         Aux_Error( ERROR_INFO, "x (%14.7e) > 1.0 !!\n", x );
//    interpolate
      for  (int i=0; i<NProf; i++)
      {
         mean_inter[i] = prof_init[i]->Data      [bin_index]*(real)(1.-x) + prof_init[i]->Data      [bin_index_right]*(real)x;
         std_inter[i]  = prof_init[i]->Data_Sigma[bin_index]*(real)(1.-x) + prof_init[i]->Data_Sigma[bin_index_right]*(real)x;
      }
   }
   else
   {
//    no left hand side bin, no interpolation
      if (bin_index==0)
      {
         for (int i=0; i<NProf; i++)
         {
            mean_inter[i] = prof_init[i]->Data      [bin_index];
            std_inter[i]  = prof_init[i]->Data_Sigma[bin_index];
         }
      }
      else
      {
         int bin_index_left  = bin_index-1;
         delta_r             = prof_init[0]->Radius[bin_index] - prof_init[0]->Radius[bin_index_left];
         x                   = (r - prof_init[0]->Radius[bin_index_left]) / delta_r;
//       check x 
         if (x<(real)(0.0))
            Aux_Error( ERROR_INFO, "x (%14.7e) < 0.0 !!\n", x );
         else if (x>(real)(1.0))
            Aux_Error( ERROR_INFO, "x (%14.7e) > 1.0 !!\n", x );
//       interpolate
         for  (int i=0; i<NProf; i++)
         {
            mean_inter[i] = prof_init[i]->Data      [bin_index_left]*(real)(1.-x) + prof_init[i]->Data      [bin_index]*(real)x;
            std_inter[i]  = prof_init[i]->Data_Sigma[bin_index_left]*(real)(1.-x) + prof_init[i]->Data_Sigma[bin_index]*(real)x;
         }
      }
   }
}


//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_ComputeCorrelation
// Description :  Compute the average radial profile of target field(s)
//
// Note        :  1. Results will be stored in the input "Correlation" object
//                   --> Correlation->Radius[]: Radial coordinate at each bin
//                       Correlation->Data  []: Profile data at each bin
//                       Correlation->Weight[]: Total weighting at each bin
//                       Correlation->NCell []: Number of cells at each bin
//                       Correlation->NBin    : Total number of bins
//                   --> See the "Profile_t" structure defined in "include/Profile.h" for details
//                   --> These arrays will be free'd when deleting "Correlation"
//                2. Maximum radius adopted when actually computing the profile may be larger than the input "r_max"
//                   --> Because "r_max" in general does not coincide with the right edge of the maximum bin
//                3. Support hybrid OpenMP/MPI parallelization
//                   --> All ranks will share the same profile data after invoking this function
//                4. Use cell volume as the weighting of each cell
//                   --> Will support other weighting functions in the future
//                5. Support computing multiple fields
//                   --> The order of fields to be returned follows TVarBitIdx[]
//
// Parameter   :  Correlation : Profile_t object array to store the correlation function
//                prof_init   : Profile_with_Sigma_t object array for storing the mean and standard deviation quantities for calculating correlation function
//                Center      : Target center coordinates
//                r_max_input : Maximum radius for computing the profile
//                              --> See also "Note-2" above
//                dr_min      : Minimum bin size
//                              --> For linear bin, this is the size of all bins
//                                  For log    bin, this is the size of the 0th bin
//                LogBin      : true/false --> log/linear bins
//                LogBinRatio : Ratio of adjacent log bins
//                              --> Right edge of log bin n = dr_min*LogBinRatio^n
//                RemoveEmpty : true  --> remove empty bins from the data
//                              false --> these empty bins will still be in the profile arrays with
//                                        Data[empty_bin]=Weight[empty_bin]=NCell[empty_bin]=0
//                TVarBitIdx  : Bitwise indices of target variables for computing the profiles
//                              --> Supported indices (defined in Macro.h):
//                                             [, _ENPY, _EINT, _POTE]
//                                     ELBDM : _DENS, _REAL, _IMAG [, _POTE]
//                              --> For a passive scalar with an integer field index FieldIdx returned by AddField(),
//                                  one can convert it to a bitwise field index by BIDX(FieldIdx)
//                NProf       : Number of Profile_t objects in Correlation
//                Min/MaxLv   : Consider patches on levels from MinLv to MaxLv
//                PatchType   : Only consider patches of the specified type
//                              --> Supported types: PATCH_LEAF, PATCH_NONLEAF, PATCH_BOTH, PATCH_LEAF_PLUS_MAXNONLEAF
//                              --> PATCH_LEAF_PLUS_MAXNONLEAF includes leaf patches on all target levels
//                                  (i.e., MinLv ~ MaxLv) and non-leaf patches only on MaxLv
//                PrepTime    : Target physical time to prepare data
//                              --> If PrepTime<0, turn off temporal interpolation and always use the most recent data
//
// Example     :  const double      Center[3]      = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
//                const double      MaxRadius      = 0.5*amr->BoxSize[0];
//                const double      MinBinSize     = amr->dh[MAX_LEVEL];
//                const bool        LogBin         = true;
//                const double      LogBinRatio    = 1.25;
//                const bool        RemoveEmptyBin = true;
//                const long        TVar[]         = { _DENS, _PRES };
//                const int         NProf          = 2;
//                const int         MinLv          = 0;
//                const int         MaxLv          = MAX_LEVEL;
//                const PatchType_t PatchType      = PATCH_LEAF_PLUS_MAXNONLEAF;
//                const double      PrepTime       = -1.0;
//
//                Profile_t Correlation_Dens, Correlation_Pres;
//                Profile_t *Correlation[] = { &Correlation_Dens, &Correlation_Pres };
//
//                Aux_ComputeProfile( Correlation, Center, MaxRadius, MinBinSize, LogBin, LogBinRatio, RemoveEmptyBin,
//                                    TVar, NProf, MinLv, MaxLv, PatchType, PrepTime );
//
//                if ( MPI_Rank == 0 )
//                {
//                   for (int p=0; p<NProf; p++)
//                   {
//                      char Filename[MAX_STRING];
//                      sprintf( Filename, "Correlation_function%d.txt", p );
//                      FILE *File = fopen( Filename, "w" );
//                      fprintf( File, "#%19s  %21s  %21s  %10s\n", "Radius", "Data", "Weight", "Cells" );
//                      for (int b=0; b<Correlation[p]->NBin; b++)
//                         fprintf( File, "%20.14e  %21.14e  %21.14e  %10ld\n",
//                                  Correlation[p]->Radius[b], Correlation[p]->Data[b], Correlation[p]->Weight[b], Correlation[p]->NCell[b] );
//                      fclose( File );
//                   }
//                }
//
// Return      :  Correlation
//-------------------------------------------------------------------------------------------------------
//void Aux_ComputeCorrelation( Profile_t *Correlation[], FieldIdx_t *Passive_idx[],  const Profile_t *prof_init[], const double Center[], 
void Aux_ComputeCorrelation( Profile_t *Correlation[], const Profile_with_Sigma_t *prof_init[], const double Center[], 
                             const double r_max_input, const double dr_min, const bool LogBin, const double LogBinRatio,
                             const bool RemoveEmpty, const long TVarBitIdx[], const int NProf, const int MinLv, const int MaxLv, 
                             const PatchType_t PatchType, const double PrepTime, const double dr_min_prof)
{

// check
#  ifdef GAMER_DEBUG
   if ( r_max_input <= 0.0 )
      Aux_Error( ERROR_INFO, "r_max_input (%14.7e) <= 0.0 !!\n", r_max_input );

   if ( dr_min <= 0.0 )
      Aux_Error( ERROR_INFO, "dr_min (%14.7e) <= 0.0 !!\n", dr_min );

   if ( LogBin  &&  LogBinRatio <= 1.0 )
      Aux_Error( ERROR_INFO, "LogBinRatio (%14.7e) <= 1.0 !!\n", LogBinRatio );

   if ( MinLv < 0  ||  MinLv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect MinLv (%d) !!\n", MinLv );

   if ( MaxLv < 0  ||  MaxLv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect MaxLv (%d) !!\n", MaxLv );

   if ( MinLv > MaxLv )
      Aux_Error( ERROR_INFO, "MinLv (%d) > MaxLv (%d) !!\n", MinLv, MaxLv );

   if ( NProf != NCOMP_PASSIVE )
      Aux_Error( ERROR_INFO, "NProf (%d) != NCOMP_PASSIVE (%d) !!\n", NProf, NCOMP_PASSIVE );
#  endif


// precompute the integer indices of intrinsic fluid fields for better performance
   const int IdxUndef = -1;
   int TFluIntIdx[NProf];

   for (int p=0; p<NProf; p++)
   {
      TFluIntIdx[p] = IdxUndef;

      for (int v=0; v<NCOMP_TOTAL; v++)
         if ( TVarBitIdx[p] & (1L<<v) )   TFluIntIdx[p] = v;
   }


// check whether _POTE is in TVarBitIdx since the potential array may have not been computed during initialization
#  ifdef GRAVITY
   bool InclPot = false;

   for (int p=0; p<NProf; p++)
      if ( TVarBitIdx[p] & _POTE )   InclPot = true;
#  endif


// initialize the profile objects
   for (int p=0; p<NProf; p++)
   {
//    get the total number of radial bins and the corresponding maximum radius
      if ( LogBin )
      {
         Correlation[p]->NBin      = int( log(r_max_input/dr_min)/log(LogBinRatio) ) + 2;
         Correlation[p]->MaxRadius = dr_min*pow( LogBinRatio, Correlation[p]->NBin-1 );
      }

      else // linear bin
      {
         Correlation[p]->NBin      = (int)ceil( r_max_input / dr_min );
         Correlation[p]->MaxRadius = dr_min*Correlation[p]->NBin;
      }


//    record profile parameters
      for (int d=0; d<3; d++)    Correlation[p]->Center[d] = Center[d];

      Correlation[p]->LogBin = LogBin;

      if ( LogBin )  Correlation[p]->LogBinRatio = LogBinRatio;


//    allocate all member arrays of Correlation
      Correlation[p]->AllocateMemory();


//    record radial coordinates
      if ( LogBin )
         for (int b=0; b<Correlation[0]->NBin; b++)    Correlation[p]->Radius[b] = dr_min*pow( LogBinRatio, b-0.5 );
      else
         for (int b=0; b<Correlation[0]->NBin; b++)    Correlation[p]->Radius[b] = (b+0.5)*dr_min;

   } // for (int p=0; p<NProf; p++)


// allocate memory for the per-thread arrays
#  ifdef OPENMP
   const int NT = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int NT = 1;
#  endif

   double ***OMP_Data=NULL, ***OMP_Weight=NULL;
   long   ***OMP_NCell=NULL;

   Aux_AllocateArray3D( OMP_Data,   NProf, NT, Correlation[0]->NBin );
   Aux_AllocateArray3D( OMP_Weight, NProf, NT, Correlation[0]->NBin );
   Aux_AllocateArray3D( OMP_NCell,  NProf, NT, Correlation[0]->NBin );


// collect profile data in this rank
   const double r_max2      = SQR( Correlation[0]->MaxRadius );
   const double HalfBox[3]  = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const bool   Periodic[3] = { OPT__BC_FLU[0] == BC_FLU_PERIODIC,
                                OPT__BC_FLU[2] == BC_FLU_PERIODIC,
                                OPT__BC_FLU[4] == BC_FLU_PERIODIC };

#  pragma omp parallel
   {
#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif

//    initialize arrays
      for (int p=0; p<NProf; p++)
      for (int b=0; b<Correlation[0]->NBin; b++)
      {
         OMP_Data  [p][TID][b] = 0.0;
         OMP_Weight[p][TID][b] = 0.0;
         OMP_NCell [p][TID][b] = 0;
      }

//    allocate passive scalar arrays
      real *Passive      = new real [NCOMP_PASSIVE];
//      real *Passive_IntT = new real [NCOMP_PASSIVE];

//    loop over all target levels
      for (int lv=MinLv; lv<=MaxLv; lv++)
      {
         const double dh = amr->dh[lv];
         const double dv = CUBE( dh );


//       determine temporal interpolation parameters
         bool FluIntTime = false;
         int  FluSg      = amr->FluSg[lv];
         int  FluSg_IntT;
         real FluWeighting, FluWeighting_IntT;

#        ifdef MHD
         bool MagIntTime = false;
         int  MagSg      = amr->MagSg[lv];
         int  MagSg_IntT;
         real MagWeighting, MagWeighting_IntT;
#        endif

#        ifdef GRAVITY
         bool PotIntTime = false;
         int  PotSg      = amr->PotSg[lv];
         int  PotSg_IntT;
         real PotWeighting, PotWeighting_IntT;
#        endif

         if ( PrepTime >= 0.0 )
         {
//          fluid
            SetTempIntPara( lv, amr->FluSg[lv], PrepTime, amr->FluSgTime[lv][0], amr->FluSgTime[lv][1],
                            FluIntTime, FluSg, FluSg_IntT, FluWeighting, FluWeighting_IntT );

//          magnetic field
#           ifdef MHD
            SetTempIntPara( lv, amr->MagSg[lv], PrepTime, amr->MagSgTime[lv][0], amr->MagSgTime[lv][1],
                            MagIntTime, MagSg, MagSg_IntT, MagWeighting, MagWeighting_IntT );
#           endif

//          potential
#           ifdef GRAVITY
            if ( InclPot )
               SetTempIntPara( lv, amr->PotSg[lv], PrepTime, amr->PotSgTime[lv][0], amr->PotSgTime[lv][1],
                               PotIntTime, PotSg, PotSg_IntT, PotWeighting, PotWeighting_IntT );
#           endif
         }


//       use the "static" schedule for reproducibility
#        pragma omp for schedule( static )
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         {
//          skip untargeted patches
            if ( amr->patch[0][lv][PID]->son != -1 )
            {
               if ( PatchType == PATCH_LEAF )                                    continue;
               if ( PatchType == PATCH_LEAF_PLUS_MAXNONLEAF  &&  lv != MaxLv )   continue;
            }

            else
            {
               if ( PatchType == PATCH_NONLEAF )                                 continue;
            }


            const real (*FluidPtr)[PS1][PS1][PS1] = amr->patch[ FluSg ][lv][PID]->fluid;
#           ifdef GRAVITY
            const real (*PotPtr  )[PS1][PS1]      = amr->patch[ PotSg ][lv][PID]->pot;
#           endif

//          pointer for temporal interpolation
            const real (*FluidPtr_IntT)[PS1][PS1][PS1] = ( FluIntTime ) ? amr->patch[ FluSg_IntT ][lv][PID]->fluid : NULL;
#           ifdef GRAVITY
            const real (*PotPtr_IntT  )[PS1][PS1]      = ( PotIntTime ) ? amr->patch[ PotSg_IntT ][lv][PID]->pot   : NULL;
#           endif


            const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh - Center[0];
            const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh - Center[1];
            const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh - Center[2];

            for (int k=0; k<PS1; k++)  {  double dz = z0 + k*dh;
                                          if ( Periodic[2] ) {
                                             if      ( dz > +HalfBox[2] )  {  dz -= amr->BoxSize[2];  }
                                             else if ( dz < -HalfBox[2] )  {  dz += amr->BoxSize[2];  }
                                          }
            for (int j=0; j<PS1; j++)  {  double dy = y0 + j*dh;
                                          if ( Periodic[1] ) {
                                             if      ( dy > +HalfBox[1] )  {  dy -= amr->BoxSize[1];  }
                                             else if ( dy < -HalfBox[1] )  {  dy += amr->BoxSize[1];  }
                                          }
            for (int i=0; i<PS1; i++)  {  double dx = x0 + i*dh;
                                          if ( Periodic[0] ) {
                                             if      ( dx > +HalfBox[0] )  {  dx -= amr->BoxSize[0];  }
                                             else if ( dx < -HalfBox[0] )  {  dx += amr->BoxSize[0];  }
                                          }

               const double r2 = SQR(dx) + SQR(dy) + SQR(dz);

               if ( r2 < r_max2 )
               {
                  const double r   = sqrt( r2 );
                  const int    bin = ( LogBin ) ? (  (r<dr_min) ? 0 : int( log(r/dr_min)/log(LogBinRatio) ) + 1  )
                                                : int( r/dr_min );
//                prevent from round-off errors
                  if ( bin >= Correlation[0]->NBin )   continue;

//                check
#                 ifdef GAMER_DEBUG
                  if ( bin < 0 )    Aux_Error( ERROR_INFO, "bin (%d) < 0 !!\n", bin );
#                 endif

//                interpolate to get mean value at r
                  real mean_value[NProf], std_value[NProf];
//                find corresponding bin index in density profile, which always uses linear bin
                  const int    bin_prof = int (r/dr_min_prof); 
                  InterpolateDensity( mean_value, std_value, prof_init, NProf, bin_prof, r );

//                prepare passive scalars (for better sustainability, always do it even when unnecessary)
                  for (int v_out=0; v_out<NCOMP_PASSIVE; v_out++)
                  {
                     const int v_in = v_out + NCOMP_FLUID;

                     Passive     [v_out] = FluidPtr     [v_in][k][j][i];
                  }

                  for (int p=0; p<NProf; p++)
                  {
//                   intrinsic fluid fields
                     if ( TFluIntIdx[p] != IdxUndef )
                     {
                        const real Weight = dv;
                        real delta  = ( FluIntTime )
                                          ? ( FluWeighting     *FluidPtr     [ TFluIntIdx[p] ][k][j][i]
                                            + FluWeighting_IntT*FluidPtr_IntT[ TFluIntIdx[p] ][k][j][i] )
                                          :                     FluidPtr     [ TFluIntIdx[p] ][k][j][i]  ;
//                        real delta_passive = amr->patch[FluSg][lv][PID]->fluid[ *(Passive_idx[p]) ][k][j][i];
                        real delta_passive = Passive[p];
                        delta         = delta/mean_value[p]         - (real)1.;
                        delta_passive = delta_passive/mean_value[p] - (real)1.;

                        OMP_Data  [p][TID][bin] += delta*delta_passive*Weight;
//                        OMP_Data  [p][TID][bin] += delta*Weight;
                        OMP_Weight[p][TID][bin] += Weight;
                        OMP_NCell [p][TID][bin] ++;
                     }

//                   other fields
                     else
                     {
                        switch ( TVarBitIdx[p] )
                        {
//                         gravitational potential
#                          ifdef GRAVITY
                           case _POTE:
                           {
                              const real Weight = ( FluIntTime )    // weighted by cell mass
                                                ? ( FluWeighting     *FluidPtr     [DENS][k][j][i]
                                                  + FluWeighting_IntT*FluidPtr_IntT[DENS][k][j][i] )*dv
                                                :                     FluidPtr     [DENS][k][j][i]  *dv;
                              

                              real delta  = ( PotIntTime )
                                                ? ( PotWeighting     *PotPtr     [k][j][i]
                                                  + PotWeighting_IntT*PotPtr_IntT[k][j][i] )
                                                       :                     PotPtr     [k][j][i]  ;
//                              real delta_passive = amr->patch[0][lv][PID]->fluid[ *(Passive_idx[p]) ][k][j][i];
                              real delta_passive = Passive[p];
                              delta         = delta/mean_value[p]         - (real)1.;
                              delta_passive = delta_passive/mean_value[p] - (real)1.;

                              OMP_Data  [p][TID][bin] += delta*delta_passive*Weight;
                              OMP_Weight[p][TID][bin] += Weight;
                              OMP_NCell [p][TID][bin] ++;
                           }
                           break;
#                          endif

/*
//                         derived fields
#                          if ( MODEL == HYDRO )
                           case _VELR:
                           {
                              const real Weight = ( FluIntTime )    // weighted by cell mass
                                                ? ( FluWeighting     *FluidPtr     [DENS][k][j][i]
                                                  + FluWeighting_IntT*FluidPtr_IntT[DENS][k][j][i] )*dv
                                                :                     FluidPtr     [DENS][k][j][i]  *dv;

                              const real MomR   = ( FluIntTime )
                                                ? ( FluWeighting     *( FluidPtr     [MOMX][k][j][i]*dx +
                                                                        FluidPtr     [MOMY][k][j][i]*dy +
                                                                        FluidPtr     [MOMZ][k][j][i]*dz )
                                                  + FluWeighting_IntT*( FluidPtr_IntT[MOMX][k][j][i]*dx +
                                                                        FluidPtr_IntT[MOMY][k][j][i]*dy +
                                                                        FluidPtr_IntT[MOMZ][k][j][i]*dz ) ) / r
                                                :                     ( FluidPtr     [MOMX][k][j][i]*dx +
                                                                        FluidPtr     [MOMY][k][j][i]*dy +
                                                                        FluidPtr     [MOMZ][k][j][i]*dz )   / r;

                              OMP_Data  [p][TID][bin] += MomR*dv;    // vr*(rho*dv)
                              OMP_Weight[p][TID][bin] += Weight;
                              OMP_NCell [p][TID][bin] ++;
                           }
                           break;

                           case _PRES:
                           {
                              const bool CheckMinPres_No = false;
                              const real Weight          = dv;
#                             ifdef MHD
                              const real Emag            = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg      );
                              const real Emag_IntT       = ( MagIntTime )
                                                         ? MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg_IntT )
                                                         : NULL_REAL;
#                             else
                              const real Emag            = NULL_REAL;
                              const real Emag_IntT       = NULL_REAL;
#                             endif
                              const real Pres = ( FluIntTime )
                                              ?   FluWeighting     *Hydro_Con2Pres( FluidPtr     [DENS][k][j][i],
                                                                                    FluidPtr     [MOMX][k][j][i],
                                                                                    FluidPtr     [MOMY][k][j][i],
                                                                                    FluidPtr     [MOMZ][k][j][i],
                                                                                    FluidPtr     [ENGY][k][j][i],
                                                                                    Passive,
                                                                                    CheckMinPres_No, NULL_REAL, Emag,
                                                                                    EoS_DensEint2Pres_CPUPtr, EoS_AuxArray_Flt,
                                                                                    EoS_AuxArray_Int, h_EoS_Table, NULL )
                                                + FluWeighting_IntT*Hydro_Con2Pres( FluidPtr_IntT[DENS][k][j][i],
                                                                                    FluidPtr_IntT[MOMX][k][j][i],
                                                                                    FluidPtr_IntT[MOMY][k][j][i],
                                                                                    FluidPtr_IntT[MOMZ][k][j][i],
                                                                                    FluidPtr_IntT[ENGY][k][j][i],
                                                                                    Passive_IntT,
                                                                                    CheckMinPres_No, NULL_REAL, Emag_IntT,
                                                                                    EoS_DensEint2Pres_CPUPtr, EoS_AuxArray_Flt,
                                                                                    EoS_AuxArray_Int, h_EoS_Table, NULL )
                                              :                     Hydro_Con2Pres( FluidPtr     [DENS][k][j][i],
                                                                                    FluidPtr     [MOMX][k][j][i],
                                                                                    FluidPtr     [MOMY][k][j][i],
                                                                                    FluidPtr     [MOMZ][k][j][i],
                                                                                    FluidPtr     [ENGY][k][j][i],
                                                                                    Passive,
                                                                                    CheckMinPres_No, NULL_REAL, Emag,
                                                                                    EoS_DensEint2Pres_CPUPtr, EoS_AuxArray_Flt,
                                                                                    EoS_AuxArray_Int, h_EoS_Table, NULL );

                              OMP_Data  [p][TID][bin] += Pres*Weight;
                              OMP_Weight[p][TID][bin] += Weight;
                              OMP_NCell [p][TID][bin] ++;
                           }
                           break;

                           case _EINT_DER:
                           {
                              const real Weight = dv;
                              const real Dens   = FluidPtr[DENS][k][j][i];

//                            use the dual-energy variable to calculate the internal energy directly, if applicable
#                             ifdef DUAL_ENERGY

#                             if   ( DUAL_ENERGY == DE_ENPY )
                              const bool CheckMinPres_No = false;
                              const real Enpy = FluidPtr[ENPY][k][j][i];
                              const real Pres = Hydro_DensEntropy2Pres( Dens, Enpy, EoS_AuxArray_Flt[1],
                                                                        CheckMinPres_No, NULL_REAL );
                              const real Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, Passive, EoS_AuxArray_Flt,
                                                                          EoS_AuxArray_Int, h_EoS_Table, NULL );
#                             elif ( DUAL_ENERGY == DE_EINT )
#                             error : DE_EINT is NOT supported yet !!
#                             endif

#                             else // #ifdef DUAL_ENERGY

                              const bool CheckMinEint_No = false;
                              const real MomX            = FluidPtr[MOMX][k][j][i];
                              const real MomY            = FluidPtr[MOMY][k][j][i];
                              const real MomZ            = FluidPtr[MOMZ][k][j][i];
                              const real Etot            = FluidPtr[ENGY][k][j][i];
#                             ifdef MHD
                              const real Emag            = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg      );
                              const real Emag_IntT       = ( MagIntTime )
                                                         ? MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg_IntT )
                                                         : NULL_REAL;
#                             else
                              const real Emag            = NULL_REAL;
                              const real Emag_IntT       = NULL_REAL;
#                             endif
                              const real Eint = ( FluIntTime )
                                              ?   FluWeighting     *Hydro_Con2Eint( FluidPtr     [DENS][k][j][i],
                                                                                    FluidPtr     [MOMX][k][j][i],
                                                                                    FluidPtr     [MOMY][k][j][i],
                                                                                    FluidPtr     [MOMZ][k][j][i],
                                                                                    FluidPtr     [ENGY][k][j][i],
                                                                                    CheckMinEint_No, NULL_REAL, Emag )
                                                + FluWeighting_IntT*Hydro_Con2Eint( FluidPtr_IntT[DENS][k][j][i],
                                                                                    FluidPtr_IntT[MOMX][k][j][i],
                                                                                    FluidPtr_IntT[MOMY][k][j][i],
                                                                                    FluidPtr_IntT[MOMZ][k][j][i],
                                                                                    FluidPtr_IntT[ENGY][k][j][i],
                                                                                    CheckMinEint_No, NULL_REAL, Emag_IntT )
                                              :                     Hydro_Con2Eint( FluidPtr     [DENS][k][j][i],
                                                                                    FluidPtr     [MOMX][k][j][i],
                                                                                    FluidPtr     [MOMY][k][j][i],
                                                                                    FluidPtr     [MOMZ][k][j][i],
                                                                                    FluidPtr     [ENGY][k][j][i],
                                                                                    CheckMinEint_No, NULL_REAL, Emag );
#                             endif // #ifdef DUAL_ENERGY ... else

                              OMP_Data  [p][TID][bin] += Eint*Weight;
                              OMP_Weight[p][TID][bin] += Weight;
                              OMP_NCell [p][TID][bin] ++;
                           }
                           break;
#                          endif // HYDRO
*/

                           default:
                              Aux_Error( ERROR_INFO, "unsupported field (%ld) !!\n", TVarBitIdx[p] );
                              exit( 1 );
                        } // switch ( TVarBitIdx[p] )
                     } // if ( TFluIntIdx[p] != IdxUndef ) ... else ...
                  } // for (int p=0; p<NProf; p++)
               } // if ( r2 < r_max2 )
            }}} // i,j,k
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      } // for (int lv=lv_min; lv<=lv_max; lv++)

      delete [] Passive;         Passive      = NULL;

   } // OpenMP parallel region


// sum over all OpenMP threads
   for (int p=0; p<NProf; p++)
   {
      for (int b=0; b<Correlation[0]->NBin; b++)
      {
         Correlation[p]->Data  [b]  = OMP_Data  [p][0][b];
         Correlation[p]->Weight[b]  = OMP_Weight[p][0][b];
         Correlation[p]->NCell [b]  = OMP_NCell [p][0][b];
      }

      for (int t=1; t<NT; t++)
      for (int b=0; b<Correlation[0]->NBin; b++)
      {
         Correlation[p]->Data  [b] += OMP_Data  [p][t][b];
         Correlation[p]->Weight[b] += OMP_Weight[p][t][b];
         Correlation[p]->NCell [b] += OMP_NCell [p][t][b];
      }
   }

// free per-thread arrays
   Aux_DeallocateArray3D( OMP_Data );
   Aux_DeallocateArray3D( OMP_Weight );
   Aux_DeallocateArray3D( OMP_NCell );


// collect data from all ranks (in-place reduction)
#  ifndef SERIAL
   for (int p=0; p<NProf; p++)
   {
      if ( MPI_Rank == 0 )
      {
         MPI_Reduce( MPI_IN_PLACE,    Correlation[p]->Data,   Correlation[p]->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
         MPI_Reduce( MPI_IN_PLACE,    Correlation[p]->Weight, Correlation[p]->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
         MPI_Reduce( MPI_IN_PLACE,    Correlation[p]->NCell , Correlation[p]->NBin, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );
      }

      else
      {
         MPI_Reduce( Correlation[p]->Data,   NULL,            Correlation[p]->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
         MPI_Reduce( Correlation[p]->Weight, NULL,            Correlation[p]->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
         MPI_Reduce( Correlation[p]->NCell,  NULL,            Correlation[p]->NBin, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );
      }
   }
#  endif


// compute profile by the root rank
   if ( MPI_Rank == 0 )
   {
      for (int p=0; p<NProf; p++)
      for (int b=0; b<Correlation[0]->NBin; b++)
      {
//       skip empty bins since both their data and weight are zero
         if ( Correlation[p]->NCell[b] > 0L )    Correlation[p]->Data[b] /= Correlation[p]->Weight[b];
      }
   }


// broadcast data to all ranks
   for (int p=0; p<NProf; p++)
   {
      MPI_Bcast( Correlation[p]->Data,   Correlation[p]->NBin, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Bcast( Correlation[p]->Weight, Correlation[p]->NBin, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Bcast( Correlation[p]->NCell,  Correlation[p]->NBin, MPI_LONG,   0, MPI_COMM_WORLD );
   }


// remove the empty bins
// --> all ranks do the same work so that no data broadcast is required
   if ( RemoveEmpty )
   {
      for (int b=0; b<Correlation[0]->NBin; b++)
      {
         if ( Correlation[0]->NCell[b] != 0L )   continue;

//       remove consecutive empty bins at the same time for better performance
         int b_up;
         for (b_up=b+1; b_up<Correlation[0]->NBin; b_up++)
            if ( Correlation[0]->NCell[b_up] != 0L )   break;

         const int stride = b_up - b;

         for (b_up=b+stride; b_up<Correlation[0]->NBin; b_up++)
         {
            const int b_up_ms = b_up - stride;

            for (int p=0; p<NProf; p++)
            {
               Correlation[p]->Radius[b_up_ms] = Correlation[p]->Radius[b_up];
               Correlation[p]->Data  [b_up_ms] = Correlation[p]->Data  [b_up];
               Correlation[p]->Weight[b_up_ms] = Correlation[p]->Weight[b_up];
               Correlation[p]->NCell [b_up_ms] = Correlation[p]->NCell [b_up];
            }
         }

//       reset the total number of bins
         for (int p=0; p<NProf; p++)
            Correlation[p]->NBin -= stride;
      } // for (int b=0; b<Correlation->NBin; b++)

//    update the maximum radius since the last bin may have been removed
      for (int p=0; p<NProf; p++)
      {
         const int LastBin = Correlation[p]->NBin-1;

         Correlation[p]->MaxRadius = ( LogBin ) ? Correlation[p]->Radius[LastBin]*sqrt( LogBinRatio )
                                         : Correlation[p]->Radius[LastBin] + 0.5*dr_min;
      }
   } // if ( RemoveEmpty )

} // FUNCTION : Aux_ComputeCorrelation
