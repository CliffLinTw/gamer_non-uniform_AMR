#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Load_StepTable
// Description :  Load the dump table from the file "Input__StepTable"
//-------------------------------------------------------------------------------------------------------
void Init_Load_StepTable()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_Load_StepTable ...\n" );


   const char FileName[] = "Input__StepTable";

   if ( !Aux_CheckFileExist(FileName) )   Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );

   FILE *File = fopen( FileName, "r" );

   const int MaxLine = 10000;
   char *input_line = NULL;
   size_t len = 0;
   int Trash, line, n;


// allocate the step table
   DumpTable = new double [MaxLine];


// skip the header
   getline( &input_line, &len, File );

// begin to read
   for (line=0; line<MaxLine; line++)
   {
      n = getline( &input_line, &len, File );

//    check
      if ( n <= 1 )
         Aux_Error( ERROR_INFO, "incorrect reading at line %d of the file <%s> !!\n", line+2, FileName );

      sscanf( input_line, "%d%lf", &Trash, &DumpTable[line] );

//    stop the reading
      if ( input_line[0] == 42 )                   // '*' == 42
      {

//       ensure that at least one step index is loaded
         if ( line == 0 )
            Aux_Error( ERROR_INFO, "please provide at least one step index in the step table !!\n" );

         DumpTable_NDump   = line;                 // record the number of step indices
         DumpTable[line]   = __INT_MAX__;          // set the next step as an extremely large number

         if ( DumpTable[line-1] < END_STEP )
         {
            END_T          = DumpTable[line-1];    // reset the ending time as the time of the last step

            if ( MPI_Rank == 0 )
               Aux_Message( stdout, "NOTE : the END_STEP is reset to the time of the last step index = %de\n",
                            END_STEP );
         }


//       verify the loaded dump table
         for (int t=1; t<=line; t++)
         {
            if ( DumpTable[t] < DumpTable[t-1] )
               Aux_Error( ERROR_INFO, "values recorded in \"%s\" must be monotonically increasing !!\n",
                          FileName );
         }

         break;

      } // if ( input_line[0] == 42 )
   } // for (line=0; line<MaxLine; line++)


   if ( line == MaxLine )
      Aux_Error( ERROR_INFO, "please prepare a symbol * in the end of the file <%s> !!\n", FileName );


   fclose( File );

   if ( input_line != NULL )     free( input_line );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_Load_StepTable ... done\n" );

} // FUNCTION : Init_Load_DumpTable
