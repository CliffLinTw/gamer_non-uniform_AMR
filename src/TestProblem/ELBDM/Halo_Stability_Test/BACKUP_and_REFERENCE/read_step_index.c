#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAX_STRING 1024

int main(void)
{
    int row_counter = 0;
    char input_filename[MAX_STRING], check_buff_1[MAX_STRING], check_buff_2[MAX_STRING]; 
    sprintf(input_filename, "./Input__StepTable");

    FILE *file_step_idx = fopen("./Input__StepTable", "r");
    if (!file_step_idx)
    {
       printf("Cannot find input file %s !! Exit!!\n", input_filename);
       exit(1);
    }
    else
       printf("Input file %s found.\n", input_filename);
   
   fscanf(file_step_idx, "%s %s", check_buff_1, check_buff_2);
   if ( strcmp(check_buff_1, "Output_ID")!=0 ) 
   {
      printf("First row Check failed!! First column %s != Output_ID !! Exit!!\n");
      exit(1);
   }
   else if ( strcmp(check_buff_2,"Step_Index")!=0 )
   {
      printf("First row Check failed!! Second column %s != Step_Index !! Exit!!\n");
      exit(1);
   }
   else
      printf("First column is %s ; Second column is %s . Check pass.\n", check_buff_1, check_buff_2);
   while (!feof(file_step_idx))
   {
      fgets(check_buff_1, MAX_STRING, file_step_idx);
      printf("%s", check_buff_1);
   }   

   fclose(file_step_idx);
   return EXIT_SUCCESS;
}
