//gcc -o ~/bin/eps_average.x eps_average.c -lm -I./include/
/*///////////////////////////////////////////////////////////////////////////////////////////////////
PARAMETERS:	1- filename #1		default: none
		2- filename #2
		3- ...
/////////////////////////////////////////////////////////////////////////////////////////////////// */
#include <stdio.h>
#include <string.h>
#include <qepp/qepp_main.h>		//~~~

int main(int argc,char *argv[])
{
#ifdef __MPI
	initialize_mpi_data();
#endif
	OPEN_IO_ENV( 0, 0, 5);

	int n_files=argc-1;
	char outname[128]="averaged.dat";

	opt_data * data[n_files];
	
//Read the files
	for(int i=0; i<n_files; i++)
	{
		READ (argv[i+1], &data[i], "#");
		if(data[i] == NULL)
			return -1;
	}
//Check for incongruencies (diff number of points, different x values, different numb of columns)
	for( int i=0; i<n_files-1; i++)
		if( data[i]->n_pt != data[i+1]->n_pt)
		{
			QEPP_PRINT("ERROR: The files contains a different number of points\nExiting.....\n");
			return -1;
		}
	for( int i=0; i<n_files-1; i++)
		if( data[i]->n_col != data[i+1]->n_col)
		{
			QEPP_PRINT("ERROR: The files contains a different number of columns\nExiting.....\n");
			return -1;
		}
	for( long int i=0; i<data[0]->n_pt; i++)
		for( int i1=0; i1<n_files-1; i1++)
			if( data[i1]->x[i] != data[i1+1]->x[i])
			{
				QEPP_PRINT("ERROR: The files dont list the same x\nExiting.....\n");
				return -1;
			}

//Calculate average
	/*pass_data pass=MAKE_PASS(data[0]->n_pt,data[0]->n_col);
	opt_data * averaged = INIT( &pass, ID_OPT_DATA);*/
	opt_data * averaged = initialize_opt_data( data[0]->n_pt, data[0]->n_col);

	for( long int i=0; i<data[0]->n_pt; i++)
	{
		averaged->x[i]=data[0]->x[i];
		for( int i1=0; i1<data[0]->n_col; i1++)
			averaged->values[i][i1]=0;
		for( int i1=0; i1<n_files; i1++)
			for( int i2=0; i2<data[0]->n_col; i2++)
				averaged->values[i][i2]+=data[i1]->values[i][i2];
		for( int i1=0; i1<data[0]->n_col; i1++)
			averaged->values[i][i1]/=n_files;
	}

	FILE * out = fopen( outname, "w");
	QEPP_OUT( out, "#File generated from:\n#");
	for( int i=0; i<argc; i++)
		QEPP_OUT( out, "%s ", argv[i]);
	QEPP_OUT( out, "\n");
	PRINT_DATA( averaged, out);
	fclose( out);

	for( int i=0; i<n_files; i++)
		FREE(data[i]);
	FREE(averaged);

#ifdef __MPI
	FREE( mpi);
#endif
	CLOSE_IO_ENV();
	return 0;
}










