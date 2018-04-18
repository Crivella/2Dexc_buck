/*///////////////////////////////////////////////////////////////////////////////////////////////////
PARAMETERS:	1- input filename	default: none
		2- output filename	default: eps_real.dat
		3- mode(define opertion to perform on data)
			1) eps_real         --> eps1 = 1 + ....
			3) restrict phase in (-PI/2,PI/2]
/////////////////////////////////////////////////////////////////////////////////////////////////// */
#define _BSD_SOURCE
#include <stdio.h>
#include <string.h>
#include <qepp/qepp_main.h>		//~~~

#define UP_FOLDER_LVL 1

int main(int argc,char *argv[])
{
#ifdef __MPI
	initialize_mpi_data();
#endif
	OPEN_IO_ENV( 0, 0, 5);


	char filename[128];
	char outname[128]="eps_real.dat";
	FILE *read;
	if(argc==1)
	{
		QEPP_PRINT("No input file specified. Exiting...\n");
		return -1;
	}
	strcpy(filename,argv[1]);
	read = fopen(filename,"r");
	if( read == NULL)
	{
		QEPP_PRINT("Cannot open file \"%s\"\nExiting...\n",filename);
		return -1;
	}
	printf("filename: %s\n",filename);
	if(argc>2)
		strcpy(outname,argv[2]);

	int mode=0;
	if(argc>3)
		mode=atoi(argv[3]);

	QEPP_PRINT("Outname: %s\n",outname);

	opt_data * data;
	READ( filename, &data, "#");
	opt_data * app = NULL;
	long int n_pt = data->n_pt;
	int n_col = data->n_col;
	int ar=0;
	//*
	if( data->x[0] > 0 && data->x[0] < 1.E-2)
	{
		QEPP_PRINT("Generating antiresonant part...\n");
		ar = 1;
		app = data;
		long int lim = n_pt*2;
		data = initialize_opt_data( lim, n_col);
		for( long int i=0; i<n_pt; i++)
		{
			data->x[i] = app->x[n_pt-i-1] * (-1);
			for( int j=0; j<n_col; j++)
				data->values[i][j] = app->values[n_pt-i-1][j] * (-1);
		}
		for( long int i=0; i<n_pt; i++)
		{
			data->x[i+n_pt] = app->x[i];
			for( int j=0; j<n_col; j++)
				data->values[i+n_pt][j] = app->values[i][j];
		}
		n_pt *= 2;
		FREE( app);
	}// */

	switch( mode)
	{
	case 3:
		QEPP_PRINT("CASE 3---------------\n");
		QEPP_PRINT("Restricting phase in (-PI/2,PI/2]...\n");
		for( long int i=0; i<n_pt; i++)
			for( int j=0; j<n_col; j++)
			{
				while( data->values[i][j] > PI/2)
					data->values[i][j] -= PI;
				while( data->values[i][j] <= -PI/2)
					data->values[i][j] += PI;
			}
		break;
	}
	opt_data * new = apply_kk_im(data);
	switch( mode)
	{
	case 1:
		QEPP_PRINT("CASE 1---------------\n");
		QEPP_PRINT("Adding 1 to result to obtain eps1 from eps2...\n");
		for( long int i=0; i<n_pt; i++)
			for( int j=0; j<n_col; j++)
				new->values[i][j] += 1;
		break;
	}
	if( ar == 1)
	{
		n_pt /= 2;
		app = new;
		new = initialize_opt_data( n_pt, n_col);
		for( long int i=0; i<n_pt; i++)
		{
			new->x[i] = app->x[n_pt+i];
			for( int j=0; j<n_col; j++)
				new->values[i][j] = app->values[i+n_pt][j];	
		}
		FREE( app);
	}

	QEPP_PRINT( "Saving output data in %s\n",outname);
	FILE * out = fopen( outname, "w");
	QEPP_OUT( out, "#File generated from:\n#");
	for( int i=0; i<argc; i++)
		QEPP_OUT( out, "%s ", argv[i]);
	QEPP_OUT( out, "\n");
	PRINT_DATA( new, out);
	fclose( out);

	FREE( data);
	FREE( new);
	fclose( read);

#ifdef __MPI
	FREE( mpi);
#endif
	CLOSE_IO_ENV();
	return 0;
}










