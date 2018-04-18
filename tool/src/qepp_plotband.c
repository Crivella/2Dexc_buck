//gcc -o ~/bin/my_plotband.x my_plotband.c -lm -I./include/
/*///////////////////////////////////////////////////////////////////////////////////////////////////
PARAMETERS:	1- input filename	default: bands_1.out
		2- output filename	default: plotted.dat
/////////////////////////////////////////////////////////////////////////////////////////////////// */
#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <qepp/qepp_main.h>		//~~~
//#include <gnuplot_i.h>

int main(int argc, char *argv[])
{
#ifdef __MPI
	initialize_mpi_data();
#endif
	OPEN_IO_ENV( 0, 0, 5);

	char filename[64]="bands_1.out";
	char outname[64]="plotted.dat";
	char outagr[64];
	if(argc>1)
		strcpy(filename,argv[1]);
	if(argc>2)
		strcpy(outname,argv[2]);
	else
	{
//Make outputfile name
		char * path=realpath(filename,NULL);
		char * tokens[128];
		tokens[0]=strtok(path,"/");
		int count=0;
		while(tokens[count] != NULL)
			tokens[++count]=strtok(NULL,"/");
		count--;
		while(strstr(tokens[count],"bands_") != NULL)
			count--;
		strcpy(outname,tokens[count]);
		strcpy(outagr,outname);
		strcat(outname,"_plotted.dat");
		strcat(outagr,"_plotted.pdf");
	}
//printf("Executing %s -> %s\n",filename,outname);

//Read data from nscf.out
	nscf_data * data;
	if( READ( filename, &data) != SUCCESS)
	{
		band_data * data1;
		if( READ( filename, &data1) != SUCCESS)
		{
			QEPP_PRINT( "%s in neither a pw.x or band.x output\n", filename);
			exit( 1);
		}
		QEPP_PRINT( "Printing output in \"%s\"\n", outname);
		QEPP_PRINT_BAND_STRUCTURE( data1, outname);
	}
	else
	{
		QEPP_PRINT( "Printing output in \"%s\"\n", outname);
		QEPP_PRINT_BAND_STRUCTURE( data, outname);
	}		
	if( data == NULL)	
	{
		fprintf( stderr, "Failed to load the data from %s\n",filename);
		return 1;
	}

//Plot band structure
	/*char command[256]="xmgrace -nxy ";
	strcat(command,outname);
	strcat(command," -saveall ");
	strcat(command,outagr);
	if( system(command)) {}*/
	//gnuplot_ctrl * h1 = gnuplot_init();
	//gnuplot_cmd( h1, "set term wxt enhanced");
	//gnuplot_cmd( h1, "plot for[i=2:%d] \'%s\' u 1:(column(i)) w l notitle", data->n_bnd, outname);
	//getchar();
	//gnuplot_close( h1);

	FREE(data);

#ifdef __MPI
	FREE( mpi);
#endif
	CLOSE_IO_ENV();

	return 0;
}
















