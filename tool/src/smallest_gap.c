//gcc -o ~/bin/smallest_gap.x smallest_gap.c -lm -I./include/
/*///////////////////////////////////////////////////////////////////////////////////////////////////
PARAMETERS:	1- input filename	default: nscf_1.out
		2/3/4 - center coord
		5 - radius
/////////////////////////////////////////////////////////////////////////////////////////////////// */
#include <stdio.h>
#include <string.h>
#include <qepp/qepp_main.h>		//~~~

int main(int argc, char *argv[])
{
#ifdef __MPI
	initialize_mpi_data();
#endif
	OPEN_IO_ENV( 0, 0, 1);
	char fileapp[128] = "nscf_1.out";
	double radius=0;
	double * comp_point=NULL;
	if(argc>1)
		strcpy(fileapp,argv[1]);
	char * filename = get_file( fileapp, ".out");
	if( filename == NULL)
		return 1;
	nscf_data * data;
	if(argc>2)
		if(argc<=4)
		{
			fprintf( stderr, "Not enough parameters to specify coordinates...\n");
			return 1;
		}
	if(argc>4)
	{
		if(argc<=5)
		{
			fprintf( stderr, "Radius not specified...\n");
			return 1;
		}
		comp_point=malloc( 3 * sizeof(double));
		for(int i=0; i<3; i++)
			comp_point[i]=atof(argv[2+i]);
	}
	if(argc>5)
	{
		radius=atof(argv[5]);
		QEPP_PRINT( "Reading coordinates around %g %g %g with radius %g\n",comp_point[0],comp_point[1],comp_point[2],radius);
	}

//Read data from nscf.out
	READ( filename, &data);
	if(argc>6)
		data->md->e_fermi=atof(argv[6]);

	QEPP_PRINT("E_fermi(from file):\t%f eV\n",data->md->e_fermi);
	double e_fermi;
	parse_errh( find_e_fermi( &e_fermi, data));
	QEPP_PRINT("E_fermi(from calc):\t%f eV\n", e_fermi);

	int vb;
	if( !data->md->spin_orbit)
	{
		QEPP_PRINT("No spin-orbit correction\n");
		vb = ceil(data->md->n_el/2)-1;
	}
	else
	{
		QEPP_PRINT("spin-orbit correction detected\n");
		vb = data->md->n_el-1;
	}
	int cb = vb+1;
	printf("vb = %d, cb = %d\n",vb+1,cb+1);
	
	long int mve;
	parse_errh( find_max_v_bnd( &mve, data, comp_point, radius));
	long int mce;
	parse_errh( find_min_c_bnd( &mce, data, comp_point, radius));
	if( mve >=0)
	{
		QEPP_PRINT("\nMax_vb_energy: vb= %lf cb= %lf eV\n",data->energies[mve][vb],data->energies[mve][vb+1]);
		QEPP_PRINT("\tat   %lf %lf %lf (2pi/a) (# %li)\n",data->kpt[mve][0],data->kpt[mve][1],data->kpt[mve][2],mve+1);
	}
	else
		QEPP_PRINT("Max_vb_energy: FAILED\n");
	if( mce >=0)
	{
		QEPP_PRINT("Min_cb_energy: vb= %lf cb= %lf eV\n",data->energies[mce][vb],data->energies[mce][vb+1]);
		QEPP_PRINT("\tat   %lf %lf %lf (2pi/a) (# %li)\n",data->kpt[mce][0],data->kpt[mce][1],data->kpt[mce][2],mce+1);
	}
	else
		QEPP_PRINT("Min_cb_energy: FAILED\n");

	double ef = data->md->e_fermi;
	long int mog;
	parse_errh( find_min_opt_gap( &mog, data, comp_point, radius));
	if( mog >= 0)
	{
		double min_opt_gap = data->energies[mog][vb+1] - data->energies[mog][vb];
		QEPP_PRINT("\nMin_opt_gap: %lf eV\n\tat   %lf %lf %lf (2pi/a) (# %li)\n",min_opt_gap,data->kpt[mog][0],data->kpt[mog][1],
			data->kpt[mog][2], mog+1);
		QEPP_PRINT("\t%lf -> %lf   Ef: %lf eV\n", data->energies[mog][vb], data->energies[mog][vb+1],ef);
	}
	else
		QEPP_PRINT("\nMin_opt_gap: FAILED\n");

	long int mg;
	parse_errh( find_smallest_gap( &mg, data, comp_point, radius));
	double mge = data->energies[mg][cb] - data->energies[mg][vb];
	QEPP_PRINT("\nMin gap energy: %lf eV\n",mge);
	QEPP_PRINT("\tat %f %f %f (2pi/a) (# %li)\n",data->kpt[mg][0],data->kpt[mg][1],data->kpt[mg][2],mg+1);
	QEPP_PRINT("\t%lf -> %lf   Ef: %lf eV\n", data->energies[mg][vb], data->energies[mg][vb+1],ef);

	parse_errh( find_smallest_gap2( &mg, data, comp_point, radius));
	mge = data->energies[mg][cb+1] - data->energies[mg][vb];
	QEPP_PRINT("\nMin gap energy(vb->cb+1): %lf eV\n",mge);
	QEPP_PRINT("\tat %f %f %f (2pi/a) (# %li)\n",data->kpt[mg][0],data->kpt[mg][1],data->kpt[mg][2],mg+1);
	QEPP_PRINT("\t%lf -> %lf   Ef: %lf eV\n", data->energies[mg][vb], data->energies[mg][vb+2],ef);

	parse_errh( find_smallest_gap3( &mg, data, comp_point, radius));
	mge = data->energies[mg][cb+1] - data->energies[mg][cb];
	QEPP_PRINT("\nMin gap energy(cb->cb+1): %lf eV\n",mge);
	QEPP_PRINT("\tat %f %f %f (2pi/a) (# %li)\n",data->kpt[mg][0],data->kpt[mg][1],data->kpt[mg][2],mg+1);
	QEPP_PRINT("\t%lf -> %lf   Ef: %lf eV\n", data->energies[mg][cb], data->energies[mg][vb+2],ef);

	free(comp_point);
	FREE(data);

#ifdef __MPI
	FREE( mpi);
#endif
	CLOSE_IO_ENV();

	return 0;
}
















