
#include <stdlib.h>
#include <stdio.h>
#include <check.h>
#include <FTAPI.h>



// ---------------------------
// Test for FTAPI_init
// TODO: document items tested
// and write more tests
// ---------------------------
START_TEST(test_fapi_init)
{
   	// Check some good parameters that it works
	{
		int dim = 3;
		int proc_grid[3] = {1, 1, 1};
		double  xmin   =  0.0,  ymin   =  0.0,  zmin   =  0.0;
		double  xmax   =  1.0,  ymax   =  1.0,  zmax   =  1.0;
		int     isize  =  32,   jsize  =  32,   ksize  =  32;
		int     ibuf   =  4,    jbuf   =  4,    kbuf   =  4;
		int xlbdry = FT_BOUNDARY_PERIODIC, 
			ylbdry = FT_BOUNDARY_PERIODIC,
			zlbdry = FT_BOUNDARY_PERIODIC;
		int xubdry = FT_BOUNDARY_PERIODIC,
			yubdry = FT_BOUNDARY_PERIODIC,
			zubdry = FT_BOUNDARY_PERIODIC;
		int geometry = FT_GEOMETRY_CARTESIAN;
		int restart = 0;
		char *restart_filename = NULL;
		int restart_filename_length = 0;


		FTAPI_Init(&dim, proc_grid,
		    &xmin, &xmax, &ymin, &ymax, &zmin, &zmax,
		    &isize, &jsize, &ksize, &ibuf, &jbuf, &kbuf,
		    &xlbdry, &xubdry, &ylbdry, &yubdry, &zlbdry, &zubdry,
		    NULL, NULL, NULL, NULL,
		    &geometry, &restart,
		    restart_filename, restart_filename_length);
		//-------- Check that it went okay.
	}
	// TODO: check other boundary types and geomtery types
	// TODO: error if peroiodic is not on both sides.

}
END_TEST

// ---------------------------
// Test for FTAPI_create_surface
// TODO: document items tested
// and write more tests
// ---------------------------
double test_fapi_create_surface_level_func(void *params, double *coords)
{
	//y=0.5
	return coords[0]-0.5;
}
int test_fapi_create_surface_vel_func(void *params, double *coords, double *vel)
{
	// zero vel
	int d;
	for(d=0; d<3; ++d)
		vel[d]=0.0f;
	return 1;
}

START_TEST(test_fapi_create_surface)
{
	//set up
	{
		int dim = 3;
		int proc_grid[3] = {1, 1, 1};
		double  xmin   =  0.0,  ymin   =  0.0,  zmin   =  0.0;
		double  xmax   =  1.0,  ymax   =  1.0,  zmax   =  1.0;
		int     isize  =  32,   jsize  =  32,   ksize  =  32;
		int     ibuf   =  4,    jbuf   =  4,    kbuf   =  4;
		int xlbdry = FT_BOUNDARY_PERIODIC, 
			ylbdry = FT_BOUNDARY_PERIODIC,
			zlbdry = FT_BOUNDARY_PERIODIC;
		int xubdry = FT_BOUNDARY_PERIODIC,
			yubdry = FT_BOUNDARY_PERIODIC,
			zubdry = FT_BOUNDARY_PERIODIC;
		int geometry = FT_GEOMETRY_CARTESIAN;
		int restart = 0;
		char *restart_filename = NULL;
		int restart_filename_length = 0;


		FTAPI_Init(&dim, proc_grid,
		    &xmin, &xmax, &ymin, &ymax, &zmin, &zmax,
		    &isize, &jsize, &ksize, &ibuf, &jbuf, &kbuf,
		    &xlbdry, &xubdry, &ylbdry, &yubdry, &zlbdry, &zubdry,
		    NULL, NULL, NULL, NULL,
		    &geometry, &restart,
		    restart_filename, restart_filename_length);
		//-------- Check that it went okay.
	}
	{
		//CHECK createsurface
		int pos_comp = 1;
		int neg_comp = 2;
		int surf_type = 1;


		FTAPI_CreateSurface(pos_comp, neg_comp, test_fapi_create_surface_level_func, NULL,
			surf_type,  test_fapi_create_surface_vel_func, NULL);
	}
}
END_TEST

// ---------------------------
// Test for FTAPI_propagate
// TODO: document items tested
// and write more tests
// ---------------------------
START_TEST(test_fapi_propagate)
{
	//set up
	{
		int dim = 3;
		int proc_grid[3] = {1, 1, 1};
		double  xmin   =  0.0,  ymin   =  0.0,  zmin   =  0.0;
		double  xmax   =  1.0,  ymax   =  1.0,  zmax   =  1.0;
		int     isize  =  32,   jsize  =  32,   ksize  =  32;
		int     ibuf   =  4,    jbuf   =  4,    kbuf   =  4;
		int xlbdry = FT_BOUNDARY_PERIODIC, 
			ylbdry = FT_BOUNDARY_PERIODIC,
			zlbdry = FT_BOUNDARY_PERIODIC;
		int xubdry = FT_BOUNDARY_PERIODIC,
			yubdry = FT_BOUNDARY_PERIODIC,
			zubdry = FT_BOUNDARY_PERIODIC;
		int geometry = FT_GEOMETRY_CARTESIAN;
		int restart = 0;
		char *restart_filename = NULL;
		int restart_filename_length = 0;


		FTAPI_Init(&dim, proc_grid,
		    &xmin, &xmax, &ymin, &ymax, &zmin, &zmax,
		    &isize, &jsize, &ksize, &ibuf, &jbuf, &kbuf,
		    &xlbdry, &xubdry, &ylbdry, &yubdry, &zlbdry, &zubdry,
		    NULL, NULL, NULL, NULL, //TODO: remove these.
		    &geometry, &restart,
		    restart_filename, restart_filename_length);
		//-------- Check that it went okay.
	}
	{
		//have to createsurface
		int pos_comp = 1;
		int neg_comp = 2;
		int surf_type = 1;


		FTAPI_CreateSurface(pos_comp, neg_comp, test_fapi_create_surface_level_func, NULL,
			surf_type,  test_fapi_create_surface_vel_func, NULL);
	}
	{
		//CHECK propagate
		double dt = 0.0;
		FTAPI_Propagate(&dt);
		//assert nothing moved
	}
	{
		//int i;
		//for(i=0; i<10; ++i)
		//{
		//	double dt = 1.0;
		//	FTAPI_Propagate(&dt);
		//}
	}
}
END_TEST

// ---------------------------
// Test for FTAPI_IO
// TODO: document items tested
// and write more tests
// ---------------------------
START_TEST(test_fapi_io)
{
}
END_TEST

//---------------------------------------------------
//---------------------------------------------------

// -----------------------------
// Test for FTAPI_cartesian_init
// TODO: document items tested
// and write more tests
// ---------------------------
START_TEST(test_fapi_init_cartesian)
{
	{
		int dim = 3;
		int proc_grid[3] = {1, 1, 1};
		double  xmin   =  0.0,  ymin   =  0.0,  zmin   =  0.0;
		double  xmax   =  1.0,  ymax   =  1.0,  zmax   =  1.0;
		int     isize  =  32,   jsize  =  32,   ksize  =  32;
		int     ibuf   =  4,    jbuf   =  4,    kbuf   =  4;
		int xlbdry = FT_BOUNDARY_PERIODIC, 
			ylbdry = FT_BOUNDARY_PERIODIC,
			zlbdry = FT_BOUNDARY_PERIODIC;
		int xubdry = FT_BOUNDARY_PERIODIC,
			yubdry = FT_BOUNDARY_PERIODIC,
			zubdry = FT_BOUNDARY_PERIODIC;
		int geometry = FT_GEOMETRY_CARTESIAN;
		int restart = 0;
		char *restart_filename = NULL;
		int restart_filename_length = 0;
		FNAME(ftapi_init_cartesian,FTAPI_INIT_CARTESIAN)(&dim,
			proc_grid,
        // min / max of real interior domain limits
        &xmin, &xmax,
        &ymin, &ymax,
        &zmin, &zmax,
        // number of interior cells in each direction
        &isize, &jsize, &ksize,
        // number of buffer cells in direction
        // TODO: should be 6 (lower and upper)
        &ibuf, &jbuf, &kbuf,
        // boundary type (FT_BOUNDARY TYPE enum)
        &xlbdry, &xubdry,
        &ylbdry, &yubdry,
        &zlbdry, &zubdry,
        // Velocity arrays
        //double *vx, double *vy, double *vz,
		// FIXME: NULL VELOCITY ARRAYS
		NULL, NULL, NULL,
		test_fapi_create_surface_level_func, NULL,
        // geomerty type
        &geometry,
        &restart,
		restart_filename,
		restart_filename_length);
	}
}
END_TEST

START_TEST(test_fapi_cartesian_create_surface)
{
	{
		int dim = 3;
		int proc_grid[3] = {1, 1, 1};
		double  xmin   =  0.0,  ymin   =  0.0,  zmin   =  0.0;
		double  xmax   =  1.0,  ymax   =  1.0,  zmax   =  1.0;
		int     isize  =  32,   jsize  =  32,   ksize  =  32;
		int     ibuf   =  4,    jbuf   =  4,    kbuf   =  4;
		int xlbdry = FT_BOUNDARY_PERIODIC, 
			ylbdry = FT_BOUNDARY_PERIODIC,
			zlbdry = FT_BOUNDARY_PERIODIC;
		int xubdry = FT_BOUNDARY_PERIODIC,
			yubdry = FT_BOUNDARY_PERIODIC,
			zubdry = FT_BOUNDARY_PERIODIC;
		int geometry = FT_GEOMETRY_CARTESIAN;
		int restart = 0;
		char *restart_filename = NULL;
		int restart_filename_length = 0;
		FNAME(ftapi_init_cartesian,FTAPI_INIT_CARTESIAN)(&dim,
			proc_grid,
        // min / max of real interior domain limits
        &xmin, &xmax,
        &ymin, &ymax,
        &zmin, &zmax,
        // number of interior cells in each direction
        &isize, &jsize, &ksize,
        // number of buffer cells in direction
        // TODO: should be 6 (lower and upper)
        &ibuf, &jbuf, &kbuf,
        // boundary type (FT_BOUNDARY TYPE enum)
        &xlbdry, &xubdry,
        &ylbdry, &yubdry,
        &zlbdry, &zubdry,
        // Velocity arrays
        //double *vx, double *vy, double *vz,
		// FIXME: NULL VELOCITY ARRAYS
		NULL, NULL, NULL,
		test_fapi_create_surface_level_func, NULL,
        // geomerty type
        &geometry,
        &restart,
		restart_filename,
		restart_filename_length);
	}
	{
		//CHECK createsurface
		int pos_comp = 1;
		int neg_comp = 2;
		int surf_type = 1;


		FTAPI_CreateSurface_cartesian(pos_comp, neg_comp, test_fapi_create_surface_level_func, NULL, surf_type);
	}
}
END_TEST

START_TEST(test_fapi_propagate_cartesian)
{
	int size;
	{
		int dim = 3;
		int proc_grid[3] = {1, 1, 1};
		double  xmin   =  0.0,  ymin   =  0.0,  zmin   =  0.0;
		double  xmax   =  1.0,  ymax   =  1.0,  zmax   =  1.0;
		int     isize  =  32,   jsize  =  32,   ksize  =  32;
		int     ibuf   =  4,    jbuf   =  4,    kbuf   =  4;
		size = isize*jsize*ksize;

		int xlbdry = FT_BOUNDARY_PERIODIC, 
			ylbdry = FT_BOUNDARY_PERIODIC,
			zlbdry = FT_BOUNDARY_PERIODIC;
		int xubdry = FT_BOUNDARY_PERIODIC,
			yubdry = FT_BOUNDARY_PERIODIC,
			zubdry = FT_BOUNDARY_PERIODIC;
		int geometry = FT_GEOMETRY_CARTESIAN;
		int restart = 0;
		char *restart_filename = NULL;
		int restart_filename_length = 0;
		FNAME(ftapi_init_cartesian,FTAPI_INIT_CARTESIAN)(&dim,
			proc_grid,
        // min / max of real interior domain limits
        &xmin, &xmax,
        &ymin, &ymax,
        &zmin, &zmax,
        // number of interior cells in each direction
        &isize, &jsize, &ksize,
        // number of buffer cells in direction
        // TODO: should be 6 (lower and upper)
        &ibuf, &jbuf, &kbuf,
        // boundary type (FT_BOUNDARY TYPE enum)
        &xlbdry, &xubdry,
        &ylbdry, &yubdry,
        &zlbdry, &zubdry,
        // Velocity arrays
        //double *vx, double *vy, double *vz,
		// FIXME: NULL VELOCITY ARRAYS
		NULL, NULL, NULL,
		test_fapi_create_surface_level_func, NULL,
        // geomerty type
        &geometry,
        &restart,
		restart_filename,
		restart_filename_length);
	}
	int surf_type = 1;
	{
		//createsurface
		int pos_comp = 1;
		int neg_comp = 2;

		FTAPI_CreateSurface_cartesian(pos_comp, neg_comp, test_fapi_create_surface_level_func, NULL, surf_type);
	}
	{
		// CHECK propagate
		double *xvel = malloc(sizeof(double)*size);
		double *yvel = malloc(sizeof(double)*size);
		double *zvel = malloc(sizeof(double)*size);
		int i;

		for(i=0; i<size; ++i)
		{
			xvel[i] = yvel[i] = zvel[i] = 0.0f;
		}
		double dt = 0.0;
		FNAME(ftapi_propagate_cartesian,FTAPI_PROPAGATE_CARTESIAN)(&dt, surf_type, xvel, yvel, zvel);
		//assert nothing moved
	}
}
END_TEST


Suite *ftapi_suite (void)
{
	Suite *s = suite_create ("FTAPI");

	/* Core test case */
	TCase *tc_core = tcase_create ("Core");
	tcase_add_test (tc_core, test_fapi_init);
	tcase_add_test (tc_core, test_fapi_create_surface);
	tcase_add_test (tc_core, test_fapi_propagate);
	tcase_add_test (tc_core, test_fapi_io);
	suite_add_tcase (s, tc_core);

	return s;
}

Suite *ftapi_cartesian_suite (void)
{
	Suite *s = suite_create ("FTAPI_cartesian");

	/* Core test case */
	TCase *tc_core = tcase_create ("Core");
	tcase_add_test (tc_core, test_fapi_init_cartesian);
	tcase_add_test (tc_core, test_fapi_propagate_cartesian);
	//tcase_add_test (tc_core, test_fapi_create_surface);
	//tcase_add_test (tc_core, test_fapi_io);
	suite_add_tcase (s, tc_core);

	return s;
}

int main (void)
{
	int number_failed;
	Suite *s1 = ftapi_suite ();
	Suite *s2 = ftapi_cartesian_suite ();
	SRunner *sr1 = srunner_create (s1);
	SRunner *sr2 = srunner_create (s2);
	srunner_run_all (sr1, CK_NORMAL);
	srunner_run_all (sr2, CK_NORMAL);
	number_failed = srunner_ntests_failed(sr1) +  srunner_ntests_failed(sr2);;
	srunner_free (sr1);
	srunner_free (sr2);
	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
