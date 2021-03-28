#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <FTAPI.h>
#include "fdecs.h"
#include "ftapi.h"
#undef free

#define RK_PRINT(...) //do{printf("%s:%d:",__FILE__,__LINE__);printf(__VA_ARGS__); fflush(stdout);}while(0)

//--------------------------------
// GLOBAL VARIABLES
//--------------------------------

// global front
// This api only provides interaction with a single front.
Front front; 

// Velocity function array.
// Serves as a map from surface ID to velocity function.
int (*VELOCITY_FUNC[256])(POINTER, double*, double*, int*, int*);
void* VELOCITY_PARAMS[256];

//--------------------------------
// PROTOTYPES FOR HELPER FUNCTIONS
//--------------------------------
static int intersection2d(  const double ax1, const double ay1, const double ax2, const double ay2, 
							const double bx1, const double by1, const double bx2, const double by2, 
							double *x, double *y);
static int pts_cross(const double *p1, const double *p2, const double *q1, const double *q2);
static int left(const double *p, const double *q1, const double *q2);
void	   sort_crossing(double [], int);

//--------------------------------
// PASSIVE TRACKING
//--------------------------------

static int curve_velocity_func(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
//    	assert(hs->velocity_func != NULL);
    	//assert(hs->velocity_params != NULL);
    	//return hs->velocity_func(hs->velocity_params, Coords(p), vel);
	int type=hs->surf_type;
    int neg_comp = hs->neg_comp;
    int pos_comp = hs->pos_comp;
	return VELOCITY_FUNC[type](VELOCITY_PARAMS[type], Coords(p), vel, &neg_comp, &pos_comp);
}

// JAM - Added for FTAPI -- FIXME: This should be in a header file or removed
// seems like a bad way to do this in this function
#define FT_BOUNDARY_PERIODIC -135
#define FT_BOUNDARY_REFLECTING -31
#define FT_BOUNDARY_FLOW -32
#define FT_BOUNDARY_SUBDOMAIN -10

static int get_bc(int apiBC)
{
	int ftBC;
	if (apiBC==FT_BOUNDARY_PERIODIC)
	    ftBC = PERIODIC_BOUNDARY;
	else if (apiBC==FT_BOUNDARY_REFLECTING)
	    ftBC =  REFLECTION_BOUNDARY;
	else if (apiBC == FT_BOUNDARY_FLOW)
	    ftBC = OPEN_BOUNDARY;
    else if (apiBC == FT_BOUNDARY_SUBDOMAIN)
        ftBC = SUBDOMAIN_BOUNDARY;
	else
	{
	    printf("FrontTracking: Unsupported boundary type %d\n",apiBC);
	    assert(0);
	}
	return ftBC;
}

void FTAPI_Init(int *dim,
	int * procGrid,
	// min / max of real interior domain limits
        double *xmin, double*xmax, 
	double *ymin, double *ymax, 
        double *zmin, double *zmax,
	// number of interior cells in each direction
        int *isize, int *jsize, int *ksize, 
	// number of buffer cells in direction
	// TODO: should be 6 (lower and upper)
	int *ibuf, int *jbuf, int *kbuf,
	// boundary type (FT_BOUNDARY TYPE enum)
        int *xlbdry, int *xubdry,
        int *ylbdry, int *yubdry,
        int *zlbdry, int *zubdry,
	double (*level_func)( POINTER func_params, double *coords),
	void* level_func_params,
	int (*vel_func)(POINTER, double*, double*, int*, int*),
	void* vel_params,
	// geomerty type
	int *geometry, 
	int *restart,
	char *restart_filename,
	int *restart_filename_length)
{
    int i; // counter


    // flash grid sizes are local sizes. Recreate the global grid size for FT
    //*isize*=procGrid[0];
    //*jsize*=procGrid[1];
    //*ksize*=procGrid[2];

    /// Generate the cmd line args for FT_Init() based pm 
    // do,emsop, partition, etc.
    int argc;
    char **argv;
    static LEVEL_FUNC_PACK level_func_pack;
    F_BASIC_DATA f_basic;
    RK_PRINT("dim=%d\n",*dim);

    if(*dim==3)
    {
        argc=9;
        
        argv = malloc(argc*sizeof(char*));
        for(i=0;i<argc;++i)
            argv[i] = malloc(40*sizeof(char));

        switch(*geometry)
        {
            case FT_GEOMETRY_CYLINDRICAL:
                sprintf(argv[4], "%c", 'c');
                break;
            
            case FT_GEOMETRY_CARTESIAN:
                sprintf(argv[4], "%c", 'r');
                level_func_pack.set_3d_bdry = 1;
                break;

            case FT_GEOMETRY_SPHERICAL:
                sprintf(argv[4], "%c", 's');
                break;
	    default:
	    	printf("Did not understand geometry: %d\n", *geometry);
		clean_up(-1);
        }
    }
    else if (*dim==2)
    {
        argc=8;
        argv = malloc(argc*sizeof(char*));

        for(i=0;i<argc;++i)
            argv[i] = malloc(40*sizeof(char));

        switch(*geometry)
        {
            case FT_GEOMETRY_CYLINDRICAL:
                sprintf(argv[4], "%c", 'c');
                break;
            case FT_GEOMETRY_CARTESIAN:
                sprintf(argv[4], "%c", 'r');
                break;
            case FT_GEOMETRY_SPHERICAL:
                sprintf(argv[4], "%c", 's');
                break;
        }
    }
    else if (*dim==1)
    {
        argc=7;
        argv = malloc(argc*sizeof(char*));

        for(i=0;i<argc;++i)
            argv[i] = malloc(40*sizeof(char));

        switch(*geometry)
        {
            case FT_GEOMETRY_CYLINDRICAL:
                sprintf(argv[4], "%c", 'c');
                break;
            case FT_GEOMETRY_CARTESIAN:
                sprintf(argv[4], "%c", 'r');
                break;
            case FT_GEOMETRY_SPHERICAL:
                sprintf(argv[4], "%c", 's');
                break;
        }
    }
    else
    {
	printf("Invalid Dimension: %d\n", *dim);
	clean_up(-1);
    }
   
    sprintf(argv[0], "./flash4");
    sprintf(argv[1], "-d");
    sprintf(argv[2], "%d",*dim);
    sprintf(argv[3], "-c");
    sprintf(argv[5], "-p");
    sprintf(argv[6], "%d",procGrid[0]);
    if(*dim > 1)
    sprintf(argv[7], "%d",procGrid[1]);
    if(*dim==3)
        sprintf(argv[8], "%d",procGrid[2]); // 3d
    f_basic.dim = *dim; 
    

    FT_Init(argc, argv, &f_basic);

    for(i=0;i<argc;++i)free(argv[i]);
    free(argv);

    // Flash gives local coordinates. Get the global domain limits.
    pp_global_max(xmax,1L);
    pp_global_max(ymax,1L);
    pp_global_max(zmax,1L);

    pp_global_min(xmin,1L);
    pp_global_min(ymin,1L);
    pp_global_min(zmin,1L);
    
    // Set the init parameters
    f_basic.gmax_new_pp_grid[0] = procGrid[0];
    f_basic.gmax_new_pp_grid[1] = procGrid[1];
    f_basic.gmax_new_pp_grid[2] = procGrid[2];

    f_basic.L[0]    = *xmin;     f_basic.L[1]    = *ymin;
    f_basic.U[0]    = *xmax;     f_basic.U[1]    = *ymax;
    f_basic.L[2]    = *zmin;     f_basic.U[2]    = *zmax;
    f_basic.gmax[0] = *isize;    f_basic.gmax[1] = *jsize;
    f_basic.gmax[2] = *ksize;
    f_basic.buffer_size[0] = *ibuf;
    f_basic.buffer_size[1] = *jbuf;
    f_basic.buffer_size[2] = *kbuf;
    RK_PRINT("Gmax = %d %d %d\n", 
	    f_basic.gmax[0],f_basic.gmax[1],f_basic.gmax[2]);

    // Convert the flash boundary types to FT enum types for FT init
    //FIXME: move call to other side.
    f_basic.boundary[0][0] = get_bc(*xlbdry);
    f_basic.boundary[0][1] = get_bc(*xubdry);
    if(*dim > 1)
    {
    f_basic.boundary[1][0] = get_bc(*ylbdry);
    f_basic.boundary[1][1] = get_bc(*yubdry);
    }
    // 3d
    if(*dim==3)
    {
	f_basic.boundary[2][0] = get_bc(*zlbdry);
	f_basic.boundary[2][1] = get_bc(*zubdry);
    }
    if(*restart == 1)
    {
	strncpy(f_basic.restart_name, restart_filename, *restart_filename_length);
	f_basic.restart_name[*restart_filename_length] = '\0';
	//strcat(f_basic.restart_name,"/intfc-ts0001790");
	if(pp_numnodes()>1)
	    sprintf(f_basic.restart_name,"%s-p%04d",f_basic.restart_name,pp_mynode());
	f_basic.RestartRun = YES;
    }

    f_basic.size_of_intfc_state = 0;

    FT_StartUp(&front, &f_basic);
    //----------- interface initialization --------------------------
    if(*restart==0)
    {
	level_func_pack.neg_component = 1;
	level_func_pack.pos_component = 2;
	level_func_pack.func_params = level_func_params;

	level_func_pack.func = level_func;

	level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;
	
	FT_InitIntfc(&front,&level_func_pack);
    }

    // JAM - Added for FTAPI
    //next_point(&front.intfc, NULL, NULL, NULL);
    //while(next_point(&front.intfc,&p,&hse,&hs))
    //{
    //    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
    //}

    //---------------- velocity initialization ---------------------
    VELO_FUNC_PACK velo_func_pack;

    velo_func_pack.func = curve_velocity_func;

    //velo_func_pack.point_propagate = first_order_point_propagate;
    velo_func_pack.point_propagate = second_order_point_propagate;
    //velo_func_pack.point_propagate = fourth_order_point_propagate;
    FT_InitVeloFunc(&front,&velo_func_pack);

    // JAM - Added for FTAPI FIXME: I don't understand the logic here
    // It seems to me that the api is only set up for 3D and surfaces
    // In addition, it isn't clear why it needs this to be set and why
    // this is different than the velo_func_pack set above...
    
    int surf_type = 0;
    VELOCITY_FUNC[surf_type] = vel_func;
    VELOCITY_PARAMS[surf_type] = vel_params;

    
    // Further init
    // JAM : Added 4/4/16
    // Spacing between grid points, default is 0.75
    Front_spacing(&front,GENERAL_WAVE) = 0.75;
    Front_spacing(&front,VECTOR_WAVE) = 0.75;
    
    FT_RedistMesh(&front);
    RECT_GRID *gr =  &topological_grid(front.interf);
    front.dt = 0.0f;

    FT_MakeGridIntfc(&front);
    //FT_Propagate(&front);
    FT_SetTimeStep(&front); // Set the timestep based off of FT
    Frequency_of_redistribution(&front,GENERAL_WAVE) = 20;
    Redistribution_count(&front) = front.step; 
}


void FTAPI_CreateSurface(int pos_comp, int neg_comp, double (*level_func)(void*, double*), void* level_params, int surf_type,
	int (*vel_func)(POINTER, double*, double*, int*, int*),
	void* vel_params)
{
    	int num_segs;
	int surf_type_ft;
	CURVE **c;
	c = (CURVE**)FT_CreateLevelHyperSurfs(front.rect_grid, front.interf, neg_comp, pos_comp, level_func, level_params, FIRST_PHYSICS_WAVE_TYPE, &num_segs);
	for(;c&&*c;++c)
	{
	    //(*c)->user_type = surf_type;

	    Hyper_surf(*c)->surf_type=surf_type;
	    //Hyper_surf(*c)->velocity_func = vel_func;
	    //Hyper_surf(*c)->velocity_params = vel_params;
	    VELOCITY_FUNC[surf_type] = vel_func;
	    VELOCITY_PARAMS[surf_type] = vel_params;
	}

	LEVEL_FUNC_PACK level_func_pack;
//	level_func_pack.neg_component = 1;
//	level_func_pack.pos_component = 2;
//	level_func_pack.func_params = level_func_params;

//	level_func_pack.func = level_func;
	
//	level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;
	memset(&level_func_pack, 0, sizeof(level_func_pack));
	
	FT_InitIntfc(&front,&level_func_pack);

	// Further init
	FT_RedistMesh(&front);
	RECT_GRID *gr =  &topological_grid(front.interf);
	front.dt = 0.0f;

	FT_MakeGridIntfc(&front);

	FT_Propagate(&front);
	FT_SetTimeStep(&front); // Set the timestep based off of FT

	//FIXME: add to API as parameter; pull from dan she to get this
	Frequency_of_redistribution(&front,GENERAL_WAVE) = 5; 
	RK_PRINT("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
			Frequency_of_redistribution(&front,GENERAL_WAVE));
	Redistribution_count(&front) = front.step;
}

void FTAPI_Propagate(const double *grid_dt, double *elapsed)
{
#ifndef NDEBUG
	{
		/*DBG*/
		POINT *p;
		HYPER_SURF *hs;
		HYPER_SURF_ELEMENT *hse;
		INTERFACE *intfc = front.interf;
		next_point(intfc, NULL, NULL, NULL);
		int d;
		while(next_point(intfc, &p, &hse, &hs))
			for(d=0; d<intfc->dim; ++d)
			assert(!isnan(Coords(p)[0]));
	}
#endif //NDEBUG

    // This is the amount of real time elapsed by subcycling
    // This should be compared to external dt at the end of subcycling
    int subCycCount;

    // Loop over subcycles
    // (ideally no subcycles necessary)
    RK_PRINT("Propagating: dt = %f\n",*grid_dt);
    subCycCount = 0;
    // JAM - 5/9/16 - Subcycling moved to client side
/*
    while(elapsed < *grid_dt)
    {
*/
        FT_SetTimeStep(&front); // Set the timestep based off of FT
        front.dt = min(((*grid_dt) - *elapsed), front.dt);
	// As subcyling was moved to client side, this error catching is no longer needed
        // JAM - 5/7/16  Added Error catching for a front dt of 0
/*
	if (front.dt == 0.0 && front.step > 1 && subCycCount > 15)
	{
	    fprintf(stdout,"ERROR, the front's dt was calculated to 0, this is most\n"
	   		   "likely a result of failures to untangle too many times.\n"
			   "Investigating the front and it's tangles are probably\n"
			   "needed at this point. Front plots will be dumped with\n"
			   "the name ft_plt_dt0ERROR-nd****.vtk\n");
            vtk_interface_plot_byname("ft_plt_dt0ERROR", front.interf, NO, 
				      front.time, front.step, front.coordinate);
	    clean_up(ERROR);
	}

#ifndef NDEBUG
		if(front.dt < *grid_dt)printf("FT subcycling dt=%f\n",front.dt);
#endif //NDEBUG
*/
        *elapsed += front.dt;

        FT_Propagate(&front);
        //subCycCount++;
        FT_AddTimeStepToCounter(&front);
//    }
#ifndef NDEBUG
	{
		/*DBG*/
		POINT *p;
		HYPER_SURF *hs;
		HYPER_SURF_ELEMENT *hse;
		INTERFACE *intfc = front.interf;
		next_point(intfc, NULL, NULL, NULL);
		int d;
		while(next_point(intfc, &p, &hse, &hs))
			for(d=0; d<intfc->dim; ++d)
				assert(!isnan(Coords(p)[0]));
	}
#endif //NDEBUG
}

void FTAPI_Output(const char *filename, int *length)
{
    char file[256];
    int i;
    for(i=0; i<=length; ++i)
    {
	if(filename[i]=='\0')
	{
	    file[i]='\0';
	    break;
	}
	else if(filename[i]==' ')
	{
	    file[i] = '\0';
	    break;
	}
      file[i] = filename[i];
    }
    assert(i<=length);
    vtk_interface_plot_byname(file, front.interf,
	NO, front.time, front.step, front.coordinate);
}
void FTAPI_WriteRestart(const char* filename, int *length)
{
    char file[256];
    int i;
    for(i=0; i<length; ++i)
    {
	if(filename[i]=='\0')
	{
	    file[i]='\0';
	    break;
	}
	else if(filename[i]==' ')
	{
	    file[i] = '\0';
	    break;
	}
	file[i] = filename[i];
    }
    assert(i<length);
    FT_Save(&front, file);
}

// This function is not needed, but is useful for debugging so I left it.
/* static void fr_getCellOfPoint(double *L,double *U,int *gmax,double *coords,int *gcoords)
{

        //int bufsize[3][2]={{0,0},{0,0},{0,0}};
        int bufsize[3][2]={{front.rect_grid->lbuf[0],front.rect_grid->ubuf[0]},
                     {front.rect_grid->lbuf[1],front.rect_grid->ubuf[1]},
                     {front.rect_grid->lbuf[2],front.rect_grid->ubuf[2]}};
        int i;

        gcoords[0] = floor((coords[0]-L[0])/(U[0]-L[0]) *gmax[0]);
        gcoords[1] = floor((coords[1]-L[1])/(U[1]-L[1]) *gmax[1]);
        gcoords[2] = floor((coords[2]-L[0])/(U[2]-L[2]) *gmax[2]);

        for(i=0; i < 3; ++i)
        {
                gcoords[i]+=bufsize[i][0];
                // FIXME: Bug fix hack
                if(gcoords[i] <0)
                {
//                    printf("got coords: %f %f %f : %f %f %f\n", 
//                            coords[0], coords[1],coords[2],
//                            (coords[0]-L[0])/(U[0]-L[0]) *gmax[0],
//                            (coords[1]-L[1])/(U[1]-L[1]) *gmax[1],
//                            (coords[2]-L[2])/(U[2]-L[2]) *gmax[2]);
//                    printf("Gcoords = %.16f %.16f %.16f : %d %d %d\n", 
//                        (coords[0]-L[0])/(U[0]-L[0]) *gmax[0],
//                        (coords[1]-L[1])/(U[1]-L[1]) *gmax[1],
//                        (coords[2]-L[2])/(U[2]-L[2]) *gmax[2],
//                        gcoords[0],gcoords[1],gcoords[2]);
                    assert(gcoords[i] == -1);
                    gcoords[i]=0;// BUG FIX
                    // FIXME: !!!!!
                    // Drifting out to sea.... whoaaaaaaaaa

                }
                assert(gcoords[i]>=0);

                if(gcoords[i] >= gmax[i]+bufsize[i][0]+bufsize[i][1])
                {
                        assert(gcoords[i] == gmax[i]+bufsize[i][0]+bufsize[i][1]);
                        gcoords[i] -= 1;
                }
                if(gcoords[i] < bufsize[i][0])
                    gcoords[i]+=gmax[i];
                else if (gcoords[i] >= gmax[i] + bufsize[i][0])
                    gcoords[i]-=gmax[i];
        }

        assert(gcoords[0] >= 0 && gcoords[0] < gmax[0]+8);
        assert(gcoords[1] >= 0 && gcoords[1] < gmax[1]+8);
        assert(gcoords[2] >= 0 && gcoords[2] < gmax[2]+8);
}*/

/* void debug_point(POINT *p, HYPER_SURF_ELEMENT *hse, HYPER_SURF *hs)
{
        double tol = sqr(1.e-5);
        double coords1[3] = {-0.241666916780643892,   0.00998216360650780628, -0.25};
        double coords2[3] ={-0.241666916780643892 + 0.5,   0.00998216360650780628, -0.25};
        double dist1 = sqr(Coords(p)[0] - coords1[0]) + sqr(Coords(p)[1] - coords1[1]) + sqr(Coords(p)[2] - coords1[2]);
        double dist2 = sqr(Coords(p)[0] - coords2[0]) + sqr(Coords(p)[1] - coords2[1]) + sqr(Coords(p)[2] - coords2[2]);
        
        double vel[3];
        if(dist1 < tol)
        {
            flash_vel_intrp( front.vparams, &front, p, hse, hs, vel);
            printf("P1: coords = %f %f %f : vel = %g %g %g\n",
            Coords(p)[0], Coords(p)[1], Coords(p)[2],
            vel[0], vel[1],vel[2]);
        }
        
        if(dist2 < tol)
        {
            flash_vel_intrp( front.vparams, &front, p, hse, hs, vel);
            printf("P2: coords = %f %f %f : vel = %g %g %g\n",
            Coords(p)[0], Coords(p)[1], Coords(p)[2],
            vel[0], vel[1],vel[2]);
        }
}*/                   

//--------------------------------
// ACTIVE TRACKING
//--------------------------------
int FTAPI_getComponent(double* coords, int* newcomp)
{
    int comp = component(coords, front.interf);
    *newcomp = comp;
    return comp;
}


void FTAPI_normal(const double* coords, double *nor)
{
	INTERFACE *intfc = front.interf;
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	POINT *p;
	//double nor[MAXD];
	double dist;
	double min_dist = HUGE_VAL;
	int dim = intfc->dim;

	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (p->hs->obj.c->_boundary != 0)
	    {
			p->_nor[0] = 0.0;
			p->_nor[1] = 0.0;
			p->_nor[2] = 0.0;
			continue;
	    }
		dist = distance_between_positions(Coords(p), coords, dim);
		if (dist < min_dist)
		{
			normal(p,hse,hs,nor,&front);
			min_dist = dist;
		}
	}

}

// JAM - This function is designed to get the grid front normal to be used when
//       a component changes and the velocity needs to be dotted with the normal
//       front velocity of the crossing front.  We use the method copied from 
//       cFcartsn.cpp:7775 get_normal_from_front() which basically does an avg
//       procedure for all front points in the bounding boxes which are between 
//       the current cell center and any neighboring cell center
void FTAPI_gridNormal(const double* coords, double* dh, double* nor)
{
    INTERFACE *intfc = front.interf;
    HYPER_SURF *hs;
    HYPER_SURF_ELEMENT *hse;
    POINT *p;
    double dist[MAXD];
    double norm[MAXD];
    double f;
    double mag;
    int dim = intfc->dim;
    int pass, i;

    for (i = 0; i < dim; i++)
        nor[i] = 0.0;

    (void) next_point(intfc,NULL,NULL,NULL);
    while (next_point(intfc,&p,&hse,&hs))
    {
        if (p->hs->obj.c->_boundary != 0)
            continue;
        
        pass = 1;
        for (i = 0; i < dim; i++)
        {
            dist[i] = fabs(coords[i] - Coords(p)[i])/dh[i];
            if (dist[i] > 1.+1.e-6*dh[i] || dist[i] < -1.e-6*dh[i])
                pass = 0;
        }
        
        if (pass == 0)
            continue;

        f = 1.0 - dist[0];
        for (i = 1; i < dim; i++)
            f *= (1.0 - dist[i]);

        normal(p,hse,hs,norm,&front);

        for (i = 0; i < dim; i++)
            nor[i] += norm[i]*f;
    }

    mag = 0.0;
    for (i = 0; i < dim; i++)
        mag += sqr(nor[i]);

    // JAM : 5/7/16: if mag is still 0, we didn't find any front points in the bounding box, 
    // I believe this should only happen if a small isolated curve was deleted, leaving cells
    // that were previously of one component needing to be changed.  Ideally they would take 
    // an average of their neighbors, but the Riemann solution seems to be a minor error.  We
    // will use a normal vector of even magnitude in all directions   FIXME: There may be a better
    // way to approach this
    if (mag == 0.0)
    {
        // This cannot happen in 1D, plus this API isn't really setup for 1D anyway
        if (dim == 2)
        {
            mag = 1;
            nor[0] = sqrt(2)/2.;
            nor[1] = sqrt(2)/2.;
        }
        else if (dim == 3)
        {
            mag = 1;
            nor[0] = sqrt(3)/3.;
            nor[1] = sqrt(3)/3.;
            nor[2] = sqrt(3)/3.;
        }
        else
        {
            fprintf(stdout,"This is not setup for a non 2D or 3D problem, fapi.c:664\n");
            clean_up(ERROR);
        }
    }
   
    mag = sqrt(mag);

    for (i = 0; i < dim; i++)
        nor[i] = nor[i]/mag;

    // JAM : I experienced some problems with floating point errors in the normals
    //       and with garbage in unused dimensions, so doing some housekeeping here
    for (i = 0; i < dim; i++)
        if (fabs(nor[i]) < 1.e-15)
            nor[i] = 0.0;

    for (i = dim; i < 3; i++)
        nor[i] = 0.0;
}

/*
void FTAPI_normal(int *r_size, double** r_coords, double **r_normal)
{
	INTERFACE *intfc = front.interf;
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	POINT *p;
	//double nor[MAXD];
	double dist;
	double min_dist = HUGE_VAL;
	int dim = intfc->dim;
	int num_pts = 0;

	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (Boundary_point(p))
	    {
			p->_nor[0] = 0.0;
			p->_nor[1] = 0.0;
			p->_nor[2] = 0.0;
			continue;
	    }
		++num_pts;
	}
	// allocate and fill pt array ... ???
	double *coords = malloc(num_pts*dim*sizeof(double));
	double *normals = malloc(num_pts*dim*sizeof(double));
	int i_c = 0;  // index in coords array
	int i_n = 0;  // index in normals array

	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (Boundary_point(p))
	    {
			p->_nor[0] = 0.0;
			p->_nor[1] = 0.0;
			p->_nor[2] = 0.0;
			continue;
	    }
		int d;
		for(d=0; d<dim; ++d)
			coords[i_c++] = Coords(p)[d];
		normal(p, hse, hs, &normals[i_n], &front);
		i_n += dim;
	}
	*r_coords  =  coords;
	*r_normal  =  normals;
	*r_size    =  num_pts;
}
*/

void FTAPI_cross(const double** line, double *cross_return)
{
	assert(front.interf->dim==2);
	const int MAX_CROSS=5;
	double cross[MAX_CROSS][2];//2D only
	//double cross_temp[MAX_CROSS][2];
	//double cross[2][MAX_CROSS];//2D only, hao modified
	//double cross_temp[2][MAX_CROSS];//2D only, hao modified, make sure cross keep the original array
	int    j;
	int    index;
	int    crossid = 0;
	double ab; // distance between point line[0] and line[1]
	double small_edge = HUGE_VAL; // sorting the smallest value
	double pa[MAX_CROSS];


	CURVE **c = front.interf->curves;
	BOND *b;
	for(c=front.interf->curves; c&&*c; ++c)
	{
		for(b=(*c)->first; b!=(*c)->last; b=b->next)
		{
			double *edge[] = {Coords(b->start), Coords(b->end)};
			if(pts_cross(line[0], line[1], edge[0], edge[1]))
			{
				double x,y;
				intersection2d( line[0][0], line[0][1], line[1][0], line[1][1],
								edge[0][0], edge[0][1], edge[1][0], edge[1][1],
								&x, &y);
				cross[crossid][0] = x;
				cross[crossid][1] = y;
				//cross[0][crossid] = x;
				//cross[1][crossid] = y;
				++crossid;
			}
		}
	}
	assert(crossid>0);
	if(crossid==1)
	{
		//cross_return[0] = cross[0][0];
		//cross_return[1] = cross[1][0];
		cross_return[0] = cross[0][0];
		cross_return[1] = cross[0][1];
		return;
	}
	/*
	else
	{
		cross_return[0] = cross[0][0];
		cross_return[1] = cross[0][1];
	}
	*/
	else
	{
		// this is the sorting process.
		// TODO: line[0] and line[1] which is input first? leftmost OR rightmost
		ab = distance_between_positions(line[0], line[1], 2);
		for (j = 0; j < crossid; j++)
		{
			pa[j] = distance_between_positions(line[0], cross[j], 2);
			if (pa[j] < small_edge)
			{
				small_edge = pa[j];
				index = j;
			}
		}
		cross_return[0] = cross[index][0];
		cross_return[1] = cross[index][1];
		return;
	}
	//assert(0&&"Implement multiple crosses here...");
}

// ======================================================
// helper functions
// TODO: move to internal
// ======================================================
static int left(const double *p, const double *q1, const double *q2)
{
	double v1[2];
	double v2[2];
	int i;
	for(i=0; i<2; ++i)
	{
		v1[i] = p[i]-q1[i];
		v2[i] = p[i]-q2[i];
	}
	/*
	printf("-----------------------------haozhang formula v1 = p -q1 and v2 = p - q2---------------------------------------------\n");
	printf("haozhang in function %s v1 = {%f %f}\n", __func__, v1[0], v1[1]);
	printf("haozhang in function %s v2 = {%f %f}\n", __func__, v2[0], v2[1]);
	printf("haozhang in function %s q1 = {%f %f}\n", __func__, q1[0], q1[1]);
	printf("haozhang in function %s q2 = {%f %f}\n", __func__, q2[0], q2[1]);
	printf("haozhang in function %s p  = {%f %f}\n", __func__, p[0], p[1]);
	*/
	double crx = v1[0]*v2[1] - v1[1]*v2[0];
	//printf("haozhang in function %s crx  = v1[0]*v2[1] - v1[1]*v2[0] = %f\n", __func__, crx);
	return crx > 0;
}

static int pts_cross(const double *p1, const double *p2, const double *q1, const double *q2)
{
	return (left(p1, q1, q2) != left(p2, q1, q2)) && 
		(left(q1, p1, p2) != left(q2, p1, p2));
}

static int intersection2d(  const double ax1, const double ay1,
			    const double ax2, const double ay2, 
			    const double bx1, const double by1,
			    const double bx2, const double by2, 
			    double *x, double *y)
{
	double dya = ay2-ay1;
	double dxa = ax2-ax1;
	double dyb = by2-by1, dxb = bx2-bx1;
	double ma = dya/dxa;
	double mb = dyb/dxb;
	const double EPSILON = 1.e-15; // FIXME: tolerance

	// We can solve the intersection simply with algebra, as long as no line is vertical
	// unfortunately, this happens quite often, since lines are often grid lines.
	// So, we have catches for this case on each solve

	// solve for x:
	if (dxa == 0.0)
	{
		assert(fabs(dxb)>0.0);
		assert(!isfinite(ma));

		*x = ax1;

		assert(*x==ax2);

		*y = mb* (*x - bx1) + by1;

		assert(fabs(mb * (*x - bx2) + by2 - *y) < EPSILON);
		return 0;
	}
	else if (dxb == 0.0)
	{
		assert(fabs(dxa)>0.0);
		assert(!isfinite(mb));

		*x = bx1;

		assert(*x==bx2);

		*y = ma* (*x - ax1) + ay1;

		assert(fabs(ma * (*x - ax2) + ay2 - *y) < EPSILON);
		return 0;
	}
	else
	{
		assert(isfinite(ma));
		assert(isfinite(mb));

		*x = (ma*ax1  - mb*bx1 - ay1 + by1)/(ma-mb);
	}

	// Solve for y:
	if(dya == 0.0)
	{
		assert(fabs(dyb)>0.0);

		*y = ay1;

		assert(*y==ay2);
	}
	else if(dyb == 0.0)
	{
		assert(fabs(dya)>0.0);

		*y = by1;

		assert(*y==by2);
	}
	else
	{
		assert(isfinite(ma));
		assert(isfinite(mb));
		assert(!isnan(ma));
		assert(!isnan(mb));
		assert(!isinf(ma));
		assert(!isinf(mb));

		*y =  ma * ((*x) - ax1) + ay1;

		assert(*y == mb * ((*x) - bx1) + by1);
	}
	return 0;
}

void FTI_get_bonds(struct FTI_BOND **i_bonds, int *num_bonds)
{
	INTERFACE *intfc = front.interf;
	CURVE **c;
	int i;

	// count the bonds
	int nb = 0;
	for (c = intfc->curves; c && *c; ++c)
	{
	    if(is_bdry_hs((*c)->hs))
		continue;
            nb += NumOfCurveBonds(*c);
	}

	FT_VectorMemoryAlloc((POINTER*)i_bonds,nb,sizeof(struct FTI_BOND));
	struct FTI_BOND *bonds = *i_bonds;

	int n=0;
	for (c = intfc->curves; c && *c; ++c)
	{
		if(is_bdry_hs((*c)->hs))
		    continue;
		BOND *b;
		for(b = (*c)->first; b != NULL; b = b->next)
		{
			for(i=0; i<3; ++i)
			{
				bonds[n].start[i] = Coords(b->start)[i];
				bonds[n].end[i] = Coords(b->end)[i];
			}
			++n;
		}
	}
	assert(n==nb);
	*num_bonds = nb;
}

void FTI_free(void *ptr)
{
	f_ree(ptr, "hello?");
}
