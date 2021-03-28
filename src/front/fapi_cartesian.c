
#include <stdlib.h>
#include <assert.h>
#include "fdecs.h"
#include "FTAPI.h"
#include "ftapi.h"
#include <stdio.h>

extern Front front;

struct _FT_CARTESIAN_VEL_PARAMS
{
    double *xvel;
    double *yvel;
    double *zvel;
};
typedef struct _FT_CARTESIAN_VEL_PARAMS FT_CARTESIAN_VEL_PARAMS;


static double getDouble(POINTER value)
{
    assert(0);
    return *((double*)(value));
}

static int FT_cartesian_vel_intrp(POINTER params, double *coords, double *vel);
static int FT_polar_vel_intrp(POINTER params, double *coords, double *vel);

void FNAME(ftapi_init_cartesian,FTAPI_INIT_CARTESIAN)(int *dim,
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
        // Velocity arrays
        double *vx, double *vy, double *vz,
	double (*level_func)( POINTER func_params, double *coords),
        void* level_func_params,
        // geomerty type
        int *geometry,
        int *restart,
	char *restart_filename,
	int restart_filename_length)
{
    	static FT_CARTESIAN_VEL_PARAMS vel_params;
	vel_params.xvel = vx;
	vel_params.yvel = vy;
	vel_params.zvel = vz;

	if(*geometry == FT_GEOMETRY_CARTESIAN)
	{
	    FNAME(ftapi_init,FTAPI_INIT)(dim, procGrid, 
		xmin, xmax, ymin, ymax, zmin, zmax,
		isize, jsize, ksize,
		ibuf, jbuf, kbuf,
		xlbdry, xubdry, ylbdry, yubdry, zlbdry, zubdry,
		level_func, level_func_params,
		FT_cartesian_vel_intrp, &vel_params,
		geometry, restart, restart_filename, restart_filename_length);
	}
	else if(*geometry == FT_GEOMETRY_POLAR)
	{
	    assert(*dim==2);
	    FNAME(ftapi_init,FTAPI_INIT)(dim, procGrid, 
		xmin, xmax, ymin, ymax, zmin, zmax,
		isize, jsize, ksize,
		ibuf, jbuf, kbuf,
		xlbdry, xubdry, ylbdry, yubdry, zlbdry, zubdry,
		level_func, level_func_params,
		FT_polar_vel_intrp, &vel_params,
		geometry, restart, restart_filename, restart_filename_length);
	}
}

void FTAPI_CreateSurface_cartesian(int pos_comp, int neg_comp, double (*level_func)(void*, double*), void* level_params, int surf_type)
{
	static FT_CARTESIAN_VEL_PARAMS vel_params;
	FTAPI_CreateSurface(pos_comp, neg_comp, level_func, level_params, surf_type,
		FT_cartesian_vel_intrp, &vel_params);
}

// JAM - Removed surf_type form argument list as it doesn't seem to be used for anything...?
void FNAME(ftapi_propagate_cartesian, FTAPI_PROPAGATE_CARTESIAN)(double *dt, double *xvel, double *yvel, double *zvel, double *fr_dt)
{
// JAM - Removed for FTAPI -- Don't understand what the point of this is?
//	extern void* VELOCITY_PARAMS[256];

//    FT_CARTESIAN_VEL_PARAMS *vparams=
//    	(FT_CARTESIAN_VEL_PARAMS*)(VELOCITY_PARAMS[surf_type]);
//    vparams->xvel = xvel;
//    vparams->yvel = yvel;
//    vparams->zvel = zvel;
    FTAPI_Propagate(dt, fr_dt);
}

// Interpolate velocity at point based on flash grid (polar)
static int FT_polar_vel_intrp(
        POINTER params,
        double *coords,
        double *vel)
{
    FT_CARTESIAN_VEL_PARAMS* veldata = (FT_CARTESIAN_VEL_PARAMS*)params;
    if(veldata->xvel != NULL) 
    {
        double coords[2];
        double x,y,r,theta;
        
        x = coords[0];
        y = coords[1];
        r = sqrt(x*x+y*y);
        theta = atan2(y,x);
        coords[0] = r;
        coords[1] = theta;

        double vr = vel[0];
        double vtheta = vel[1];
            
        double st = sin(theta);
        double ct = cos(theta);
            
        double vx = r*vtheta*st + vr*ct;
        double vy = r*vtheta*ct + vr*st;
        vel[0] = vx;
        vel[1] = vy;
     }

    else
    {
        vel[0] = vel[1] = vel[2]= 0.0f;
#ifndef NDEBUG 
        if(front.step != 0)
        {
            printf("ERROR: Velocity not set after step 0!\n");
            clean_up(-1);
        }   
#endif  
    }

    return 0;

}

/*
void FT_cartesian_vel_set( double *xvel, double *yvel, double *zvel)
{
    FT_CARTESIAN_VEL_PARAMS *vparams=(FT_CARTESIAN_VEL_PARAMS*)(front.vparams);
    vparams->xvel = xvel;
    vparams->yvel = yvel;
    vparams->zvel = zvel;
}
*/


// Interpolate velocity at point based on flash grid
static int FT_cartesian_vel_intrp(
        POINTER params,
        double *coords,
        double *vel)
{
    FT_CARTESIAN_VEL_PARAMS* veldata = (FT_CARTESIAN_VEL_PARAMS*)params;
    if(veldata->xvel != NULL)
    {
	// TODO: set component here for active tracking
	//	 i.e., use component, not NO_COMP
		
	if(!(FT_IntrpStateVarAtCoords(&front, NO_COMP, coords, veldata->xvel, 
                    getDouble, &(vel[0])) &&
            FT_IntrpStateVarAtCoords(&front, NO_COMP, coords, veldata->yvel, 
                    getDouble, &(vel[1])) &&
            FT_IntrpStateVarAtCoords(&front, NO_COMP, coords, veldata->zvel, 
                    getDouble, &(vel[2]))))
        {
	    // interface boundary points outside the comp domain
		// should not be propagated. Thus, we set the velocity
		// to zero.
	
		memset(vel, 0.0f, 3*sizeof(double));
		return 0;
        }
    }
    else
    {
        vel[0] = vel[1] = vel[2]= 0.0f;
#ifndef NDEBUG
        if(front.step != 0)
        {
            printf("ERROR: Velocity not set after step 0!\n");
            clean_up(-1);
        }
#endif
    }
    return 0;

}
