#ifndef _FTAPI_H
#define _FTAPI_H
#define FNAME(x,X) x##_

#define FTAPI_Init FNAME(ftapi_init, FTAPI_INIT)
#define FTAPI_Propagate FNAME(ftapi_propagate, FTAPI_PROPAGATE)
#define FTAPI_Output FNAME(ftapi_output, FTAPI_OUTPUT)
#define FTAPI_WriteRestart FNAME(ftapi_writerestart, FTAPI_WRITERESTART)
#define FTAPI_getComponent FNAME(ftapi_getcomponent, FTAPI_GETCOMPONENT)
#define FTAPI_normal FNAME(ftapi_normal, FTAPI_NORMAL)
#define FTAPI_gridNormal FNAME(ftapi_gridnormal, FTAPI_GRIDNORMAL)

#ifdef __cplusplus
extern "C"
{
#endif //__cplusplus

struct FTI_BOND
{
    double start[3];
    double end[3];
};
void FTI_get_bonds(struct FTI_BOND **bonds, int *num_bonds);
void FTI_free(void *ptr);


//TODO: some const correctness on pointers.
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
        double (*level_func)( void* func_params, double *coords),
        void* level_func_params,
        int (*vel_func)(void*, double*, double*, int*, int*),
        void* vel_params,
        // geomerty type
        int *geometry,
        int *restart,
	char *restart_filename,
	int *restart_filename_length);

void FTAPI_Output(const char *filename, int *length);
void FTAPI_WriteRestart(const char *filename, int *length);
void FTAPI_Propagate(const double *dt, double *elapsed);

void FTAPI_CreateSurface(
	int pos_comp, 
	int neg_comp, 
	double (*level_func)(void*, double*), 
	void* level_params, 
	int surf_type, 
	int (*vel_func)(void*, double*, double*, int*, int*),
	void* vel_params);


enum FT_BOUNDARY_TYPE
{
    FT_BOUNDARY_PERIODIC,
    FT_BOUNDARY_REFLECTING,
    FT_BOUNDARY_FLOW
};

enum FT_GEOMETRY_TYPE
{
    FT_GEOMETRY_SPHERICAL,
    FT_GEOMETRY_CARTESIAN,
    FT_GEOMETRY_CYLINDRICAL,
    FT_GEOMETRY_POLAR
};
///---------------------------------------------------
///--ACTIVE TRACKING----------------------------------
///---------------------------------------------------
int FTAPI_getComponent(double* coords, int* newcomp);

//----------------------------------------------------
//-Cartesian Interface--------------------------------
//----------------------------------------------------
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
	double (*level_func)(void* func_params, double *coords),
        void* level_func_params,
        // geomerty type
        int *geometry,
        int *restart,
	char *restart_filename,
	int restart_filename_length);

void FTAPI_CreateSurface_cartesian(int pos_comp, int neg_comp,
	    double (*level_func)(void*, double*), void* level_params,
	    int surf_type);

// Removed surf_type (2nd) from argument list
void FNAME(ftapi_propagate_cartesian, FTAPI_PROPAGATE_CARTESIAN)(double *dt, double *xvel, double *yvel, double *zvel, double *elapsed);

//--------------------- ACTIVE TRACKING ----------------

//void FTAPI_normal(int *size, double** coords, double **nor);
void FTAPI_normal(const double* coords, double *nor);
void FTAPI_cross(const double** line, double *cross);
void FTAPI_gridNormal(const double* coords, double* dh, double* nor);

#ifdef __cplusplus
}
#endif //__cplusplus

#endif //_FTAPI_H
