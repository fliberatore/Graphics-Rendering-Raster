/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include	"float.h"

#define M_PI 3.14159265358979323846f 

#ifndef DDA
#define DDA
typedef struct{
	GzCoord start, end, current;
	float slopeX, slopeZ;
} GzDDA;
#endif

#ifndef SPAN
#define SPAN
typedef struct{
	GzCoord start, end, current;
	float slopeZ;
} GzSpan;
#endif

#define NUM_DIMENSIONS 3
#define NUM_VERTICES 3
#define NUM_EDGES 3

//Functions prototypes external to API
//Rasterizer's functions
int GzFrustumCull(GzRender *render, GzCoord vertexList[NUM_VERTICES], bool &culled);
int GzScanLineRasterizer(GzRender *render, GzCoord vertexList[NUM_VERTICES]);
int GzPixelToDisplay(GzSpan *span, GzDisplay *display, GzColor *flatcolor);

//Transforms' functions
int GzSetupXsp(GzRender *render);
int GzSetupXpi(GzRender *render);
int GzSetupXiw(GzRender *render);

//Algebra functions
float GzVectorNorm(GzCoord v);
float GzDotProduct(GzCoord a, GzCoord b);

//Usefull functions
short ctoi(float color);
int GzInitDefaultCamera(GzRender *render);
int GzXformTriangle(GzRender *render, GzCoord *vertexList);

int GzRotXMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value
	float theta = degree * M_PI / 180.0f;
	mat[0][0] = 1.0f; mat[0][1] = 0.0f; mat[0][2] = 0.0f; mat[0][3] = 0.0f;
	mat[1][0] = 0.0f; mat[1][1] = cos(theta); mat[1][2] = - sin(theta); mat[1][3] = 0.0f;
	mat[2][0] = 0.0f; mat[2][1] = sin(theta); mat[2][2] = cos(theta); mat[2][3] = 0.0f;
	mat[3][0] = 0.0f; mat[3][1] = 0.0f; mat[3][2] = 0.0f; mat[3][3] = 1.0f;

	return GZ_SUCCESS;
}


int GzRotYMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value
	float theta = degree * M_PI / 180.0f;
	mat[0][0] = cos(theta); mat[0][1] = 0.0f; mat[0][2] = sin(theta); mat[0][3] = 0.0f;
	mat[1][0] = 0.0f; mat[1][1] = 1.0f; mat[1][2] = 0.0f; mat[1][3] = 0.0f;
	mat[2][0] = -sin(theta); mat[2][1] = 0.0f; mat[2][2] = cos(theta); mat[2][3] = 0.0f;
	mat[3][0] = 0.0f; mat[3][1] = 0.0f; mat[3][2] = 0.0f; mat[3][3] = 1.0f;

	return GZ_SUCCESS;
}


int GzRotZMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value
	float theta = degree * M_PI / 180.0f;
	mat[0][0] = cos(theta); mat[0][1] = -sin(theta); mat[0][2] = 0.0f; mat[0][3] = 0.0f;
	mat[1][0] = sin(theta); mat[1][1] = cos(theta); mat[1][2] = 0.0f; mat[1][3] = 0.0f;
	mat[2][0] = 0.0f; mat[2][1] = 0.0f; mat[2][2] = 1.0f; mat[2][3] = 0.0f;
	mat[3][0] = 0.0f; mat[3][1] = 0.0f; mat[3][2] = 0.0f; mat[3][3] = 1.0f;

	return GZ_SUCCESS;
}


int GzTrxMat(GzCoord translate, GzMatrix mat)
{
// Create translation matrix
// Pass back the matrix using mat value
	mat[0][0] = 1.0f; mat[0][1] = 0.0f; mat[0][2] = 0.0f; mat[0][3] = translate[0];
	mat[1][0] = 0.0f; mat[1][1] = 1.0f; mat[1][2] = 0.0f; mat[1][3] = translate[1];
	mat[2][0] = 0.0f; mat[2][1] = 0.0f; mat[2][2] = 1.0f; mat[2][3] = translate[2];
	mat[3][0] = 0.0f; mat[3][1] = 0.0f; mat[3][2] = 0.0f; mat[3][3] = 1.0f;
	
	return GZ_SUCCESS;
}


int GzScaleMat(GzCoord scale, GzMatrix mat)
{
// Create scaling matrix
// Pass back the matrix using mat value
	mat[0][0] = scale[0]; mat[0][1] = 0.0f; mat[0][2] = 0.0f; mat[0][3] = 0.0f;
	mat[1][0] = 0.0f; mat[1][1] = scale[1]; mat[1][2] = 0.0f; mat[1][3] = 0.0f;
	mat[2][0] = 0.0f; mat[2][1] = 0.0f; mat[2][2] = scale[2]; mat[2][3] = 0.0f;
	mat[3][0] = 0.0f; mat[3][1] = 0.0f; mat[3][2] = 0.0f; mat[3][3] = 1.0f;
	
	return GZ_SUCCESS;
}


//----------------------------------------------------------
// Begin main functions

int GzNewRender(GzRender **render, GzDisplay	*display)
{
	if(!render || !display)
		return GZ_FAILURE;
	
	//malloc a renderer struct
	GzRender* pRender =(GzRender*) malloc(sizeof(GzRender));
	if(!pRender)
		return GZ_FAILURE;

	//span interpolator needs pointer to display for pixel writes
	pRender->display = display;

	//setup Xsp
	if(GzSetupXsp(pRender))
		return GZ_FAILURE;

	//init default camera
	if(GzInitDefaultCamera(pRender))
		return GZ_FAILURE;

	*render = pRender;
	return GZ_SUCCESS;
}


int GzFreeRender(GzRender *render)
{
	//free all renderer resources
	if(render) free(render);
	return GZ_SUCCESS;
}


int GzBeginRender(GzRender *render)
{

	//setup for start of each frame - init frame buffer color,alpha,z
	if(GzInitDisplay(render->display))
		return GZ_FAILURE;  /* init for new frame */

	//compute Xiw and projection xform Xpi from camera definition
	if(GzSetupXpi(render) || GzSetupXiw(render))
		return GZ_FAILURE;

	//init Ximage - put Xsp at base of stack, push on Xpi and Xiw 
	render->matlevel = -1;
	if(GzPushMatrix(render, render->Xsp) || GzPushMatrix(render, render->camera.Xpi) || GzPushMatrix(render, render->camera.Xiw))
		return GZ_FAILURE;

	//now stack contains Xsw and app can push model Xforms when needed 
	return GZ_SUCCESS;
}

int GzPutCamera(GzRender *render, GzCamera *camera)
{
	if(!render || !camera)
		return GZ_FAILURE;
/*
- overwrite renderer camera structure with new camera definition
*/
	memcpy(render->camera.position, camera->position, sizeof(GzCoord));
	memcpy(render->camera.lookat, camera->lookat, sizeof(GzCoord));
	memcpy(render->camera.worldup, camera->worldup, sizeof(GzCoord));
	render->camera.FOV = camera->FOV;

	return GZ_SUCCESS;	
}

//Push a matrix into the Ximage stack
int GzPushMatrix(GzRender *render, GzMatrix	matrix)
{

	//Check for errors and stack overflow
	if(!render || render->matlevel >= MATLEVELS)
		return GZ_FAILURE;

	if(render->matlevel == -1)
	{
		render->matlevel = 0;
		memcpy(render->Ximage[render->matlevel], matrix, sizeof(GzMatrix));
	} else
	{
		//Matrix multiplication
		GzMatrix tmpMat;
		for(unsigned int i = 0; i < 4; ++i){
			for(unsigned int j = 0; j < 4; ++j){
				tmpMat[i][j] = 0.0f;
				for(unsigned int pos = 0; pos < 4; ++pos){
					tmpMat[i][j] += render->Ximage[render->matlevel][i][pos] * matrix[pos][j];
				}
			}
		}

		//Increase Top of Stack and copy tmpMat into top position
		render->matlevel++;
		memcpy(render->Ximage[render->matlevel], tmpMat, sizeof(GzMatrix));
	}
	return GZ_SUCCESS;
}

//pop a matrix off the Ximage stack
int GzPopMatrix(GzRender *render)
{
	//Check for errors and stack underflow
	if(!render || render->matlevel < 0)
		return GZ_FAILURE;
	
	//Decrease Top of Stack
	render->matlevel--;
	return GZ_SUCCESS;
}


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer	*valueList) /* void** valuelist */
{
	if(!render || numAttributes <= 0 || !nameList || !valueList)
		return GZ_FAILURE;

	int status = 0;
	for(unsigned int i = 0; i < numAttributes; ++i){
		switch(nameList[i]){
		case GZ_NULL_TOKEN: break;
		case GZ_RGB_COLOR:
			memcpy(render->flatcolor, valueList[i], sizeof(GzColor));
			break;
		default: status = GZ_FAILURE;
		}
		if(status) return GZ_FAILURE;
	}

	return GZ_SUCCESS;
}

int GzPutTriangle(GzRender	*render, int numParts, GzToken *nameList, GzPointer	*valueList)
{
/* numParts : how many names and values */

/* 
pass in a triangle description with tokens and values corresponding to
      GZ_NULL_TOKEN:		do nothing - no values
      GZ_POSITION:		3 vert positions 
*/
	if(!render || numParts <= 0 || !nameList || !valueList)
		return GZ_FAILURE;

	int status = GZ_SUCCESS;
	//Loop on every part in the lists
	for(unsigned int p = 0; p < numParts; ++p){
		switch(nameList[p]){
		case GZ_NULL_TOKEN: break;
		case GZ_POSITION:
			GzCoord vertexList[NUM_VERTICES];
			memcpy(vertexList, valueList[p], NUM_VERTICES*sizeof(GzCoord));
			for(unsigned int v = 0; v < NUM_VERTICES; ++v){
				//Xform positions of verts using matrix on top of stack 
				if(GzXformTriangle(render, &(vertexList[v])))
					return GZ_FAILURE;
				// Clip - just discard any triangle with any vert(s) behind view plane 
				if(vertexList[v][Z] < 0.0f)
					continue;
			}
			//Test for triangles with all three verts off-screen (trivial frustum cull)
			bool culled;
			if(GzFrustumCull(render, vertexList, culled))
				return GZ_FAILURE;
			if(culled)
				continue;
			//Invoke the scan converter and return an error code
			status = GzScanLineRasterizer(render, vertexList);
			break;
		default: status = GZ_FAILURE;
		}
		if(status) return status;
	}

	return GZ_SUCCESS;
}

/* NOT part of API - just for general assistance */
short	ctoi(float color)		/* convert float color to GzIntensity short */
{
  return(short)((int)(color * ((1 << 12) - 1)));
}

int GzXformTriangle(GzRender *render, GzCoord *vertex)
{
	if(!render || !vertex)
		return GZ_FAILURE;

	//Getting Xsm
	GzMatrix Xsm;
	memcpy(Xsm, render->Ximage[render->matlevel], sizeof(GzMatrix));

	//Building 4D vector
	float modVert[4];
	memcpy(modVert, *vertex, sizeof(GzCoord));
	modVert[3] = 1.0f;

	//Applying transform
	float scrVert[4];
	for(unsigned int i = 0; i < 4; ++i){
		scrVert[i] = 0.0f;
		for(unsigned int j = 0; j < 4; ++j){
			scrVert[i] += Xsm[i][j] * modVert[j];
		}
	}

	//Converting 4D vector to 3D
	for(unsigned int i = 0; i < 3; ++i){
		(*vertex)[i] = scrVert[i] / scrVert[3];
	}

	return GZ_SUCCESS;
}

/////////
// RASTERIZER'S FUNCTIONS
/////////
/* Scan Line Rasterizer algorithm */
int GzScanLineRasterizer(GzRender *render, GzCoord vertexList[3]){
	if(!render)
		return GZ_FAILURE;

	// 1 - Sorting vertices
	GzCoord tmpVertex;
	for(unsigned int i = 0; i < NUM_VERTICES-1; ++i){
		for(unsigned int j = i+1; j < NUM_VERTICES; ++j){
			if(vertexList[i][Y] > vertexList[j][Y]){
				memcpy(tmpVertex, vertexList[i], sizeof(GzCoord));
				memcpy(vertexList[i], vertexList[j], sizeof(GzCoord));
				memcpy(vertexList[j], tmpVertex, sizeof(GzCoord));
			}
		}
	}
	float Ymax = vertexList[NUM_VERTICES-1][Y];
	float Ymin = vertexList[0][Y];

	// 2 - Setup edge DDAs
	GzDDA DDAList[NUM_EDGES];
	unsigned int num_h = 0, num_v = 0;
	GzDDA* DDA_H = NULL; //Horizontal DDA
	bool DDA_H_isTop; //Horizontal DDA is a top edge
	GzDDA* DDA_NC[NUM_EDGES]; //Unclassified DDA
	unsigned int DDA_NC_size = 0;
	{
		//Edge 0
		memcpy(DDAList[0].start, vertexList[0], sizeof(GzCoord));
		memcpy(DDAList[0].current, vertexList[0], sizeof(GzCoord));
		memcpy(DDAList[0].end, vertexList[1], sizeof(GzCoord));
		if(DDAList[0].end[Y] != DDAList[0].start[Y]){
			DDAList[0].slopeX = (DDAList[0].end[X] - DDAList[0].start[X])/(DDAList[0].end[Y] - DDAList[0].start[Y]);
			DDAList[0].slopeZ = (DDAList[0].end[Z] - DDAList[0].start[Z])/(DDAList[0].end[Y] - DDAList[0].start[Y]);
			DDA_NC[DDA_NC_size++] = &(DDAList[0]);
		}else{
			DDAList[0].slopeX = DDAList[0].slopeZ = 0;
			DDA_H = &DDAList[0];
			++num_h;
			DDA_H_isTop = true;
		}
		if(DDAList[0].end[X] == DDAList[0].start[X]) ++num_v;

		//Edge 1
		memcpy(DDAList[1].start, vertexList[0], sizeof(GzCoord));
		memcpy(DDAList[1].current, vertexList[0], sizeof(GzCoord));
		memcpy(DDAList[1].end, vertexList[2], sizeof(GzCoord));
		if(DDAList[1].end[Y] != DDAList[1].start[Y]){
			DDAList[1].slopeX = (DDAList[1].end[X] - DDAList[1].start[X])/(DDAList[1].end[Y] - DDAList[1].start[Y]);
			DDAList[1].slopeZ = (DDAList[1].end[Z] - DDAList[1].start[Z])/(DDAList[1].end[Y] - DDAList[1].start[Y]);
			DDA_NC[DDA_NC_size++] = &(DDAList[1]);
		}else{
			DDAList[1].slopeX = DDAList[1].slopeZ = 0;
			DDA_H = &DDAList[1];
			++num_h;
		}
		if(DDAList[1].end[X] == DDAList[1].start[X]) ++num_v;

		//Edge 2
		memcpy(DDAList[2].start, vertexList[1], sizeof(GzCoord));
		memcpy(DDAList[2].current, vertexList[1], sizeof(GzCoord));
		memcpy(DDAList[2].end, vertexList[2], sizeof(GzCoord));
		if(DDAList[2].end[Y] != DDAList[2].start[Y]){
			DDAList[2].slopeX = (DDAList[2].end[X] - DDAList[2].start[X])/(DDAList[2].end[Y] - DDAList[2].start[Y]);
			DDAList[2].slopeZ = (DDAList[2].end[Z] - DDAList[2].start[Z])/(DDAList[2].end[Y] - DDAList[2].start[Y]);
			DDA_NC[DDA_NC_size++] = &(DDAList[2]);
		}else{
			DDAList[2].slopeX = DDAList[2].slopeZ = 0;
			DDA_H = &DDAList[2];
			++num_h;
			DDA_H_isTop = false;
		}
		if(DDAList[2].end[X] == DDAList[2].start[X]) ++num_v;
	}
	//Check for ill-defined triangles
	if(num_v > 1 || num_h > 1)
		return GZ_SUCCESS;

	
	// 3 - Sort edges by L or R
	GzDDA* DDA_L[NUM_EDGES-1]; //Left DDA
	GzDDA* DDA_R[NUM_EDGES-1]; //Right DDA
	if(num_h == 0){
		if(DDA_NC[0]->slopeX < DDA_NC[1]->slopeX){
			DDA_L[0] = DDA_NC[0];
			DDA_R[0] = DDA_NC[1];
		}else{
			DDA_L[0] = DDA_NC[1];
			DDA_R[0] = DDA_NC[0];
		}

		if(DDA_L[0]->end[Y] < DDA_R[0]->end[Y]){
			DDA_L[1] = DDA_NC[2];
			DDA_R[1] = NULL;
		}else{
			DDA_R[1] = DDA_NC[2];
			DDA_L[1] = NULL;
		}
	}else{
		//Horizontal edge
		DDA_L[1] = NULL;
		DDA_R[1] = NULL;
		if(DDA_H_isTop){
			if(DDA_NC[0]->slopeX > DDA_NC[1]->slopeX){
				DDA_L[0] = DDA_NC[0];
				DDA_R[0] = DDA_NC[1];
			}else{
				DDA_L[0] = DDA_NC[1];
				DDA_R[0] = DDA_NC[0];
			}
		}else{
			if(DDA_NC[0]->slopeX < DDA_NC[1]->slopeX){
				DDA_L[0] = DDA_NC[0];
				DDA_R[0] = DDA_NC[1];
			}else{
				DDA_L[0] = DDA_NC[1];
				DDA_R[0] = DDA_NC[0];
			}
		}
	}

	//Rasterizing loop
	unsigned int pos_L = 0, pos_R = 0;
	float Ypos = ceil(Ymin);
	while(Ypos < Ymax){
		//Avoiding useless iterations
		if(Ypos < 0.0f)
			Ypos = 0.0f;
		if(Ypos > render->display->yres)
			break;

		//Jumping to the next edge when required
		if(Ypos > DDA_L[pos_L]->end[Y]){
			++pos_L;
			continue;
		}
		if(Ypos > DDA_R[pos_R]->end[Y]){
			++pos_R;
			continue;
		}

		// 5 - Advance DDAs current positions to y-scan line
		float deltaY;
		if(DDA_L[pos_L]->current[Y] != Ypos){
			deltaY = Ypos - DDA_L[pos_L]->current[Y];
			DDA_L[pos_L]->current[X] += deltaY * DDA_L[pos_L]->slopeX;
			DDA_L[pos_L]->current[Y] = Ypos;
			DDA_L[pos_L]->current[Z] += deltaY * DDA_L[pos_L]->slopeZ;
		}
		if(DDA_R[pos_R]->current[Y] != Ypos){
			deltaY = Ypos - DDA_R[pos_R]->current[Y];
			DDA_R[pos_R]->current[X] += deltaY * DDA_R[pos_R]->slopeX;
			DDA_R[pos_R]->current[Y] = Ypos;
			DDA_R[pos_R]->current[Z] += deltaY * DDA_R[pos_R]->slopeZ;
		}

		// 7 - Set span DDA current and end positions to right and left edge values
		GzSpan span;
		memcpy(span.start, DDA_L[pos_L]->current, sizeof(GzCoord));
		memcpy(span.end, DDA_R[pos_R]->current, sizeof(GzCoord));
		memcpy(span.current, DDA_L[pos_L]->current, sizeof(GzCoord));
		span.slopeZ = (DDA_R[pos_R]->current[Z] - DDA_L[pos_L]->current[Z])/(DDA_R[pos_R]->current[X] - DDA_L[pos_L]->current[X]);

		// 8 - Advance span current position to left-most covered pixel (ceiling)
		span.current[X] = (ceil(span.start[X]) < 0.0f)?0.0f:ceil(span.start[X]);
		float deltaX = span.current[X] - span.start[X];
		span.current[Z] += deltaX * span.slopeZ;

		//NOT DRAWING pixels covered by bottom horizontal edges
		if(!(num_h == 1 && !DDA_H_isTop && Ypos == Ymax)){
			// 9 - Interpolate span position and parameters (Z) until current position > end
			//NOT DRAWING pixels covered by right edges
			while(span.current[X] < span.end[X] && span.current[X] <= render->display->xres){
				if(GzPixelToDisplay(&span, render->display, &(render->flatcolor)))
					return GZ_FAILURE;
				span.current[X] += 1.0f;
				span.current[Z] += span.slopeZ;		
			}
		}
		//Advance to the next scan line
		Ypos += 1.0f;
	} //END (Ypos < Ymax)

	return GZ_SUCCESS;
}

int GzPixelToDisplay(GzSpan *span, GzDisplay *display, GzColor *flatcolor){
	if(!span || !display || !flatcolor)
		return GZ_FAILURE;

	if(span->current[Z] <= 0 || 
		span->current[X] < 0 || span->current[X] >= display->xres ||
		span->current[Y] < 0 || span->current[Y] >= display->yres)
		return GZ_SUCCESS; //clipping pixel

	int pos = ARRAY((int) span->current[X], (int) span->current[Y]);

	//Test interpolated-Z against FB-Z for each pixel - low Z wins
	if(span->current[Z] < display->fbuf[pos].z)
		//Write color value into FB pixel  (default or computed color)
		GzPutDisplay(display, (int) span->current[X], (int) span->current[Y], 
			ctoi((*flatcolor)[RED]), ctoi((*flatcolor)[GREEN]), ctoi((*flatcolor)[BLUE]), 
			MAX_INTENSITY, (int) span->current[Z]);

	return GZ_SUCCESS;
}

int GzFrustumCull(GzRender *render, GzCoord vertexList[NUM_VERTICES], bool &culled){
	if(!render)
		return GZ_FAILURE;
	
	culled = false;
	
	if(vertexList[0][X] < 0.0f && vertexList[1][X] < 0.0f && vertexList[2][X] < 0.0F){
		culled = true;
		return GZ_SUCCESS;
	}

	if(vertexList[0][X] > render->display->xres && vertexList[1][X] > render->display->xres && vertexList[2][X] > render->display->xres){
		culled = true;
		return GZ_SUCCESS;
	}

	if(vertexList[0][Y] < 0.0f && vertexList[1][Y] < 0.0f && vertexList[2][Y] < 0.0F){
		culled = true;
		return GZ_SUCCESS;
	}

	if(vertexList[0][Y] > render->display->yres && vertexList[1][Y] > render->display->yres && vertexList[2][Y] > render->display->yres){
		culled = true;
		return GZ_SUCCESS;
	}
	
	return GZ_SUCCESS;
}

/////////
//TRANSFORMS' FUNCTIONS
/////////

//Setup Xsp transform
int GzSetupXsp(GzRender *render)
{
	if(!render || !render->display)
		return GZ_FAILURE;
	
	//Populating by row
	render->Xsp[0][0] = render->display->xres/2.0f;
	render->Xsp[0][1] = 0.0f;
	render->Xsp[0][2] = 0.0f;
	render->Xsp[0][3] = render->display->xres/2.0f;
	
	render->Xsp[1][0] = 0.0f;
	render->Xsp[1][1] = -render->display->yres/2.0f;
	render->Xsp[1][2] = 0.0f;
	render->Xsp[1][3] = render->display->yres/2.0f;

	render->Xsp[2][0] = 0.0f;
	render->Xsp[2][1] = 0.0f;
	render->Xsp[2][2] =(float) MAXINT;
	render->Xsp[2][3] = 0.0f;
	
	render->Xsp[3][0] = 0.0f;
	render->Xsp[3][1] = 0.0f;
	render->Xsp[3][2] = 0.0f;
	render->Xsp[3][3] = 1.0f;

	return GZ_SUCCESS;
}

int GzSetupXpi(GzRender *render)
{
	if(!render)
		return GZ_FAILURE;
	GzCamera *pCam = &(render->camera);
	float dInv = tan((pCam->FOV/2.0f) * M_PI / 180.0f);
	
	//Populating the matrix by row
	pCam->Xpi[0][0] = 1.0f;
	pCam->Xpi[0][1] = 0.0f;
	pCam->Xpi[0][2] = 0.0f;
	pCam->Xpi[0][3] = 0.0f;

	pCam->Xpi[1][0] = 0.0f;
	pCam->Xpi[1][1] = 1.0f;
	pCam->Xpi[1][2] = 0.0f;
	pCam->Xpi[1][3] = 0.0f;
	
	pCam->Xpi[2][0] = 0.0f;
	pCam->Xpi[2][1] = 0.0f;
	pCam->Xpi[2][2] = dInv;
	pCam->Xpi[2][3] = 0.0f;

	pCam->Xpi[3][0] = 0.0f;
	pCam->Xpi[3][1] = 0.0f;
	pCam->Xpi[3][2] = dInv;
	pCam->Xpi[3][3] = 1.0f;

	return GZ_SUCCESS;
}

int GzSetupXiw(GzRender *render)
{
	if(!render)
		return GZ_FAILURE;

	GzCamera *pCam = &(render->camera);
	GzCoord vX, vY, vZ;

	//Defining Z-axis
	for(unsigned int i = 0; i < NUM_DIMENSIONS; ++i){
		vZ[i] = pCam->lookat[i] - pCam->position[i];
	}
	float tmpFloat = GzVectorNorm(vZ);
	if(tmpFloat == 0.0f)
		return GZ_FAILURE;
	for(unsigned int i = 0; i < NUM_DIMENSIONS; ++i){
		vZ[i] = vZ[i] / tmpFloat; 
	}

	//Defining Y-axis
	tmpFloat = GzDotProduct(pCam->worldup, vZ);
	for(unsigned int i = 0; i < NUM_DIMENSIONS; ++i){
		vY[i] = pCam->worldup[i] - tmpFloat * vZ[i];
	}
	tmpFloat = GzVectorNorm(vY);
	if(tmpFloat == 0.0f)
		return GZ_FAILURE;
	for(unsigned int i = 0; i < NUM_DIMENSIONS; ++i){
		vY[i] = vY[i] / tmpFloat; 
	}

	//Defining X-axis
	vX[X] = vY[Y] * vZ[Z] - vY[Z] * vZ[Y];
	vX[Y] = vY[Z] * vZ[X] - vY[X] * vZ[Z];
	vX[Z] = vY[X] * vZ[Y] - vY[Y] * vZ[X];

	//Populating Xiw by row
	memcpy(pCam->Xiw[0], vX, sizeof(GzCoord));
	pCam->Xiw[0][3] = - GzDotProduct(vX, pCam->position);
	
	memcpy(pCam->Xiw[1], vY, sizeof(GzCoord));
	pCam->Xiw[1][3] = - GzDotProduct(vY, pCam->position);
	
	memcpy(pCam->Xiw[2], vZ, sizeof(GzCoord));
	pCam->Xiw[2][3] = - GzDotProduct(vZ, pCam->position);

	pCam->Xiw[3][0] = 0.0f;
	pCam->Xiw[3][1] = 0.0f;
	pCam->Xiw[3][2] = 0.0f;
	pCam->Xiw[3][3] = 1.0f;
	
	return GZ_SUCCESS;
}

int GzInitDefaultCamera(GzRender *render)
{
	if(!render)
		return GZ_FAILURE;

	render->camera.position[X] = DEFAULT_IM_X;      
  	render->camera.position[Y] = DEFAULT_IM_Y;
  	render->camera.position[Z] = DEFAULT_IM_Z;

  	render->camera.lookat[X] = 0.0f;
  	render->camera.lookat[Y] = 0.0f;
  	render->camera.lookat[Z] = 0.0f;

  	render->camera.worldup[X] = 0.0f;
  	render->camera.worldup[Y] = 1.0f;
  	render->camera.worldup[Z] = 0.0f;

	render->camera.FOV = DEFAULT_FOV;              /* degrees */

	return GZ_SUCCESS;
}

/////////
//ALGEBRA FUNCTIONS
/////////

float GzVectorNorm(GzCoord v){
	return sqrt(pow(v[X],2) + pow(v[Y],2) + pow(v[Z],2));
}
float GzDotProduct(GzCoord a, GzCoord b){
	return a[X]*b[X] + a[Y]*b[Y] + a[Z]*b[Z];
}
