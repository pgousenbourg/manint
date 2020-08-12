////////////////////////////////////////////////////////////
//
// Name:   serialize.cpp
//
// Author: Steven Michael
//
// Date:   6/25/05
//
// Description:
//
//   Allows for stored serialized version of base "Model3D" class
//   in MATLAB.  A "model3d" object can be constructed by calling
//   the "model3d" constructor with the output of this function as
//   the argument.
//
// Copyright (C) 2005 MIT Lincoln Laboratory
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, 
// Boston, MA  02110-1301  USA
//
////////////////////////////////////////////////////////////


#include <mex.h>
#include <model3d.h>

#ifndef WIN32
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#else
#include <fstream>
#endif

#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  if(nrhs < 2) {
    mexPrintf("Must pass in a model & a filename\n");
    return;
  }
  
  if(mxIsStruct(prhs[0])==0) {
    mexPrintf("First argument must be of type Model3D\n");
    return;
  }
  if(mxIsChar(prhs[1])==0) {
    mexPrintf("Second argument must be a string\n");
    return;
  }
  
  Model3D *model = new Model3D;
  // Populate the class from the MATLAB variable
  mxArray *mxLayer = mxGetField(prhs[0],0,"layers");
  if(!mxLayer) {
    mexPrintf("Cannot find layer member\n");
    return;
  }
  model->nLayers = (int)mxGetNumberOfElements(mxLayer);
  model->layers = (Layer **)malloc(sizeof(Layer *)*model->nLayers);
  
  for(int i=0;i<model->nLayers;i++) {
    Layer *layer = (Layer *)malloc(sizeof(Layer));
    mxArray *mxFacet = mxGetField(mxLayer,i,"facets");
    mxArray *mxFacetIdx = mxGetField(mxLayer,i,"facetidx");
    mxArray *mxVertices = mxGetField(mxLayer,i,"vertices");
    
    if(mxFacetIdx) 
      layer->nFacetIdx = (int) mxGetNumberOfElements(mxFacetIdx)/4;
    else {
      layer->nFacetIdx = 0;
      layer->facetIdx = (FacetIdx *)0;
    }
    if(mxFacet)
      layer->nFacets = (int) mxGetNumberOfElements(mxFacet)/3;
    else {
      layer->nFacets = 0;
      layer->facets = (Facet *)0;
    }
    if(mxVertices)
      layer->nVertices = (int) mxGetNumberOfElements(mxVertices)/3;
    else {
      layer->nVertices = 0;
      layer->vertices = (Vector *)0;
    }
    if(layer->nFacets) {
      layer->facets = (Facet *)malloc(sizeof(Facet)*layer->nFacets);
      memcpy(layer->facets,mxGetPr(mxFacet),sizeof(Facet)*layer->nFacets);
    }
     if(layer->nFacetIdx) {
      layer->facetIdx = (FacetIdx *)malloc(sizeof(FacetIdx)*layer->nFacetIdx);
      memcpy(layer->facetIdx,mxGetPr(mxFacetIdx),sizeof(FacetIdx)*layer->nFacetIdx);
    }
    if(layer->nVertices) {
      layer->vertices =(Vector *)malloc(sizeof(Vector)*layer->nVertices);
      memcpy(layer->vertices,mxGetPr(mxVertices),sizeof(Vector)*layer->nVertices);
    }
    memcpy(layer->ambient,
           mxGetPr(mxGetField(mxLayer,i,"ambient")),sizeof(Vector));
    memcpy(layer->diffuse,
           mxGetPr(mxGetField(mxLayer,i,"diffuse")),sizeof(Vector));
    memcpy(layer->specular,
           mxGetPr(mxGetField(mxLayer,i,"specular")),sizeof(Vector));
    mxGetString(mxGetField(mxLayer,i,"name"),
                layer->name,64);
    layer->shininess = (float) mxGetScalar(mxGetField(mxLayer,i,"shininess"));
    layer->shinystrength = (float) mxGetScalar(mxGetField(mxLayer,i,"shinystrength"));
    layer->transparency = (float) mxGetScalar(mxGetField(mxLayer,i,"transparency")); 
    model->layers[i] = layer;
  }
  mxArray *mxLight = mxGetField(prhs[0],0,"lights");
  if(mxLight)
    model->nLights = (int)mxGetNumberOfElements(mxLight);
  else
    model->nLights = 0;
  for(int i=0;i<model->nLights;i++) {
    Light *light = (Light *)malloc(sizeof(Light));
    memcpy(light->pos,mxGetPr(mxGetField(mxLight,i,"pos")),
           sizeof(Vector));
    memcpy(light->dir,mxGetPr(mxGetField(mxLight,i,"dir")),
           sizeof(Vector));
    model->lights[i] = light;    
  } 
  mxArray *mxCamPos = mxGetField(prhs[0],0,"campos");
  mxArray *mxCamTar = mxGetField(prhs[0],0,"camtarget");
  if(mxCamPos && mxCamTar) {
    model->useCam = 1;
    memcpy(model->camPos,mxGetPr(mxCamPos),sizeof(Vector));
    memcpy(model->camTarget,mxGetPr(mxCamTar),sizeof(Vector));
  }
  else
    model->useCam = 0; 
  
  // Do the serialization
  void *mem;
  int len;  
  model->setlog(mexPrintf);
  model->serialize(&mem,len);
 
  char filename[256];
  mxGetString(prhs[1],filename,256);
#if !defined(WIN32)
  int fd = open(filename,O_WRONLY|O_CREAT|O_TRUNC,0644);
  if(fd == -1) {
    mexPrintf("Cannot open output file\n");
    delete model;
    free(mem);
    return;
  }
  write(fd,mem,len);
  close(fd);
#else
	std::ofstream of(filename,std::ios::binary | std::ios::out);
  if(of.is_open()==0) {
    mexPrintf("Cannot open output file: \"%s\"\n",filename);
    delete model;
    free(mem);
    return;
  }
  of.write((const char *)mem,len);
  of.close();
#endif
    
  delete model;
  free(mem);  
  return;
} // end of mexFunction
