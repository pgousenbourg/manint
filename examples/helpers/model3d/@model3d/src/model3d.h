////////////////////////////////////////////////////////////
//
// Name:   model3d.h
//
// Author: Steven Michael
//
// Date:   6/25/05
//
// Description:
//
//   Header of base class for storing 3D scene data.
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

#ifndef _MODEL3D_H_
#define _MODEL3D_H_

typedef float Vector[3];
typedef int FacetIdx[4];
typedef Vector Facet[4];
typedef float RGB[3];

typedef struct {
  FacetIdx *facetIdx;
  int nFacetIdx;
  Facet *facets;
  int nFacets;
  Vector *vertices;
  int nVertices;
  RGB diffuse;
  RGB ambient;
  RGB specular;
  float transparency;
  float shininess;
  float shinystrength;
	int useTexture;
  char name[64];
	char textureFile[64];
} Layer;    

typedef struct {
  Vector pos;
  Vector dir;
}Light;

#define MAX_LIGHTS 10

class Model3D 
{
  public:
    Model3D();
    virtual ~Model3D();
    
    Layer **layers;
    int    nLayers;

    virtual void setlog(int (*setlog)(const char *,...));    
  
    // This must be overloaded with a function
    // that returns the model type
    virtual const char *model_type();
  
    int nLights;
    Light *lights[MAX_LIGHTS];
        
    // camera position & direction
    // if defined in the file.
    Vector camPos;
    Vector camTarget;
    bool   useCam;
  
    static const char *header();
    int serialize(void **mem, int &len);
    static Model3D *unserialize(void *mem);
    
  protected:
    
    int (*logmsg)(const char *,...);
}; // end of Model3D

#endif
