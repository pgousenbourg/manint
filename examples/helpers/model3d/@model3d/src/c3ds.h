////////////////////////////////////////////////////////////
//
// Name:   c3ds.h
//
// Author: Steven Michael
//
// Date:   6/25/05
//
// Description:
//
//   Definition of class for reading ".3ds" files.  
//   Class is derived from base "Model3D" class for storing
//   of 3-D data.
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

#ifndef _C3DS_H_
#define _C3DS_H_

#include <model3d.h>

typedef unsigned char byte;

typedef struct  {
  char name[64];
  RGB  ambient;
  RGB  diffuse;
  RGB  specular;
  float shininess;
  float shinystrength;
  float transparency;
	char textureFile[64];
	int useTexture;
} Material3DS;

// Forward declaration 
class C3DS;

typedef int (C3DS::*ChunkHandler)(byte *, int);

typedef struct {
  unsigned short ID;
  char name[64];
  ChunkHandler handler;
} ChunkStruct;

class C3DS : public Model3D {
  public:
    C3DS(const char *filename = (const char *)0);
    virtual ~C3DS();

    int read_file(const char *filename);

    const char *model_type() {
      static char nm[] = "3DS";
      return nm;
    }
    
  protected:
    int read_main_chunk(byte *s, int len);
    int read_0x4000_chunk(byte *s, int len);
    int read_0x4120_chunk(byte *s, int len);
    int read_0x4130_chunk(byte *s, int len);
    int read_0xb000_chunk(byte *s, int len);
    int read_0xb005_chunk(byte *s, int len);
    int read_0xa000_chunk(byte *s, int len);
    int read_0xa010_chunk(byte *s, int len);
    int read_0xa020_chunk(byte *s, int len);
    int read_0xa030_chunk(byte *s, int len);
    int read_0xa040_chunk(byte *s, int len);
    int read_0xa041_chunk(byte *s, int len);
    int read_0xa050_chunk(byte *s, int len);
		int read_0xa200_chunk(byte *s, int len);
		int read_0xa300_chunk(byte *s, int len);
    int read_0x4110_chunk(byte *s, int len);
    int read_0x4160_chunk(byte *s, int len);
    int read_0x4700_chunk(byte *s, int len);
    int read_0x4600_chunk(byte *s, int len);
    int read_0x4610_chunk(byte *s, int len);
    int read_color(byte *,RGB);
    int read_number(byte *, float &);
    
    static ChunkStruct chunkStruct[];
    
    int (*logmsg)(const char *,...); 
    int add_facetidx(FacetIdx *, int);
    int add_vertex(Vector *, int);
    
    // This creates a current layer and adds it
    // to the end of the list
    Layer *delLayer;
    int delete_layer();
    int create_layer();
    Layer *curLayer;
    
    Light *curLight;
  
    int apply_material(const char *matname);
    int create_material();
    Material3DS **materials;
    int           nMaterials;
    Material3DS  *curMaterial;
};

#endif
