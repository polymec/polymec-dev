// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// This implements polymesher's capability for writing VTK plots of meshes.

#include <string.h>
#include "core/polymec.h"
#include "core/interpreter.h"
#include "io/generate_face_node_conn.h"
#include "io/generate_cell_node_conn.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// write_vtk_plot(args) -- This function writes a given mesh to a file 
// on disk. Arguments (passed in a table according to Chapter 5.3 of the 
// Lua reference manual) are:
//
// mesh -> mesh object 
// filename -> name of the file to write (1 file only)
int write_t2v_mesh(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (((num_args == 2) && (!lua_ismesh(lua, 1) || !lua_isstring(lua, 2))) || 
      ((num_args == 3) && (!lua_ismesh(lua, 1) || !lua_istable(lua, 2) || !lua_isstring(lua, 3))) || 
      ((num_args != 2) && (num_args != 3)))
  {
    return luaL_error(lua, "write_t2v_mesh: invalid arguments. Usage:\n"
                      "write_t2v_mesh(mesh, filename)");
  }

  // Get the argument(s).
  mesh_t* mesh = lua_tomesh(lua, 1);
  char* filename = strdup(lua_tostring(lua, 2));

  log_info("Writing TOUGH2Viewer mesh with filename '%s'...", filename);

  // Figure out the face-node connectivity.
  int num_faces = mesh->num_faces;
  int face_node_offsets[num_faces+1];
  int *face_nodes;
  generate_face_node_conn(mesh, &face_nodes, face_node_offsets);

  // From that, figure out the cell-node connectivity.
  int num_cells = mesh->num_cells;
  int cell_node_offsets[num_cells+1];
  int *cell_nodes;
  generate_cell_node_conn(mesh, face_nodes, face_node_offsets,
                          &cell_nodes, cell_node_offsets);

  // The file format is simple enough that we don't need an I/O interface
  // to support it.
  FILE* fd = fopen(filename, "w");
  for (int c = 0; c < num_cells; ++c)
  {
    // Construct a string containing the nodes for the cell.
    int num_cell_nodes = cell_node_offsets[c+1] - cell_node_offsets[c];
    char cell_nodes_string[16384];
    int offset = 0;
    for (int n = 0; n < num_cell_nodes; ++n)
    {
      char coords[64];
      int nid = cell_nodes[cell_node_offsets[c]+n];
      snprintf(coords, 64, "(%g,%g,%g) ", mesh->nodes[nid].x, mesh->nodes[nid].y, mesh->nodes[nid].z);
      int len = strlen(coords);
      memcpy(&cell_nodes_string[offset], coords, len*sizeof(char));
      offset += len;
      ASSERT(offset < 16384);
    }
    cell_nodes_string[offset] = '\0';

    // Create a list of node indices for nodes belonging to faces.
    offset = 0;
    char face_nodes_string[16384];
    for (int f = 0; f < mesh->cells[c].num_faces; ++f)
    {
      face_t* face = mesh->cells[c].faces[f];
      int fid = face - &mesh->faces[0];
      int num_nodes = face_node_offsets[fid+1] - face_node_offsets[fid];
      face_nodes_string[offset++] = '(';
      for (int n = 0; n < num_nodes; ++n)
      {
        int nid = face_nodes[face_node_offsets[fid] + n];

        // Search for this node in the list of cell nodes.
        int nn;
        for (nn = 0; nn < num_cell_nodes; ++nn)
        {
          int cnid = cell_nodes[cell_node_offsets[c]+nn];
          if (nid == cnid) break;
        }
        char cn_str[10];
        if (n < num_nodes-1)
          snprintf(cn_str, 10, "%d, ", nn);
        else
          snprintf(cn_str, 10, "%d) ", nn);
        int len = strlen(cn_str);
        memcpy(&face_nodes_string[offset], cn_str, len*sizeof(char));
        offset += len;
      }
      ASSERT(offset < 16384);
    }
    face_nodes_string[offset] = '\0';

    // Create a list of normal vectors for faces.
    offset = 0;
    char face_normals_string[16384];
    for (int f = 0; f < mesh->cells[c].num_faces; ++f)
    {
      face_t* face = mesh->cells[c].faces[f];

      // Normal points from cell 1 to 2.
      vector_t n;
      char n_str[64];
      if (face->cell1 == &mesh->cells[c])
        snprintf(n_str, 64, "(%g,%g,%g) ", face->normal.x, face->normal.y, face->normal.z);
      else
        snprintf(n_str, 64, "(%g,%g,%g) ", -face->normal.x, -face->normal.y, -face->normal.z);
      int len = strlen(n_str);
      memcpy(&face_normals_string[offset], n_str, len*sizeof(char));
      offset += len;
      ASSERT(offset < 16384);
    }
    face_normals_string[offset] = '\0';

    fprintf(fd, "%d %g %g %g %d %s %d %s %s\n",
      c, mesh->cells[c].center.x, mesh->cells[c].center.y, mesh->cells[c].center.z,
      num_cell_nodes, cell_nodes_string, mesh->cells[c].num_faces, face_nodes_string, 
      face_normals_string);
  }
  fclose(fd);

  // Clean up.
  free(cell_nodes);
  free(face_nodes);
  free(filename);
  return 0;
}

