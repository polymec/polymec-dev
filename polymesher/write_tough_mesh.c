// This implements polymesher's capability for writing meshes that can 
// be used by TOUGH2 and TOUGH+.

#include <string.h>
#include "core/interpreter.h"
#include "core/point.h"

#ifdef __cplusplus
extern "C" {
#endif

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static void make_elem_name(int elem_name_len, 
                           int elem_number,
                           char* elem_name)
{
  if (elem_name_len == 5)
  {
    // Name is 1 letter + 4 numbers.
    elem_name[0] = 'A' + elem_number / 10000;
    elem_number -= (elem_name[0] - 'A') * 10000;
    elem_name[1] = '0' + elem_number/1000;
    elem_number -= (elem_number/1000) * 1000;
    elem_name[2] = '0' + elem_number/100;
    elem_number -= (elem_number/100) * 100;
    elem_name[3] = '0' + elem_number/10;
    elem_number -= (elem_number/10) * 10;
    elem_name[4] = '0' + elem_number;
    elem_name[5] = ' ';
    elem_name[6] = ' ';
    elem_name[7] = ' ';
    elem_name[8] = '\0';
  }
  else
  {
    // Name is 4 letters + 4 numbers.
    elem_name[0] = 'A' + elem_number / (10000 * 26 * 26 * 26);
    elem_number -= (elem_name[0] - 'A') * (10000 * 26 * 26 * 26);
    elem_name[1] = 'A' + elem_number / (10000 * 26 * 26);
    elem_number -= (elem_name[1] - 'A') * (10000 * 26 * 26);
    elem_name[2] = 'A' + elem_number / (10000 * 26);
    elem_number -= (elem_name[1] - 'A') * (10000 * 26);
    elem_name[3] = 'A' + elem_number / 10000;
    elem_number -= (elem_name[1] - 'A') * 10000;
    elem_name[4] = '0' + elem_number/1000;
    elem_number -= (elem_number/1000) * 1000;
    elem_name[5] = '0' + elem_number/100;
    elem_number -= (elem_number/100) * 100;
    elem_name[6] = '0' + elem_number/10;
    elem_number -= (elem_number/10) * 10;
    elem_name[7] = '0' + elem_number;
    elem_name[8] = '\0';
  }
}

static void make_inactive_flag(int* tag,
                               int tag_len,
                               int elem_number,
                               char* inactive)
{
  *inactive = ' ';
  if (tag == NULL)
    return;
  else
  {
    for (int i = 0; i < tag_len; ++i)
    {
      if (elem_number == tag[i])
      {
        *inactive = 'I';
        break;
      }
    }
  }
}

static void write_tough2_mesh(mesh_t* mesh, 
                              const char* filename, 
                              const char* inactive_tag,
                              int elem_name_len)
{
  FILE* file = fopen(filename, "w");
  if (file == NULL)
    polymec_error("write_tough2_mesh: Could not open file '%s' for writing.", filename);

  fprintf(file, "ELEME\n");
  char inactive;
  int tag_len = 0;
  int* tag = NULL;
  if (inactive_tag != NULL)
    tag = mesh_tag(mesh->cell_tags, inactive_tag, &tag_len);

  char** elem_names = malloc(sizeof(char*) * mesh->num_cells);
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    cell_t* cell = &mesh->cells[c];
    // Figure out the element name, active/inactive flag.
    elem_names[c] = malloc(sizeof(char)*9);
    make_elem_name(elem_name_len, c, elem_names[c]);
    make_inactive_flag(tag, tag_len, c, &inactive);
    fprintf(file, "%s           %.4e%.4e%.4e%.4e%.4e%.4e %c\n", elem_names[c], cell->volume, 0.0, 0.0, cell->center.x, cell->center.y, cell->center.z, inactive);
  }

  fprintf(file, "\n");
  fprintf(file, "CONNE\n");
  char conn_name[17];
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    face_t* face = &mesh->faces[f];

    // Figure out the connection name.
    int cell1 = face->cell1 - &mesh->cells[0];
    if (face->cell2 != NULL)
    {
      int cell2 = face->cell2 - &mesh->cells[0];
      strcpy(&conn_name[0], elem_names[cell1]); 
      strcpy(&conn_name[elem_name_len], elem_names[cell2]); 
      for (int i = 2*elem_name_len; i < 16; ++i)
        conn_name[i] = ' ';
      conn_name[16] = '\0';
      double d1 = point_distance(&face->cell1->center, &face->center);
      double d2 = point_distance(&face->cell2->center, &face->center);
      double A = face->area;
      double beta = 0.0; // FIXME
      fprintf(file, "%s             %d%.4e%.4e%.4e%.4e%.4e\n", conn_name, 3, d1, d2, A, beta, 0.0);
    }
  }

  fclose(file);

  // Clean up.
  for (int c = 0; c < mesh->num_cells; ++c)
    free(elem_names[c]);
  free(elem_names);
}

static void write_tough_plus_mesh(mesh_t* mesh, 
                                  const char* filename, 
                                  const char* inactive_tag,
                                  int elem_name_len)
{
  FILE* file = fopen(filename, "w");
  if (file == NULL)
    polymec_error("write_tough_plus_mesh: Could not open file '%s' for writing.", filename);

  fclose(file);
}

// write_tough_mesh(args) -- This function writes a given mesh to a file 
// on disk. Arguments (passed in a table according to Chapter 5.3 of the 
// Lua reference manual) are:
//
// filename -> name of the file to write (1 file only)
// format -> 'T2', 'T+' (mesh format to use)
// mesh -> mesh object 
// inactive_tag -> the tag within the mesh object that denotes inactive elements.
int write_tough_mesh(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_istable(lua, 1))
  {
    lua_pushstring(lua, "write_tough_mesh: invalid arguments. Usage:\n"
                        "write_tough_mesh{filename [= 'MESH'], format [= 'T2'/'T+'], mesh [= mesh], inactive_tag [= 'inactive'], elem_name_len [= 5]}).");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the argument(s).
  mesh_t* mesh = NULL;
  lua_getfield(lua, 1, "mesh"); // required!
  if (!lua_isnoneornil(lua, 2))
    mesh = lua_tomesh(lua, 2);
  else
  {
    lua_pushstring(lua, "write_tough_mesh: mesh argument is required!");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  lua_pop(lua, 1);

  char* filename = NULL;
  lua_getfield(lua, 1, "filename");
  if (!lua_isnoneornil(lua, 2))
    filename = strdup(lua_tostring(lua, 2));
  lua_pop(lua, 1);

  char* format = NULL;
  lua_getfield(lua, 1, "format");
  if (!lua_isnoneornil(lua, 2))
    format = strdup(lua_tostring(lua, 2));
  lua_pop(lua, 1);

  char* inactive_tag = NULL;
  lua_getfield(lua, 1, "inactive_tag");
  if (!lua_isnoneornil(lua, 2))
    inactive_tag = strdup(lua_tostring(lua, 2));
  lua_pop(lua, 1);

  int elem_name_len = 5;
  lua_getfield(lua, 1, "elem_name_len");
  if (lua_isnumber(lua, 2))
    elem_name_len = (int)lua_tonumber(lua, 2);
  lua_pop(lua, 1);
  if (!lua_isnoneornil(lua, 2) && ((elem_name_len != 5) && (elem_name_len != 8)))
  {
    lua_pushstring(lua, "write_tough_mesh: elem_name_len must be 5 or 8.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Provide defaults.
  if (filename == NULL)
    filename = strdup("MESH");
  if (format == NULL)
    format = strdup("T2");

  // Check our arguments.
  if ((strcasecmp(format, "t2") != 0) && (strcasecmp(format, "t+") != 0))
  {
    char err[128];
    snprintf(err, 128, "write_tough_mesh: unrecognized format: '%s'", format);
    lua_pushstring(lua, err);
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Write the mesh to a file.
  if (!strcasecmp(format, "t2"))
    write_tough2_mesh(mesh, filename, inactive_tag, elem_name_len);
  else
    write_tough_plus_mesh(mesh, filename, inactive_tag, elem_name_len);

  return 1;
}

#ifdef __cplusplus
}
#endif


