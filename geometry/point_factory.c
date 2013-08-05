// point_factory.c - Implementations of interpreter functions for generating
// sets of points.

#include "core/polymec.h"
#include "core/point.h"
#include "core/interpreter.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

int point_factory_cubic_lattice(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 5)
  {
    lua_pushstring(lua, "Invalid arguments. Usage:\n"
                  "points = point_factory.cubic_lattice(nx, ny, nz, bounding_box, ng)");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the arguments.
  int nx = (int)lua_tonumber(lua, 1);
  int ny = (int)lua_tonumber(lua, 2);
  int nz = (int)lua_tonumber(lua, 3);
  if ((nx <= 0) || (ny <= 0) || (nz <= 0))
  {
    lua_pushstring(lua, "nx, ny, and nz must all be positive.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Number of ghost points.
  if (!lua_isnumber(lua, 5))
  {
    lua_pushstring(lua, "ng must be an integer.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  int ng = (int)lua_tonumber(lua, 5);
  if (ng < 0)
  {
    lua_pushstring(lua, "ng must be non-negative.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Bounding box? 
  if (!lua_isboundingbox(lua, 4))
  {
    lua_pushstring(lua, "bounding_box must be a bounding box.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  bbox_t* bbox = lua_toboundingbox(lua, 4);

  // Now check the bounds.
  if (bbox->x1 >= bbox->x2)
  {
    lua_pushstring(lua, "x1 must be less than x2.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  if (bbox->y1 >= bbox->y2)
  {
    lua_pushstring(lua, "y1 must be less than y2.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  if (bbox->z1 >= bbox->z2)
  {
    lua_pushstring(lua, "z1 must be less than z2.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Pop all the previous arguments off the stack.
  lua_pop(lua, lua_gettop(lua));

  // Create the lattice of points.
  int num_points = (nx + 2*ng) * (ny + 2*ng) * (nz + 2*ng), offset = 0;
  point_t* points = malloc(sizeof(point_t) * num_points);
  double dx = (bbox->x2 - bbox->x1) / nx;
  double dy = (bbox->y2 - bbox->y1) / ny;
  double dz = (bbox->z2 - bbox->z1) / nz;
  for (int i = -ng; i < nx+ng; ++i)
  {
    double xi = (i+0.5) * dx;
    for (int j = -ng; j < ny+ng; ++j)
    {
      double yj = (j+0.5) * dy;
      for (int k = -ng; k < nz+ng; ++k, ++offset)
      {
        double zk = (k+0.5) * dz;
        points[offset].x = xi;
        points[offset].y = yj;
        points[offset].z = zk;
      }
    }
  }

  // Push the points onto the stack.
  lua_pushpointlist(lua, points, num_points);
  return 1;
}

int point_factory_cylinder(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 7)
  {
    lua_pushstring(lua, "Invalid arguments. Usage:\n"
                  "points = point_factory.cylinder(nr, ntheta, nz, ng, r, x0, Z)");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the arguments.
  int nr = (int)lua_tonumber(lua, 1);
  int ntheta = (int)lua_tonumber(lua, 2);
  int nz = (int)lua_tonumber(lua, 3);
  if ((nr <= 0) || (ntheta <= 0) || (nz <= 0))
  {
    lua_pushstring(lua, "nr, ntheta, and nz must all be positive.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  // Number of ghost points.
  if (!lua_isnumber(lua, 4))
  {
    lua_pushstring(lua, "ng must be an integer.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  int ng = (int)lua_tonumber(lua, 4);
  if (ng < 0)
  {
    lua_pushstring(lua, "ng must be non-negative.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Radius.
  if (!lua_isnumber(lua, 5))
  {
    lua_pushstring(lua, "r must be a positive radius.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  double r = lua_tonumber(lua, 5);
  if (r <= 0)
  {
    lua_pushstring(lua, "r must be positive.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Axial point.
  if (!lua_ispoint(lua, 6))
  {
    lua_pushstring(lua, "x0 must be an axial point.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  point_t* x0 = lua_topoint(lua, 6);

  // Axial vector.
  if (!lua_isvector(lua, 7))
  {
    lua_pushstring(lua, "Z must be an axial vector.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  vector_t* Z = lua_tovector(lua, 7);
  double Zmag = vector_mag(Z);
  if (Zmag == 0)
  {
    lua_pushstring(lua, "Z must not be the zero vector.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Pop all the previous arguments off the stack.
  lua_pop(lua, lua_gettop(lua));

  // Create the points. They should be a set of concentric points with equal
  // angular spacing, with a single point at the center.
  int num_points_in_disk = 1 + (nr - 1) * ntheta;
  int num_disks = nz;
  int num_points = num_points_in_disk * nz;
  point_t* points = malloc(sizeof(point_t) * num_points);
  double dr = r / (1.0*nr - 0.5);
  double dtheta = 2.0 * M_PI / ntheta;
  double dz = Zmag / nz;

  // Compute an orthonormal basis for constructing disks.
  vector_t ex, ey;
  compute_orthonormal_basis(Z, &ex, &ey);

  // Find the point at the "bottom" of the cylinder.
  point_t x_bottom = {.x = x0->x - 0.5*Z->x,
                      .y = x0->y - 0.5*Z->y,
                      .z = x0->z - 0.5*Z->z};
  int offset = 0;
  for (int i = 0; i < num_disks; ++i)
  {
    // Find the location of the center of the disk.
    point_t x_center = {.x = x_bottom.x + i * dz * Z->x,
                        .y = x_bottom.y + i * dz * Z->y,
                        .z = x_bottom.z + i * dz * Z->z};

    // Plant a point in the center of the disk.
    points[offset++] = x_center;

    // Construct the other points in the disk.
    for (int j = 0; j < nr-1; ++j)
    {
      double rj = j * dr;
      for (int k = 0; k < ntheta; ++k, ++offset)
      {
        double thetak = k * dtheta;
        points[offset].x = x_center.x + rj * (cos(thetak) * ex.x + sin(thetak) * ey.x);
        points[offset].y = x_center.y + rj * (cos(thetak) * ex.y + sin(thetak) * ey.y);
        points[offset].z = x_center.z + rj * (cos(thetak) * ex.z + sin(thetak) * ey.z);
      }
    }
  }

  // Push the points onto the stack.
  lua_pushpointlist(lua, points, num_points);
  return 1;
}

