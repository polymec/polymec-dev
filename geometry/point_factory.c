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
  // The argument should be a single table of named values.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || (!lua_istable(lua, 1)))
  {
    lua_pushstring(lua, "Invalid arguments. Usage:\n"
                  "points = point_factory.cubic_lattice{nx, ny, nz, num_ghost, bounding_box}");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Extract arguments.
  const char* entries[] = {"nx", "ny", "nz", "num_ghost", "bounding_box"};
  int nx, ny, nz, ng;
  bbox_t* bbox = NULL;
  for (int i = 0; i < 5; ++i)
  {
    lua_pushstring(lua, entries[i]);
    lua_gettable(lua, 1);
    if (i < 4)
    {
      if (!lua_isnumber(lua, -1))
      {
        char err[1024];
        snprintf(err, 1024, "Missing integer argument: %s", entries[i]);
        lua_pushstring(lua, err);
        lua_error(lua);
        return LUA_ERRRUN;
      }
      switch(i) 
      {
        case 0: nx = (int)lua_tonumber(lua, -1); break;
        case 1: ny = (int)lua_tonumber(lua, -1); break;
        case 2: nz = (int)lua_tonumber(lua, -1); break;
        case 3: ng = (int)lua_tonumber(lua, -1); break;
        default: break;
      }
    }
    else 
    {
      if (!lua_isboundingbox(lua, -1))
      {
        lua_pushstring(lua, "bounding_box must be a bounding box.");
        lua_error(lua);
        return LUA_ERRRUN;
      }
      bbox = lua_toboundingbox(lua, -1);
    }
  }

  // Validate arguments.
  if ((nx <= 0) || (ny <= 0) || (nz <= 0))
  {
    lua_pushstring(lua, "nx, ny, and nz must all be positive.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  if (ng < 0)
  {
    lua_pushstring(lua, "ng must be non-negative.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  if (bbox->x1 >= bbox->x2)
  {
    lua_pushstring(lua, "In bounding_box: x1 must be less than x2.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  if (bbox->y1 >= bbox->y2)
  {
    lua_pushstring(lua, "In bounding_box: y1 must be less than y2.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  if (bbox->z1 >= bbox->z2)
  {
    lua_pushstring(lua, "In bounding_box: z1 must be less than z2.");
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
  // The argument should be a single table of named values.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || (!lua_istable(lua, 1)))
  {
    lua_pushstring(lua, "Invalid arguments. Usage:\n"
                  "points = point_factory.cylinder{nr, ntheta, nz, num_ghost, radius, axial_point, axial_vector}");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Extract arguments.
  const char* entries[] = {"nr", "ntheta", "nz", "num_ghost", 
                           "radius", "axial_point", "axial_vector"};
  int nr, ntheta, nz, ng;
  double r;
  point_t* x0 = NULL;
  vector_t* Z = NULL;
  for (int i = 0; i < 7; ++i)
  {
    lua_pushstring(lua, entries[i]);
    lua_gettable(lua, 1);
    if (i < 4)
    {
      if (!lua_isnumber(lua, -1))
      {
        char err[1024];
        snprintf(err, 1024, "Missing integer argument: %s", entries[i]);
        lua_pushstring(lua, err);
        lua_error(lua);
        return LUA_ERRRUN;
      }
      switch(i) 
      {
        case 0: nr = (int)lua_tonumber(lua, -1); break;
        case 1: ntheta = (int)lua_tonumber(lua, -1); break;
        case 2: nz = (int)lua_tonumber(lua, -1); break;
        case 3: ng = (int)lua_tonumber(lua, -1); break;
        default: break;
      }
    }
    else if (i == 4)
    {
      if (!lua_isnumber(lua, -1))
      {
        lua_pushstring(lua, "radius must be a positive number.");
        lua_error(lua);
        return LUA_ERRRUN;
      }
      r = lua_tonumber(lua, -1);
    }
    else if (i == 5)
    {
      if (!lua_ispoint(lua, -1))
      {
        lua_pushstring(lua, "axial_point must be a point.");
        lua_error(lua);
        return LUA_ERRRUN;
      }
      x0 = lua_topoint(lua, -1);
    }
    else // (i == 6)
    {
      if (!lua_isvector(lua, -1))
      {
        lua_pushstring(lua, "axial_vector must be a vector.");
        lua_error(lua);
        return LUA_ERRRUN;
      }
      Z = lua_tovector(lua, -1);
    }
  }

  // Validate inputs.
  if ((nr <= 0) || (ntheta <= 0) || (nz <= 0))
  {
    lua_pushstring(lua, "nr, ntheta, and nz must all be positive.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  if (ng < 0)
  {
    lua_pushstring(lua, "num_ghost must be non-negative.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  if (r <= 0)
  {
    lua_pushstring(lua, "radius must be positive.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  double Zmag = vector_mag(Z);
  if (Zmag == 0)
  {
    lua_pushstring(lua, "axial_vector must not be the zero vector.");
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
    point_t x_center = {.x = x_bottom.x + (i+0.5) * dz * Z->x,
                        .y = x_bottom.y + (i+0.5) * dz * Z->y,
                        .z = x_bottom.z + (i+0.5) * dz * Z->z};

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

