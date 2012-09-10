#include <stdlib.h>
#include <string.h>
#include "silo.h"
#include "core/silo_io.h"
#include "math.h"

#ifdef USE_MPI
#include <mpi.h>
#include "pmpio.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

static void* silo_create_file(void* context, 
                              const char* filename,
                              const char* dirname)
{
  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
  DBMkDir(file, dirname);
  DBSetDir(file, dirname);
  return (void*)file;
}

static void* silo_open_file(void* context, 
                            const char* filename, 
                            const char* dirname,
                            io_mode_t mode)
{
  int driver = DB_HDF5;
  DBfile* file;
  if (mode == IO_WRITE)
  { 
    file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
    DBMkDir(file, dirname);
    DBSetDir(file, dirname);
  }
  else
  {
    file = DBOpen(filename, driver, DB_READ);
    DBSetDir(file, dirname);
  }
  return (void*)file;
}

static void silo_close_file(void* context, void* file)
{
  // In this method, we do all the writing.
  DBClose(file);
}

static int silo_get_num_datasets(void* context, void* file, int* num_datasets)
{
  return 1;
}

static void silo_read_datasets(void* context, void* f, io_dataset_t** datasets, int num_datasets)
{
  DBfile* file = (DBfile*)f;

#if 0
  for (int d = 0; d < num_datasets; ++d)
  {
    const char* name = "default"; // FIXME 
    DBtoc* contents = DBGetToc(file);
    int num_fields = contents->nucdvar;
    int num_sources = 0; // FIXME

    io_dataset_t* dataset = io_dataset_new(name, num_fields, num_sources); 
    datasets[d] = dataset;

    // Retrieve the mesh. Note that we must deallocate the storage 
    // for this object after we're through!
    DBucdmesh* dbmesh = DBGetUcdmesh(file, "mesh");

    // Read the elements from the dbmesh.
    int num_nodes = dbmesh->nnodes;
    int num_faces = dbmesh->nfaces;
    int num_cells = dbmesh->nzones;
    int num_ghost_cells = 0; // FIXME

    // The dbmesh doesn't contain edge information, so we have to 
    // count them up here. The number of edges is the number of distinct 
    // pairs of nodes that appear consecutively in a face.
    int num_edges = 0;
    int noffset = 0;
    for (int f = 0; f < num_faces; ++f)
    {

    mesh.faces[f].resize(dbmesh->faces->shapesize[f]);
    for (int n = 0; n < mesh.faces[f].size(); ++n, ++noffset)
      mesh.faces[f][n] = dbmesh->faces->nodelist[noffset];
    }

    // Finally, we set up our proper mesh.
    mesh_t* mesh = mesh_new(num_cells, num_ghost_cells, num_faces,
                            num_edges, num_nodes);
    dataset->mesh = mesh;

    // Node coordinates.
    for (int n = 0; n < num_nodes; ++n)
    {
      mesh->nodes[i].x  = ((double*)(dbmesh->coords[i]))[0];
      mesh->nodes[i].y  = ((double*)(dbmesh->coords[i]))[1];
      mesh->nodes[i].z  = ((double*)(dbmesh->coords[i]))[2];
    }

    // Reconstruct the faces.
    int noffset = 0;
    for (int f = 0; f < mesh.faces.size(); ++f)
    {
      mesh_
    mesh.faces[f].resize(dbmesh->faces->shapesize[f]);
    for (int n = 0; n < mesh.faces[f].size(); ++n, ++noffset)
      mesh.faces[f][n] = dbmesh->faces->nodelist[noffset];
  }

  // Reconstruct the cell-face connectivity.
  mesh.cells.resize(dbmesh->zones->nzones);
  DBcompoundarray* conn = DBGetCompoundarray(file, "cell-face-conn");
  if (conn == 0)
  {
    DBClose(file);
    char err[1024];
    snprintf(err, 1024, "Could not find cell-face connectivity in file %s.", filename);
    error(err);
  }
  // First element is the number of faces in each zone.
  // Second element is the list of face indices in each zone.
  // Third element is a pair of cells for each face.
  if ((conn->nelems != 3) or 
      (conn->elem_lengths[0] != dbmesh->zones->nzones) or 
      (conn->elem_lengths[2] != 2*dbmesh->faces->nfaces))
  {
    DBClose(file);
    char err[1024];
    snprintf(err, 1024, "Found invalid cell-face connectivity in file %s.", filename);
    error(err);
  }
  int* connData = (int*)conn->values;
  int foffset = dbmesh->zones->nzones;
  for (int c = 0; c < dbmesh->zones->nzones; ++c)
  {
    int nfaces = connData[c];
    mesh.cells[c].resize(nfaces);
    copy(connData + foffset, connData + foffset + nfaces, mesh.cells[c].begin());
    foffset += nfaces;
  }
  mesh.faceCells.resize(mesh.faces.size());
  for (size_t f = 0; f < mesh.faceCells.size(); ++f)
  {
    mesh.faceCells[f].resize(2);
    mesh.faceCells[f][0] = connData[foffset];
    mesh.faceCells[f][1] = connData[foffset+1];
    foffset += 2;
  }
  DBFreeUcdmesh(dbmesh);
  DBFreeCompoundarray(conn);

  // Reconstruct the face-edge connectivity.
  
  // Check for convex hull data.
  // First element is the number of facets.
  // Second element is the array of numbers of nodes per facet.
  // Third element is the array of node indices for the facets.
  DBcompoundarray* hull = DBGetCompoundarray(file, "convexhull");
  if (hull != 0)
  {
    if ((hull->nelems != 3) or (hull->elem_lengths[0] != 1))
    {
      DBClose(file);
      char err[1024];
      snprintf(err, 1024, "Found invalid convex hull data in file %s.", filename);
      error(err);
    }
    int* hullData = (int*)conn->values;
    int nfacets = hullData[0];
    mesh.convexHull.facets.resize(nfacets);
    int foffset = 1;
    for (int f = 0; f < nfacets; ++f, ++foffset)
    {
      int nnodes = hullData[foffset];
      mesh.convexHull.facets[f].resize(nnodes);
    }
    for (int f = 0; f < nfacets; ++f, ++foffset)
    {
      for (int n = 0; n < mesh.convexHull.facets[f].size(); ++n)
        mesh.convexHull.facets[f][n] = hullData[foffset];
    }
    DBFreeCompoundarray(hull);
  }

  // Retrieve the fields.
  for (int f = 0; f < num_fields; ++f)
  {
    DBucdvar* dbvar = DBGetUcdvar(file, contents->ucdvar_names[f]);
    ASSERT(dbvar != NULL);
    int ncomp = 1; // FIXME
    double* field = malloc(ncomp*dbvar->nels*sizeof(double));
    memcpy(field, &dbvar->vals[0], sizeof(double)*dbvar->nels*ncomp);
    mesh_centering_t centering = MESH_CELL;
    io_dataset_write_field(dataset, contents->ucdvar_names[f], field, ncomp, centering);

    // Clean up.
    DBFreeUcdvar(dbvar);
  }
#endif
}

static void silo_write_datasets(void* context, void* file, io_dataset_t** datasets, int num_datasets, int rank_in_group)
{
}

static void silo_dtor(void* context)
{
  free(context);
}

io_interface_t* silo_io_new(MPI_Comm comm,
                            int num_files,
                            int mpi_tag)
{
  io_vtable vtable = {.create_file = &silo_create_file,
                      .open_file = &silo_open_file, 
                      .close_file = &silo_close_file,
                      .get_num_datasets = &silo_get_num_datasets,
                      .read_datasets = &silo_read_datasets,
                      .write_datasets = &silo_write_datasets,
                      .dtor = &silo_dtor};
  return io_interface_new(NULL, "Silo", vtable, comm, num_files, mpi_tag);
}

// This is used to generate .silo plot files.
static void silo_plot_write_datasets(void* context, void* file, io_dataset_t** datasets, int num_datasets, int rank_in_group)
{
#if 0
  io_dataset_t* dataset = datasets[0];

  // This is optional for now, but we'll give it anyway.
  char *coordnames[3];
  coordnames[0] = (char*)"xcoords";
  coordnames[1] = (char*)"ycoords";
  coordnames[2] = (char*)"zcoords";

  // Node coordinates.
  int num_nodes = dataset->mesh->num_nodes;
  double x[num_nodes], y[num_nodes], z[num_nodes];
  for (int i = 0; i < num_nodes; ++i)
  {
    x[i] = mesh->nodes[i].x;
    y[i] = mesh->nodes[i].y;
    z[i] = mesh->nodes[i].z;
  }
  double* coords[3];
  coords[0] = &(x[0]);
  coords[1] = &(y[0]);
  coords[2] = &(z[0]);

  // Figure out face-node connectivity. We do this by computing centers
  // for all the cells and then using them to define face normals, 
  // from which node orderings can be determining using a convex hull 
  // determination algorithm (gift wrapping).
  int num_cells = mesh.cells.size();
  int numFaces = mesh.faces.size();
  cell_centers[num_cells];
  for (int c = 0; c < num_cells; ++c)
  {
    const vector<int>& cellFaces = mesh.cells[c];
    int num_nodes = 0;
    for (int f = 0; f < cellFaces.size(); ++f)
    {
      const vector<unsigned>& face_nodes = mesh.faces[cellFaces[f]];
      for (int n = 0; n < face_nodes.size(); ++n)
      {
        cell_centers[3*c]   += mesh.nodes[3*face_nodes[n]];
        cell_centers[3*c+1] += mesh.nodes[3*face_nodes[n]+1];
        cell_centers[3*c+2] += mesh.nodes[3*face_nodes[n]+2];
        ++num_nodes;
      }
    }
    cell_centers[3*c+0] /= num_nodes;
    cell_centers[3*c+1] /= num_nodes;
    cell_centers[3*c+2] /= num_nodes;
  }
  vector<int> face_node_counts(numFaces), 
    all_face_nodes;
  for (int f = 0; f < mesh.faces.size(); ++f)
  {
    const vector<unsigned>& face_nodes = mesh.faces[f];

    // Compute the normal vector for the face, pointing outward from 
    // its first cell.
    ASSERT(face_nodes.size() >= 3);
    double face_center[3];
    for (int n = 0; n < face_nodes.size(); ++n)
    {
      face_center[0] += mesh.nodes[3*face_nodes[n]];
      face_center[1] += mesh.nodes[3*face_nodes[n]+1];
      face_center[2] += mesh.nodes[3*face_nodes[n]+2];
    }
    face_center[0] /= face_nodes.size();
    face_center[1] /= face_nodes.size();
    face_center[2] /= face_nodes.size();

    // Construct vectors v1, v2, and v3, where v1 is the vector pointing from the 
    // face center to the first face node, v2 is a vector pointing from the face 
    // center to any other face, node, and v3 is their cross product.
    double v1[3];
    for (int d = 0; d < 3; ++d)
      v1[d] = mesh.nodes[3*face_nodes[0]+d] - face_center[d];

    double v2[3], normal[3], normalMag;
    for (int n = 1; n < face_nodes.size(); ++n)
    {
      for (int d = 0; d < 3; ++d)
        v2[d] = mesh.nodes[3*face_nodes[n]+d] - face_center[d];

      // normal = v1 x v2.
      normal[0] = v1[1]*v2[2] - v1[2]*v2[1];
      normal[1] = v1[2]*v2[0] - v1[0]*v2[2];
      normal[1] = v1[0]*v2[1] - v1[1]*v2[0];
      normalMag = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
      if (normalMag > 1e-14) break;
    }
    normal[0] /= normalMag; normal[1] /= normalMag; normal[2] /= normalMag;

    double v3[3], cell_center[3];
    for (int d = 0; d < 3; ++d)
    {
      cell_center[d] = cell_centers[3*mesh.faceCells[f][0]+d];
      v3[d] = face_center[d] - cell_center[d];
    }
    if ((normal[0]*v3[0] + normal[1]*v3[1] + normal[2]*v3[2]) < 0.0)
    {
      normal[0] *= -1.0; normal[1] *= -1.0; normal[2] *= -1.0;
    }

    // Now project the coordinates of the face's nodes to the plane
    // with the given normal and centered about the face center.
    vector<RealType> points(2*face_nodes.size()); // NOTE: planar coordinates (2D)
    double e1[3], e2[3]; // Basis vectors in the plane.
    double v1Mag = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    for (int d = 0; d < 3; ++d)
      e1[d] = v1[d] / v1Mag;

    // e2 = normal x e1.
    e2[0] = normal[1]*e1[2] - normal[2]*e1[1];
    e2[1] = normal[2]*e1[0] - normal[0]*e1[2];
    e2[1] = normal[0]*e1[1] - normal[1]*e1[0];
    for (int p = 0; p < points.size(); ++p)
    {
      // v = node center - cell center.
      double v[3];
      for (int d = 0; d < 3; ++d)
        v[d] = mesh.nodes[3*face_nodes[p]+d] - cell_center[d];

      // Compute the perpendicular component of the point
      // with location v:
      // vPerp = v - (n o v)n.
      double vPerp[3];
      for (int d = 0; d < 3; ++d)
        vPerp[d] = v[d] - (normal[0]*v[0] + normal[1]*v[1] + normal[2]*v[2]) * normal[d];

      // Project it to the plane.
      points[2*p]   = vPerp[0]*e1[0] + vPerp[1]*e1[1] + vPerp[2]*e1[2];
      points[2*p+1] = vPerp[0]*e2[0] + vPerp[1]*e2[1] + vPerp[2]*e2[2];
    }

    // Find the node order by traversing the convex hull of 
    // the points within the plane, appending them to all_face_nodes.
    vector<int> indices;
    traverseConvexHull(points, indices);
    face_node_counts[f] = indices.size();
    for (int n = 0; n < indices.size(); ++n)
      all_face_nodes.push_back(face_nodes[indices[n]]);
  }

  // Figure out cell-face connectivity.
  vector<int> cellFaceCounts(num_cells), 
    all_cell_faces;
  for (int c = 0; c < num_cells; ++c)
  {
    const vector<int>& cellFaces = mesh.cells[c];
    cellFaceCounts[c] = cellFaces.size();
    all_cell_faces.insert(all_cell_faces.end(), 
        cellFaces.begin(), cellFaces.end());
  }

  // The polyhedral zone list is referred to in the options list.
  DBAddOption(optlist, DBOPT_PHZONELIST, (char*)"mesh_zonelist");

  // Write out the 3D polyhedral mesh.
  DBPutUcdmesh(file, (char*)"mesh", 3, coordnames, coords,
      num_nodes, num_cells, 0, 0,
      DB_DOUBLE, optlist); 

  // Write the connectivity information.
  DBPutPHZonelist(file, (char*)"mesh_zonelist", 
      face_node_counts.size(), &face_node_counts[0], 
      all_face_nodes.size(), &all_face_nodes[0], 0, 
      cellFaceCounts.size(), &cellFaceCounts[0],
      all_cell_faces.size(), &all_cell_faces[0], 
      0, 0, num_cells-1, optlist);

  // Write out the cell-face connectivity data.
  vector<int> conn(num_cells);
  int elem_lengths[3];
  char* elem_names[3];
  for (int c = 0; c < num_cells; ++c)
    conn[c] = mesh.cells[c].size();
  for (int c = 0; c < num_cells; ++c)
  {
    for (int f = 0; f < mesh.cells[c].size(); ++f)
      conn.push_back(mesh.cells[c][f]);
  }
  for (int f = 0; f < mesh.faceCells.size(); ++f)
  {
    conn.push_back(mesh.faceCells[f][0]);
    conn.push_back(mesh.faceCells[f][0]);
  }
  elem_names[0] = strdup("ncellfaces");
  elem_lengths[0] = num_cells;
  elem_names[2] = strdup("facecells");
  elem_lengths[2] = conn.size() - 2*mesh.faces.size();
  elem_names[1] = strdup("cellfaces");
  elem_lengths[1] = conn.size() - elem_lengths[2] - elem_lengths[0];
  DBPutCompoundarray(file, "conn", elem_names, elem_lengths, 3, 
      (void*)&conn[0], conn.size(), DB_INT, 0);
  free(elem_names[0]);
  free(elem_names[1]);
  free(elem_names[2]);

  // Write out convex hull data.
  vector<int> hull(1+mesh.convexHull.facets.size());
  hull[0] = mesh.convexHull.facets.size();
  for (int f = 0; f < mesh.convexHull.facets.size(); ++f)
    hull[1+f] = mesh.convexHull.facets[f].size();
  for (int f = 0; f < mesh.convexHull.facets.size(); ++f)
    for (int n = 0; n < mesh.convexHull.facets[f].size(); ++n)
      hull.push_back(mesh.convexHull.facets[f][n]);
  elem_names[0] = strdup("nfacets");
  elem_lengths[0] = 1;
  elem_names[1] = strdup("nfacetnodes");
  elem_lengths[1] = mesh.convexHull.facets.size();
  elem_names[2] = strdup("facetnodes");
  elem_lengths[2] = hull.size() - elem_lengths[0] - elem_lengths[1];
  DBPutCompoundarray(file, "convexhull", elem_names, elem_lengths, 3, 
      (void*)&hull[0], hull.size(), DB_INT, 0);
  free(elem_names[0]);
  free(elem_names[1]);
  free(elem_names[2]);

  // Write out the cell-centered mesh data.

  // Scalar fields.
  for (typename map<string, RealType*>::const_iterator iter = fields.begin();
      iter != fields.end(); ++iter)
  {
    DBPutUcdvar1(file, (char*)iter->first.c_str(), (char*)"mesh",
        (void*)iter->second, num_cells, 0, 0,
        DB_DOUBLE, DB_ZONECENT, optlist);
  }

  // Clean up.
  DBFreeOptlist(optlist);

#ifdef HAVE_MPI
  // Write the multi-block objects to the file if needed.
  int num_chunks = nproc / num_files;
  if (rankInGroup == 0)
  {
    char* mesh_names[num_chunks];
    int mesh_types[num_chunks];
    char* var_names[num_fields];
    for (int i = 0; i < num_chunks; ++i)
    {
      mesh_types[i] = DB_UCDMESH;
      var_types[i] = DB_UCDVAR;
    }
    vector<vector<char*> > var_names(fields.size());
    vector<int> var_types(num_chunks, DB_UCDVAR);
    for (int i = 0; i < num_chunks; ++i)
    {
      // Mesh.
      char mesh_name[1024];
      snprintf(mesh_name, 1024, "domain_%d/mesh", i);
      mesh_names[i] = strdup(mesh_name);

      // Field data.
      int field_index = 0;
      for (typename map<string, RealType*>::const_iterator iter = fields.begin();
          iter != fields.end(); ++iter, ++field_index)
      {
        char var_name[1024];
        snprintf(var_name, 1024, "domain_%d/%s", i, iter->first.c_str());
        var_names[field_index].push_back(strdup(var_name));
      }
    }
#endif

    // Write the mesh and variable data.
    DBSetDir(file, "/");
    DBPutMultimesh(file, "mesh", num_chunks, &mesh_names[0], 
        &mesh_types[0], optlist);
    int field_index = 0;
    for (typename map<string, RealType*>::const_iterator iter = fields.begin();
        iter != fields.end(); ++iter, ++field_index)
    {
      DBPutMultivar(file, iter->first.c_str(), num_chunks, 
          &var_names[field_index][0], &var_types[0], optlist);
    }

    // Clean up.
    DBFreeOptlist(optlist);
    for (int i = 0; i < num_chunks; ++i)
      free(mesh_names[i]);
    for (int f = 0; f < var_names.size(); ++f)
      for (int i = 0; i < num_chunks; ++i)
        free(var_names[f][i]);
  }
#endif
}

static void silo_plot_write_master(void* context, void* file, const char* prefix, io_dataset_t** datasets, int num_datasets, int num_files, int procs_per_file)
{
  io_dataset_t* dataset = datasets[0];

  char* mesh_names[num_files*procs_per_file];
  int mesh_types[num_files*procs_per_file];
  int var_types[dataset->num_fields][num_files*procs_per_file];
  for (int i = 0; i < num_files*procs_per_file; ++i)
  {
    mesh_types[i] = DB_UCDMESH;
    for (int f = 0; f < dataset->num_fields; ++f)
      var_types[f][i] = DB_UCDVAR;
  }
  char* var_names[dataset->num_fields][num_files*procs_per_file];

  for (int i = 0; i < num_files; ++i)
  {
    for (int c = 0; c < procs_per_file; ++c)
    {
      // Mesh.
      char mesh_name[1024];
      snprintf(mesh_name, 1024, "%d/%s.silo:/domain_%d/mesh", i, prefix, c);
      mesh_names[i*procs_per_file+c] = strdup(mesh_name);

      // Field data.
      for (int f = 0; f < dataset->num_fields; ++f)
      {
        char var_name[1024];
        snprintf(var_name, 1024, "%d/%s.silo:/domain_%d/%s", i, prefix, c, dataset->field_names[f]);
        var_names[f][i] = strdup(var_name);
      }
    }
  }

  DBoptlist* optlist = DBMakeOptlist(10);

  // Write the multimesh and variable data, and close the file.
  DBPutMultimesh(file, "mesh", num_files*procs_per_file, &mesh_names[0], 
      &mesh_types[0], optlist);
  for (int f = 0; f < dataset->num_fields; ++f)
  {
    DBPutMultivar(file, dataset->field_names[f], num_files*procs_per_file, 
      var_names[f], var_types[f], optlist);
  }
  DBClose(file);

  // Clean up.
  DBFreeOptlist(optlist);
  for (int i = 0; i < num_files*procs_per_file; ++i)
    free(mesh_names[i]);
  for (int f = 0; f < dataset->num_fields; ++f)
    for (int i = 0; i < num_files*procs_per_file; ++i)
      free(var_names[f][i]);
}

io_interface_t* silo_plot_io_new(MPI_Comm comm,
                                 int num_files,
                                 int mpi_tag)
{
  io_vtable vtable = {.create_file = &silo_create_file,
                      .open_file = &silo_open_file, 
                      .close_file = &silo_close_file,
                      .get_num_datasets = &silo_get_num_datasets,
                      .write_datasets = &silo_plot_write_datasets,
                      .write_master = &silo_plot_write_master};
  return io_interface_new(NULL, "Silo-plot", vtable, comm, num_files, mpi_tag);
}

#ifdef __cplusplus
}
#endif

