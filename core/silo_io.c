#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "silo.h"
#include "core/silo_io.h"
#include "core/point.h"
#include "core/slist.h"
#include "core/avl_tree.h"

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
#ifdef HAVE_HDF5
  int driver = DB_HDF5;
#else
  int driver = DB_PDB;
#endif
  DBfile* file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
  if (strcmp(dirname, "/"))
    DBMkDir(file, dirname);
  DBSetDir(file, dirname);
  return (void*)file;
}

static void* silo_open_file(void* context, 
                            const char* filename, 
                            const char* dirname,
                            io_mode_t mode)
{
#ifdef HAVE_HDF5
  int driver = DB_HDF5;
#else
  int driver = DB_PDB;
#endif
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
  // Fiddle with the table of contents to retrieve our datasets and their 
  // metadata. All our data is stored in arrays, so we just have to figure 
  // out the names of the datasets.
  DBtoc* contents = DBGetToc(file);
  int N = 0;
  for (int a = 0; a < contents->narrays; ++a)
  {
    // Look for arrays named "xyz_mesh".
    char* token = strstr(contents->array_names[a], "_mesh");
    if (token != NULL) N++;
  }
  return N;
}

static void silo_read_datasets(void* context, void* f, io_dataset_t** datasets, int num_datasets)
{
#if 0
  DBfile* file = (DBfile*)f;

  // Fiddle with the table of contents to retrieve our datasets and their 
  // metadata. All our data is stored in arrays, so we just have to figure 
  // out the names of the datasets.
  DBtoc* contents = DBGetToc(file);
  int dset = 0;
  for (int a = 0; a < contents->narrays; ++a)
  {
    // Look for arrays named "xyz_mesh".
    char* token = strstr(contents->array_names[a], "_mesh");
    if (token != NULL)
    {
      // Figure out the name of the dataset.
      int len = token - contents->array_names[a];
      char name[len];
      memcpy(name, contents->array_names[a], len);

      // Find all the fields and source codes in this dataset.
      int num_fields = 0, num_codes = 0;
      for (int f = 0; f < contents->narrays; ++f)
      {
        char* nametok = strstr(contents->array_names[f], name);
        if (nametok != NULL)
        {
          char* fieldtok = strstr(contents->array_names[f], "_field");
          if (fieldtok != NULL)
            ++num_fields;
          else
          {
            char* codetok = strstr(contents->array_names[f], "_code");
            if (codetok != NULL)
              ++num_codes;
          }
        }
      }

      io_dataset_t* dataset = io_dataset_new(name, num_fields, num_codes); 
      datasets[dset] = dataset;
      dset++;
    }
  }
  ASSERT(dset == num_datasets);

  for (int d = 0; d < num_datasets; ++d)
  {

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
  }
#endif
}

static void silo_write_datasets(void* context, void* f, io_dataset_t** datasets, 
                                int num_datasets, int rank_in_group, int procs_per_file)
{
  DBfile* file = (DBfile*)f;

  for (int d = 0; d < num_datasets; ++d)
  {
    io_dataset_t* dataset = datasets[d];
    char* name = dataset->name;
    mesh_t* mesh = dataset->mesh;
    ASSERT(mesh != NULL);

    int num_cells = mesh->num_cells;
    int num_faces = mesh->num_faces;
    int num_edges = mesh->num_edges;
    int num_nodes = mesh->num_nodes;

    {
      // Figure out the cell-face connectivity data.
      slist_t* cf_conn_list = slist_new(NULL);
      for (int c = 0; c < num_cells; ++c)
        slist_append(cf_conn_list, (void*)mesh->cells[c].num_faces);
      for (int c = 0; c < num_cells; ++c)
      {
        for (int f = 0; f < mesh->cells[c].num_faces; ++f)
        {
          int face_id = mesh->cells[c].faces[f] - &mesh->faces[0];
          slist_append(cf_conn_list, (void*)face_id);
        }
      }
      for (int f = 0; f < mesh->num_faces; ++f)
      {
        int cell1_id = mesh->faces[f].cell1 - &mesh->cells[0];
        int cell2_id = (mesh->faces[f].cell2 != NULL) ? mesh->faces[f].cell2 - &mesh->cells[0] : -1;
        slist_append(cf_conn_list, (void*)cell1_id);
        slist_append(cf_conn_list, (void*)cell2_id);
      }

      // Figure out the face-edge connectivity data.
      slist_t* fe_conn_list = slist_new(NULL);
      for (int f = 0; f < num_faces; ++f)
        slist_append(fe_conn_list, (void*)mesh->faces[f].num_edges);
      for (int f = 0; f < num_faces; ++f)
      {
        for (int e = 0; e < mesh->faces[f].num_edges; ++e)
        {
          int edge_id = mesh->faces[f].edges[e] - &mesh->edges[0];
          slist_append(fe_conn_list, (void*)edge_id);
        }
      }

      // Assemble all the connectivity data into a mesh connectivity array.
      int cf_conn_size = slist_size(cf_conn_list);
      int fe_conn_size = slist_size(fe_conn_list);
      int ne_conn_size = 2*num_edges;
      int conn[1 + cf_conn_size + fe_conn_size + ne_conn_size];
      int counter = 0;
      // NOTE: The first value in conn is the number of ghost cells!
      conn[counter++] = mesh->num_ghost_cells;
      for (int i = 0; i < cf_conn_size; ++i, ++counter)
        conn[counter] = (int)slist_pop(cf_conn_list);
      for (int i = 0; i < fe_conn_size; ++i, ++counter)
        conn[counter] = (int)slist_pop(fe_conn_list);
      for (int e = 0; e < num_edges; ++e)
      {
        int node1_id = mesh->edges[e].node1 - &mesh->nodes[0];
        int node2_id = mesh->edges[e].node2 - &mesh->nodes[0];
        conn[counter++] = node1_id;
        conn[counter++] = node2_id;
      }

      char conn_name[1024];
      snprintf(conn_name, 1024, "%s_conn", name);
      int conn_lengths[6];
      char* conn_names[6];
      conn_names[0] = strdup("ncell_faces");
      conn_lengths[0] = num_cells;
      conn_names[2] = strdup("face_cells");
      conn_lengths[2] = cf_conn_size - 2*mesh->num_faces;
      conn_names[1] = strdup("cell_faces");
      conn_lengths[1] = cf_conn_size - conn_lengths[2] - conn_lengths[0];
      conn_names[3] = strdup("nface_edges");
      conn_lengths[3] = num_faces;
      conn_names[4] = strdup("face_edges");
      conn_lengths[4] = fe_conn_size - num_faces;
      conn_names[5] = strdup("edge_nodes");
      conn_lengths[5] = ne_conn_size;

      // Write it out.
      DBPutCompoundarray(file, conn_name, conn_names, conn_lengths, 6, 
          (void*)&conn[0], cf_conn_size + fe_conn_size, DB_INT, 0);

      // Clean up.
      slist_free(fe_conn_list);
      slist_free(cf_conn_list);
      free(conn_names[0]);
      free(conn_names[1]);
      free(conn_names[2]);
      free(conn_names[3]);
      free(conn_names[4]);
      free(conn_names[5]);
    }

    {
      // Now write out node positions.
      double nodes[3*num_nodes];
      for (int n = 0; n < num_nodes; ++n)
      {
        nodes[3*n]   = mesh->nodes[n].x;
        nodes[3*n+1] = mesh->nodes[n].y;
        nodes[3*n+2] = mesh->nodes[n].z;
      }
      char nodes_name[1024];
      snprintf(nodes_name, 1024, "%s_nodes", name);
      char* nodes_names[1];
      int nodes_lengths[1];
      nodes_names[0] = strdup("positions");
      nodes_lengths[0] = num_nodes;
      DBPutCompoundarray(file, nodes_name, nodes_names, nodes_lengths, 1, 
          (void*)&nodes[0], 3*num_nodes, DB_DOUBLE, 0);

      // Clean up.
      free(nodes_names[0]);
    }

    {
      // Write fields.
      for (int f = 0; f < dataset->num_fields; ++f)
      {
        char* fname = dataset->field_names[f];
        double* field = dataset->fields[f];
        mesh_centering_t centering = dataset->field_centerings[f];
        int  num_comp = dataset->field_num_comps[f];

        char field_name[1024];
        char cstr[5];
        int len;
        if (centering == MESH_CELL)
        {
          sprintf(cstr, "cell");
          len = num_comp * num_cells;
        }
        else if (centering == MESH_FACE)
        {
          sprintf(cstr, "face");
          len = num_comp * num_faces;
        }
        else if (centering == MESH_EDGE)
        {
          sprintf(cstr, "edge");
          len = num_comp * num_edges;
        }
        else if (centering == MESH_NODE)
        {
          sprintf(cstr, "node");
          len = num_comp * num_nodes;
        }
        snprintf(field_name, 1024, "%s_%s_field_%s", name, cstr, fname);
        char* f_names[1];
        f_names[0] = strdup("data");
        int f_lengths[1];
        f_lengths[0] = len;

        DBPutCompoundarray(file, field_name, f_names, f_lengths, 1, 
            (void*)&field[0], len, DB_DOUBLE, 0);

        // Clean up.
        free(f_names[0]);
      }
    }

    {
      // Write source codes.
      for (int c = 0; c < dataset->num_codes; ++c)
      {
        char* cname = dataset->code_names[c];
        char* code  = dataset->codes[c];
        int len     = dataset->code_lengths[c];

        char code_name[1024];
        snprintf(code_name, 1024, "%s_code_%s", name, cname);
        char* c_names[1];
        c_names[0] = strdup("code");
        int c_lengths[1];
        c_lengths[0] = len;

        DBPutCompoundarray(file, code_name, c_names, c_lengths, 1, 
            (void*)&code[0], len, DB_CHAR, 0);

        // Clean up.
        free(c_names[0]);
      }
    }
  }
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

// Traverses the given points of a polygonal facet along their convex
// hull, writing their indices to indices in order.
static void traverse_convex_hull(double* points, int num_points, int* indices, int* count)
{
  *count = 0;

  // Find the "lowest" point in the set.
  double ymin = FLT_MAX;
  int index0 = -1;
  for (int p = 0; p < num_points; ++p)
  {
    if (ymin > points[2*p+1])
    {
      ymin = points[2*p+1];
      index0 = p;
    }
  }

  // We start with this point and a horizontal angle.
  double theta_prev = 0.0;
  indices[(*count)++] = index0;

  // Now start gift wrapping.
  int i = index0;
  do 
  {
    double dtheta_min = 2.0*M_PI;
    int j_min = -1;
    for (int j = 0; j < num_points; ++j)
    {
      if (j != i)
      {
        double dx = points[2*j] - points[2*i],
               dy = points[2*j+1] - points[2*i+1];
        double theta = atan2(dy, dx);
        double dtheta = theta - theta_prev;
        if (dtheta < 0.0)
          dtheta += 2.0*M_PI;
        if (dtheta_min > dtheta)
        {
          dtheta_min = dtheta;
          j_min = j;
        }
      }
    }
    if (j_min != index0)
      indices[(*count)++] = j_min;
    theta_prev += dtheta_min;
    i = j_min;
  }
  while (i != index0);

  // The convex hull should be a polygon unless the input points 
  // don't form a polygon.
  ASSERT((num_points <= 2) || 
         ((num_points > 2) && (*count > 2)));
}

// This is used to generate .silo plot files.
static void silo_plot_write_datasets(void* context, void* f, io_dataset_t** datasets, int num_datasets, int rank_in_group, int procs_per_file)
{
  DBfile* file = (DBfile*)f;
  io_dataset_t* dataset = datasets[0];
  mesh_t* mesh = dataset->mesh;
  ASSERT(mesh != NULL);

  // This is optional for now, but we'll give it anyway.
  char *coordnames[3];
  coordnames[0] = (char*)"xcoords";
  coordnames[1] = (char*)"ycoords";
  coordnames[2] = (char*)"zcoords";

  // Node coordinates.
  int num_nodes = mesh->num_nodes;
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

  // Make a list of all nodes attached to faces (unordered).
  int num_cells = mesh->num_cells;
  int num_faces = mesh->num_faces;
  int face_node_counts[num_faces];
  int* face_nodes[num_faces];
  for (int f = 0; f < num_faces; ++f)
  {
    int rays = 0, ne = mesh->faces[f].num_edges;
    for (int e = 0; e < ne; ++e)
    {
      if (mesh->faces[f].edges[e]->node2 == NULL)
        ++rays;
    }
    face_node_counts[f] = ne - rays;
    ASSERT(face_node_counts[f] >= 3);
    face_nodes[f] = malloc(face_node_counts[f]*sizeof(int));
  }

  avl_tree_t* fnodes = int_avl_tree_new();
  for (int f = 0; f < num_faces; ++f)
  {
    int counter = 0, ne = mesh->faces[f].num_edges;
    for (int e = 0; e < ne; ++e)
    {
      edge_t* edge = mesh->faces[f].edges[e];
      int node1_id = edge->node1 - &mesh->nodes[0];
      if (avl_tree_find(fnodes, (void*)node1_id) == NULL)
      {
        face_nodes[f][counter++] = node1_id;
        avl_tree_insert(fnodes, (void*)node1_id);
      }
      if (edge->node2 != NULL)
      {
        int node2_id = edge->node2 - &mesh->nodes[0];
        if (avl_tree_find(fnodes, (void*)node2_id) == NULL)
        {
          face_nodes[f][counter++] = node2_id;
          avl_tree_insert(fnodes, (void*)node2_id);
        }
      }
    }
    ASSERT(counter == face_node_counts[f]);
    avl_tree_clear(fnodes);
  }
  avl_tree_free(fnodes);

  // Compute cell centers from face nodes.
  point_t cell_centers[num_cells];
  memset(cell_centers, 0, num_cells*sizeof(point_t));
  avl_tree_t* cell_nodes = int_avl_tree_new();
  for (int c = 0; c < num_cells; ++c)
  {
    int num_nodes = 0;
    for (int f = 0; f < mesh->cells[c].num_faces; ++f)
    {
      for (int n = 0; n < face_node_counts[f]; ++n)
      {
        int node_id = face_nodes[f][n];
        if (avl_tree_find(cell_nodes, (void*)node_id) == NULL)
        {
          avl_tree_insert(cell_nodes, (void*)node_id);
          node_t* node = &mesh->nodes[face_nodes[f][n]];
          cell_centers[c].x += node->x;
          cell_centers[c].y += node->y;
          cell_centers[c].z += node->z;
          ++num_nodes;
        }
      }
    }
    cell_centers[c].x /= num_nodes;
    cell_centers[c].y /= num_nodes;
    cell_centers[c].z /= num_nodes;
    avl_tree_clear(cell_nodes);
  }
  avl_tree_free(cell_nodes);

  slist_t* all_face_nodes_list = slist_new(NULL);
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    // Compute the normal vector for the face, pointing outward from 
    // its first cell.
    int nn = face_node_counts[f];
    int* nodes = face_nodes[f];
    ASSERT(nn >= 3);
    point_t face_center = {.x = 0.0, .y = 0.0, .z = 0.0};
    for (int n = 0; n < nn; ++n)
    {
      node_t* node = &mesh->nodes[nodes[n]];
      face_center.x += node->x;
      face_center.y += node->y;
      face_center.z += node->z;
    }
    face_center.x /= nn;
    face_center.y /= nn;
    face_center.z /= nn;

    // Construct vectors v1, v2, and v3, where v1 is the vector pointing from the 
    // face center to the first face node, v2 is a vector pointing from the face 
    // center to any other face, node, and v3 is their cross product.
    vector_t v1;
    v1.x = mesh->nodes[nodes[0]].x - face_center.x;
    v1.y = mesh->nodes[nodes[0]].y - face_center.y;
    v1.z = mesh->nodes[nodes[0]].z - face_center.z;

    vector_t v2, normal;
    double normal_mag;
    for (int n = 1; n < nn; ++n)
    {
      v2.x = mesh->nodes[nodes[n]].x - face_center.x;
      v2.y = mesh->nodes[nodes[n]].y - face_center.y;
      v2.z = mesh->nodes[nodes[n]].z - face_center.z;

      // normal = v1 x v2.
      vector_cross(v1, v2, &normal);
      normal_mag = sqrt(vector_dot(normal, normal));
      if (normal_mag > 1e-14) break;
    }
    ASSERT(normal_mag > 1e-14);
    normal.x /= normal_mag; normal.y /= normal_mag; normal.z /= normal_mag;

    vector_t v3;
    point_t cell_center;
    int cell1 = mesh->faces[f].cell1 - &mesh->cells[0];
    cell_center.x = cell_centers[cell1].x;
    cell_center.y = cell_centers[cell1].y;
    cell_center.z = cell_centers[cell1].z;
    v3.x = face_center.x - cell_center.x;
    v3.y = face_center.y - cell_center.y;
    v3.z = face_center.z - cell_center.z;
    if (vector_dot(normal, v3) < 0.0)
    {
      normal.x *= -1.0; normal.y *= -1.0; normal.z *= -1.0;
    }

    // Now project the coordinates of the face's nodes to the plane
    // with the given normal and centered about the face center.
    double points[2*nn]; // NOTE: planar coordinates (2D)
    vector_t e1, e2; // Basis vectors in the plane.
    double v1_mag = sqrt(vector_dot(v1, v1));
    e1.x = v1.x / v1_mag;
    e1.y = v1.y / v1_mag;
    e1.z = v1.z / v1_mag;

    // e2 = normal x e1.
    vector_cross(normal, e1, &e2);
    for (int p = 0; p < nn; ++p)
    {
      // v = node center - cell center.
      vector_t v;
      v.x = mesh->nodes[nodes[p]].x - cell_center.x;
      v.y = mesh->nodes[nodes[p]].y - cell_center.y;
      v.z = mesh->nodes[nodes[p]].z - cell_center.z;

      // Compute the perpendicular component of the point
      // with location v:
      // v_perp = v - (n o v)n.
      vector_t v_perp;
      double nov = vector_dot(normal, v);
      v_perp.x = v.x - nov * normal.x;
      v_perp.y = v.y - nov * normal.y;
      v_perp.z = v.z - nov * normal.z;

      // Project it to the plane.
      points[2*p]   = vector_dot(v_perp, e1);
      points[2*p+1] = vector_dot(v_perp, e2);
    }

    // Find the node order by traversing the convex hull of 
    // the points within the plane, appending them to all_face_nodes.
    int indices[nn], count;
    traverse_convex_hull(points, nn, indices, &count);
    face_node_counts[f] = nn;
    for (int n = 0; n < count; ++n)
      slist_append(all_face_nodes_list, (void*)face_nodes[f][indices[n]]);
  }

  // Figure out cell-face connectivity.
  int cell_face_counts[num_cells];
  slist_t* all_cell_faces_list = slist_new(NULL);
  for (int c = 0; c < num_cells; ++c)
  {
    cell_face_counts[c] = mesh->cells[c].num_faces;
    for (int f = 0; f < mesh->cells[c].num_faces; ++f)
    {
      int face_id = mesh->cells[c].faces[f] - &mesh->faces[0];
      slist_append(all_cell_faces_list, (void*)face_id);
    }
  }

  // The polyhedral zone list is referred to in the options list.
  DBoptlist* optlist = DBMakeOptlist(10);
  DBAddOption(optlist, DBOPT_PHZONELIST, (char*)"mesh_zonelist");

  // Write out the 3D polyhedral mesh.
  DBPutUcdmesh(file, (char*)"mesh", 3, coordnames, coords,
      num_nodes, num_cells, 0, 0,
      DB_DOUBLE, optlist); 

  // Write the connectivity information.
  int all_face_nodes_len = slist_size(all_face_nodes_list);
  int all_face_nodes[all_face_nodes_len];
  for (int i = 0; i < all_face_nodes_len; ++i)
    all_face_nodes[i] = (int)slist_pop(all_face_nodes_list);
  int all_cell_faces_len = slist_size(all_cell_faces_list);
  int all_cell_faces[all_cell_faces_len];
  for (int i = 0; i < all_cell_faces_len; ++i)
    all_cell_faces[i] = (int)slist_pop(all_cell_faces_list);
  DBPutPHZonelist(file, (char*)"mesh_zonelist", 
      num_faces, &face_node_counts[0], 
      all_face_nodes_len, &all_face_nodes[0], 0, 
      num_cells, &cell_face_counts[0],
      all_cell_faces_len, &all_cell_faces[0], 
      0, 0, num_cells-1, optlist);

  // Write out the cell-centered mesh data.

  // Scalar fields.
  for (int i = 0; i < dataset->num_fields; ++i)
  {
    ASSERT(dataset->field_centerings[i] == MESH_CELL); // FIXME
    if (dataset->field_centerings[i] == MESH_CELL)
    {
      DBPutUcdvar1(file, dataset->field_names[i], (char*)"mesh",
          dataset->fields[i], num_cells, 0, 0,
          DB_DOUBLE, DB_ZONECENT, optlist);
    }
  }

#ifdef HAVE_MPI
  // Write the multi-block objects to the file if needed.
  if (rank_in_group == 0)
  {
    char* mesh_names[procs_per_file];
    int mesh_types[procs_per_file];
    int var_types[procs_per_file];
    for (int i = 0; i < procs_per_file; ++i)
    {
      mesh_types[i] = DB_UCDMESH;
      var_types[i] = DB_UCDVAR;
    }
    char** var_names[num_fields];
    for (int f = 0; f < num_fields; ++f)
      var_names[f] = malloc(sizeof(char*)*procs_per_file);
    for (int i = 0; i < procs_per_file; ++i)
    {
      // Mesh.
      char mesh_name[1024];
      snprintf(mesh_name, 1024, "domain_%d/mesh", i);
      mesh_names[i] = strdup(mesh_name);

      // Field data.
      int field_index = 0;
      for (int f = 0; f < num_fields; ++f)
      {
        char var_name[1024];
        snprintf(var_name, 1024, "domain_%d/%s", i, dataset->field_names[f]);
        var_names[f][i] = strdup(var_name);
      }
    }

    // Write the mesh and variable data.
    DBSetDir(file, "/");
    DBPutMultimesh(file, "mesh", procs_per_file, &mesh_names[0], 
        &mesh_types[0], optlist);

    for (int f = 0; f < num_fields; ++f)
    {
      DBPutMultivar(file, dataset->field_names[f], procs_per_file, 
          &var_names[f][0], &var_types[0], optlist);
    }

    for (int i = 0; i < procs_per_file; ++i)
      free(mesh_names[i]);
    for (int f = 0; f < var_names.size(); ++f)
      for (int i = 0; i < procs_per_file; ++i)
        free(var_names[f][i]);
  }

#endif

  // Clean up.
  slist_free(all_cell_faces_list);
  slist_free(all_face_nodes_list);
  DBFreeOptlist(optlist);
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

