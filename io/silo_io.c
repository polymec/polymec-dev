#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "silo.h"
#include "io/silo_io.h"
#include "core/edit_mesh.h"
#include "core/point.h"
#include "core/slist.h"
#include "io/generate_face_node_conn.h"

#ifdef USE_MPI
#include <mpi.h>
#include "pmpio.h"
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
  DBfile* file = (DBfile*)f;

  // Fiddle with the table of contents to retrieve our datasets and their 
  // metadata. All our data is stored in arrays, so we just have to figure 
  // out the names of the datasets.
  DBtoc* contents = DBGetToc(file);
  int dset = 0;
  string_slist_t* field_names = string_slist_new();
  string_slist_t* string_names = string_slist_new();
  int num_fields = 0, num_strings = 0;
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
      for (int f = 0; f < contents->narrays; ++f)
      {
        char* nametok = strstr(contents->array_names[f], name);
        if (nametok != NULL)
        {
          char* fieldtok = strstr(contents->array_names[f], "_field");
          if (fieldtok != NULL)
            string_slist_append(field_names, contents->array_names[f]);
          else
          {
            char* codetok = strstr(contents->array_names[f], "_string");
            if (codetok != NULL)
              string_slist_append(string_names, contents->array_names[f]);
          }
        }
      }
      num_fields = field_names->size;
      num_strings = string_names->size;

      io_dataset_t* dataset = io_dataset_new(name);
      datasets[dset] = dataset;
      dset++;
    }
  }
  ASSERT(dset == num_datasets);

  for (int d = 0; d < num_datasets; ++d)
  {
    io_dataset_t* dataset = datasets[dset];

    // Reconstruct the connectivity and node positions.
    char conn_name[1024];
    snprintf(conn_name, 1024, "%s_conn", io_dataset_name(dataset));
    DBcompoundarray* conn = DBGetCompoundarray(file, conn_name);
    if (conn == NULL)
    {
      DBClose(file);
      char err[1024];
      snprintf(err, 1024, "Could not find mesh connectivity in file.");
      polymec_error(err);
    }
    // Element 0 is the number of ghost cells and the number of faces in each cell.
    // Element 1 is the list of face indices for each cell.
    // Element 2 is a pair of cells for each face.
    // Element 3 is the number of edges for each face.
    // Element 4 is the list of edge indices for each face.
    // Element 5 is the pair of nodes for each edge.
    if ((conn->nelems != 6) || 
        strcmp(conn->elemnames[0], "ncell_faces") ||
        strcmp(conn->elemnames[1], "cell_faces") ||
        strcmp(conn->elemnames[2], "face_cells") ||
        strcmp(conn->elemnames[3], "nface_edges") ||
        strcmp(conn->elemnames[4], "face_edges") ||
        strcmp(conn->elemnames[5], "edge_nodes"))
    {
      DBClose(file);
      char err[1024];
      snprintf(err, 1024, "Found invalid mesh connectivity.");
      polymec_error(err);
    }

    char pos_name[1024];
    snprintf(pos_name, 1024, "%s_nodes", io_dataset_name(datasets[dset]));
    DBcompoundarray* pos = DBGetCompoundarray(file, pos_name);
    if (pos == NULL)
    {
      DBClose(file);
      char err[1024];
      snprintf(err, 1024, "Could not find mesh node positions.");
      polymec_error(err);
    }
    if ((pos->nelems != 1) || strcmp(pos->elemnames[0], "positions"))
    {
      DBClose(file);
      char err[1024];
      snprintf(err, 1024, "Found invalid mesh node positions.");
      polymec_error(err);
    }
    int* conndata = (int*)conn->values;
    double* posdata = (double*)pos->values;
    int num_cells = conn->elemlengths[0] - 1;
    int num_ghost_cells = conndata[0];
    ASSERT((conn->elemlengths[3] % 2) == 0);
    int num_faces = conn->elemlengths[3] / 2;
    ASSERT((conn->elemlengths[5] % 2) == 0);
    int num_edges = conn->elemlengths[5] / 2;
    int num_nodes = pos->elemlengths[0];

    // Build a mesh.
    mesh_t* mesh = mesh_new(num_cells, num_ghost_cells, num_faces, num_edges, num_nodes);

    // Connect everything.
    int offset = 1+num_cells;
    for (int c = 0; c < num_cells; ++c)
    {
      int nf = conndata[1+c];
      for (int f = 0; f < nf; ++f)
      {
        int fi = conndata[offset++];
        mesh_attach_face_to_cell(mesh, &mesh->faces[fi], &mesh->cells[c]);
      }
    }
    for (int f = 0; f < num_faces; ++f)
    {
      int c1 = conndata[offset++], c2 = conndata[offset++];
      mesh->faces[f].cell1 = &mesh->cells[c1];
      mesh->faces[f].cell2 = &mesh->cells[c2];
    }
    int eoffset = offset + num_faces;
    for (int f = 0; f < num_faces; ++f)
    {
      int ne = conndata[offset+f];
      for (int e = 0; e < ne; ++e)
      {
        int ei = conndata[eoffset++];
        mesh_attach_edge_to_face(mesh, &mesh->edges[ei], &mesh->faces[f]);
      }
    }
    for (int e = 0; e < num_edges; ++e)
    {
      int n1 = conndata[eoffset++], n2 = conndata[eoffset++];
      mesh->edges[e].node1 = &mesh->nodes[n1];
      mesh->edges[e].node2 = &mesh->nodes[n2];
    }

    // Read in node positions.
    for (int n = 0; n < num_nodes; ++n)
    {
      mesh->nodes[n].x = posdata[3*n];
      mesh->nodes[n].y = posdata[3*n+1];
      mesh->nodes[n].z = posdata[3*n+2];
    }
    io_dataset_put_mesh(dataset, mesh);

    // Read in fields.
    for (int f = 0; f < num_fields; ++f)
    {
      char* field_name = string_slist_pop(field_names, NULL);
      DBcompoundarray* field = DBGetCompoundarray(file, field_name);
      if (field == NULL)
      {
        DBClose(file);
        char err[1024];
        snprintf(err, 1024, "Could not find field %s", field_name);
        polymec_error(err);
      }
      if ((field->nelems != 1) || strcmp(field->elemnames[0], "data"))
      {
        DBClose(file);
        char err[1024];
        snprintf(err, 1024, "Found invalid field %s.", field_name);
        polymec_error(err);
      }

      // Find the centering and number of components of the field.
      char* cstr = strstr(field_name, "_field") - 4;
      mesh_centering_t centering;
      int num_comps;
      if (!strncmp(cstr, "cell", 4))
      {
        centering = MESH_CELL;
        num_comps = field->elemlengths[0]/mesh->num_cells;
      }
      else if (!strncmp(cstr, "face", 4))
      {
        centering = MESH_FACE;
        num_comps = field->elemlengths[0]/mesh->num_faces;
      }
      else if (!strncmp(cstr, "edge", 4))
      {
        centering = MESH_EDGE;
        num_comps = field->elemlengths[0]/mesh->num_edges;
      }
      else 
      {
        ASSERT(!strncmp(cstr, "node", 4));
        centering = MESH_NODE;
        num_comps = field->elemlengths[0]/mesh->num_nodes;
      }

      // Extract the name of the field.
      char* fname = strstr(field_name, "_field_") + 7;

      // Copy the field to the dataset.
      io_dataset_put_field(dataset, fname, field->values, num_comps, centering, true);
      DBFreeCompoundarray(field);
    }

#if 0
    // Read in strings
    for (int c = 0; c < num_strings; ++c)
    {
      char* code_name = string_slist_pop(code_names);
      DBcompoundarray* code = DBGetCompoundarray(file, code_name);
      if (code == NULL)
      {
        DBClose(file);
        char err[1024];
        snprintf(err, 1024, "Could not find source code %s", code_name);
        polymec_error(err);
      }
      if ((code->nelems != 1) || strcmp(code->elemnames[0], "code"))
      {
        DBClose(file);
        char err[1024];
        snprintf(err, 1024, "Found invalid source code %s.", code_name);
        polymec_error(err);
      }

      // Extract the name of the source code.
      char* cname = strstr(code_name, "_code_") + 6;
      dataset->code_lengths[d] = code->elemlengths[0];
      memcpy(dataset->codes[d], code->values, code->elemlengths[0]);
      dataset->code_names[d] = strdup(cname);
      DBFreeCompoundarray(code);
    }
#endif

    // Clean up.
    DBFreeCompoundarray(pos);
    DBFreeCompoundarray(conn);
    string_slist_free(string_names);
    string_slist_free(field_names);
  }
}

static void silo_write_datasets(void* context, void* f, io_dataset_t** datasets, 
                                int num_datasets, int rank_in_group, int procs_per_file)
{
  DBfile* file = (DBfile*)f;

  for (int d = 0; d < num_datasets; ++d)
  {
    io_dataset_t* dataset = datasets[d];
    const char* name = io_dataset_name(dataset);
    mesh_t* mesh = io_dataset_get_mesh(dataset);
    ASSERT(mesh != NULL);

    int num_cells = mesh->num_cells;
    int num_faces = mesh->num_faces;
    int num_edges = mesh->num_edges;
    int num_nodes = mesh->num_nodes;

    {
      // Figure out the cell-face connectivity data.
      int_slist_t* cf_conn_list = int_slist_new();
      for (int c = 0; c < num_cells; ++c)
        int_slist_append(cf_conn_list, mesh->cells[c].num_faces);
      for (int c = 0; c < num_cells; ++c)
      {
        for (int f = 0; f < mesh->cells[c].num_faces; ++f)
        {
          int face_id = mesh->cells[c].faces[f] - &mesh->faces[0];
          int_slist_append(cf_conn_list, face_id);
        }
      }
      for (int f = 0; f < mesh->num_faces; ++f)
      {
        int cell1_id = mesh->faces[f].cell1 - &mesh->cells[0];
        int cell2_id = (mesh->faces[f].cell2 != NULL) ? mesh->faces[f].cell2 - &mesh->cells[0] : -1;
        int_slist_append(cf_conn_list, cell1_id);
        int_slist_append(cf_conn_list, cell2_id);
      }

      // Figure out the face-edge connectivity data.
      int_slist_t* fe_conn_list = int_slist_new();
      for (int f = 0; f < num_faces; ++f)
        int_slist_append(fe_conn_list, mesh->faces[f].num_edges);
      for (int f = 0; f < num_faces; ++f)
      {
        for (int e = 0; e < mesh->faces[f].num_edges; ++e)
        {
          int edge_id = mesh->faces[f].edges[e] - &mesh->edges[0];
          int_slist_append(fe_conn_list, edge_id);
        }
      }

      // Assemble all the connectivity data into a mesh connectivity array.
      int cf_conn_size = cf_conn_list->size;
      int fe_conn_size = fe_conn_list->size;
      int ne_conn_size = 2*num_edges;
      int conn[1 + cf_conn_size + fe_conn_size + ne_conn_size];
      int counter = 0;
      // NOTE: The first value in conn is the number of ghost cells!
      conn[counter++] = mesh->num_ghost_cells;
      for (int i = 0; i < cf_conn_size; ++i, ++counter)
        conn[counter] = int_slist_pop(cf_conn_list, NULL);
      for (int i = 0; i < fe_conn_size; ++i, ++counter)
        conn[counter] = int_slist_pop(fe_conn_list, NULL);
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
      conn_lengths[0] = 1 + num_cells;
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
      int_slist_free(fe_conn_list);
      int_slist_free(cf_conn_list);
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
      int pos = 0, num_comps;
      char* fname;
      double *field;
      mesh_centering_t centering;
      while (io_dataset_next_field(dataset, &pos, &fname, &field, &num_comps, &centering))
      {
        char field_name[1024];
        char cstr[5];
        int len;
        if (centering == MESH_CELL)
        {
          sprintf(cstr, "cell");
          len = num_comps * num_cells;
        }
        else if (centering == MESH_FACE)
        {
          sprintf(cstr, "face");
          len = num_comps * num_faces;
        }
        else if (centering == MESH_EDGE)
        {
          sprintf(cstr, "edge");
          len = num_comps * num_edges;
        }
        else if (centering == MESH_NODE)
        {
          sprintf(cstr, "node");
          len = num_comps * num_nodes;
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

#if 0
    {
      // Write strings.
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
#endif
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
  return io_interface_new(NULL, "Silo", "silo", "silo", vtable, comm, num_files, mpi_tag);
}

// This is used to generate .silo plot files.
static void silo_plot_write_datasets(void* context, void* f, io_dataset_t** datasets, int num_datasets, int rank_in_group, int procs_per_file)
{
  DBfile* file = (DBfile*)f;
  io_dataset_t* dataset = datasets[0];
  mesh_t* mesh = io_dataset_get_mesh(dataset);
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

  // Figure out the face-node connectivity.
  int num_faces = mesh->num_faces;
  int face_node_offsets[num_faces+1];
  int *face_nodes;
  generate_face_node_conn(mesh, &face_nodes, face_node_offsets);
  int face_node_counts[num_faces];
  for (int f = 0; f < num_faces; ++f)
    face_node_counts[f] = face_node_offsets[f+1] - face_node_offsets[f];

  // Figure out cell-face connectivity.
  int num_cells = mesh->num_cells;
  int cell_face_counts[num_cells];
  int cell_faces_len = 0;
  for (int c = 0; c < num_cells; ++c)
  {
    cell_face_counts[c] = mesh->cells[c].num_faces;
    cell_faces_len += mesh->cells[c].num_faces;
  }
  int offset = 0;
  int cell_faces[cell_faces_len];
  for (int c = 0; c < num_cells; ++c)
  {
    for (int f = 0; f < mesh->cells[c].num_faces; ++f, ++offset)
    {
      int face_id = mesh->cells[c].faces[f] - &mesh->faces[0];
      cell_faces[offset] = face_id;
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
  int all_face_nodes_len = 0;
  for (int i = 0; i < num_faces; ++i)
    all_face_nodes_len += face_node_counts[i];
  DBPutPHZonelist(file, (char*)"mesh_zonelist", 
      num_faces, &face_node_counts[0], 
      face_node_offsets[num_faces], &face_nodes[0], 0, 
      num_cells, &cell_face_counts[0],
      cell_faces_len, &cell_faces[0], 
      0, 0, num_cells-1, optlist);

  // Write out the cell-centered mesh data.

  // Scalar fields.
  int pos = 0, num_comps, num_fields = 0;
  char* field_name;
  double* field;
  mesh_centering_t centering;
  while (io_dataset_next_field(dataset, &pos, &field_name, &field, &num_comps, &centering))
  {
    ASSERT(centering == MESH_CELL); // FIXME
    if (centering == MESH_CELL)
    {
      DBPutUcdvar1(file, field_name, (char*)"mesh",
          field, num_cells, 0, 0,
          DB_DOUBLE, DB_ZONECENT, optlist);
      ++num_fields;
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
    char* field_names[num_fields];
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
      int pos = 0, field_index = 0, num_comps;
      char* field_name;
      double* field;
      mesh_centering_t centering;
      while (io_dataset_next_field(dataset, &pos, &field_name, &field, &num_comps, &centering))
      {
        ASSERT(centering == MESH_CELL); // FIXME
        if (centering == MESH_CELL)
        {
          if (i == 0)
            field_names[i] = strdup(field_name);
          char var_name[1024];
          snprintf(var_name, 1024, "domain_%d/%s", i, field_name);
          var_names[field_index][i] = strdup(var_name);
          ++field_index;
        }
      }
    }

    // Write the mesh and variable data.
    DBSetDir(file, "/");
    DBPutMultimesh(file, "mesh", procs_per_file, &mesh_names[0], 
        &mesh_types[0], optlist);

    for (int f = 0; f < num_fields; ++f)
    {
      DBPutMultivar(file, field_names[f], procs_per_file, 
          &var_names[f][0], &var_types[0], optlist);
    }

    for (int i = 0; i < procs_per_file; ++i)
      free(mesh_names[i]);
    for (int f = 0; f < num_fields; ++f)
    {
      free(field_names[f]);
      for (int i = 0; i < procs_per_file; ++i)
        free(var_names[f][i]);
    }
  }

#endif

  // Clean up.
  free(face_nodes);
  DBFreeOptlist(optlist);
}

static void silo_plot_write_master(void* context, void* file, const char* prefix, io_dataset_t** datasets, int num_datasets, int num_files, int procs_per_file)
{
  io_dataset_t* dataset = datasets[0];

  char* mesh_names[num_files*procs_per_file];
  int mesh_types[num_files*procs_per_file];
  int num_fields = io_dataset_num_fields(dataset);
  int var_types[num_fields][num_files*procs_per_file];
  for (int i = 0; i < num_files*procs_per_file; ++i)
  {
    mesh_types[i] = DB_UCDMESH;
    for (int f = 0; f < num_fields; ++f)
      var_types[f][i] = DB_UCDVAR;
  }
  char* var_names[num_fields][num_files*procs_per_file];

  for (int i = 0; i < num_files; ++i)
  {
    for (int c = 0; c < procs_per_file; ++c)
    {
      // Mesh.
      char mesh_name[1024];
      snprintf(mesh_name, 1024, "%d/%s.silo:/domain_%d/mesh", i, prefix, c);
      mesh_names[i*procs_per_file+c] = strdup(mesh_name);

      // Field data.
      int pos = 0, f = 0, num_comps;
      char* field_name;
      double* field;
      mesh_centering_t centering;
      while (io_dataset_next_field(dataset, &pos, &field_name, &field, &num_comps, &centering))
      {
        char var_name[1024];
        snprintf(var_name, 1024, "%d/%s.silo:/domain_%d/%s", i, prefix, c, field_name);
        var_names[f][i] = strdup(var_name);
        ++f;
      }
    }
  }

  DBoptlist* optlist = DBMakeOptlist(10);

  // Write the multimesh and variable data, and close the file.
  {
    DBPutMultimesh(file, "mesh", num_files*procs_per_file, &mesh_names[0], 
        &mesh_types[0], optlist);
    int pos = 0, f = 0, num_comps;
    char* field_name;
    double* field;
    mesh_centering_t centering;
    while (io_dataset_next_field(dataset, &pos, &field_name, &field, &num_comps, &centering))
    {
      DBPutMultivar(file, field_name, num_files*procs_per_file, 
          var_names[f], var_types[f], optlist);
      ++f;
    }
  }
  DBClose(file);

  // Clean up.
  DBFreeOptlist(optlist);
  for (int i = 0; i < num_files*procs_per_file; ++i)
    free(mesh_names[i]);
  for (int f = 0; f < num_fields; ++f)
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
  return io_interface_new(NULL, "Silo-plot", "silo", "silo", vtable, comm, num_files, mpi_tag);
}

