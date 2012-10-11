#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "silo.h"
#include "io/silo_io.h"
#include "core/edit_mesh.h"
#include "core/point.h"
#include "core/slist.h"
#include "io/generate_cell_face_node_connectivity.h"

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
  DBfile* file = (DBfile*)f;

  // Fiddle with the table of contents to retrieve our datasets and their 
  // metadata. All our data is stored in arrays, so we just have to figure 
  // out the names of the datasets.
  DBtoc* contents = DBGetToc(file);
  int dset = 0;
  slist_t* field_names = slist_new(NULL);
  slist_t* code_names = slist_new(NULL);
  int num_fields = 0, num_codes = 0;
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
            slist_append(field_names, contents->array_names[f]);
          else
          {
            char* codetok = strstr(contents->array_names[f], "_code");
            if (codetok != NULL)
              slist_append(code_names, contents->array_names[f]);
          }
        }
      }
      num_fields = slist_size(field_names);
      num_codes = slist_size(code_names);

      io_dataset_t* dataset = io_dataset_new(name, num_fields, num_codes); 
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
    snprintf(conn_name, 1024, "%s_conn", dataset->name);
    DBcompoundarray* conn = DBGetCompoundarray(file, conn_name);
    if (conn == NULL)
    {
      DBClose(file);
      char err[1024];
      snprintf(err, 1024, "Could not find mesh connectivity in file.");
      arbi_error(err);
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
      arbi_error(err);
    }

    char pos_name[1024];
    snprintf(pos_name, 1024, "%s_nodes", datasets[dset]->name);
    DBcompoundarray* pos = DBGetCompoundarray(file, pos_name);
    if (pos == NULL)
    {
      DBClose(file);
      char err[1024];
      snprintf(err, 1024, "Could not find mesh node positions.");
      arbi_error(err);
    }
    if ((pos->nelems != 1) || strcmp(pos->elemnames[0], "positions"))
    {
      DBClose(file);
      char err[1024];
      snprintf(err, 1024, "Found invalid mesh node positions.");
      arbi_error(err);
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
        mesh_add_face_to_cell(mesh, &mesh->faces[fi], &mesh->cells[c]);
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
        mesh_add_edge_to_face(mesh, &mesh->edges[ei], &mesh->faces[f]);
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
    dataset->mesh = mesh;

    // Read in fields.
    for (int f = 0; f < num_fields; ++f)
    {
      char* field_name = (char*)slist_pop(field_names);
      DBcompoundarray* field = DBGetCompoundarray(file, field_name);
      if (field == NULL)
      {
        DBClose(file);
        char err[1024];
        snprintf(err, 1024, "Could not find field %s", field_name);
        arbi_error(err);
      }
      if ((field->nelems != 1) || strcmp(field->elemnames[0], "data"))
      {
        DBClose(file);
        char err[1024];
        snprintf(err, 1024, "Found invalid field %s.", field_name);
        arbi_error(err);
      }

      // Find the centering and number of components of the field.
      char* cstr = strstr(field_name, "_field") - 4;
      int num_comps;
      if (!strncmp(cstr, "cell", 4))
      {
        dataset->field_centerings[d] = MESH_CELL;
        num_comps = field->elemlengths[0]/mesh->num_cells;
      }
      else if (!strncmp(cstr, "face", 4))
      {
        dataset->field_centerings[d] = MESH_FACE;
        num_comps = field->elemlengths[0]/mesh->num_faces;
      }
      else if (!strncmp(cstr, "edge", 4))
      {
        dataset->field_centerings[d] = MESH_EDGE;
        num_comps = field->elemlengths[0]/mesh->num_edges;
      }
      else 
      {
        ASSERT(!strncmp(cstr, "node", 4));
        dataset->field_centerings[d] = MESH_NODE;
        num_comps = field->elemlengths[0]/mesh->num_nodes;
      }
      dataset->field_num_comps[d] = num_comps;

      // Extract the name of the field.
      char* fname = strstr(field_name, "_field_") + 7;
      memcpy(dataset->fields[d], field->values, field->elemlengths[0]);
      dataset->field_names[d] = strdup(fname);
      DBFreeCompoundarray(field);
    }

    // Read in source codes.
    for (int c = 0; c < num_codes; ++c)
    {
      char* code_name = (char*)slist_pop(code_names);
      DBcompoundarray* code = DBGetCompoundarray(file, code_name);
      if (code == NULL)
      {
        DBClose(file);
        char err[1024];
        snprintf(err, 1024, "Could not find source code %s", code_name);
        arbi_error(err);
      }
      if ((code->nelems != 1) || strcmp(code->elemnames[0], "code"))
      {
        DBClose(file);
        char err[1024];
        snprintf(err, 1024, "Found invalid source code %s.", code_name);
        arbi_error(err);
      }

      // Extract the name of the source code.
      char* cname = strstr(code_name, "_code_") + 6;
      dataset->code_lengths[d] = code->elemlengths[0];
      memcpy(dataset->codes[d], code->values, code->elemlengths[0]);
      dataset->code_names[d] = strdup(cname);
      DBFreeCompoundarray(code);
    }

    // Clean up.
    DBFreeCompoundarray(pos);
    DBFreeCompoundarray(conn);
    slist_free(code_names);
    slist_free(field_names);
  }
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

  // Figure out the cell-face-node connectivity.
  int num_cells = mesh->num_cells;
  int cell_face_counts[num_cells];
  int num_faces = mesh->num_faces;
  int face_node_counts[num_faces];
  int *all_face_nodes, *all_cell_faces;
  generate_cell_face_node_connectivity(mesh, face_node_counts, 
                                       &all_face_nodes, cell_face_counts,
                                       &all_cell_faces);

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
  int all_cell_faces_len = 0;
  for (int i = 0; i < num_cells; ++i)
    all_cell_faces_len += cell_face_counts[i];
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
  free(all_cell_faces);
  free(all_face_nodes);
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

