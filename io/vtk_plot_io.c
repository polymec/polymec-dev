#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <libxml/encoding.h>
#include "libxml/xmlwriter.h"
#include "io/vtk_plot_io.h"
#include "core/edit_mesh.h"
#include "core/point.h"
#include "core/slist.h"
#include "io/generate_face_node_conn.h"
#include "io/generate_cell_node_conn.h"

#ifdef USE_MPI
#include <mpi.h>
#include "pmpio.h"
#endif

#define VTK_ENCODING "ISO-8859-1"
#define VTK_POLYHEDRON 42  // From VTK source -- make sure this is right!

#ifdef __cplusplus
extern "C" {
#endif

static void* vtk_create_file(void* context, 
                             const char* filename,
                             const char* dirname)
{

  // Initialize the library and check its version.
  LIBXML_TEST_VERSION

  // Open an XML file for writing.
  xmlTextWriterPtr file = xmlNewTextWriterFilename(filename, 0);

  return (void*)file;
}

static void* vtk_open_file(void* context, 
                           const char* filename, 
                           const char* dirname,
                           io_mode_t mode)
{
  ASSERT(false); // Shouldn't get here!
  return NULL;
}

static void vtk_close_file(void* context, void* file)
{
  // In this method, we do all the writing.
  xmlTextWriterEndDocument(file);
  xmlFreeTextWriter(file);
}

static int vtk_get_num_datasets(void* context, void* file, int* num_datasets)
{
  return 1;
}

// libxml shortcuts.
static inline void start_element(xmlTextWriterPtr writer, const char* element)
{
  int status = xmlTextWriterStartElement(writer, (xmlChar*)element);
  if (status < 0)
    arbi_error("Could not start element '%s' in VTK file.", element);
}

static inline void write_attribute(xmlTextWriterPtr writer, const char* attr, const char* value)
{
  int status = xmlTextWriterWriteAttribute(writer, (xmlChar*)attr, (xmlChar*)value);
  if (status < 0)
    arbi_error("Could not write attribute '%s' in VTK file.", attr);
}

static inline void write_format_attribute(xmlTextWriterPtr writer, const char* attr, const char* format, ...)
{
  va_list args;
  va_start(args, format);
  char str[1024]; // WARNING!
  vsnprintf(str, 1024, format, args);
  va_end(args);
  write_attribute(writer, attr, str);
}

static inline void write_string(xmlTextWriterPtr writer, const char* string)
{
  int status = xmlTextWriterWriteString(writer, (xmlChar*)string);
  if (status < 0)
    arbi_error("Could not write string in VTK file.");
}

static inline void end_element(xmlTextWriterPtr writer, const char* element)
{
  int status = xmlTextWriterEndElement(writer);
  if (status < 0)
    arbi_error("Could not end element '%s' in VTK file.", element);
  ASSERT(status >= 0);
}

// ASCI version
static void vtk_plot_write_asci_datasets(void* context, void* f, io_dataset_t** datasets, int num_datasets, int rank_in_group, int procs_per_file)
{
  ASSERT(procs_per_file == 1);

  xmlTextWriterPtr writer = (xmlTextWriterPtr)f;
  io_dataset_t* dataset = datasets[0];
  mesh_t* mesh = dataset->mesh;
  ASSERT(mesh != NULL);
  int num_cells = mesh->num_cells;
  int num_faces = mesh->num_faces;
  int num_nodes = mesh->num_nodes;

  // Figure out the face-node connectivity.
  int face_node_offsets[num_faces+1];
  int *face_nodes;
  generate_face_node_conn(mesh, &face_nodes, face_node_offsets);

  // From that, figure out the cell-node connectivity.
  int cell_node_offsets[num_cells+1];
  int *cell_nodes;
  generate_cell_node_conn(mesh, face_nodes, face_node_offsets,
                          &cell_nodes, cell_node_offsets);

  // Start the VTKFile element.
  {
    start_element(writer, "VTKFile");
    write_attribute(writer, "type", "UnstructuredGrid");
    write_attribute(writer, "version", "0.1");
    write_attribute(writer, "byte_order", "LittleEndian");
  }

  // Start the UnstructuredGrid element.
  {
    start_element(writer, "UnstructuredGrid");
  }

  // Start the Piece element.
  {
    start_element(writer, "Piece");
    write_format_attribute(writer, "NumberOfPoints", "%d", mesh->num_nodes);
    write_format_attribute(writer, "NumberOfCells", "%d", mesh->num_cells);
  }

  // Write node-centered field data.
  {
    bool has_nc_fields = false;
    char first_nc_field[1024];
    for (int f = 0; f < dataset->num_fields; ++f)
    {
      if (dataset->field_centerings[f] == MESH_NODE)
      {
        strncpy(first_nc_field, dataset->field_names[f], 1024);
        has_nc_fields = true;
      }
      break;
    }

    if (has_nc_fields)
    {
      start_element(writer, "PointData");
      write_attribute(writer, "Scalars", first_nc_field);
    }

    for (int f = 0; f < dataset->num_fields; ++f)
    {
      if (dataset->field_centerings[f] == MESH_NODE)
      {
        start_element(writer, "DataArray");
        write_attribute(writer, "type", "Float32");
        write_attribute(writer, "Name", dataset->field_names[f]);
        write_attribute(writer, "format", "ascii");
        if (dataset->field_num_comps[f] > 1)
        {
          write_format_attribute(writer, "NumberOfComponents", "%d", dataset->field_num_comps[f]);
        }

        // Write data.
        char values[16*num_nodes];
        char value[16];
        values[0] = '\0';
        for (int i = 0; i < num_cells; ++i)
        {
          snprintf(value, 16, "%g ", dataset->fields[f][i]);
          strcat(values, value);
        }
        write_string(writer, values);

        end_element(writer, "DataArray");
      }
    }

    if (has_nc_fields)
      end_element(writer, "PointData");
  }

  // Write cell-centered field data.
  {
    bool has_cc_fields = false;
    char first_cc_field[1024];
    for (int f = 0; f < dataset->num_fields; ++f)
    {
      if (dataset->field_centerings[f] == MESH_CELL)
      {
        strncpy(first_cc_field, dataset->field_names[f], 1024);
        has_cc_fields = true;
      }
      break;
    }

    if (has_cc_fields)
    {
      start_element(writer, "CellData");
      write_attribute(writer, "Scalars", first_cc_field);
    }

    for (int f = 0; f < dataset->num_fields; ++f)
    {
      if (dataset->field_centerings[f] == MESH_CELL)
      {
        start_element(writer, "DataArray");
        write_attribute(writer, "type", "Float32");
        write_attribute(writer, "Name", dataset->field_names[f]);
        write_attribute(writer, "format", "ascii");
        if (dataset->field_num_comps[f] > 1)
        {
          write_format_attribute(writer, "NumberOfComponents", "%d", dataset->field_num_comps[f]);
        }

        // Write data.
        char values[16*num_cells];
        char value[16];
        values[0] = '\0';
        for (int i = 0; i < num_cells; ++i)
        {
          snprintf(value, 16, "%g ", dataset->fields[f][i]);
          strcat(values, value);
        }
        write_string(writer, values);

        end_element(writer, "DataArray");
      }
    }

    if (has_cc_fields)
      end_element(writer, "CellData");
  }

  // Write out the node positions.
  {
    start_element(writer, "Points");
    start_element(writer, "DataArray");
    write_attribute(writer, "type", "Float32");
    write_attribute(writer, "NumberOfComponents", "3");
    write_attribute(writer, "format", "ascii");

    // Write data.
    char* positions = malloc(16*3*num_nodes*sizeof(char));
    char position[3*16];
    positions[0] = '\0';
    for (int i = 0; i < num_nodes; ++i)
    {
      snprintf(position, 3*16, "%g %g %g ", mesh->nodes[i].x, mesh->nodes[i].y, mesh->nodes[i].z);
      strcat(positions, position);
    }
    write_string(writer, positions);
    free(positions);

    end_element(writer, "DataArray");
    end_element(writer, "Points");
  }

  // Write out the cells.
  {
    start_element(writer, "Cells");

    start_element(writer, "DataArray");
    write_attribute(writer, "type", "Int32");
    write_attribute(writer, "Name", "connectivity");
    write_attribute(writer, "format", "ascii");
    char* cnodes = malloc(10*cell_node_offsets[mesh->num_cells]*sizeof(char));
    cnodes[0] = '\0';
    for (int i = 0; i < cell_node_offsets[mesh->num_cells]; ++i)
    {
      char node[10];
      snprintf(node, 10, "%d ", cell_nodes[i]);
      strcat(cnodes, node);
    }
    write_string(writer, cnodes);
    free(cnodes);
    end_element(writer, "DataArray");

    start_element(writer, "DataArray");
    write_attribute(writer, "type", "Int32");
    write_attribute(writer, "Name", "offsets");
    write_attribute(writer, "format", "ascii");
    char* coffsets = malloc(10*(mesh->num_cells+1)*sizeof(char));
    coffsets[0] = '\0';
    for (int i = 0; i <= mesh->num_cells; ++i)
    {
      char offset[10];
      snprintf(offset, 10, "%d ", cell_node_offsets[i]+1);
      strcat(coffsets, offset);
    }
    write_string(writer, coffsets);
    free(coffsets);
    end_element(writer, "DataArray");

    start_element(writer, "DataArray");
    write_attribute(writer, "type", "UInt8");
    write_attribute(writer, "Name", "types");
    write_attribute(writer, "format", "ascii");

    // Write VTK_POLYHEDRON data
    char types[10*num_cells];
    char type[10];
    snprintf(type, 10, "%d ", VTK_POLYHEDRON);
    types[0] = '\0';
    for (int i = 0; i < num_cells; ++i)
      strcat(types, type);
    write_string(writer, types);

    end_element(writer, "DataArray");

    start_element(writer, "DataArray");
    write_attribute(writer, "type", "UInt32");
    write_attribute(writer, "Name", "faces");
    write_attribute(writer, "format", "ascii");

    // Write data.
    int faceoffsets[num_cells];
    for (int c = 0; c < num_cells; ++c)
    {
      faceoffsets[c] = (c > 0) ? faceoffsets[c-1] : 0;
      faceoffsets[c] += 1 + mesh->cells[c].num_faces;
      for (int f = 0; f < mesh->cells[c].num_faces; ++f)
        faceoffsets[c] += face_node_offsets[f+1] - face_node_offsets[f];
    }
    int faces_data_len = faceoffsets[num_cells-1] - faceoffsets[0];
    char* data = malloc(16*faces_data_len*sizeof(char));
    char datum[16];
    data[0] = '\0';
    for (int c = 0; c < num_cells; ++c)
    {
      snprintf(datum, 16, "%d ", mesh->cells[c].num_faces);
      strcat(data, datum);
      for (int f = 0; f < mesh->cells[c].num_faces; ++f)
      {
        int face_id = mesh->cells[c].faces[f] - &mesh->faces[0];
        int nnodes = face_node_offsets[face_id+1] - face_node_offsets[face_id];
        snprintf(datum, 16, "%d ", nnodes);
        strcat(data, datum);
        for (int n = 0; n < nnodes; ++n)
        {
          snprintf(datum, 16, "%d ", face_nodes[face_node_offsets[face_id]+n]);
          strcat(data, datum);
        }
      }
    }
    write_string(writer, data);

    end_element(writer, "DataArray");

    start_element(writer, "DataArray");
    write_attribute(writer, "type", "UInt32");
    write_attribute(writer, "Name", "faceoffsets");
    write_attribute(writer, "format", "ascii");

    // Write data.
    data[0] = '\0';
    for (int c = 0; c < num_cells; ++c)
    {
      snprintf(datum, 16, "%d ", faceoffsets[c]);
      strcat(data, datum);
    }
    write_string(writer, data);
    free(data);

    end_element(writer, "DataArray");

    end_element(writer, "Cells");
  }

  // End the document.
  {
    xmlTextWriterEndDocument(writer);
  }

  // Clean up.
  free(cell_nodes);
  free(face_nodes);
}

static void vtk_plot_write_binary_datasets(void* context, void* f, io_dataset_t** datasets, int num_datasets, int rank_in_group, int procs_per_file)
{
  arbi_error("vtk_plot_write_binary_datasets: Not yet implemented!");
}

static void vtk_plot_write_master(void* context, void* file, const char* prefix, io_dataset_t** datasets, int num_datasets, int num_files, int procs_per_file)
{
  ASSERT(procs_per_file == 1);

  xmlTextWriterPtr writer = (xmlTextWriterPtr)file;
  io_dataset_t* dataset = datasets[0];

  // Start the VTKFile element.
  {
    start_element(writer, (char*)"VTKFile");
    write_attribute(writer, "type", "PUnstructuredGrid");
    write_attribute(writer, "version", "0.1");
    write_attribute(writer, "byte_order", "LittleEndian");
  }

  // Start the PUnstructuredGrid element.
  {
    start_element(writer, (char*)"PUnstructuredGrid");
    write_attribute(writer, "GhostLevel", "0");
  }

  // Write node-centered fields.
  {
    start_element(writer, (char*)"PPointData");
    write_attribute(writer, "Scalars", "scalars");

    for (int f = 0; f < dataset->num_fields; ++f)
    {
      if (dataset->field_centerings[f] == MESH_NODE)
      {
        start_element(writer, (char*)"PDataArray");
        write_attribute(writer, "type", "Float32");
        write_attribute(writer, "Name", dataset->field_names[f]);
        if (dataset->field_num_comps[f] > 1)
        {
          write_format_attribute(writer, "NumberOfComponents", "%d", dataset->field_num_comps[f]);
        }
        end_element(writer, (char*)"PDataArray");
      }
    }

    end_element(writer, "PPointData");
  }

  // Write cell-centered fields.
  {
    start_element(writer, (char*)"PCellData");
    write_attribute(writer, "Scalars", "scalars");

    for (int f = 0; f < dataset->num_fields; ++f)
    {
      if (dataset->field_centerings[f] == MESH_CELL)
      {
        start_element(writer, (char*)"PDataArray");
        write_attribute(writer, "type", "Float32");
        write_attribute(writer, "Name", dataset->field_names[f]);
        if (dataset->field_num_comps[f] > 1)
        {
          write_format_attribute(writer, "NumberOfComponents", "%d", dataset->field_num_comps[f]);
        }
        end_element(writer, (char*)"PDataArray");
      }
    }

    end_element(writer, "PCellData");
  }

  // We use 3D points.
  {
    start_element(writer, (char*)"PPoints");
    start_element(writer, (char*)"PDataArray");
    write_attribute(writer, "type", "Float32");
    write_attribute(writer, "NumberOfComponents", "3");
    end_element(writer, "PDataArray");
    end_element(writer, "PPoints");
  }

  // Write out the pieces.
  for (int i = 0; i < num_files; ++i)
  {
    start_element(writer, (char*)"Piece");
    write_format_attribute(writer, "Source", "%s_%d.vtu", prefix, i);
  }

  // End the document.
  {
    xmlTextWriterEndDocument(writer);
  }
}

io_interface_t* vtk_plot_io_new(MPI_Comm comm,
                                int mpi_tag,
                                bool binary)
{
  io_vtable vtable = {.create_file = &vtk_create_file,
                      .open_file = &vtk_open_file, 
                      .close_file = &vtk_close_file,
                      .get_num_datasets = &vtk_get_num_datasets,
                      .write_master = &vtk_plot_write_master};
  if (binary)
    vtable.write_datasets = &vtk_plot_write_binary_datasets;
  else
    vtable.write_datasets = &vtk_plot_write_asci_datasets;
  int num_files;
  MPI_Comm_size(comm, &num_files);
  return io_interface_new(NULL, "VTK-plot", "vtu", "pvtu", vtable, comm, num_files, mpi_tag);
}

#ifdef __cplusplus
}
#endif

