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
#include "io/generate_cell_face_node_connectivity.h"

#ifdef USE_MPI
#include <mpi.h>
#include "pmpio.h"
#endif

#define VTK_ENCODING "ISO-8859-1"

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

static void vtk_close_file(void* context, void* file)
{
  // In this method, we do all the writing.
  xmlTextWriterEndDocument(file);
}

static int vtk_get_num_datasets(void* context, void* file, int* num_datasets)
{
  return 1;
}

// libxml shortcuts.
static inline void start_element(xmlTextWriterPtr writer, const char* element)
{
  int status = xmlTextWriterStartElement(writer, (xmlChar*)element);
  ASSERT(status == 0);
}

static inline void write_attribute(xmlTextWriterPtr writer, const char* attr, const char* value)
{
  int status = xmlTextWriterWriteAttribute(writer, (xmlChar*)attr, (xmlChar*)value);
  ASSERT(status == 0);
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

static inline void end_element(xmlTextWriterPtr writer, const char* element)
{
  int status = xmlTextWriterEndElement(writer);
  ASSERT(status == 0);
}

// ASCI version
static void vtk_plot_write_asci_datasets(void* context, void* f, io_dataset_t** datasets, int num_datasets, int rank_in_group, int procs_per_file)
{
  ASSERT(procs_per_file == 1);

  xmlTextWriterPtr writer = (xmlTextWriterPtr)f;
  io_dataset_t* dataset = datasets[0];
  mesh_t* mesh = dataset->mesh;
  ASSERT(mesh != NULL);

  // Figure out the cell-face-node connectivity.
  int num_cells = mesh->num_cells;
  int cell_face_counts[num_cells];
  int num_faces = mesh->num_faces;
  int face_node_counts[num_faces];
  int *all_face_nodes, *all_cell_faces;
  generate_cell_face_node_connectivity(mesh, face_node_counts, 
                                       &all_face_nodes, cell_face_counts,
                                       &all_cell_faces);

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
    start_element(writer, "PPointData");
    write_attribute(writer, "Scalars", "scalars");

    for (int f = 0; f < dataset->num_fields; ++f)
    {
      if (dataset->field_centerings[f] == MESH_NODE)
      {
        start_element(writer, "DataArray");
        write_attribute(writer, "type", "Float32");
        write_attribute(writer, "Name", dataset->field_names[f]);
        if (dataset->field_num_comps[f] > 1)
        {
          write_format_attribute(writer, "NumberOfComponents", "%d", dataset->field_num_comps[f]);
        }

        // FIXME: Write values

        end_element(writer, "DataArray");
      }
    }

    end_element(writer, "PointData");
  }

  // Write cell-centered field data.
  {
    start_element(writer, "CellData");
    write_attribute(writer, "Scalars", "scalars");

    for (int f = 0; f < dataset->num_fields; ++f)
    {
      if (dataset->field_centerings[f] == MESH_CELL)
      {
        start_element(writer, "DataArray");
        write_attribute(writer, "type", "Float32");
        write_attribute(writer, "Name", dataset->field_names[f]);
        if (dataset->field_num_comps[f] > 1)
        {
          write_format_attribute(writer, "NumberOfComponents", "%d", dataset->field_num_comps[f]);
        }

        // FIXME: Write values

        end_element(writer, "DataArray");
      }
    }

    end_element(writer, "CellData");
  }

  // Write out the node positions.
  {
    start_element(writer, "Points");
    start_element(writer, "DataArray");
    write_attribute(writer, "type", "Float32");
    write_attribute(writer, "NumberOfComponents", "3");
    write_attribute(writer, "format", "ascii");

    // FIXME: Write data

    end_element(writer, "DataArray");
    end_element(writer, "Points");
  }

  // Write out the cells.
  // FIXME: This needs more study!
  {
    start_element(writer, "Cells");

#if 0
    start_element(writer, "DataArray");
    write_attribute(writer, "type", "Int32");
    write_attribute(writer, "Name", "connectivity");
    write_attribute(writer, "format", "ascii");

    // FIXME: Write data

    end_element(writer, "DataArray");

    start_element(writer, "DataArray");
    write_attribute(writer, "type", "Int32");
    write_attribute(writer, "Name", "offsets");
    write_attribute(writer, "format", "ascii");

    // FIXME: Write data

    end_element(writer, "DataArray");
#endif

    start_element(writer, "DataArray");
    write_attribute(writer, "type", "UInt8");
    write_attribute(writer, "Name", "types");
    write_attribute(writer, "format", "ascii");

    // FIXME: Write VTK_POLYHEDRON data

    end_element(writer, "DataArray");

    start_element(writer, "DataArray");
    write_attribute(writer, "type", "UInt8");
    write_attribute(writer, "Name", "faces");
    write_attribute(writer, "format", "ascii");

    // FIXME: Write faces data

    end_element(writer, "DataArray");

    start_element(writer, "DataArray");
    write_attribute(writer, "type", "UInt8");
    write_attribute(writer, "Name", "faceoffsets");
    write_attribute(writer, "format", "ascii");

    // FIXME: Write faces data

    end_element(writer, "DataArray");

    end_element(writer, "Cells");
  }

  // End the document.
  {
    xmlTextWriterEndDocument(writer);
  }

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
                      .close_file = &vtk_close_file,
                      .get_num_datasets = &vtk_get_num_datasets,
                      .write_master = &vtk_plot_write_master};
  if (binary)
    vtable.write_datasets = &vtk_plot_write_binary_datasets;
  else
    vtable.write_datasets = &vtk_plot_write_asci_datasets;
  int num_files;
  MPI_Comm_size(comm, &num_files);
  return io_interface_new(NULL, "VTK-plot", vtable, comm, num_files, mpi_tag);
}

#ifdef __cplusplus
}
#endif

