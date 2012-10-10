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
#include "core/avl_tree.h"

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

// ASCI version
static void vtk_plot_write_asci_datasets(void* context, void* f, io_dataset_t** datasets, int num_datasets, int rank_in_group, int procs_per_file)
{
  ASSERT(procs_per_file == 1);

  xmlTextWriterPtr writer = (xmlTextWriterPtr)f;
  io_dataset_t* dataset = datasets[0];
  mesh_t* mesh = dataset->mesh;
  ASSERT(mesh != NULL);

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

    start_element(writer, "DataArray");
    write_attribute(writer, "type", "UInt8");
    write_attribute(writer, "Name", "types");
    write_attribute(writer, "format", "ascii");

    // FIXME: Write VTK_POLYHEDRON data

    end_element(writer, "DataArray");

    end_element(writer, "Cells");
  }

#if 0
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

  // Clean up.
  slist_free(all_cell_faces_list);
  slist_free(all_face_nodes_list);
#endif

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

