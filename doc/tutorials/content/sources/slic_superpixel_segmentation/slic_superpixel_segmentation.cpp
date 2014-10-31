#include <iostream>
#include <vtkImageData.h>
#include <vtkPNGReader.h>
#include <vtkPNGWriter.h>
#include <pcl/io/pcd_io.h>
#include <pcl/segmentation/slic_superpixel_segmentation.h>
#include <pcl/visualization/pcl_visualizer.h>

typedef pcl::PointXYZRGBA PointT;
typedef typename pcl::PointCloud<PointT> PointCloud;
typedef typename PointCloud::Ptr PointCloudPtr;
typedef typename PointCloud::ConstPtr PointCloudConstPtr;
typedef pcl::Label PointL;
typedef typename pcl::PointCloud<PointL> PointCloudL;
typedef typename pcl::PointCloud<PointL>::Ptr PointCloudLPtr;
typedef typename pcl::PointCloud<PointL>::ConstPtr PointCloudLConstPtr;

int
main(int argc, char **argv)
{
  // Load pcd file
  PointCloudPtr cloud (new PointCloud);
  if (pcl::io::loadPCDFile (argv[1], *cloud) != 0)
  {
    std::cerr << "can't find pcd file" << std::endl;
    return -1;
  }
  std::cout << "PCD file loaded!" << std::endl;

  // Slic superpixel segmentation
  pcl::SLICSuperpixelSegmentation<PointT, pcl::Label> slic;
  PointCloudL labels;
  std::vector<pcl::PointIndices> label_indices;
  pcl::PointIndices boundary_indices;
  slic.setInputCloud (cloud);
  slic.setNumberOfSuperpixels (100);
  slic.setMaximumIteration (1);
  slic.setRefineSeeds (false);
  slic.setEnforceConnectivity(false);
  std::cout << "begin segment" << std::endl;
  slic.segment (labels, label_indices);
  std::cout << "end segment" << std::endl;
  std::cout << "begin find boundary" << std::endl;
  slic.findBoundary (labels, label_indices, boundary_indices);
  std::cout << "end find boundary" << std::endl;
  std::cout << "pixels count: " << label_indices.size () << std::endl;
  for (size_t i = 0; i < label_indices.size (); ++i)
  {
    std::cout << label_indices[i].indices.size () << " ";
  }
  std::cout << std::endl;

  // Write boundary to PNG file
  vtkSmartPointer<vtkImageData> image_data;
  vtkSmartPointer<vtkPNGReader> image_reader = vtkSmartPointer<vtkPNGReader>::New ();
  vtkSmartPointer<vtkPNGWriter> image_writer = vtkSmartPointer<vtkPNGWriter>::New ();
  int width, height;
  int offset, x, y;
  image_reader->SetFileName (argv[2]);
  image_data = image_reader->GetOutput ();
  image_data->Update ();
  width = image_data->GetDimensions()[0];
  height = image_data->GetDimensions()[1];
  std::cout << "width: " << width << " height: " << height << std::endl;
  for (size_t i = 0; i < boundary_indices.indices.size (); ++i)
  {
    offset = boundary_indices.indices[i];
    x = offset % width;
    y = offset / width;
    image_data->SetScalarComponentFromFloat (x, y, 0, 0, 255.0f);
  }
  image_writer->SetInput (image_data);
  image_writer->SetFileName ("a.png");
  image_writer->Write ();

  // Visualization
  pcl::visualization::PCLVisualizer viewer ("SLIC Superpixel Segmentation");
  PointCloudPtr boundary (new PointCloud);
  pcl::copyPointCloud (*cloud, boundary_indices, *boundary);
  pcl::visualization::PointCloudColorHandlerCustom<PointT> boundary_color (boundary, 255, 0, 0);
  viewer.addPointCloud (cloud, "cloud");
  viewer.addPointCloud (boundary, boundary_color, "boundary");
  viewer.setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "boundary");
  //viewer.setCameraPosition (0, 0, 0, 0, 0, 1, 0, -1, 0);
  viewer.setSize (cloud->width, cloud->height);
  viewer.resetCamera ();
  while (!viewer.wasStopped ())
  {
    viewer.spinOnce (10);
  }

  return 0;
}
