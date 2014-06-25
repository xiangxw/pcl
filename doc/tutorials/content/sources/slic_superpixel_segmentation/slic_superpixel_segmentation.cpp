#include <iostream>
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

  // Slic superpixel segmentation
  pcl::SLICSuperpixelSegmentation<PointT, pcl::Label> slic;
  PointCloudL labels;
  std::vector<pcl::PointIndices> label_indices;
  pcl::PointIndices boundary_indices;
  slic.setInputCloud (cloud);
  slic.setNumberOfSuperpixels (100);
  //slic.setMaximumIteration (1);
  //slic.setRefineSeeds (false);
  slic.segment (labels, label_indices);
  slic.findBoundary (labels, label_indices, boundary_indices);
  std::cout << "pixels count: " << label_indices.size () << std::endl;
  for (size_t i = 0; i < label_indices.size (); ++i)
  {
    std::cout << label_indices[i].indices.size () << " ";
  }
  std::cout << std::endl;

  // Visualization
  pcl::visualization::PCLVisualizer viewer ("SLIC Superpixel Segmentation");
  PointCloudPtr boundary (new PointCloud);
  pcl::copyPointCloud (*cloud, boundary_indices, *boundary);
  pcl::visualization::PointCloudColorHandlerCustom<PointT> boundary_color (boundary, 255, 0, 0);
  viewer.addPointCloud (cloud, "cloud");
  viewer.addPointCloud (boundary, boundary_color, "boundary");
  viewer.setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "boundary");
  viewer.setCameraPosition (0, 0, 0, 0, 0, 1, 0, -1, 0);
  viewer.setSize (cloud->height, cloud->width);
  while (!viewer.wasStopped ())
  {
    viewer.spinOnce (10);
  }

  return 0;
}
