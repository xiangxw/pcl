/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2010-2012, Willow Garage, Inc.
 *  Copyright (c) 2012-, Open Perception, Inc.
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the copyright holder(s) nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef PCL_SEGMENTATION_IMPL_SLIC_SUPERPIXEL_SEGMENTATION_H_
#define PCL_SEGMENTATION_IMPL_SLIC_SUPERPIXEL_SEGMENTATION_H_

#include <pcl/segmentation/slic_superpixel_segmentation.h>
#include <pcl/console/print.h>
#include <vtkMath.h>

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT> void
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::segment (PointCloudL &labels, std::vector<pcl::PointIndices> &label_indices)
{
  std::vector<int> seeds;

  // Initialize compute
  if (!initCompute ())
  {
    PCL_ERROR ("[segment] Initialize compute failed!\n");
    return;
  }

  // Seeding
  seeding (seeds);

  // Refine seeds
  if (refine_seeds_)
  {
    refineSeeds (seeds);
  }

  // Iterative clustering
  iterativeCluster (seeds, labels, label_indices);

  // Enfore connectivity
  if (enfore_connectivity_)
  {
    enforeConnectivity (seeds, labels, label_indices);
  }

  // Deinit compute
  deinitCompute ();
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT>  bool
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::initCompute ()
{
  double tmp;

  if (!pcl::PCLBase<PointT>::initCompute ())
  {
    return (false);
  }

  // Input point cloud should be organized
  if (input_->height == 1)
  {
    PCL_ERROR ("[initCompute] Input point cloud should be organized!\n");
    return (false);
  }

  // Convert RGB to Lab color space
  labs_.resize (input_->size ());
  for (size_t i = 0; i < input_->size (); ++i)
  {
    vtkMath::RGBToLab (input_->points[i].r, input_->points[i].g, input_->points[i].b,
                       &(labs_[i].l), &(labs_[i].a), &(labs_[i].b));
  }

  // Calculate step and offset
  tmp = sqrt (static_cast<double> (input_->size ()) / static_cast<double> (num_superpixels_));
  offset_ = static_cast<int> (tmp / 2.0);
  step_ = offset_ + offset_ + 1;

  // Distance cache
  double *cache;
  int rows = step_ + step_ + 1;
  int cols = rows;
  int size = rows * cols;
  size_t index;
  size_t x;
  double sum;
  size_t count;

  lab_dist_.resize (input_->size ());
  sum = 0.0;
  count = 0;
  for (size_t i = 0; i < lab_dist_.size (); ++i)
  {
    cache = new double[size];
    lab_dist_[i] = cache;

    x = i % input_->width;
    for (int row = 0; row < rows; ++row)
    {
      index = i + (row - step_) * input_->width - x;
      for (int col = 0; col < cols; ++col)
      {
        if (index >= 0 && index < input_->size ())
        {
          tmp = calculateColorDistance (i, index);
          *cache = tmp;
          sum += tmp;
          ++count;
        }
        ++cache;
        ++index;
      }
    }
  }
  mean_lab_dist_ = sum / count;

  xyz_dist_.resize (input_->size ());
  sum = 0.0;
  count = 0;
  for (size_t i = 0; i < xyz_dist_.size (); ++i)
  {
    cache = new double[size];
    xyz_dist_[i] = cache;

    x = i % input_->width;
    for (int row = 0; row < rows; ++row)
    {
      index = i + (row - step_) * input_->width - x;
      for (int col = 0; col < cols; ++col)
      {
        if (index >= 0 && index < input_->size () && pcl::isFinite (input_->points[i]))
        {
          tmp = calculateSpatialDistance (i, index);
          *cache = tmp;
          sum += tmp;
          ++count;
        }
        ++cache;
        ++index;
      }
    }
  }
  mean_xyz_dist_ = sum / count;

  return (true);
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT> void
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::seeding (std::vector<int> &seeds)
{
  int x, y;

  // Perform hex grid seeding
  seeds.clear ();
  for (int col = 0, magic = 0; ; ++col, ++magic)
  {
    y = col * step_ + offset_;
    if (y > input_->height - 1)
    {
      break;
    }
    for (int row = 0; ; ++row)
    {
      x = row * step_ + (offset_ << (magic & 0x1)); // hex grid
      if (x > input_->width - 1)
      {
        break;
      }
      seeds.push_back (y * input_->width + x);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT> void
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::refineSeeds (std::vector<int> &seeds)
{
  int old, index;
  double min;
  int offset[8];
  double tmp;

  // 8-connected neightorhoods
  offset[0] = -input_->width - 1; // up-left
  offset[1] = -input_->width;     // up
  offset[2] = -input_->width + 1; // up-right
  offset[3] = -1;                 // left
  offset[4] = 1;                  // right
  offset[5] = input_->width - 1;  // down-left
  offset[6] = input_->width;      // down
  offset[7] = input_->width + 1;  // down-right

  // Find the point with minimum gradient in 8-connected neighborhoods
  for (std::vector<int>::iterator it = seeds.begin (); it != seeds.end (); ++it)
  {
    old = *it;
    min = calculateGradient (old);

    for (int i = 0; i < 8; ++i)
    {
      index = old + offset[i];
      if (index >= 0 && index < input_->size ())
      {
        tmp = calculateGradient (index);
        if (tmp < min)
        {
          min = tmp;
          *it = index;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT> void
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::iterativeCluster (const std::vector<int> &seeds,
                                                                    PointCloudL &labels,
                                                                    std::vector<pcl::PointIndices> &label_indices)
{
  labels.resize (input_->size ());
  label_indices.resize (seeds.size ());
  (void)seeds;
  (void)labels;
  (void)label_indices;

  for (int i = 0; i < max_iteration_; ++i)
  {

  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT> void
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::enforeConnectivity (const std::vector<int> &seeds,
                                                                      PointCloudL &labels,
                                                                      std::vector<pcl::PointIndices> &label_indices)
{
  (void)seeds;
  (void)labels;
  (void)label_indices;
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT> inline double
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::calculateGradient (int index) const
{
  int offset[4];
  double sum = 0.0;
  int count = 0;
  int tmp;

  offset[0] = -input_->width; // up
  offset[1] = input_->width;  // down
  offset[2] = -1;             // left
  offset[3] = 1;              // right

  for (int i = 0; i < 4; ++i)
  {
    tmp = index + offset[i];
    if (tmp >= 0 && tmp < input_->size () && pcl::isFinite (input_->points[tmp]))
    {
      sum += calculateDistance (index, tmp);
      ++count;
    }
  }

  if (count == 0)
  {
    return std::numeric_limits<double>::max ();
  }
  return sum / count;
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT> inline double
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::calculateDistance (int index1, int index2) const
{
  double sum;

  // Color distance
  sum = calculateColorDistance (index1, index2) / mean_lab_dist_ / mean_lab_dist_;

  // Spatial distance
  if (pcl::isFinite (input_->points[index1]) && pcl::isFinite (input_->points[index2]))
  {
    sum += calculateSpatialDistance (index1, index2) / mean_xyz_dist_ / mean_xyz_dist_;
  }
  else
  {
    sum += sum;
  }

  return sum;
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT> inline double
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::calculateColorDistance (int index1, int index2) const
{
  double sum;
  double tmp;

  tmp = labs_[index1].l - labs_[index2].l;
  sum = tmp * tmp;
  tmp = labs_[index1].a - labs_[index2].a;
  sum += tmp * tmp;
  tmp = labs_[index1].b - labs_[index2].b;
  sum += tmp * tmp;

  return sum;
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT> inline double
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::calculateSpatialDistance (int index1, int index2) const
{
  double sum;
  double tmp;

  tmp = input_->points[index1].x - input_->points[index2].x;
  sum = tmp * tmp;
  tmp = input_->points[index1].y - input_->points[index2].y;
  sum += tmp * tmp;
  tmp = input_->points[index1].z - input_->points[index2].z;
  sum += tmp * tmp;

  return sum;
}

#define PCL_INSTANTIATE_SLICSuperpixelSegmentation(T,LT) template class PCL_EXPORTS pcl::SLICSuperpixelSegmentation<T,LT>;

#endif // end PCL_SEGMENTATION_IMPL_SLIC_SUPERPIXEL_SEGMENTATION_H_
