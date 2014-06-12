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
  // Initialize compute
  if (!initCompute ())
  {
    PCL_ERROR ("[segment] Initialize compute failed!\n");
    return;
  }

  // Seeding
  seeding ();

  // Refine seeds
  if (refine_seeds_)
  {
    refineSeeds ();
  }

  // Iterative clustering
  iterativeCluster (labels, label_indices);

  // Enfore connectivity
  if (enfore_connectivity_)
  {
    enforeConnectivity (labels, label_indices);
  }

  // Deinit compute
  deinitCompute ();
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT>  bool
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::initCompute ()
{
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
  double tmp;
  tmp = sqrt (static_cast<double> (input_->size ()) / static_cast<double> (num_superpixels_));
  offset_ = static_cast<int> (tmp / 2.0);
  step_ = offset_ + offset_ + 1;

  // Calculate mean distance of lab color space and xyz spatial space
  int rows = step_ + step_ + 1;
  int cols = rows;
  size_t index;
  size_t x;
  double sum;
  size_t count;

  sum = 0.0;
  count = 0;
  for (size_t i = 0; i < input_->size (); ++i)
  {
    x = i % input_->width;
    // For each point, calculate distance with points within a 2*step_ * 2*step_ sized rect
    for (int row = 0; row < rows; ++row)
    {
      index = i + (row - step_) * input_->width - x;
      for (int col = 0; col < cols; ++col)
      {
        if (index >= 0 && index < input_->size ())
        {
          sum += calculateColorDistance (i, index);
          ++count;
        }
        ++index;
      }
    }
  }
  mean_lab_dist_ = sum / count;

  sum = 0.0;
  count = 0;
  for (size_t i = 0; i < input_->size (); ++i)
  {
    x = i % input_->width;
    // For each point, calculate distance with points within a 2*step_ * 2*step_ sized rect
    for (int row = 0; row < rows; ++row)
    {
      index = i + (row - step_) * input_->width - x;
      for (int col = 0; col < cols; ++col)
      {
        if (index >= 0 && index < input_->size () && pcl::isFinite (input_->points[i]))
        {
          sum += calculateSpatialDistance (i, index);
          ++count;
        }
        ++index;
      }
    }
  }
  mean_xyz_dist_ = sum / count;

  return (true);
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT> void
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::seeding ()
{
  Seed seed;
  int x, y;

  // Perform hex grid seeding
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
      seed.index = y * input_->width + x;
      seed.l = labs_[seed.index].l;
      seed.a = labs_[seed.index].a;
      seed.b = labs_[seed.index].b;
      seed.x = input_->points[seed.index].x;
      seed.y = input_->points[seed.index].y;
      seed.z = input_->points[seed.index].z;
      seeds_.push_back (seed);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT> void
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::refineSeeds ()
{
  size_t old, index;
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
  for (std::vector<int>::iterator it = seeds_.begin (); it != seeds_.end (); ++it)
  {
    old = it->index;
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
          it->index = index;
        }
      }
    }

    if (old != it->index)
    {
      it->l = labs_[it->index].l;
      it->a = labs_[it->index].a;
      it->b = labs_[it->index].b;
      it->x = input_->points[it->index].x;
      it->y = input_->points[it->index].y;
      it->z = input_->points[it->index].z;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT> void
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::iterativeCluster (PointCloudL &labels,
                                                                    std::vector<pcl::PointIndices> &label_indices)
{
  std::vector<double> dist_vec ( // Distance to the closest seed
                                input_->size (), std::numeric_limits<double>::max ());
  int rows = step_ + step_ + 1;
  int cols = rows;
  size_t index, count;
  size_t x;
  double tmp;

  labels.resize (input_->size ());
  label_indices.resize (seeds_.size ());

  for (unsigned int i = 0; i < max_iteration_; ++i)
  {
    // Find the closest seed and calculate distance to the seed
    for (size_t i = 0; i < seeds_.size (); ++i)
    {
      x = seeds_[i].index % input_->width;
      count = 0;
      // For each point, calculate distance with points within a 2*step_ * 2*step_ sized rect
      for (int row = 0; row < rows; ++row)
      {
        index = seeds_[i].index + (row - step_) * input_->width - x;
        for (int col = 0; col < cols; ++col)
        {
          if (index >= 0 && index < input_->size ())
          {
            tmp = calculateDistance (seeds_[i].index, index);
            if (dist_vec[index] > tmp)
            {
              dist_vec[index] = tmp;
            }
          }
          ++index;
          ++count;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT> void
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::enforeConnectivity (PointCloudL &labels,
                                                                      std::vector<pcl::PointIndices> &label_indices)
{
  (void)seeds;
  (void)labels;
  (void)label_indices;
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT> inline double
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::calculateGradient (size_t index) const
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
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::calculateDistance (size_t index1, size_t index2) const
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
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::calculateColorDistance (size_t index1, size_t index2) const
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
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::calculateSpatialDistance (size_t index1, size_t index2) const
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
