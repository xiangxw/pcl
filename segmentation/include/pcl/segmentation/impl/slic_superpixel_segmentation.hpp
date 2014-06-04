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
#include <vtkMath.h>

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename PointLT> void
pcl::SLICSuperpixelSegmentation<PointT, PointLT>::segment (PointCloudL &labels, std::vector<pcl::PointIndices> &label_indices) const
{
  Lab lab;

  // Init compute
  if (!initCompute ())
  {
    PCL_ERROR ("[segment] Input point cloud not valid!\n");
    return;
  }
  if (input_->height == 1)
  {
    PCL_ERROR ("[segment] Input point cloud should be organized!\n");
    return;
  }

  // Convert RGB to Lab color space
  labs_.clear ();
  for (size_t i = 0; i < input_->points.size (); ++i)
  {
    vtkMath::RGBToLab (input_[i].r, input_[i].g, input_[i].b, &(lab.l), &(lab.a), &(lab.b));
    labs_.push_back (lab);
  }

  // Seeding
  seeding ();

  // Deinit compute
  deinitCompute ();
}

#define PCL_INSTANTIATE_SLICSuperpixelSegmentation(T,LT) template class PCL_EXPORTS pcl::SLICSLICSuperpixelSegmentation<T,LT>;

#endif // end PCL_SEGMENTATION_IMPL_SLIC_SUPERPIXEL_SEGMENTATION_H_
