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

#ifndef PCL_SEGMENTATION_SLIC_SUPERPIXEL_SEGMENTATION_H_
#define PCL_SEGMENTATION_SLIC_SUPERPIXEL_SEGMENTATION_H_

#include <pcl/pcl_base.h>

namespace pcl
{
  /** \brief Simple linear iterative clustering (SLIC) superpixel segmentation.
    *
    * \note The input point cloud should be organized
    * \note For more information please see
    * <b>Radhakrishna Achanta, Appu Shaji, Kevin Smith, Aurelien Lucchi, Pascal Fua, and Sabine Süsstrunk,
    * SLIC Superpixels Compared to State-of-the-art Superpixel Methods,
    * IEEE Transactions on Pattern Analysis and Machine Intelligence,
    * vol. 34, num. 11, p. 2274 - 2282, May 2012.</b>
    *
    * \author xiangxw
    */
  template <typename PointT, typename PointLT>
  class SLICSuperpixelSegmentation : public PCLBase<PointT>
  {
    using PCLBase<PointT>::input_;
    using PCLBase<PointT>::initCompute;
    using PCLBase<PointT>::deinitCompute;

    public:
      typedef typename pcl::PointCloud<PointT> PointCloud;
      typedef typename PointCloud::Ptr PointCloudPtr;
      typedef typename PointCloud::ConstPtr PointCloudConstPtr;

      typedef typename pcl::PointCloud<PointLT> PointCloudL;
      typedef typename PointCloudL::Ptr PointCloudLPtr;
      typedef typename PointCloudL::ConstPtr PointCloudLConstPtr;

      /** \brief Constructor for SLICSuperpixelSegmentation. */
      SLICSuperpixelSegmentation ()
        : labs_ (), num_superpixels_ (500), refine_seeds_ (true)
      {
      }

      /** \brief Destructor for SLICSuperpixelSegmentation. */
      virtual
      ~SLICSuperpixelSegmentation ()
      {
      }

      /** \brief Set number of superpixels to segment.
        * \param[in] superpixels number of superpixels
        */
      inline void
      setNumberOfSuperpixels (unsigned int superpixels)
      {
        num_superpixels_ = superpixels;
      }
      /** \brief Get number of superpixels. */
      inline unsigned int
      getNumberOfSuperpixels () const
      {
        return (num_superpixels_);
      }

      /** \brief Set whether or not to refine seeds.
        * \sa refine_seeds_
        */
      inline void
      setRefineSeeds (bool refine_seeds)
      {
        refine_seeds_ = refine_seeds;
      }
      /** \brief Get whether or not to refine seeds.
        * \sa refine_seeds_
        */
      inline bool
      getRefineSeeds () const
      {
        return (refine_seeds_);
      }

      /** \brief Perform SLIC superpixel segmentation.
        * \param[out] labels a PointCloud of labels: each superpixel will have a unique id
        * \param[out] label_indices a vector of PointIndices corresponding to each label
        */
      void
      segment (PointCloudL &labels, std::vector<pcl::PointIndices> &label_indices) const;

    protected:
      /** \brief CIELAB color space. */
      struct Lab
      {
        double l, a, b;
      };

      /** \brief Seeding.
        * \param[out] seeds index of seed points
        */
      void
      seeding (std::vector<int> &seeds) const;

      /** \brief Refine seeds. 
        * \param[out] seeds index of seed points
        */
      void
      refineSeeds (std::vector<int> &seeds) const;

      /** \brief Interative cluster.
        * \param[in] seeds index of seed points
        * \param[out] labels a PointCloud of labels: each superpixel will have a unique id
        * \param[out] label_indices a vector of PointIndices corresponding to each label
        */
      void
      interativeCluster (const std::vector<int> &seeds, PointCloudL &labels, std::vector<pcl::PointIndices> &label_indices) const;

      /** \brief Enfore connectivity.
        * \param[in] seeds index of seed points
        * \param[out] labels a PointCloud of labels: each superpixel will have a unique id
        * \param[out] label_indices a vector of PointIndices corresponding to each label
        */
      void
      enforeConnectivity (const std::vector<int> &seeds, PointCloudL &labels, std::vector<pcl::PointIndices> &label_indices) const;

      /** \brief Values of CIELAB color space for all pixels. */
      std::vector<CIELab> labs_;
      /** \brief Number of superpixels to segment. */
      unsigned int num_superpixels_;
      /** \brief Refine seeds or not.
        * If true, seeds will be moved to the lowest gradient position in a 3 x 3 neighborhood.
        * This is done to avoid centering a superpixel on an edge, and to reduce the chance of
        * seeding a superpixel with a noisy pixel. For more information please see the paper.
        * Default value is true.
        */
      bool refine_seeds_;
  };
}

#ifdef PCL_NO_PRECOMPILE
#include <pcl/segmentation/impl/slic_superpixel_segmentation.hpp>
#endif

#endif // end PCL_SEGMENTATION_SLIC_SUPERPIXEL_SEGMENTATION_H_
