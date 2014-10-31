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
    * \note The input point cloud should be organized.
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
        : num_superpixels_ (500), refine_seeds_ (true), max_iteration_ (10), enfore_connectivity_ (true)
        , labs_ (), seeds_ (), mean_lab_dist_ (), mean_xyz_dist_ (), step_ (), offset_ ()
      {
      }

      /** \brief Destructor for SLICSuperpixelSegmentation. */
      virtual
      ~SLICSuperpixelSegmentation ()
      {
      }

      /** \brief Set number of superpixels to segment.
        * \param[in] superpixels number of superpixels
        * \sa num_superpixels_
        */
      inline void
      setNumberOfSuperpixels (unsigned int superpixels)
      {
        num_superpixels_ = superpixels;
      }
      /** \brief Get number of superpixels.
        * \sa num_superpixels_
        */
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

      /** \brief Set maximum iteration count.
        * \sa max_iteration_
        */
      inline void
      setMaximumIteration (unsigned int max_iteration)
      {
        max_iteration_ = max_iteration;
      }
      /** \brief Get maximum iteration count.
        * \sa max_iteration_
        */
      inline unsigned int
      getMaximumIteration () const
      {
        return (max_iteration_);
      }

      /** \brief Set whether or not to enfore connectivity.
        * \sa enfore_connectivity_
        */
      inline void
      setEnforceConnectivity (bool enfore_connectivity)
      {
        enfore_connectivity_ = enfore_connectivity;
      }
      /** \brief Get whether or not to enfore connectivity.
        * \sa enfore_connectivity_
        */
      inline bool
      getEnforceConnectivity () const
      {
        return enfore_connectivity_;
      }

      /** \brief Perform SLIC superpixel segmentation.
        * \param[out] labels a PointCloud of labels. Each superpixel will have a unique id. Id of unlabeled points is 0
        * \param[out] label_indices a vector of PointIndices corresponding to each label
        */
      void
      segment (PointCloudL &labels, std::vector<pcl::PointIndices> &label_indices);

      /** \brief Find boundary.
        * \param[in] labels a PointCloud of labels. Each superpixel will have a unique id. Id of unlabeled points is 0
        * \param[in] label_indices a vector of PointIndices corresponding to each label
        * \param[out] boundary_indices boundary indices
        */
      static void
      findBoundary (const PointCloudL &labels, const std::vector<pcl::PointIndices> &label_indices, pcl::PointIndices &boundary_indices);

    protected:
      /** Initialize compute. */
      bool
      initCompute ();

    private:
      /** \brief CIELAB color space. */
      struct Lab
      {
        double l, a, b;
      };

      struct Seed
      {
        double l, a, b;
        double x, y, z; // 3d pos
        size_t xx, yy;  // 2d pos in input_->width * input_->height
        size_t index;   // index in input_
      };

      /** \brief Uniform spatial seeding. */
      void
      seeding ();

      /** \brief Refine seeds. */
      void
      refineSeeds ();

      /** \brief Iterative cluster.
        * \param[out] labels a PointCloud of labels: each superpixel will have a unique id
        * \param[out] label_indices a vector of PointIndices corresponding to each label
        */
      void
      iterativeCluster (PointCloudL &labels, std::vector<pcl::PointIndices> &label_indices);

      /** \brief Enfore connectivity.
        * 1. Find an adjacent label for each pixel at the start.
        * 2. If a certain pixel is too small, assign the previously found
        *    adjacent label to this pixel, and not increase the label.
        * \param[out] labels a PointCloud of labels: each superpixel will have a unique id
        * \param[out] label_indices a vector of PointIndices corresponding to each label
        */
      void
      enforeConnectivity (PointCloudL &labels, std::vector<pcl::PointIndices> &label_indices);

      /** \brief Calculate gradient of a point with 4-connected neighborhoods.
        * For more information please see the paper.
        * \param[in] index index of the point 
        */
      double
      calculateGradient (size_t index) const;
      /** \brief Calculate distance of two points.
        * For more information please see the paper.
        * \param[in] index1 index of point 1
        * \param[in] index2 index of point 2
        */
      double
      calculateDistance (size_t index1, size_t index2) const;
      /** \brief Calculate color distance of two points.
        * \param[in] index1 index of point 1
        * \param[in] index2 index of point 2
        */
      double
      calculateColorDistance (size_t index1, size_t index2) const;
      /** \brief Calculate spatial distance of two points.
        * \param[in] index1 index of point 1
        * \param[in] index2 index of point 2
        */
      double
      calculateSpatialDistance (size_t index1, size_t index2) const;

      /** \brief Number of superpixels to segment.
        * \note number of segmented superpixels may not be exactly equal to this value
        */
      unsigned int num_superpixels_;
      /** \brief Refine seeds or not.
        * If true, seeds will be moved to the lowest gradient position in a 3 x 3 neighborhood.
        * This is done to avoid centering a superpixel on an edge, and to reduce the chance of
        * seeding a superpixel with a noisy pixel. For more information please see the paper.
        * Default value is true.
        */
      bool refine_seeds_;
      /** \brief Maximum iteration count. Default value is 10. */
      unsigned int max_iteration_;
      /** \brief Whether or not to enfore connectivity.
        * If true, we will enforce connectivity by re-assigning disjoint pixels to nearby superpixels. For more information please see the paper.
        * Default value is true.
        */
      bool enfore_connectivity_;

      /** \brief Values of Lab color space for all pixels. */
      std::vector<Lab> labs_;
      /** \brief Seeds. */
      std::vector<Seed> seeds_;
      /** Mean lab color distance. */
      double mean_lab_dist_;
      /** Mean xyz spatial distance. */
      double mean_xyz_dist_;
      /** \brief Step size of the grid interval. Equal to :math:`sqrt (N / num_superpixels_)`. */
      int step_;
      /** \brief Offset of the grid interval. Equal to :math:`sqrt (N / num_superpixels_) / 2`. */
      int offset_;
  };
}

#ifdef PCL_NO_PRECOMPILE
#include <pcl/segmentation/impl/slic_superpixel_segmentation.hpp>
#endif

#endif // end PCL_SEGMENTATION_SLIC_SUPERPIXEL_SEGMENTATION_H_
