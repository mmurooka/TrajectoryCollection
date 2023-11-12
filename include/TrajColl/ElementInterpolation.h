#pragma once

#include <SpaceVecAlg/SpaceVecAlg>

namespace TrajColl
{
/** \brief Calculate the value interpolating from start to end.
    \tparam T value type
    \param start start value
    \param end end value
    \param ratio interpolation ratio
*/
template<class T>
inline T interpolate(const T & start, const T & end, double ratio)
{
  return (1 - ratio) * start + ratio * end;
}

/** \brief Calculate the Quaternion interpolating from start to end.
    \param start start value
    \param end end value
    \param ratio interpolation ratio
*/
template<>
inline Eigen::Quaterniond interpolate(const Eigen::Quaterniond & start, const Eigen::Quaterniond & end, double ratio)
{
  return start.slerp(ratio, end);
}

/** \brief Calculate the 3D matrix interpolating from start to end.
    \param start start value
    \param end end value
    \param ratio interpolation ratio
*/
template<>
inline Eigen::Matrix3d interpolate(const Eigen::Matrix3d & start, const Eigen::Matrix3d & end, double ratio)
{
  return interpolate<Eigen::Quaterniond>(Eigen::Quaterniond(start), Eigen::Quaterniond(end), ratio).toRotationMatrix();
}

/** \brief Calculate the pose interpolating from start to end.
    \param start start value
    \param end end value
    \param ratio interpolation ratio
*/
template<>
inline sva::PTransformd interpolate(const sva::PTransformd & start, const sva::PTransformd & end, double ratio)
{
  Eigen::Vector3d startPos = start.translation();
  Eigen::Vector3d endPos = end.translation();
  Eigen::Vector3d startNormal = (start.rotation().row(0).cross(Eigen::Vector3d::UnitZ())).normalized();
  Eigen::Vector3d endNormal = (end.rotation().row(0).cross(Eigen::Vector3d::UnitZ())).normalized();

  constexpr double parallelDotThre = 1e-10;
  if(std::abs(startNormal.dot(endNormal)) > 1.0 - parallelDotThre)
  {
    return sva::interpolate(start, end, ratio);
  }

  // Ref: http://www.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/tech0044.html
  Eigen::Vector3d centerPos =
      0.5
      * (startPos + endPos
         + (((startNormal - (startNormal.dot(endNormal) * endNormal)).dot(endPos - startPos)) * startNormal
            - ((endNormal - (startNormal.dot(endNormal) * startNormal)).dot(endPos - startPos)) * endNormal)
               / (1.0 - std::pow(startNormal.dot(endNormal), 2)));
  sva::PTransformd relTrans = start * sva::PTransformd(start.rotation(), centerPos).inv();
  return relTrans
         * sva::PTransformd(
             interpolate<Eigen::Matrix3d>(start.rotation().transpose(), end.rotation().transpose(), ratio).transpose(),
             centerPos);
}

/** \brief Calculate the derivative value interpolating from start to end.
    \tparam T value type
    \tparam U derivative type
    \param start start value
    \param end end value
    \param ratio interpolation ratio
    \param order derivative order
*/
template<class T, class U = T>
inline U interpolateDerivative(const T & start,
                               const T & end,
                               double, // ratio
                               int order = 1)
{
  if(order == 1)
  {
    return end - start;
  }
  else
  {
    T ret = start; // Dummy initialization for dynamic size class
    ret.setZero();
    return ret;
  }
}

/** \brief Calculate the derivative of scalar value interpolating from start to end.
    \param start start value
    \param end end value
    \param ratio interpolation ratio
    \param order derivative order
*/
template<>
inline double interpolateDerivative(const double & start,
                                    const double & end,
                                    double, // ratio
                                    int order)
{
  if(order == 1)
  {
    return end - start;
  }
  else
  {
    return 0.0;
  }
}

/** \brief Calculate the derivative of Quaternion interpolating from start to end.
    \param start start value
    \param end end value
    \param ratio interpolation ratio
    \param order derivative order
*/
template<>
inline Eigen::Vector3d interpolateDerivative(const Eigen::Quaterniond & start,
                                             const Eigen::Quaterniond & end,
                                             double, // ratio
                                             int order)
{
  if(order == 1)
  {
    Eigen::AngleAxisd aa(start.inverse() * end);
    return aa.angle() * aa.axis();
  }
  else
  {
    return Eigen::Vector3d::Zero();
  }
}

/** \brief Calculate the derivative of 3D matrix interpolating from start to end.
    \param start start value
    \param end end value
    \param ratio interpolation ratio
    \param order derivative order
*/
template<>
inline Eigen::Vector3d interpolateDerivative(const Eigen::Matrix3d & start,
                                             const Eigen::Matrix3d & end,
                                             double, // ratio
                                             int order)
{
  if(order == 1)
  {
    Eigen::AngleAxisd aa(start.transpose() * end);
    return aa.angle() * aa.axis();
  }
  else
  {
    return Eigen::Vector3d::Zero();
  }
}

/** \brief Calculate the derivative of pose interpolating from start to end.
    \param start start value
    \param end end value
    \param ratio interpolation ratio
    \param order derivative order
*/
template<>
inline sva::MotionVecd interpolateDerivative(const sva::PTransformd & start,
                                             const sva::PTransformd & end,
                                             double, // ratio
                                             int order)
{
  if(order == 1)
  {
    return sva::transformError(start, end);
  }
  else
  {
    return sva::MotionVecd::Zero();
  }
}

/** \brief Calculate the derivative of wrench interpolating from start to end.
    \param start start value
    \param end end value
    \param ratio interpolation ratio
    \param order derivative order
*/
template<>
inline sva::ForceVecd interpolateDerivative(const sva::ForceVecd & start,
                                            const sva::ForceVecd & end,
                                            double, // ratio
                                            int order)
{
  if(order == 1)
  {
    return end - start;
  }
  else
  {
    return sva::ForceVecd::Zero();
  }
}
} // namespace TrajColl
