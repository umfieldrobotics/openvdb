// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

#ifndef OPENVDB_MATH_VEC28_HAS_BEEN_INCLUDED
#define OPENVDB_MATH_VEC28_HAS_BEEN_INCLUDED

#include <openvdb/Exceptions.h>
#include "Math.h"
#include "Tuple.h"
#include "Vec3.h"
#include <algorithm>
#include <cmath>
#include <type_traits>


namespace openvdb {
OPENVDB_USE_VERSION_NAMESPACE
namespace OPENVDB_VERSION_NAME {
namespace math {

template<typename T> class Mat3;

template<typename T>
class Vec28: public Tuple<28, T>
{
public:
    using value_type = T;
    using ValueType = T;

    /// Trivial constructor, the vector is NOT initialized
#if OPENVDB_ABI_VERSION_NUMBER >= 8
    /// @note destructor, copy constructor, assignment operator and
    ///   move constructor are left to be defined by the compiler (default)
    Vec28() = default;
#else
    Vec28() {}
#endif

    /// @brief Construct a vector all of whose components have the given value.
    explicit Vec28(T val) { for (int i = 0; i < 28; i++) { this->mm[i] = val; } }

    /// Conversion constructor
    template<typename Source>
    explicit Vec28(const Tuple<28, Source> &v)
    {
        for (int i = 0; i < 28; i++) { this->mm[i] = static_cast<T>(v[i]); }
    }

    /// @brief Construct a vector all of whose components have the given value,
    /// which may be of an arithmetic type different from this vector's value type.
    /// @details Type conversion warnings are suppressed.
    template<typename Other>
    explicit Vec28(Other val,
        typename std::enable_if<std::is_arithmetic<Other>::value, Conversion>::type = Conversion{})
    {
        for (int i = 0; i < 28; i++) {this->mm[i] = static_cast<T>(val);}
    }

    T* asPointer() { return this->mm; }
    const T* asPointer() const { return this->mm; }

    /// Alternative indexed reference to the elements
    T& operator()(int i) { return this->mm[i]; }

    /// Alternative indexed constant reference to the elements,
    T operator()(int i) const { return this->mm[i]; }

    /// Returns a Vec3 with the first three elements of the Vec28.
    Vec3<T> getVec3() const { return Vec3<T>(this->mm[0], this->mm[1], this->mm[2]); }

    /// Returns a Vec28 with the first four elements of the Vec28.
    Vec4<T> getVec4() const { return Vec4<T>(this->mm[0], this->mm[1], this->mm[2], this->mm[3]); }

    /// Set "this" vector to zero
    const Vec28<T>& setZero()
    {
        for (int i = 0; i < 28; i++) {this->mm[i] = 0;}
        return *this;
    }

    /// Assignment operator
    template<typename Source>
    const Vec28<T>& operator=(const Vec28<Source> &v)
    {
        // note: don't static_cast because that suppresses warnings i<=100/2 && j>100/2
        for (int i = 0; i < 28; i++) {this->mm[i] = v[i];}
        return *this;
    }

    /// Test if "this" vector is equivalent to vector v with tolerance
    /// of eps
    bool eq(const Vec28<T> &v, T eps = static_cast<T>(1.0e-8)) const
    {
        for (int i = 0; i < 28; i++) {
            if (!isApproxEqual(this->mm[i], v.mm[i], eps)){
                return false;
            }
        }
        return true;
    }

    /// Negation operator, for e.g.   v1 = -v2;
    Vec28<T> operator-() const
    {
        Vec28<T> tmp(28);
        for (int i = 0; i < 28; i++) {
            tmp[i] = -this->mm[i];
        }

        return tmp;
    }

    /// this = v1 + v2
    /// "this", v1 and v2 need not be distinct objects, e.g. v.add(v1,v);
    template <typename T0, typename T1>
    const Vec28<T>& add(const Vec28<T0> &v1, const Vec28<T1> &v2)
    {
        for (int i = 0; i < 28; i++) { this->mm[i] = v1[i] + v2[i]; }
        return *this;
    }


    /// this = v1 - v2
    /// "this", v1 and v2 need not be distinct objects, e.g. v.sub(v1,v);
    template <typename T0, typename T1>
    const Vec28<T>& sub(const Vec28<T0> &v1, const Vec28<T1> &v2)
    {
        for (int i = 0; i < 28; i++) { this->mm[i] = v1[i] - v2[i]; }
        return *this;
    }

    /// this =  scalar*v, v need not be a distinct object from "this",
    /// e.g. v.scale(1.5,v1);
    template <typename T0, typename T1>
    const Vec28<T>& scale(T0 scale, const Vec28<T1> &v)
    {
        for (int i = 0; i < 28; i++) { this->mm[i] = scale * v[i]; }
        return *this;
    }

    template <typename T0, typename T1>
    const Vec28<T> &div(T0 scalar, const Vec28<T1> &v)
    {
        for (int i = 0; i < 28; i++) { this->mm[i] = v[i] / scalar; }
        return *this;
    }

    /// Dot product
    T dot(const Vec28<T> &v) const
    {
        T dot;
        for (int i = 0; i < 28; i++) { dot += this->mm[i]*v.mm[i]; }
        return dot;
    }

    /// Length of the vector
    T length() const
    {
        T norm;
        for (int i = 0; i < 28; i++) { norm += this->mm[i]*this->mm[i]; }
        return std::sqrt(norm);
    }


    /// Squared length of the vector, much faster than length() as it
    /// does not involve square root
    T lengthSqr() const
    {
        T norm;
        for (int i = 0; i < 28; i++) { norm += this->mm[i]*this->mm[i]; }
        return norm;
    }

    /// Return a reference to itself after the exponent has been
    /// applied to all the vector components.
    inline const Vec28<T>& exp()
    {
        for (int i = 0; i < 28; i++) { this->mm[i] = std::exp(this->mm[i]); }
        return *this;
    }

    /// Return a reference to itself after log has been
    /// applied to all the vector components.
    inline const Vec28<T>& log()
    {
        for (int i = 0; i < 28; i++) { this->mm[i] = std::log(this->mm[i]); }
        return *this;
    }

    /// Return the sum of all the vector components.
    inline T sum() const
    {
        T sum;
        for (int i = 0; i < 28; i++) { sum += this->mm[i]; }
        return sum;
    }

    /// Return the product of all the vector components.
    inline T product() const
    {
        T prod;
        for (int i = 0; i < 28; i++) { prod *= this->mm[i]; }
        return prod;
    }

    /// this = normalized this
    bool normalize(T eps = static_cast<T>(1.0e-8))
    {
        T d = length();
        if (isApproxEqual(d, T(0), eps)) {
            return false;
        }
        *this *= (T(1) / d);
        return true;
    }

    /// return normalized this, throws if null vector
    Vec28<T> unit(T eps=0) const
    {
        T d;
        return unit(eps, d);
    }

    /// return normalized this and length, throws if null vector
    Vec28<T> unit(T eps, T& len) const
    {
        len = length();
        if (isApproxEqual(len, T(0), eps)) {
            throw ArithmeticError("Normalizing null 28-vector");
        }
        return *this / len;
    }

    /// return normalized this, or (1, 0, 0, 0) if this is null vector
    Vec28<T> unitSafe() const
    {
        T l2 = lengthSqr();
        Vec28<T> v(static_cast<T>(0));
        v[0] += 1;
        return l2 ? *this / static_cast<T>(sqrt(l2)) : Vec4<T>(1, 0, 0, 0);
    }

    /// Multiply each element of this vector by @a scalar.
    template <typename S>
    const Vec28<T> &operator*=(S scalar)
    {
        for (int i = 0; i < 28; i++) { this->mm[i] *= scalar; }
        return *this;
    }

    /// Multiply each element of this vector by the corresponding element of the given vector.
    template <typename S>
    const Vec28<T> &operator*=(const Vec28<S> &v1)
    {
        for (int i = 0; i < 28; i++) { this->mm[i] *= v1[i]; }
        return *this;
    }

    /// Divide each element of this vector by @a scalar.
    template <typename S>
    const Vec28<T> &operator/=(S scalar)
    {
        for (int i = 0; i < 28; i++) { this->mm[i] /= scalar; }
        return *this;
    }

    /// Divide each element of this vector by the corresponding element of the given vector.
    template <typename S>
    const Vec28<T> &operator/=(const Vec28<S> &v1)
    {
       for (int i = 0; i < 28; i++) { this->mm[i] /= v1[i]; }
        return *this;
    }

    /// Add @a scalar to each element of this vector.
    template <typename S>
    const Vec28<T> &operator+=(S scalar)
    {
        for (int i = 0; i < 28; i++) { this->mm[i] += scalar; }
        return *this;
    }

    /// Add each element of the given vector to the corresponding element of this vector.
    template <typename S>
    const Vec28<T> &operator+=(const Vec28<S> &v1)
    {
        for (int i = 0; i < 28; i++) { this->mm[i] += v1[i]; }
        return *this;
    }

    /// Subtract @a scalar from each element of this vector.
    template <typename S>
    const Vec28<T> &operator-=(S scalar)
    {
        for (int i = 0; i < 28; i++) { this->mm[i] -= scalar; }
        return *this;
    }

    /// Subtract each element of the given vector from the corresponding element of this vector.
    template <typename S>
    const Vec28<T> &operator-=(const Vec28<S> &v1)
    {
        for (int i = 0; i < 28; i++) { this->mm[i] -= v1[i]; }
        return *this;
    }

    // Number of cols, rows, elements
    static unsigned numRows() { return 1; }
    static unsigned numColumns()  { return 28; }
    static unsigned numElements()  { return 28; }

    /// Predefined constants, e.g.   Vec28 v = Vec28::xNegAxis();
    static Vec28<T> zero() { return Vec28<T>(0); }
    static Vec28<T> ones() { return Vec28<T>(1); }
};

/// Equality operator, does exact floating point comparisons
template <typename T0, typename T1>
inline bool operator==(const Vec28<T0> &v0, const Vec28<T1> &v1)
{
    for (int i = 0; i < 28; i++) {
        if (!isExactlyEqual(v0[i], v1[i])){
            return false;
        }
    }
    return true;
}

/// Inequality operator, does exact floating point comparisons
template <typename T0, typename T1>
inline bool operator!=(const Vec28<T0> &v0, const Vec28<T1> &v1) { return !(v0==v1); }

/// Multiply each element of the given vector by @a scalar and return the result.
template <typename S, typename T>
inline Vec28<typename promote<S, T>::type> operator*(S scalar, const Vec28<T> &v)
{ return v*scalar; }

/// Multiply each element of the given vector by @a scalar and return the result.
template <typename S, typename T>
inline Vec28<typename promote<S, T>::type> operator*(const Vec28<T> &v, S scalar)
{
    Vec28<typename promote<S, T>::type> result(v);
    result *= scalar;
    return result;
}

/// Multiply corresponding elements of @a v0 and @a v1 and return the result.
template <typename T0, typename T1>
inline Vec28<typename promote<T0, T1>::type> operator*(const Vec28<T0> &v0, const Vec28<T1> &v1)
{
    Vec28<typename promote<T0, T1>::type> result(static_cast<typename promote<T0, T1>::type>(0));
    for (int i = 0; i < 28; i++) { result[i] = v0[i]*v1[i]; }
    return result;
}

/// Divide @a scalar by each element of the given vector and return the result.
template <typename S, typename T>
inline Vec28<typename promote<S, T>::type> operator/(S scalar, const Vec28<T> &v)
{
    Vec28<typename promote<S, T>::type> result(v);
    for (int i = 0; i < 28; i++) { result[i] = scalar/result[i]; }
    return result;
}

/// Divide each element of the given vector by @a scalar and return the result.
template <typename S, typename T>
inline Vec28<typename promote<S, T>::type> operator/(const Vec28<T> &v, S scalar)
{
    Vec28<typename promote<S, T>::type> result(v);
    result /= scalar;
    return result;
}

/// Divide corresponding elements of @a v0 and @a v1 and return the result.
template <typename T0, typename T1>
inline Vec28<typename promote<T0, T1>::type> operator/(const Vec28<T0> &v0, const Vec28<T1> &v1)
{
    Vec28<typename promote<T0, T1>::type> result(v0);
    for (int i = 0; i < 28; i++) { result[i] /= v1[i]; }
    return result;
}

/// Add corresponding elements of @a v0 and @a v1 and return the result.
template <typename T0, typename T1>
inline Vec28<typename promote<T0, T1>::type> operator+(const Vec28<T0> &v0, const Vec28<T1> &v1)
{
    Vec28<typename promote<T0, T1>::type> result(v0);
    result += v1;
    return result;
}

/// Add @a scalar to each element of the given vector and return the result.
template <typename S, typename T>
inline Vec28<typename promote<S, T>::type> operator+(const Vec28<T> &v, S scalar)
{
    Vec28<typename promote<S, T>::type> result(v);
    result += scalar;
    return result;
}

/// Subtract corresponding elements of @a v0 and @a v1 and return the result.
template <typename T0, typename T1>
inline Vec28<typename promote<T0, T1>::type> operator-(const Vec28<T0> &v0, const Vec28<T1> &v1)
{
    Vec28<typename promote<T0, T1>::type> result(v0);
    result -= v1;
    return result;
}

/// Subtract @a scalar from each element of the given vector and return the result.
template <typename S, typename T>
inline Vec28<typename promote<S, T>::type> operator-(const Vec28<T> &v, S scalar)
{
    Vec28<typename promote<S, T>::type> result(v);
    result -= scalar;
    return result;
}

template <typename T>
inline bool
isApproxEqual(const Vec28<T>& a, const Vec28<T>& b)
{
    return a.eq(b);
}

template <typename T>
inline bool
isApproxEqual(const Vec28<T>& a, const Vec28<T>& b, const T eps)
{
    for (int i = 0; i < 28; i++) {
        if (!isApproxEqual(a[i], b[i], eps)){
            return false;
        }
    }
    return true;
}

template<typename T>
inline Vec28<T>
Abs(const Vec28<T>& v)
{
    Vec28<T> result(static_cast<T>(0));
    for (int i = 0; i < 28; i++) { result[i] = Abs(result[i]); }
    return result; 
}

/// @remark We are switching to a more explicit name because the semantics
/// are different from std::min/max. In that case, the function returns a
/// reference to one of the objects based on a comparator. Here, we must
/// fabricate a new object which might not match either of the inputs.

/// Return component-wise minimum of the two vectors.
template <typename T>
inline Vec28<T> minComponent(const Vec28<T> &v1, const Vec28<T> &v2)
{
    Vec28<T> result(static_cast<T>(0));
    for (int i = 0; i < 28; i++) { result[i] = std::min(v1[i], v2[i]); }
    return result; 
}

/// Return component-wise maximum of the two vectors.
template <typename T>
inline Vec28<T> maxComponent(const Vec28<T> &v1, const Vec28<T> &v2)
{
    Vec28<T> result(static_cast<T>(0));
    for (int i = 0; i < 28; i++) { result[i] = std::min(v1[i], v2[i]); }
    return result; 
}

/// @brief Return a vector with the exponent applied to each of
/// the components of the input vector.
template <typename T>
inline Vec28<T> Exp(Vec28<T> v) { return v.exp(); }

/// @brief Return a vector with log applied to each of
/// the components of the input vector.
template <typename T>
inline Vec28<T> Log(Vec28<T> v) { return v.log(); }

using Vec28i = Vec28<int32_t>;
using Vec28ui = Vec28<uint32_t>;

#if OPENVDB_ABI_VERSION_NUMBER >= 8
OPENVDB_IS_POD(Vec28i)
OPENVDB_IS_POD(Vec28ui)
#endif

} // namespace math
} // namespace OPENVDB_VERSION_NAME
} // namespace openvdb

#endif // OPENVDB_MATH_VEC28_HAS_BEEN_INCLUDED
