// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

#ifndef OPENVDB_MATH_VECX_HAS_BEEN_INCLUDED
#define OPENVDB_MATH_VECX_HAS_BEEN_INCLUDED

#include <openvdb/Exceptions.h>
#include "Math.h"
#include "Tuple.h"
#include <algorithm>
#include <cmath>
#include <type_traits>


namespace openvdb {
OPENVDB_USE_VERSION_NAMESPACE
namespace OPENVDB_VERSION_NAME {
namespace math {

template<int SIZE, typename T>
class VecX: public Tuple<SIZE, T>
{
public:
    using value_type = T;
    using ValueType = T;

    /// Trivial constructor, the vector is NOT initialized
    /// @note destructor, copy constructor, assignment operator and
    ///   move constructor are left to be defined by the compiler (default)
    VecX() = default;

    /// @brief Construct a vector all of whose components have the given value.
    explicit VecX(T val) { for (int i = 0; i < SIZE; i++) { this->mm[i] = val; } }

    /// Conversion constructor
    template<int SRC_SIZE, typename Source>
    explicit VecX(const Tuple<SRC_SIZE, Source> &v)
    {
        for (int i = 0; i < SRC_SIZE; i++) { this->mm[i] = static_cast<T>(v[i]); }
    }

    /// @brief Construct a vector all of whose components have the given value,
    /// which may be of an arithmetic type different from this vector's value type.
    /// @details Type conversion warnings are suppressed.
    template<typename Other>
    explicit VecX(Other val,
        typename std::enable_if<std::is_arithmetic<Other>::value, Conversion>::type = Conversion{})
    {
        for (int i = 0; i < SIZE; i++) {this->mm[i] = static_cast<T>(val);}
    }

    T* asPointer() { return this->mm; }
    const T* asPointer() const { return this->mm; }

    /// Alternative indexed reference to the elements
    T& operator()(int i) { return this->mm[i]; }

    /// Alternative indexed constant reference to the elements,
    T operator()(int i) const { return this->mm[i]; }

    /// Cast to double
    operator double() const {
        int max_i = 0;
        for (int i = 0; i < SIZE; ++i) {
            if (this->mm[i] > this->mm[max_i]) max_i = i;
        }
        return static_cast<double>(this->mm[max_i]);
    }

    /// Returns size of VecX
    const int getSize() const { return this->size; }

    /// Returns a Vec3 with the first three elements of the VecX.
    Vec2<T> getVec2() const { return Vec2<T>(this->mm[0], this->mm[1]); }

    /// Returns a Vec3 with the first three elements of the VecX.
    Vec3<T> getVec3() const { return Vec3<T>(this->mm[0], this->mm[1], this->mm[2]); }

    /// Returns a Vec4 with the first four elements of the VecX.
    Vec4<T> getVec4() const { return Vec4<T>(this->mm[0], this->mm[1], this->mm[2], this->mm[3]); }

    /// Set "this" vector to zero
    const VecX<SIZE, T>& setZero()
    {
        for (int i = 0; i < SIZE; i++) {this->mm[i] = 0;}
        return *this;
    }

    /// Assignment operator
    template<int SRC_SIZE, typename Source>
    const VecX<SIZE, T>& operator=(const VecX<SRC_SIZE, Source> &v)
    {
        // note: don't static_cast because that suppresses warnings i<=100/2 && j>100/2
        for (int i = 0; i < getSize() && i < v.getSize(); i++) {this->mm[i] = v[i];}
        return *this;
    }

    /// Test if "this" vector is equivalent to vector v with int
    template<int SRC_SIZE>
    bool eq(const VecX<SRC_SIZE, T> &v, T eps = static_cast<T>(1.0e-8)) const
    {
        if (SRC_SIZE != getSize()) { return false; }
        for (int i = 0; i < SRC_SIZE; i++) {
            if (!isApproxEqual(this->mm[i], v.mm[i], eps)){
                return false;
            }
        }
        return true;
    }

    /// Negation operator, for e.g.   v1 = -v2;
    VecX operator-() const
    {   
        VecX<SIZE, T> tmp(SIZE);     
        for (int i = 0; i < SIZE; i++) {
            tmp[i] = -this->mm[i];
        }
        return tmp;
    }

    /// this = v1 + v2
    /// "this", v1 and v2 need not be distinct objects, e.g. v.add(v1,v);
    template <typename T0, typename T1, int SIZE_0, int SIZE_1>
    const VecX<SIZE, T>& add(const VecX<SIZE_0, T0> &v1, const VecX<SIZE_1, T1> &v2)
    {
        if (SIZE_0 == SIZE_1)
        {
            for (int i = 0; i < SIZE_0; i++) { this->mm[i] = v1[i] + v2[i]; }
            return *this;
        } else if (SIZE_0 > SIZE_1){
            for (int i = 0; i < SIZE_1; i++) { this->mm[i] = v1[i] + v2[i]; }
            for (int i = SIZE_1; i < SIZE_0; i++) { this->mm[i] = v1[i]; }
            return *this;
        } else if (SIZE_1 > SIZE_0){
            for (int i = 0; i < SIZE_0; i++) { this->mm[i] = v1[i] + v2[i]; }
            for (int i = SIZE_0; i < SIZE_1; i++) { this->mm[i] = v2[i]; }
            return *this;
        }
    }


    /// this = v1 - v2
    /// "this", v1 and v2 need not be distinct objects, e.g. v.sub(v1,v);
    template <typename T0, typename T1, int SIZE_0, int SIZE_1>
    const VecX<SIZE, T>& sub(const VecX<SIZE_0, T0> &v1, const VecX<SIZE_1, T1> &v2)
    {
        if (SIZE_0 == SIZE_1)
        {
            for (int i = 0; i < SIZE_0; i++) { this->mm[i] = v1[i] - v2[i]; }
            return *this;
        } else if (SIZE_0 > SIZE_1){
            for (int i = 0; i < SIZE_1; i++) { this->mm[i] = v1[i] - v2[i]; }
            for (int i = SIZE_1; i < SIZE_0; i++) { this->mm[i] = v1[i]; }
            return *this;
        } else if (SIZE_1 > SIZE_0){
            for (int i = 0; i < SIZE_0; i++) { this->mm[i] = v1[i] - v2[i]; }
            for (int i = SIZE_0; i < SIZE_1; i++) { this->mm[i] = -v2[i]; }
            return *this;
        }
    }

    /// this =  scalar*v, v need not be a distinct object from "this",
    /// e.g. v.scale(1.5,v1);
    template <typename T0, typename T1>
    const VecX<SIZE, T>& scale(T0 scale, const VecX<SIZE, T1> &v)
    {
        for (int i = 0; i < SIZE; i++) { this->mm[i] = scale * v[i]; }
        return *this;
    }

    template <typename T0, typename T1>
    const VecX<SIZE, T> &div(T0 scalar, const VecX<SIZE, T1> &v)
    {
        for (int i = 0; i < SIZE; i++) { this->mm[i] = v[i] / scalar; }
        return *this;
    }

    /// Dot product
    T dot(const VecX<SIZE, T> &v) const
    {   
        T dot;
        for (int i = 0; i < SIZE; i++) { dot += this->mm[i]*v.mm[i]; }
        return dot;
    }

    /// Length of the vector
    T length() const
    {
        T norm;
        for (int i = 0; i < SIZE; i++) { norm += this->mm[i]*this->mm[i]; }
        return std::sqrt(norm);
    }

    /// Squared length of the vector, much faster than length() as it
    /// does not involve square root
    T lengthSqr() const
    {
        T norm;
        for (int i = 0; i < SIZE; i++) { norm += this->mm[i]*this->mm[i]; }
        return norm;
    }

    /// Return a reference to itself after the exponent has been
    /// applied to all the vector components.
    inline const VecX<SIZE, T>& exp()
    {
        for (int i = 0; i < SIZE; i++) { this->mm[i] = std::exp(this->mm[i]); }
        return *this;
    }

    /// Return a reference to itself after log has been
    /// applied to all the vector components.
    inline const VecX<SIZE, T>& log()
    {
        for (int i = 0; i < SIZE; i++) { this->mm[i] = std::log(this->mm[i]); }
        return *this;
    }

    /// Return the sum of all the vector components.
    inline T sum() const
    {
        T sum;
        for (int i = 0; i < SIZE; i++) { sum += this->mm[i]; }
        return sum;
    }

    /// Return the product of all the vector components.
    inline T product() const
    {
        T prod;
        for (int i = 0; i < SIZE; i++) { prod *= this->mm[i]; }
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
    VecX<SIZE, T> unit(T eps=0) const
    {
        T d;
        return unit(eps, d);
    }

    /// return normalized this and length, throws if null vector
    VecX<SIZE, T> unit(T eps, T& len) const
    {
        len = length();
        if (isApproxEqual(len, T(0), eps)) {
            throw ArithmeticError("Normalizing null Xd-vector");
        }
        return *this / len;
    }

    /// return normalized this, or (1, zeros(N-1)) if this is null vector
    VecX<SIZE, T> unitSafe() const
    {
        T l2 = lengthSqr();
        VecX<SIZE,T> v(static_cast<T>(0));
        v[0] += 1;
        return l2 ? *this / static_cast<T>(sqrt(l2)) : Vec4<T>(1, 0, 0, 0);
    }

    /// Multiply each element of this vector by @a scalar.
    template <typename S>
    const VecX<SIZE, T> &operator*=(S scalar)
    {
        for (int i = 0; i < SIZE; i++) { this->mm[i] *= scalar; }
        return *this;
    }

    /// Multiply each element of this vector by the corresponding element of the given vector.
    template <typename S>
    const VecX<SIZE, T> &operator*=(const VecX<SIZE, S> &v1)
    {
        for (int i = 0; i < SIZE; i++) { this->mm[i] *= v1[i]; }
        return *this;
    }

    /// Divide each element of this vector by @a scalar.
    template <typename S>
    const VecX<SIZE, T> &operator/=(S scalar)
    {
        for (int i = 0; i < SIZE; i++) { this->mm[i] /= scalar; }
        return *this;
    }

    /// Divide each element of this vector by the corresponding element of the given vector.
    template <typename S>
    const VecX<SIZE, T> &operator/=(const VecX<SIZE, S> &v1)
    {
        for (int i = 0; i < SIZE; i++) { this->mm[i] /= v1[i]; }
        return *this;
    }

    /// Add @a scalar to each element of this vector.
    template <typename S>
    const VecX<SIZE, T> &operator+=(S scalar)
    {
        for (int i = 0; i < SIZE; i++) { this->mm[i] += scalar; }
        return *this;
    }

    /// Add each element of the given vector to the corresponding element of this vector.
    template <typename S>
    const VecX<SIZE, T> &operator+=(const VecX<SIZE, S> &v1)
    {
        for (int i = 0; i < SIZE; i++) { this->mm[i] += v1[i]; }
        return *this;
    }

    /// Subtract @a scalar from each element of this vector.
    template <typename S>
    const VecX<SIZE, T> &operator-=(S scalar)
    {
        for (int i = 0; i < SIZE; i++) { this->mm[i] -= scalar; }
        return *this;
    }

    /// Subtract each element of the given vector from the corresponding element of this vector.
    template <typename S>
    const VecX<SIZE, T> &operator-=(const VecX<SIZE, S> &v1)
    {
        for (int i = 0; i < SIZE; i++) { this->mm[i] -= v1[i]; }
        return *this;
    }

    // Number of cols, rows, elements
    static unsigned numRows() { return 1; }
    static unsigned numColumns()  { return getSize(); }
    static unsigned numElements()  { return getSize(); }

    /// Predefined constants, e.g.   Vec4f v = Vec4f::xNegAxis();
    static VecX<SIZE, T> zero() { return VecX<SIZE, T>(0); }
    static VecX<SIZE, T> ones() { return VecX<SIZE, T>(1); }
};

/// Equality operator, does exact floating point comparisons
template <typename T0, typename T1, int SIZE_0, int SIZE_1>
inline bool operator==(const VecX<SIZE_0, T0> &v0, const VecX<SIZE_1, T1> &v1)
{
    if (SIZE_0 != SIZE_0) { return false; }
        for (int i = 0; i < SIZE_0; i++) {
            if (!isExactlyEqual(v0[i], v1[i])){
                return false;
            }
        }
        return true;
}

/// Inequality operator, does exact floating point comparisons
template <typename T0, typename T1, int SIZE_0, int SIZE_1>
inline bool operator!=(const VecX<SIZE_0, T0> &v0, const VecX<SIZE_1, T1> &v1) { return !(v0==v1); }

/// Multiply each element of the given vector by @a scalar and return the result.
template <int SIZE, typename S, typename T>
inline VecX<SIZE, typename promote<S, T>::type> operator*(S scalar, const VecX<SIZE, T> &v)
{ return v*scalar; }

/// Multiply each element of the given vector by @a scalar and return the result.
template <int SIZE, typename S, typename T>
inline VecX<SIZE, typename promote<S, T>::type> operator*(const VecX<SIZE, T> &v, S scalar)
{
    VecX<SIZE, typename promote<S, T>::type> result(v);
    result *= scalar;
    return result;
}

/// Multiply corresponding elements of @a v0 and @a v1 and return the result.
template <int SIZE, typename T0, typename T1>
inline VecX<SIZE, typename promote<T0, T1>::type> operator*(const VecX<SIZE, T0> &v0, const VecX<SIZE, T1> &v1)
{
    VecX<SIZE, typename promote<T0, T1>::type> result(static_cast<typename promote<T0, T1>::type>(0));
    for (int i = 0; i < SIZE; i++) { result[i] = v0[i]*v1[i]; }
    return result;
}

/// Divide @a scalar by each element of the given vector and return the result.
template <int SIZE, typename S, typename T>
inline VecX<SIZE, typename promote<S, T>::type> operator/(S scalar, const VecX<SIZE, T> &v)
{
    VecX<SIZE, typename promote<S, T>::type> result(v);
    for (int i = 0; i < SIZE; i++) { result[i] = scalar/result[i]; }
    return result;
}

/// Divide each element of the given vector by @a scalar and return the result.
template <int SIZE, typename S, typename T>
inline VecX<SIZE, typename promote<S, T>::type> operator/(const VecX<SIZE, T> &v, S scalar)
{
    VecX<SIZE, typename promote<S, T>::type> result(v);
    result /= scalar;
    return result;
}

/// Divide corresponding elements of @a v0 and @a v1 and return the result.
template <int SIZE, typename T0, typename T1>
inline VecX<SIZE, typename promote<T0, T1>::type> operator/(const VecX<SIZE, T0> &v0, const VecX<SIZE, T1> &v1)
{
    VecX<SIZE, typename promote<T0, T1>::type> result(v0);
    for (int i = 0; i < SIZE; i++) { result[i] /= v1[i]; }
    return result;
}

/// Add corresponding elements of @a v0 and @a v1 and return the result.
template <int SIZE, typename T0, typename T1>
inline VecX<SIZE, typename promote<T0, T1>::type> operator+(const VecX<SIZE, T0> &v0, const VecX<SIZE, T1> &v1)
{
    VecX<SIZE, typename promote<T0, T1>::type> result(v0);
    result += v1;
    return result;
}

/// Add @a scalar to each element of the given vector and return the result.
template <int SIZE, typename S, typename T>
inline VecX<SIZE, typename promote<S, T>::type> operator+(const VecX<SIZE, T> &v, S scalar)
{
    VecX<SIZE, typename promote<S, T>::type> result(v);
    result += scalar;
    return result;
}

/// Subtract corresponding elements of @a v0 and @a v1 and return the result.
template <int SIZE, typename T0, typename T1>
inline VecX<SIZE, typename promote<T0, T1>::type> operator-(const VecX<SIZE, T0> &v0, const VecX<SIZE, T1> &v1)
{
    VecX<SIZE, typename promote<T0, T1>::type> result(v0);
    result -= v1;
    return result;
}

/// Subtract @a scalar from each element of the given vector and return the result.
template <int SIZE, typename S, typename T>
inline VecX<SIZE, typename promote<S, T>::type> operator-(const VecX<SIZE, T> &v, S scalar)
{
    VecX<SIZE, typename promote<S, T>::type> result(v);
    result -= scalar;
    return result;
}

template <int SIZE, typename T>
inline bool
isApproxEqual(const VecX<SIZE, T>& a, const VecX<SIZE, T>& b)
{
    return a.eq(b);
}

template <typename T, int SIZE_0, int SIZE_1>
inline bool
isApproxEqual(const VecX<SIZE_0, T>& a, const VecX<SIZE_1, T>& b, const T eps)
{
    if (SIZE_0 != SIZE_1) { return false; }
        for (int i = 0; i < SIZE_0; i++) {
            if (!isApproxEqual(a[i], b[i], eps)){
                return false;
            }
        }
        return true;
}

template<int SIZE, typename T>
inline VecX<SIZE, T>
Abs(const VecX<SIZE, T>& v)
{
    VecX<SIZE, T> result(static_cast<T>(0));
    for (int i = 0; i < SIZE; i++) { result[i] = Abs(result[i]); }
    return result; 
}

/// @remark We are switching to a more explicit name because the semantics
/// are different from std::min/max. In that case, the function returns a
/// reference to one of the objects based on a comparator. Here, we must
/// fabricate a new object which might not match either of the inputs.

/// Return component-wise minimum of the two vectors.
template <int SIZE, typename T>
inline VecX<SIZE, T> minComponent(const VecX<SIZE, T> &v1, const VecX<SIZE, T> &v2)
{
    VecX<SIZE, T> result(static_cast<T>(0));
    for (int i = 0; i < SIZE; i++) { result[i] = std::min(v1[i], v2[i]); }
    return result; 
}

/// Return component-wise maximum of the two vectors.
template <int SIZE, typename T>
inline VecX<SIZE, T> maxComponent(const VecX<SIZE, T> &v1, const VecX<SIZE, T> &v2)
{
    VecX<SIZE, T> result(static_cast<T>(0));
    for (int i = 0; i < SIZE; i++) { result[i] = std::min(v1[i], v2[i]); }
    return result; 
}

/// @brief Return a vector with the exponent applied to each of
/// the components of the input vector.
template <int SIZE, typename T>
inline VecX<SIZE, T> Exp(VecX<SIZE, T> v) { return v.exp(); }

/// @brief Return a vector with log applied to each of
/// the components of the input vector.
template <int SIZE, typename T>
inline VecX<SIZE, T> Log(VecX<SIZE, T> v) { return v.log(); }

} // namespace math
} // namespace OPENVDB_VERSION_NAME
} // namespace openvdb

#endif // OPENVDB_MATH_VECX_HAS_BEEN_INCLUDED