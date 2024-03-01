use std::ops::*;

use crate::traits::GeometricMul;
use crate::util::EqualsEps;
use crate::*;

macro_rules! vec3s {
    ($(($v2t:ident, $n:ident, $bn:ident, $rn:ident, $v4t:ident) => $t:ident),+) => {
        $(/// A set of three coordinates which may be interpreted as a point or vector in 3d space,
        /// or as a homogeneous 2d vector or point.
        ///
        /// Generally this distinction between a point and vector is more of a pain than it is worth
        /// to distinguish on a type level, however when converting to and from homogeneous
        /// coordinates it is quite important.
        #[derive(Clone, Copy, Debug, Default, PartialEq)]
        #[repr(C)]
        pub struct $n {
            pub x: $t,
            pub y: $t,
            pub z: $t,
        }

        impl $n {
            #[inline]
            pub const fn new(x: $t, y: $t, z: $t) -> Self {
                $n { x, y, z }
            }

            #[inline]
            pub const fn broadcast(val: $t) -> Self {
                Self::new(val, val, val)
            }

            #[inline]
            pub fn unit_x() -> Self {
                $n{ x: $t::splat(1.0), y: $t::splat(0.0), z: $t::splat(0.0) }
            }

            #[inline]
            pub fn unit_y() -> Self {
                $n{ x: $t::splat(0.0), y: $t::splat(1.0), z: $t::splat(0.0) }
            }

            #[inline]
            pub fn unit_z() -> Self {
                $n{ x: $t::splat(0.0), y: $t::splat(0.0), z: $t::splat(1.0) }
            }

            /// Create a homogeneous 3d *point* from this vector interpreted as a point,
            /// meaning the homogeneous component will start with a value of 1.0.
            #[inline]
            pub fn into_homogeneous_point(self) -> $v4t {
                $v4t { x: self.x, y: self.y, z: self.z, w: $t::splat(1.0) }
            }

            /// Create a homogeneous 3d *vector* from this vector,
            /// meaning the homogeneous component will always have a value of 0.0.
            #[inline]
            pub fn into_homogeneous_vector(self) -> $v4t {
                $v4t { x: self.x, y: self.y, z: self.z, w: $t::splat(0.0) }
            }

            /// Create a 3d point from a homogeneous 3d *point*, performing
            /// division by the homogeneous component. This should not be used
            /// for homogeneous 3d *vectors*, which will have 0 as their
            /// homogeneous component.
            #[inline]
            pub fn from_homogeneous_point(v: $v4t) -> Self {
                Self { x: v.x / v.w, y: v.y / v.w, z: v.z / v.w }
            }

            /// Create a 3d vector from homogeneous 2d *vector*, which simply
            /// discards the homogeneous component.
            #[inline]
            pub fn from_homogeneous_vector(v: $v4t) -> Self {
                v.into()
            }

            #[inline]
            pub fn dot(&self, other: $n) -> $t {
                GeometricMul::dot(self, &other)
            }

            /// The wedge (aka exterior) product of two vectors.
            ///
            /// This operation results in a bivector, which represents
            /// the plane parallel to the two vectors, and which has a
            /// 'oriented area' equal to the parallelogram created by extending
            /// the two vectors, oriented such that the positive direction is the
            /// one which would move `self` closer to `other`.
            #[inline]
            pub fn wedge(&self, other: $n) -> $bn {
                GeometricMul::wedge(self, &other)
            }

            /// The geometric product of this and another vector, which
            /// is defined as the sum of the dot product and the wedge product.
            ///
            /// This operation results in a 'rotor', named as such as it may define
            /// a rotation. The rotor which results from the geometric product
            /// will rotate in the plane parallel to the two vectors, by twice the angle between
            /// them and in the opposite direction (i.e. it will rotate in the direction that would
            /// bring `other` towards `self`, and rotate in that direction by twice the angle between them).
            #[inline]
            pub fn geom(&self, other: $n) -> $rn {
                GeometricMul::gmul(self, &other)
            }

            #[inline]
            pub fn rotate_by(&mut self, rotor: $rn) {
                rotor.rotate_vec(self);
            }

            #[inline]
            pub fn rotated_by(mut self, rotor: $rn) -> Self {
                rotor.rotate_vec(&mut self);
                self
            }

            #[inline]
            pub fn cross(&self, other: $n) -> Self {
                $n::new(
                    (self.y * other.z) + (-self.z * other.y),
                    (self.z * other.x) + (-self.x * other.z),
                    (self.x * other.y) + (-self.y * other.x),
                )
            }

            #[inline]
            pub fn reflect(&mut self, normal: $n) {
                *self -= $t::splat(2.0) * self.dot(normal) * normal;
            }

            #[inline]
            pub fn reflected(&self, normal: $n) -> Self {
                let mut a = *self;
                a.reflect(normal);
                a
            }

            #[inline]
            pub fn mag_sq(&self) -> $t {
                (self.x * self.x) + (self.y * self.y) + (self.z * self.z)
            }

            #[inline]
            pub fn mag(&self) -> $t {
                self.mag_sq().sqrt()
            }

            #[inline]
            pub fn normalize(&mut self) {
                let r_mag = $t::splat(1.0) / self.mag();
                self.x *= r_mag;
                self.y *= r_mag;
                self.z *= r_mag;
            }

            #[inline]
            #[must_use = "Did you mean to use `.normalize()` to normalize `self` in place?"]
            pub fn normalized(&self) -> Self {
                let mut r = self.clone();
                r.normalize();
                r
            }

            /// Normalize `self` in-place by interpreting it as a homogeneous point, i.e.
            /// scaling the vector to ensure the homogeneous component has length 1.
            #[inline]
            pub fn normalize_homogeneous_point(&mut self) {
                let recip_z = $t::splat(1.0) / self.z;
                self.x *= recip_z;
                self.y *= recip_z;
                self.z = $t::splat(1.0);
            }

            /// Normalize `self` by interpreting it as a homogeneous point, i.e.
            /// scaling the vector to ensure the homogeneous component has length 1.
            #[inline]
            #[must_use = "Did you mean to use `.normalize_homogeneous_point()` to normalize `self` in place?"]
            pub fn normalized_homogeneous_point(&self) -> Self {
                let mut r = self.clone();
                r.normalize_homogeneous_point();
                r
            }

            /// Convert `self` into a Vec2 by simply removing its `z` component.
            #[inline]
            pub fn truncated(&self) -> $v2t {
                $v2t::new(
                    self.x,
                    self.y
                )
            }

            #[inline]
            pub fn mul_add(&self, mul: $n, add: $n) -> Self {
                $n::new(
                    self.x.mul_add(mul.x, add.x),
                    self.y.mul_add(mul.y, add.y),
                    self.z.mul_add(mul.z, add.z),
                )
            }

            #[inline]
            pub fn abs(&self) -> Self {
                Self::new(self.x.abs(), self.y.abs(), self.z.abs())
            }

            #[inline]
            pub fn clamp(&mut self, min: Self, max: Self) {
                self.x = self.x.max(min.x).min(max.x);
                self.y = self.y.max(min.y).min(max.y);
                self.z = self.z.max(min.z).min(max.z);
            }

            #[inline]
            pub fn clamped(mut self, min: Self, max: Self) -> Self {
                self.clamp(min, max);
                self
            }

            #[inline]
            pub fn map<F>(&self, mut f: F) -> Self
                where F: FnMut($t) -> $t
            {
                $n::new(
                    f(self.x),
                    f(self.y),
                    f(self.z)
                )
            }

            #[inline]
            pub fn apply<F>(&mut self, mut f: F)
                where F: FnMut($t) -> $t
            {
                self.x = f(self.x);
                self.y = f(self.y);
                self.z = f(self.z);
            }

            #[inline]
            pub fn max_by_component(mut self, other: Self) -> Self {
                self.x = self.x.max(other.x);
                self.y = self.y.max(other.y);
                self.z = self.z.max(other.z);
                self
            }

            #[inline]
            pub fn min_by_component(mut self, other: Self) -> Self {
                self.x = self.x.min(other.x);
                self.y = self.y.min(other.y);
                self.z = self.z.min(other.z);
                self
            }

            #[inline]
            pub fn component_max(&self) -> $t {
                self.x.max(self.y).max(self.z)
            }

            #[inline]
            pub fn component_min(&self) -> $t {
                self.x.min(self.y).min(self.z)
            }

            #[inline]
            pub fn zero() -> Self {
                Self::broadcast($t::splat(0.0))
            }

            #[inline]
            pub fn one() -> Self {
                Self::broadcast($t::splat(1.0))
            }

            #[inline]
            pub const fn xy(&self) -> $v2t {
                $v2t::new(self.x, self.y)
            }

            #[inline]
            pub fn xyzw(&self) -> $v4t {
                $v4t::new(self.x, self.y, self.z, $t::splat(0.0))
            }

            /// Get the [`core::alloc::Layout`] of `Self`
            #[inline]
            pub fn layout() -> alloc::alloc::Layout {
                alloc::alloc::Layout::from_size_align(std::mem::size_of::<Self>(), std::mem::align_of::<$t>()).unwrap()
            }

            /// Interpret `self` as a statically-sized array of its base numeric type
            #[inline]
            pub fn as_array(&self) -> &[$t; 3] {
                let ptr = self as *const $n as *const [$t; 3];
                unsafe { &*ptr }
            }

            /// Interpret `self` as a statically-sized array of its base numeric type
            #[inline]
            pub fn as_mut_array(&mut self) -> &mut [$t; 3] {
                let ptr = self as *mut $n as *mut [$t; 3];
                unsafe { &mut *ptr }
            }

            /// Interpret `self` as a slice of its base numeric type
            #[inline]
            pub fn as_slice(&self) -> &[$t] {
                // This is safe because we are statically bounding our slices to the size of these
                // vectors
                unsafe {
                    std::slice::from_raw_parts(self as *const $n as *const $t, 3)
                }
            }

            /// Interpret `self` as a slice of its base numeric type
            #[inline]
            pub fn as_mut_slice(&mut self) -> &mut [$t] {
                // This is safe because we are statically bounding our slices to the size of these
                // vectors
                unsafe {
                    std::slice::from_raw_parts_mut(self as *mut $n as *mut $t, 3)
                }
            }

            #[inline]
            pub fn as_byte_slice(&self) -> &[u8] {
                // This is safe because we are statically bounding our slices to the size of these
                // vectors
                unsafe {
                    std::slice::from_raw_parts(self as *const $n as *const u8, 3 * std::mem::size_of::<$t>())
                }
            }

            #[inline]
            pub fn as_mut_byte_slice(&mut self) -> &mut [u8] {
                // This is safe because we are statically bounding our slices to the size of these
                // vectors
                unsafe {
                    std::slice::from_raw_parts_mut(self as *mut $n as *mut u8, 3 * std::mem::size_of::<$t>())
                }
            }

            /// Returns a constant unsafe pointer to the underlying data in the underlying type.
            /// This function is safe because all types here are repr(C) and can be represented
            /// as their underlying type.
            ///
            /// # Safety
            ///
            /// It is up to the caller to correctly use this pointer and its bounds.
            #[inline]
            pub const fn as_ptr(&self) -> *const $t {
                self as *const $n as *const $t
            }

            /// Returns a mutable unsafe pointer to the underlying data in the underlying type.
            /// This function is safe because all types here are repr(C) and can be represented
            /// as their underlying type.
            ///
            /// # Safety
            ///
            /// It is up to the caller to correctly use this pointer and its bounds.
            #[inline]
            pub fn as_mut_ptr(&mut self) -> *mut $t {
                self as *mut $n as *mut $t
            }
        }

        impl GeometricMul<$n> for $n {
            type Lower = $t;
            type Upper = $bn;
            type Full = $rn;

            #[inline]
            fn dot(&self, other: &$n) -> $t {
                (self.x * other.x) + (self.y * other.y) + (self.z * other.z)
            }

            /// The wedge (aka exterior) product of two vectors.
            ///
            /// This operation results in a bivector, which represents
            /// the plane parallel to the two vectors, and which has a
            /// 'oriented area' equal to the parallelogram created by extending
            /// the two vectors, oriented such that the positive direction is the
            /// one which would move `self` closer to `other`.
            #[inline]
            fn wedge(&self, other: &$n) -> $bn {
                $bn::new(
                    (self.x * other.y) - (self.y * other.x),
                    (self.x * other.z) - (self.z * other.x),
                    (self.y * other.z) - (self.z * other.y),
                )
            }
        }

        impl EqualsEps for $n {
            fn eq_eps(self, other: Self) -> bool {
                self.x.eq_eps(other.x) && self.y.eq_eps(other.y) && self.z.eq_eps(other.z)
            }
        }

        impl From<$n> for [$t; 3] {
            #[inline]
            fn from(v: $n) -> Self {
                [v.x, v.y, v.z]
            }
        }

        impl From<[$t; 3]> for $n {
            #[inline]
            fn from(comps: [$t; 3]) -> Self {
                Self::new(comps[0], comps[1], comps[2])
            }
        }

        impl From<&[$t; 3]> for $n {
            #[inline]
            fn from(comps: &[$t; 3]) -> Self {
                Self::from(*comps)
            }
        }

        impl From<&mut [$t; 3]> for $n {
            #[inline]
            fn from(comps: &mut [$t; 3]) -> Self {
                Self::from(*comps)
            }
        }

        impl From<($t, $t, $t)> for $n {
            #[inline]
            fn from(comps: ($t, $t, $t)) -> Self {
                Self::new(comps.0, comps.1, comps.2)
            }
        }

        impl From<&($t, $t, $t)> for $n {
            #[inline]
            fn from(comps: &($t, $t, $t)) -> Self {
                Self::from(*comps)
            }
        }

        impl From<$n> for ($t, $t, $t) {
            #[inline]
            fn from(v: $n) -> Self {
                (v.x, v.y, v.z)
            }
        }

        impl Add for $n {
            type Output = Self;
            #[inline]
            fn add(self, rhs: $n) -> Self {
                $n::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
            }
        }

        impl AddAssign for $n {
            #[inline]
            fn add_assign(&mut self, rhs: $n) {
                self.x += rhs.x;
                self.y += rhs.y;
                self.z += rhs.z;
            }
        }

        impl Sub for $n {
            type Output = Self;
            #[inline]
            fn sub(self, rhs: $n) -> Self {
                $n::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
            }
        }

        impl SubAssign for $n {
            #[inline]
            fn sub_assign(&mut self, rhs: $n) {
                self.x -= rhs.x;
                self.y -= rhs.y;
                self.z -= rhs.z;
            }
        }

        impl Mul for $n {
            type Output = Self;
            #[inline]
            fn mul(self, rhs: $n) -> Self {
                $n::new(self.x * rhs.x, self.y * rhs.y, self.z * rhs.z)
            }
        }

        impl Mul<$n> for $t {
            type Output = $n;
            #[inline]
            fn mul(self, rhs: $n) -> $n {
                $n::new(self * rhs.x, self * rhs.y, self * rhs.z)
            }
        }

        impl Mul<$t> for $n {
            type Output = $n;
            #[inline]
            fn mul(self, rhs: $t) -> $n {
                $n::new(self.x * rhs, self.y * rhs, self.z * rhs)
            }
        }

        impl MulAssign for $n {
            #[inline]
            fn mul_assign(&mut self, rhs: $n) {
                self.x *= rhs.x;
                self.y *= rhs.y;
                self.z *= rhs.z;
            }
        }

        impl MulAssign<$t> for $n {
            #[inline]
            fn mul_assign(&mut self, rhs: $t) {
                self.x *= rhs;
                self.y *= rhs;
                self.z *= rhs;
            }
        }

        impl Div for $n {
            type Output = Self;
            #[inline]
            fn div(self, rhs: $n) -> Self {
                $n::new(self.x / rhs.x, self.y / rhs.y, self.z / rhs.z)
            }
        }

        impl Div<$t> for $n {
            type Output = $n;
            #[inline]
            fn div(self, rhs: $t) -> $n {
                $n::new(self.x / rhs, self.y / rhs, self.z / rhs)
            }
        }

        impl DivAssign for $n {
            #[inline]
            fn div_assign(&mut self, rhs: $n) {
                self.x /= rhs.x;
                self.y /= rhs.y;
                self.z /= rhs.z;
            }
        }

        impl DivAssign<$t> for $n {
            #[inline]
            fn div_assign(&mut self, rhs: $t) {
                self.x /= rhs;
                self.y /= rhs;
                self.z /= rhs;
            }
        }

        impl Neg for $n {
            type Output = $n;
            #[inline]
            fn neg(self) -> $n {
                self * $t::splat(-1.0)
            }
        }

        impl Index<usize> for $n {
            type Output = $t;

            fn index(&self, index: usize) -> &Self::Output {
                match index {
                    0 => &self.x,
                    1 => &self.y,
                    2 => &self.z,
                    _ => panic!("Invalid for vector of type: {}", std::any::type_name::<$n>()),
                }
            }
        }

        impl IndexMut<usize> for $n {
            fn index_mut(&mut self, index: usize) -> &mut Self::Output {
                match index {
                    0 => &mut self.x,
                    1 => &mut self.y,
                    2 => &mut self.z,
                    _ => panic!("Invalid for vector of type: {}", std::any::type_name::<$n>()),
                }
            }
        }

        impl std::iter::Sum<$n> for $n {
            fn sum<I>(iter: I) -> Self where I: Iterator<Item = Self> {
                // Kahan summation algorithm
                // https://en.wikipedia.org/wiki/Kahan_summation_algorithm
                let mut sum = $n::zero();
                let mut c = $n::zero();
                for v in iter {
                    let y = v - c;
                    let t = sum + y;
                    c = (t - sum) - y;
                    sum = t;
                }
                sum
            }
        }
        )+
    }
}

// SCALAR VEC3 IMPLS

macro_rules! impl_scalar_vec3s {
    ($(($vt:ident, $v2t:ident, $v4t:ident) => $t:ident),+) => {
        $(impl $vt {
            #[inline]
            pub fn refract(&mut self, normal: Self, eta: $t) {
                *self = self.refracted(normal, eta);
            }

            #[inline]
            pub fn refracted(&self, normal: Self, eta: $t) -> Self {
                let n = normal;
                let i = *self;
                let ndi = n.dot(i);
                let k = 1.0 - eta * eta * (1.0 - ndi * ndi);
                if k < 0.0 {
                    Self::zero()
                } else {
                    i * eta - (eta * ndi + k.sqrt()) * n
                }
            }
        }

        impl From<$v2t> for $vt {
            #[inline]
            fn from(vec: $v2t) -> Self {
                Self {
                    x: vec.x,
                    y: vec.y,
                    z: 0.0,
                }
            }
        }

        impl From<$v4t> for $vt {
            #[inline]
            fn from(vec: $v4t) -> Self {
                Self {
                    x: vec.x,
                    y: vec.y,
                    z: vec.z,
                }
            }
        })+
    };
}

// WIDE VEC3 IMPLS

macro_rules! impl_wide_vec3s {
    ($($vt:ident => $tt:ident, $t:ident, $maskt:ident, $nonwidet:ident, $v2t:ident, $v4t:ident),+) => {
        $(impl $vt {
            #[inline]
            pub fn new_splat(x: $tt, y: $tt, z: $tt) -> Self {
                Self {
                    x: $t::splat(x),
                    y: $t::splat(y),
                    z: $t::splat(z),
                }
            }

            #[inline]
            pub fn splat(vec: $nonwidet) -> Self {
                Self {
                    x: $t::splat(vec.x),
                    y: $t::splat(vec.y),
                    z: $t::splat(vec.z),
                }
            }

            /// Blend two vectors together lanewise using `mask` as a mask.
            ///
            /// This is essentially a bitwise blend operation, such that any point where
            /// there is a 1 bit in `mask`, the output will put the bit from `tru`, while
            /// where there is a 0 bit in `mask`, the output will put the bit from `fals`
            #[inline]
            pub fn blend(mask: $maskt, tru: Self, fals: Self) -> Self {
                Self {
                    x: mask.blend(tru.x, fals.x),
                    y: mask.blend(tru.y, fals.y),
                    z: mask.blend(tru.z, fals.z),
                }
            }

            #[inline]
            pub fn refract(&mut self, normal: Self, eta: $t) {
                *self = self.refracted(normal, eta);
            }

            #[inline]
            pub fn refracted(&self, normal: Self, eta: $t) -> Self {
                let n = normal;
                let i = *self;
                let one = $t::splat(1.0);
                let ndi = n.dot(i);

                let k = one - eta * eta * (one - ndi * ndi);
                let mask = k.cmp_lt($t::splat(0.0));

                let out = i.mul_add(Self::broadcast(eta), -(eta * ndi + k.sqrt()) * n);

                Self::blend(mask, Self::zero(), out)
            }
        }

        impl From<$v2t> for $vt {
            #[inline]
            fn from(vec: $v2t) -> Self {
                Self {
                    x: vec.x,
                    y: vec.y,
                    z: $t::splat(0.0),
                }
            }
        }

        impl From<$v4t> for $vt {
            #[inline]
            fn from(vec: $v4t) -> Self {
                Self {
                    x: vec.x,
                    y: vec.y,
                    z: vec.z,
                }
            }
        })+
    }
}

impl From<Vec3x4> for [Vec3; 4] {
    #[inline]
    fn from(v: Vec3x4) -> Self {
        let xs: [f32; 4] = v.x.into();
        let ys: [f32; 4] = v.y.into();
        let zs: [f32; 4] = v.z.into();
        [
            Vec3::new(xs[0], ys[0], zs[0]),
            Vec3::new(xs[1], ys[1], zs[1]),
            Vec3::new(xs[2], ys[2], zs[2]),
            Vec3::new(xs[3], ys[3], zs[3]),
        ]
    }
}

impl From<[Vec3; 4]> for Vec3x4 {
    #[inline]
    fn from(vecs: [Vec3; 4]) -> Self {
        Self {
            x: f32x4::from([vecs[0].x, vecs[1].x, vecs[2].x, vecs[3].x]),
            y: f32x4::from([vecs[0].y, vecs[1].y, vecs[2].y, vecs[3].y]),
            z: f32x4::from([vecs[0].z, vecs[1].z, vecs[2].z, vecs[3].z]),
        }
    }
}

impl From<Vec3x8> for [Vec3; 8] {
    #[inline]
    fn from(v: Vec3x8) -> [Vec3; 8] {
        let xs: [f32; 8] = v.x.into();
        let ys: [f32; 8] = v.y.into();
        let zs: [f32; 8] = v.z.into();
        [
            Vec3::new(xs[0], ys[0], zs[0]),
            Vec3::new(xs[1], ys[1], zs[1]),
            Vec3::new(xs[2], ys[2], zs[2]),
            Vec3::new(xs[3], ys[3], zs[3]),
            Vec3::new(xs[4], ys[4], zs[4]),
            Vec3::new(xs[5], ys[5], zs[5]),
            Vec3::new(xs[6], ys[6], zs[6]),
            Vec3::new(xs[7], ys[7], zs[7]),
        ]
    }
}

impl From<[Vec3; 8]> for Vec3x8 {
    #[inline]
    fn from(vecs: [Vec3; 8]) -> Self {
        Self {
            x: f32x8::from([
                vecs[0].x, vecs[1].x, vecs[2].x, vecs[3].x, vecs[4].x, vecs[5].x, vecs[6].x,
                vecs[7].x,
            ]),
            y: f32x8::from([
                vecs[0].y, vecs[1].y, vecs[2].y, vecs[3].y, vecs[4].y, vecs[5].y, vecs[6].y,
                vecs[7].y,
            ]),
            z: f32x8::from([
                vecs[0].z, vecs[1].z, vecs[2].z, vecs[3].z, vecs[4].z, vecs[5].z, vecs[6].z,
                vecs[7].z,
            ]),
        }
    }
}

#[cfg(feature = "f64")]
impl From<DVec3x2> for [DVec3; 2] {
    #[inline]
    fn from(v: DVec3x2) -> Self {
        let xs: [f64; 2] = v.x.into();
        let ys: [f64; 2] = v.y.into();
        let zs: [f64; 2] = v.z.into();
        [
            DVec3::new(xs[0], ys[0], zs[0]),
            DVec3::new(xs[1], ys[1], zs[1]),
        ]
    }
}

#[cfg(feature = "f64")]
impl From<[DVec3; 2]> for DVec3x2 {
    #[inline]
    fn from(vecs: [DVec3; 2]) -> Self {
        Self {
            x: f64x2::from([vecs[0].x, vecs[1].x]),
            y: f64x2::from([vecs[0].y, vecs[1].y]),
            z: f64x2::from([vecs[0].z, vecs[1].z]),
        }
    }
}

#[cfg(feature = "f64")]
impl From<DVec3x4> for [DVec3; 4] {
    fn from(v: DVec3x4) -> Self {
        let xs: [f64; 4] = v.x.into();
        let ys: [f64; 4] = v.y.into();
        let zs: [f64; 4] = v.z.into();
        [
            DVec3::new(xs[0], ys[0], zs[0]),
            DVec3::new(xs[1], ys[1], zs[1]),
            DVec3::new(xs[2], ys[2], zs[2]),
            DVec3::new(xs[3], ys[3], zs[3]),
        ]
    }
}

#[cfg(feature = "f64")]
impl From<[DVec3; 4]> for DVec3x4 {
    #[inline]
    fn from(vecs: [DVec3; 4]) -> Self {
        Self {
            x: f64x4::from([vecs[0].x, vecs[1].x, vecs[2].x, vecs[3].x]),
            y: f64x4::from([vecs[0].y, vecs[1].y, vecs[2].y, vecs[3].y]),
            z: f64x4::from([vecs[0].z, vecs[1].z, vecs[2].z, vecs[3].z]),
        }
    }
}

vec3s!(
    (Vec2, Vec3, Bivec3, Rotor3, Vec4) => f32,
    (Vec2x4, Vec3x4, Bivec3x4, Rotor3x4, Vec4x4) => f32x4,
    (Vec2x8, Vec3x8, Bivec3x8, Rotor3x8, Vec4x8) => f32x8
);

#[cfg(feature = "f64")]
vec3s!(
    (DVec2, DVec3, DBivec3, DRotor3, DVec4) => f64,
    (DVec2x2, DVec3x2, DBivec3x2, DRotor3x2, DVec4x2) => f64x2,
    (DVec2x4, DVec3x4, DBivec3x4, DRotor3x4, DVec4x4) => f64x4
);

impl_scalar_vec3s!(
    (Vec3, Vec2, Vec4) => f32
);

#[cfg(feature = "f64")]
impl_scalar_vec3s!(
    (DVec3, DVec2, DVec4) => f64
);

impl_wide_vec3s!(
    Vec3x4 => f32, f32x4, m32x4, Vec3, Vec2x4, Vec4x4,
    Vec3x8 => f32, f32x8, m32x8, Vec3, Vec2x8, Vec4x8
);

#[cfg(feature = "f64")]
impl_wide_vec3s!(
    DVec3x2 => f64, f64x2, m64x2, DVec3, DVec2x2, DVec4x2,
    DVec3x4 => f64, f64x4, m64x4, DVec3, DVec2x4, DVec4x4
);
