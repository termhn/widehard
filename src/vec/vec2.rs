use std::ops::*;

use crate::util::*;
use crate::*;

macro_rules! vec2s {
    ($(($n:ident, $bn:ident, $rn:ident, $v3t:ident, $v4t:ident) => $t:ident),+) => {
        $(
        /// A set of two coordinates which may be interpreted as a vector or point in 2d space.
        ///
        /// Generally this distinction between a point and vector is more of a pain than it is worth
        /// to distinguish on a type level, however when converting to and from homogeneous
        /// coordinates it is quite important.
        #[derive(Clone, Copy, Debug, Default)]
        #[repr(C)]
        pub struct $n {
            pub x: $t,
            pub y: $t,
        }

        impl $n {
            #[inline]
            pub fn new(x: $t, y: $t) -> Self {
                $n { x, y }
            }

            #[inline]
            pub fn broadcast(val: $t) -> Self {
                Self::new(val, val)
            }

            #[inline]
            pub fn unit_x() -> Self {
                $n{ x: $t::splat(1.0), y: $t::splat(0.0) }
            }

            #[inline]
            pub fn unit_y() -> Self {
                $n{ x: $t::splat(0.0), y: $t::splat(1.0) }
            }

            /// Create a homogeneous 2d *point* from this vector interpreted as a point,
            /// meaning the homogeneous component will start with a value of 1.0.
            #[inline]
            pub fn into_homogeneous_point(self) -> $v3t {
                $v3t { x: self.x, y: self.y, z: $t::splat(1.0) }
            }

            /// Create a homogeneous 2d *vector* from this vector,
            /// meaning the homogeneous component will always have a value of 0.0.
            #[inline]
            pub fn into_homogeneous_vector(self) -> $v3t {
                $v3t { x: self.x, y: self.y, z: $t::splat(0.0) }
            }

            /// Create a 2d point from a homogeneous 2d *point*, performing
            /// division by the homogeneous component. This should not be used
            /// for homogeneous 2d *vectors*, which will have 0 as their
            /// homogeneous component.
            #[inline]
            pub fn from_homogeneous_point(v: $v3t) -> Self {
                Self { x: v.x / v.z, y: v.y / v.z }
            }

            /// Create a 2d vector from homogeneous 2d *vector*, which simply
            /// discards the homogeneous component.
            #[inline]
            pub fn from_homogeneous_vector(v: $v3t) -> Self {
                v.into()
            }

            #[inline]
            pub fn dot(&self, other: $n) -> $t {
                self.x.mul_add(other.x, self.y * other.y)
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
                $bn::new(self.x.mul_add(other.y, -(other.x * self.y)))
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
                $rn::new(self.dot(other), self.wedge(other))
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
            pub fn reflected(&self, normal: $n) -> Self {
                *self - ($t::splat(2.0) * self.dot(normal) * normal)
            }


            #[inline]
            pub fn mag_sq(&self) -> $t {
                self.x.mul_add(self.x, self.y * self.y)
            }

            #[inline]
            pub fn mag(&self) -> $t {
                self.mag_sq().sqrt()
            }

            #[inline]
            pub fn normalize(&mut self) {
                let mag = self.mag();
                self.x /= mag;
                self.y /= mag;
            }

            #[inline]
            pub fn normalized(&self) -> Self {
                let mut r = self.clone();
                r.normalize();
                r
            }

            #[inline]
            pub fn mul_add(&self, mul: $n, add: $n) -> Self {
                $n::new(
                    self.x.mul_add(mul.x, add.x),
                    self.y.mul_add(mul.y, add.y),
                )
            }

            #[inline]
            pub fn abs(&self) -> Self {
                Self::new(self.x.abs(), self.y.abs())
            }

            #[inline]
            pub fn clamp(&mut self, min: Self, max: Self) {
                self.x = self.x.max(min.x).min(max.x);
                self.y = self.y.max(min.y).min(max.y);
            }

            #[inline]
            pub fn clamped(mut self, min: Self, max: Self) -> Self {
                self.clamp(min, max);
                self
            }

            #[inline]
            pub fn map<F>(&self, f: F) -> Self
                where F: Fn($t) -> $t
            {
                $n::new(
                    f(self.x),
                    f(self.y),
                )
            }

            #[inline]
            pub fn apply<F>(&mut self, f: F)
                where F: Fn($t) -> $t
            {
                self.x = f(self.x);
                self.y = f(self.y);
            }

            #[inline]
            pub fn max_by_component(mut self, other: Self) -> Self {
                self.x = self.x.max(other.x);
                self.y = self.y.max(other.y);
                self
            }

            #[inline]
            pub fn min_by_component(mut self, other: Self) -> Self {
                self.x = self.x.min(other.x);
                self.y = self.y.min(other.y);
                self
            }

            #[inline]
            pub fn component_max(&self) -> $t {
                self.x.max(self.y)
            }

            #[inline]
            pub fn component_min(&self) -> $t {
                self.x.min(self.y)
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
            pub fn xyz(&self) -> $v3t {
                $v3t::new(self.x, self.y, $t::splat(0.0))
            }

            #[inline]
            pub fn xyzw(&self) -> $v4t {
                $v4t::new(self.x, self.y, $t::splat(0.0), $t::splat(0.0))
            }

            #[inline]
            pub fn layout() -> alloc::alloc::Layout {
                alloc::alloc::Layout::from_size_align(std::mem::size_of::<Self>(), std::mem::align_of::<$t>()).unwrap()
            }

            #[inline]
            pub fn as_array(&self) -> &[$t; 2] {
                use std::convert::TryInto;
                self.as_slice().try_into().unwrap()
            }

            #[inline]
            pub fn as_slice(&self) -> &[$t] {
                // This is safe because we are statically bounding our slices to the size of these
                // vectors
                unsafe {
                    std::slice::from_raw_parts(self as *const $n as *const $t, 2)
                }
            }


            #[inline]
            pub fn as_byte_slice(&self) -> &[u8] {
                // This is safe because we are statically bounding our slices to the size of these
                // vectors
                unsafe {
                    std::slice::from_raw_parts(self as *const $n as *const u8, 2 * std::mem::size_of::<$t>())
                }
            }

            #[inline]
            pub fn as_mut_slice(&mut self) -> &mut [$t] {
                // This is safe because we are statically bounding our slices to the size of these
                // vectors
                unsafe {
                    std::slice::from_raw_parts_mut(self as *mut $n as *mut $t, 2)
                }
            }

            #[inline]
            pub fn as_mut_byte_slice(&mut self) -> &mut [u8] {
                // This is safe because we are statically bounding our slices to the size of these
                // vectors
                unsafe {
                    std::slice::from_raw_parts_mut(self as *mut $n as *mut u8, 2 * std::mem::size_of::<$t>())
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
            pub fn as_ptr(&self) -> *const $t {
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

        impl Into<[$t; 2]> for $n {
            #[inline]
            fn into(self) -> [$t; 2] {
                [self.x, self.y]
            }
        }

        impl From<[$t; 2]> for $n {
            #[inline]
            fn from(comps: [$t; 2]) -> Self {
                Self::new(comps[0], comps[1])
            }
        }

        impl From<&[$t; 2]> for $n {
            #[inline]
            fn from(comps: &[$t; 2]) -> Self {
                Self::from(*comps)
            }
        }

        impl From<&mut [$t; 2]> for $n {
            #[inline]
            fn from(comps: &mut [$t; 2]) -> Self {
                Self::from(*comps)
            }
        }

        impl From<($t, $t)> for $n {
            #[inline]
            fn from(comps: ($t, $t)) -> Self {
                Self::new(comps.0, comps.1)
            }
        }

        impl From<&($t, $t)> for $n {
            #[inline]
            fn from(comps: &($t, $t)) -> Self {
                Self::from(*comps)
            }
        }

        impl From<$n> for ($t, $t) {
            #[inline]
            fn from(v: $n) -> Self {
                (v.x, v.y)
            }
        }

        impl EqualsEps for $n {
            fn eq_eps(self, other: Self) -> bool {
                self.x.eq_eps(other.x) && self.y.eq_eps(other.y)
            }
        }

        impl Add for $n {
            type Output = Self;
            #[inline]
            fn add(self, rhs: $n) -> Self {
                $n::new(self.x + rhs.x, self.y + rhs.y)
            }
        }

        impl AddAssign for $n {
            #[inline]
            fn add_assign(&mut self, rhs: $n) {
                self.x += rhs.x;
                self.y += rhs.y;
            }
        }

        impl Sub for $n {
            type Output = Self;
            #[inline]
            fn sub(self, rhs: $n) -> Self {
                $n::new(self.x - rhs.x, self.y - rhs.y)
            }
        }

        impl SubAssign for $n {
            #[inline]
            fn sub_assign(&mut self, rhs: $n) {
                self.x -= rhs.x;
                self.y -= rhs.y;
            }
        }

        impl Mul for $n {
            type Output = Self;
            #[inline]
            fn mul(self, rhs: $n) -> Self {
                $n::new(self.x * rhs.x, self.y * rhs.y)
            }
        }

        impl Mul<$n> for $t {
            type Output = $n;
            #[inline]
            fn mul(self, rhs: $n) -> $n {
                $n::new(self * rhs.x, self * rhs.y)
            }
        }

        impl Mul<$t> for $n {
            type Output = $n;
            #[inline]
            fn mul(self, rhs: $t) -> $n {
                $n::new(self.x * rhs, self.y * rhs)
            }
        }

        impl MulAssign for $n {
            #[inline]
            fn mul_assign(&mut self, rhs: $n) {
                self.x *= rhs.x;
                self.y *= rhs.y;
            }
        }

        impl MulAssign<$t> for $n {
            #[inline]
            fn mul_assign(&mut self, rhs: $t) {
                self.x *= rhs;
                self.y *= rhs;
            }
        }

        impl Div for $n {
            type Output = Self;
            #[inline]
            fn div(self, rhs: $n) -> Self {
                $n::new(self.x / rhs.x, self.y / rhs.y)
            }
        }

        impl Div<$t> for $n {
            type Output = $n;
            #[inline]
            fn div(self, rhs: $t) -> $n {
                $n::new(self.x / rhs, self.y / rhs)
            }
        }

        impl DivAssign for $n {
            #[inline]
            fn div_assign(&mut self, rhs: $n) {
                self.x /= rhs.x;
                self.y /= rhs.y;
            }
        }

        impl DivAssign<$t> for $n {
            #[inline]
            fn div_assign(&mut self, rhs: $t) {
                self.x /= rhs;
                self.y /= rhs;
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
                    _ => panic!("Invalid for vector of type: {}", std::any::type_name::<$n>()),
                }
            }
        }

        impl IndexMut<usize> for $n {
            fn index_mut(&mut self, index: usize) -> &mut Self::Output {
                match index {
                    0 => &mut self.x,
                    1 => &mut self.y,
                    _ => panic!("Invalid for vector of type: {}", std::any::type_name::<$n>()),
                }
            }
        }
        )+
    };
}

vec2s!(
    (Vec2, Bivec2, Rotor2, Vec3, Vec4) => f32,
    (Vec2x4, Bivec2x4, Rotor2x4, Vec3x4, Vec4x4) => f32x4,
    (Vec2x8, Bivec2x8, Rotor2x8, Vec3x8, Vec4x8) => f32x8,

    (DVec2, DBivec2, DRotor2, DVec3, DVec4) => f64,
    (DVec2x2, DBivec2x2, DRotor2x2, DVec3x2, DVec4x2) => f64x2,
    (DVec2x4, DBivec2x4, DRotor2x4, DVec3x4, DVec4x4) => f64x4
);

#[cfg(feature = "nightly")]
vec2s!(
    (Vec2x16, Bivec2x16, Rotor2x16, Vec3x16, Vec4x16) => f32x16,

    (DVec2x8, DBivec2x8, DRotor2x8, DVec3x8, DVec4x8) => f64x8
);
// (Vec2x8, Bivec2x8, Rotor2x8, Vec3x8, Vec4x8) => f32x8);

// SCALAR VEC2 IMPLS

macro_rules! impl_scalar_vec2s {
    ($(($vt:ident, $v3t:ident) => $t:ident),+) => {
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

        impl From<$v3t> for $vt {
            #[inline]
            fn from(vec: $v3t) -> Self {
                Self { x: vec.x, y: vec.y }
            }
        }

        impl PartialEq for $vt {
            fn eq(&self, other: &Self) -> bool {
                self.x == other.x && self.y == other.y
            }
        })+
    };
}

impl_scalar_vec2s!(
    (Vec2, Vec3) => f32,
    (DVec2, DVec3) => f64
);

// WIDE VEC2 IMPLS

macro_rules! impl_wide_vec2s {
    ($($vt:ident => $tt:ident, $t:ident, $maskt:ident, $nonwidet:ident, $v3t:ident),+) => {
        $(impl $vt {
            #[inline]
            pub fn new_splat(x: $tt, y: $tt) -> Self {
                Self {
                    x: $t::splat(x),
                    y: $t::splat(y),
                }
            }

            #[inline]
            pub fn splat(vec: $nonwidet) -> Self {
                Self {
                    x: $t::splat(vec.x),
                    y: $t::splat(vec.y),
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

                let out = i * eta - (eta * ndi + k.sqrt()) * n;

                Self::blend(mask, Self::zero(), out)
            }
        }

        impl From<$nonwidet> for $vt {
            #[inline]
            fn from(vec: $nonwidet) -> Self {
                Self::splat(vec)
            }
        }

        impl From<$v3t> for $vt {
            #[inline]
            fn from(vec: $v3t) -> Self {
                Self { x: vec.x, y: vec.y }
            }
        })+
    }
}

impl_wide_vec2s!(
    Vec2x4 => f32, f32x4, m32x4, Vec2, Vec3x4,
    Vec2x8 => f32, f32x8, m32x8, Vec2, Vec3x8,

    DVec2x2 => f64, f64x2, m64x2, DVec2, DVec3x2,
    DVec2x4 => f64, f64x4, m64x4, DVec2, DVec3x4
);

#[cfg(feature = "nightly")]
impl_wide_vec2s!(
    Vec2x16 => f32, f32x16, m32x16, Vec2, Vec3x16,

    DVec2x8 => f64, f64x8, m64x8, DVec2, DVec3x8
);

impl Into<[Vec2; 4]> for Vec2x4 {
    #[inline]
    fn into(self) -> [Vec2; 4] {
        let xs: [f32; 4] = self.x.into();
        let ys: [f32; 4] = self.y.into();
        [
            Vec2::new(xs[0], ys[0]),
            Vec2::new(xs[1], ys[1]),
            Vec2::new(xs[2], ys[2]),
            Vec2::new(xs[3], ys[3]),
        ]
    }
}

impl From<[Vec2; 4]> for Vec2x4 {
    #[inline]
    fn from(vecs: [Vec2; 4]) -> Self {
        Self {
            x: f32x4::from([vecs[0].x, vecs[1].x, vecs[2].x, vecs[3].x]),
            y: f32x4::from([vecs[0].y, vecs[1].y, vecs[2].y, vecs[3].y]),
        }
    }
}

impl Into<[Vec2; 8]> for Vec2x8 {
    #[inline]
    fn into(self) -> [Vec2; 8] {
        let xs: [f32; 8] = self.x.into();
        let ys: [f32; 8] = self.y.into();
        [
            Vec2::new(xs[0], ys[0]),
            Vec2::new(xs[1], ys[1]),
            Vec2::new(xs[2], ys[2]),
            Vec2::new(xs[3], ys[3]),
            Vec2::new(xs[4], ys[4]),
            Vec2::new(xs[5], ys[5]),
            Vec2::new(xs[6], ys[6]),
            Vec2::new(xs[7], ys[7]),
        ]
    }
}

impl From<[Vec2; 8]> for Vec2x8 {
    #[inline]
    fn from(vecs: [Vec2; 8]) -> Self {
        Self {
            x: f32x8::from([
                vecs[0].x, vecs[1].x, vecs[2].x, vecs[3].x, vecs[4].x, vecs[5].x, vecs[6].x,
                vecs[7].x,
            ]),
            y: f32x8::from([
                vecs[0].y, vecs[1].y, vecs[2].y, vecs[3].y, vecs[4].y, vecs[5].y, vecs[6].y,
                vecs[7].y,
            ]),
        }
    }
}

#[cfg(feature = "nightly")]
impl Into<[Vec2; 16]> for Vec2x16 {
    #[inline]
    fn into(self) -> [Vec2; 16] {
        let xs: [f32; 16] = self.x.into();
        let ys: [f32; 16] = self.y.into();
        [
            Vec2::new(xs[0], ys[0]),
            Vec2::new(xs[1], ys[1]),
            Vec2::new(xs[2], ys[2]),
            Vec2::new(xs[3], ys[3]),
            Vec2::new(xs[4], ys[4]),
            Vec2::new(xs[5], ys[5]),
            Vec2::new(xs[6], ys[6]),
            Vec2::new(xs[7], ys[7]),
            Vec2::new(xs[8], ys[8]),
            Vec2::new(xs[9], ys[9]),
            Vec2::new(xs[10], ys[10]),
            Vec2::new(xs[11], ys[11]),
            Vec2::new(xs[12], ys[12]),
            Vec2::new(xs[13], ys[13]),
            Vec2::new(xs[14], ys[14]),
            Vec2::new(xs[15], ys[15]),
        ]
    }
}

#[cfg(feature = "nightly")]
impl From<[Vec2; 16]> for Vec2x16 {
    #[inline]
    fn from(vecs: [Vec2; 16]) -> Self {
        Self {
            x: f32x16::from([
                vecs[0].x, vecs[1].x, vecs[2].x, vecs[3].x, vecs[4].x, vecs[5].x, vecs[6].x,
                vecs[7].x, vecs[8].x, vecs[9].x, vecs[10].x, vecs[11].x, vecs[12].x, vecs[13].x,
                vecs[14].x, vecs[15].x,
            ]),
            y: f32x16::from([
                vecs[0].y, vecs[1].y, vecs[2].y, vecs[3].y, vecs[4].y, vecs[5].y, vecs[6].y,
                vecs[7].y, vecs[8].y, vecs[9].y, vecs[10].y, vecs[11].y, vecs[12].y, vecs[13].y,
                vecs[14].y, vecs[15].y,
            ]),
        }
    }
}

impl Into<[DVec2; 2]> for DVec2x2 {
    #[inline]
    fn into(self) -> [DVec2; 2] {
        let xs: [f64; 2] = self.x.into();
        let ys: [f64; 2] = self.y.into();
        [DVec2::new(xs[0], ys[0]), DVec2::new(xs[1], ys[1])]
    }
}

impl From<[DVec2; 2]> for DVec2x2 {
    #[inline]
    fn from(vecs: [DVec2; 2]) -> Self {
        Self {
            x: f64x2::from([vecs[0].x, vecs[1].x]),
            y: f64x2::from([vecs[0].y, vecs[1].y]),
        }
    }
}

impl Into<[DVec2; 4]> for DVec2x4 {
    #[inline]
    fn into(self) -> [DVec2; 4] {
        let xs: [f64; 4] = self.x.into();
        let ys: [f64; 4] = self.y.into();
        [
            DVec2::new(xs[0], ys[0]),
            DVec2::new(xs[1], ys[1]),
            DVec2::new(xs[2], ys[2]),
            DVec2::new(xs[3], ys[3]),
        ]
    }
}

impl From<[DVec2; 4]> for DVec2x4 {
    #[inline]
    fn from(vecs: [DVec2; 4]) -> Self {
        Self {
            x: f64x4::from([vecs[0].x, vecs[1].x, vecs[2].x, vecs[3].x]),
            y: f64x4::from([vecs[0].y, vecs[1].y, vecs[2].y, vecs[3].y]),
        }
    }
}

#[cfg(feature = "nightly")]
impl Into<[DVec2; 8]> for DVec2x8 {
    #[inline]
    fn into(self) -> [DVec2; 8] {
        let xs: [f64; 8] = self.x.into();
        let ys: [f64; 8] = self.y.into();
        [
            DVec2::new(xs[0], ys[0]),
            DVec2::new(xs[1], ys[1]),
            DVec2::new(xs[2], ys[2]),
            DVec2::new(xs[3], ys[3]),
            DVec2::new(xs[4], ys[4]),
            DVec2::new(xs[5], ys[5]),
            DVec2::new(xs[6], ys[6]),
            DVec2::new(xs[7], ys[7]),
        ]
    }
}

#[cfg(feature = "nightly")]
impl From<[DVec2; 8]> for DVec2x8 {
    #[inline]
    fn from(vecs: [DVec2; 8]) -> Self {
        Self {
            x: f64x8::from([
                vecs[0].x, vecs[1].x, vecs[2].x, vecs[3].x, vecs[4].x, vecs[5].x, vecs[6].x,
                vecs[7].x,
            ]),
            y: f64x8::from([
                vecs[0].y, vecs[1].y, vecs[2].y, vecs[3].y, vecs[4].y, vecs[5].y, vecs[6].y,
                vecs[7].y,
            ]),
        }
    }
}
