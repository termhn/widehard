//! Bivectors, i.e. oriented areas.
//!
//! A bivector is an *oriented area*, and is equivalent
//! to the result of the exterior (wedge) product of two vectors, i.e.
//! `u ∧ v`. This means it is the *oriented area* of the parallelogram
//! created by attaching two vectors and then extending them into a parallelogram.
//!
//! This may be hard to visualize at first, but bivectors are as fundamental as vectors. If vectors
//! are a representation of *lines*, then bivectors are a representation of *planes*.
//!
//! A normalized bivector can be thought of as representing a plane of rotation and the *direction of rotation*
//! inside that plane such that a *positive* rotation follows the orientation of the bivector. When
//! you obtain a bivector by taking the exterior product of two vectors, the positive direction of rotation
//! is defined as the one that *brings the first vector closer to the second*. For example, a bivector
//! created by taking the exterior product `x ∧ y` of the x and y basis vectors will create a unit
//! bivector that represents the xy plane, with orientation such that a positive rotation of `x` inside
//! the plane would bring `x` closer to `y`. This is why positive rotation is generally defined as
//! "counter clockwise" in 2d, since such a rotation brings `x` to `y`.
//!
//! Much like vectors can be represented as a linear combination of *basis vectors*, i.e.
//! a vector "component representation," bivectors can be represented as a linear combination
//! of *basis bivectors*. If the basis vectors are the unit vectors in the direction of each
//! canonical axis of a space, then the basis bivectors are the *unit area planes* in each of the
//! canonical planes.
//!
//! In 2d, there is only one basis plane, the xy plane, which represents all of 2d space. As such, in 2d
//! there is only *one* basis bivector, while there are *two* basis vectors. This means that a 2d bivector
//! has only one component.
//!
//! In 3d, there are three basis planes, the xy plane, the xz plane, and the yz plane, which are respectively
//! the planes parallel to those combinations of the x, y, and z basis vectors. Therefore, a 3d bivector has
//! three components, each of which represents the *projected area* of that bivector onto one of the three
//! basis bivectors. This is analogous to how vector components represent the *projected length* of that vector
//! onto each unit vector.
use crate::*;

use crate::traits::GeometricMul;
use crate::util::*;

use std::ops::*;

macro_rules! bivec2s {
    ($(($bn:ident) => $t:ident),+) => {
        $(
        /// A bivector in 2d space.
        ///
        /// Since in 2d there is only one plane in the whole of 2d space, a 2d bivector
        /// has only one component.
        ///
        /// Please see the module level documentation for more information on bivectors generally!
        #[derive(Clone, Copy, Debug, Default, PartialEq)]
        #[repr(C)]
        pub struct $bn {
            pub xy: $t
        }

        impl $bn {
            #[inline]
            pub const fn new(xy: $t) -> Self {
                Self {
                    xy
                }
            }

            #[inline]
            pub fn zero() -> Self {
                Self::new($t::splat(0.0))
            }

            #[inline]
            pub fn unit_xy() -> Self {
                Self::new($t::splat(1.0))
            }

            #[inline]
            pub fn mag_sq(&self) -> $t {
                self.xy * self.xy
            }

            #[inline]
            pub fn mag(&self) -> $t {
                self.mag_sq().sqrt()
            }

            #[inline]
            pub fn normalize(&mut self) {
                let mag = self.mag();
                self.xy /= mag;
            }

            #[inline]
            #[must_use = "Did you mean to use `.normalize()` to normalize `self` in place?"]
            pub fn normalized(&self) -> Self {
                let mut r = self.clone();
                r.normalize();
                r
            }

            #[inline]
            pub fn dot(&self, rhs: Self) -> $t {
                self.xy * rhs.xy
            }

            #[inline]
            pub fn layout() -> alloc::alloc::Layout {
                alloc::alloc::Layout::from_size_align(std::mem::size_of::<Self>(), std::mem::align_of::<$t>()).unwrap()
            }

            #[inline]
            pub fn as_slice(&self) -> &[$t] {
                // This is safe because we are statically bounding our slices to the size of these
                // vectors
                unsafe {
                    std::slice::from_raw_parts(self as *const $bn as *const $t, 1)
                }
            }


            #[inline]
            pub fn as_byte_slice(&self) -> &[u8] {
                // This is safe because we are statically bounding our slices to the size of these
                // vectors
                unsafe {
                    std::slice::from_raw_parts(self as *const $bn as *const u8, std::mem::size_of::<$t>())
                }
            }

            #[inline]
            pub fn as_mut_slice(&mut self) -> &mut [$t] {
                // This is safe because we are statically bounding our slices to the size of these
                // vectors
                unsafe {
                    std::slice::from_raw_parts_mut(self as *mut $bn as *mut $t, 1)
                }
            }

            #[inline]
            pub fn as_mut_byte_slice(&mut self) -> &mut [u8] {
                // This is safe because we are statically bounding our slices to the size of these
                // vectors
                unsafe {
                    std::slice::from_raw_parts_mut(self as *mut $bn as *mut u8, std::mem::size_of::<$t>())
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
                self as *const $bn as *const $t
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
                self as *mut $bn as *mut $t
            }
        }

        impl EqualsEps for $bn {
            fn eq_eps(self, other: Self) -> bool {
                self.xy.eq_eps(other.xy)
            }
        }

        impl Add for $bn {
            type Output = Self;
            #[inline]
            fn add(mut self, rhs: $bn) -> Self {
                self += rhs;
                self
            }
        }

        impl AddAssign for $bn {
            #[inline]
            fn add_assign(&mut self, rhs: $bn) {
                self.xy += rhs.xy;
            }
        }

        impl Sub for $bn {
            type Output = Self;
            #[inline]
            fn sub(mut self, rhs: $bn) -> Self {
                self -= rhs;
                self
            }
        }

        impl SubAssign for $bn {
            #[inline]
            fn sub_assign(&mut self, rhs: $bn) {
                self.xy -= rhs.xy;
            }
        }

        impl Mul for $bn {
            type Output = Self;
            #[inline]
            fn mul(mut self, rhs: $bn) -> Self {
                self *= rhs;
                self
            }
        }

        impl Mul<$bn> for $t {
            type Output = $bn;
            #[inline]
            fn mul(self, mut rhs: $bn) -> $bn {
                rhs *= self;
                rhs
            }
        }

        impl Mul<$t> for $bn {
            type Output = Self;
            #[inline]
            fn mul(mut self, rhs: $t) -> Self {
                self *= rhs;
                self
            }
        }

        impl MulAssign for $bn {
            #[inline]
            fn mul_assign(&mut self, rhs: Self) {
                self.xy *= rhs.xy;
            }
        }

        impl MulAssign<$t> for $bn {
            #[inline]
            fn mul_assign(&mut self, rhs: $t) {
                self.xy *= rhs;
            }
        }

        impl Div for $bn {
            type Output = Self;
            #[inline]
            fn div(mut self, rhs: $bn) -> Self {
                self /= rhs;
                self
            }
        }

        impl Div<$t> for $bn {
            type Output = $bn;
            #[inline]
            fn div(mut self, rhs: $t) -> $bn {
                self.xy /= rhs;
                self
            }
        }

        impl DivAssign for $bn {
            #[inline]
            fn div_assign(&mut self, rhs: $bn) {
                self.xy /= rhs.xy;
            }
        }

        impl DivAssign<$t> for $bn {
            #[inline]
            fn div_assign(&mut self, rhs: $t) {
                self.xy /= rhs;
            }
        }

        impl Neg for $bn {
            type Output = Self;
            #[inline]
            fn neg(mut self) -> Self {
                self.xy = -self.xy;
                self
            }
        }
        )+
    }
}

macro_rules! bivec3s {
    ($($bn:ident => ($vt:ident, $rt:ident, $t:ident)),+) => {
        $(
        /// A bivector in 3d space.
        ///
        /// In 3d, a bivector has 3 components, each one representing the signed *projected area* of the bivector
        /// onto one of the 3 *basis bivectors*, which can be thought of as corresponding to each of the
        /// three basis planes. This is analogous to the components of a 3d vector, which correspond to the
        /// *projected length* of the vector onto the three basis *vectors. Since in 3d, there are three
        /// components for both vectors and bivectors, 3d bivectors have been historically confused with
        /// 3d vectors quite a lot.
        ///
        /// Please see the module level documentation for more information on bivectors generally!
        #[derive(Clone, Copy, Debug, Default, PartialEq)]
        #[repr(C)]
        pub struct $bn {
            pub xy: $t,
            pub xz: $t,
            pub yz: $t,
        }

        impl EqualsEps for $bn {
            fn eq_eps(self, other: Self) -> bool {
                self.xy.eq_eps(other.xy) && self.xz.eq_eps(other.xz) && self.yz.eq_eps(other.yz)
            }
        }

        impl $bn {
            #[inline]
            pub const fn new(xy: $t, xz: $t, yz: $t) -> Self {
                Self {
                    xy, xz, yz
                }
            }

            #[inline]
            pub fn zero() -> Self {
                Self::new($t::splat(0.0), $t::splat(0.0), $t::splat(0.0))
            }

            /// Create the bivector which represents the same plane of rotation as a given
            /// normalized 'axis vector'
            #[inline]
            pub fn from_normalized_axis(v: $vt) -> Self {
                Self::new(v.z, -v.y, v.x)
            }

            #[inline]
            pub fn unit_xy() -> Self {
                Self::new($t::splat(1.0), $t::splat(0.0), $t::splat(0.0))
            }

            #[inline]
            pub fn unit_xz() -> Self {
                Self::new($t::splat(0.0), $t::splat(1.0), $t::splat(0.0))
            }

            #[inline]
            pub fn unit_yz() -> Self {
                Self::new($t::splat(0.0), $t::splat(0.0), $t::splat(1.0))
            }

            #[inline]
            pub fn mag_sq(&self) -> $t {
                (self.xy * self.xy) + (self.xz * self.xz) + (self.yz * self.yz)
            }

            #[inline]
            pub fn mag(&self) -> $t {
                self.mag_sq().sqrt()
            }

            #[inline]
            pub fn normalize(&mut self) {
                let mag = self.mag();
                self.xy /= mag;
                self.xz /= mag;
                self.yz /= mag;
            }

            #[inline]
            #[must_use = "Did you mean to use `.normalize()` to normalize `self` in place?"]
            pub fn normalized(&self) -> Self {
                let mut r = self.clone();
                r.normalize();
                r
            }

            #[inline]
            pub fn dot(&self, rhs: Self) -> $t {
                (self.xy * rhs.xy) + (self.xz * rhs.xz) + (self.yz * rhs.yz)
            }

            #[inline]
            pub fn reverse(&self) -> Self {
                -*self
            }

            #[inline]
            pub fn layout() -> alloc::alloc::Layout {
                alloc::alloc::Layout::from_size_align(std::mem::size_of::<Self>(), std::mem::align_of::<$t>()).unwrap()
            }

            #[inline]
            pub fn as_slice(&self) -> &[$t] {
                // This is safe because we are statically bounding our slices to the size of these
                // vectors
                unsafe {
                    std::slice::from_raw_parts(self as *const $bn as *const $t, 3)
                }
            }


            #[inline]
            pub fn as_byte_slice(&self) -> &[u8] {
                // This is safe because we are statically bounding our slices to the size of these
                // vectors
                unsafe {
                    std::slice::from_raw_parts(self as *const $bn as *const u8, 3 * std::mem::size_of::<$t>())
                }
            }

            #[inline]
            pub fn as_mut_slice(&mut self) -> &mut [$t] {
                // This is safe because we are statically bounding our slices to the size of these
                // vectors
                unsafe {
                    std::slice::from_raw_parts_mut(self as *mut $bn as *mut $t, 3)
                }
            }

            #[inline]
            pub fn as_mut_byte_slice(&mut self) -> &mut [u8] {
                // This is safe because we are statically bounding our slices to the size of these
                // vectors
                unsafe {
                    std::slice::from_raw_parts_mut(self as *mut $bn as *mut u8, 3 * std::mem::size_of::<$t>())
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
                self as *const $bn as *const $t
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
                self as *mut $bn as *mut $t
            }
        }

        impl Add for $bn {
            type Output = Self;
            #[inline]
            fn add(mut self, rhs: $bn) -> Self {
                self += rhs;
                self
            }
        }

        impl GeometricMul<$bn> for $bn {
            type Lower = $t;
            type Upper = $bn;
            type Full = $rt;

            //       other
            //     xy  xz  yz
            // xy  -1 -yz  xz
            // xz  yz  -1 -xy
            // yz -xz  xy  -1

            #[inline]
            fn dot(&self, other: &$bn) -> Self::Lower {
                (self.xy * -other.xy) + (self.yz * -other.yz) + (self.xz * -other.xz)
            }

            /// The wedge (aka exterior) product of two bivectors.
            #[inline]
            fn wedge(&self, other: &$bn) -> Self::Upper {
                $bn::new(
                    (self.xz * -other.yz) + -(self.yz * -other.xz),
                    (self.xy * other.yz)  + -(-self.yz * -other.xy),
                    (-self.xy * other.xz) + -(-self.xz * other.xy),
                )
            }
        }

        impl Add<$t> for $bn {
            type Output = $rt;
            #[inline]
            fn add(self, rhs: $t) -> Self::Output {
                $rt::new(rhs, self)
            }
        }

        impl AddAssign for $bn {
            #[inline]
            fn add_assign(&mut self, rhs: $bn) {
                self.xy += rhs.xy;
                self.xz += rhs.xz;
                self.yz += rhs.yz;
            }
        }

        impl Sub for $bn {
            type Output = Self;
            #[inline]
            fn sub(mut self, rhs: $bn) -> Self {
                self -= rhs;
                self
            }
        }

        impl SubAssign for $bn {
            #[inline]
            fn sub_assign(&mut self, rhs: $bn) {
                self.xy -= rhs.xy;
                self.xz -= rhs.xz;
                self.yz -= rhs.yz;
            }
        }

        impl Mul for $bn {
            type Output = Self;
            #[inline]
            fn mul(mut self, rhs: $bn) -> Self {
                self *= rhs;
                self
            }
        }

        impl Mul<$bn> for $t {
            type Output = $bn;
            #[inline]
            fn mul(self, mut rhs: $bn) -> $bn {
                rhs *= self;
                rhs
            }
        }

        impl Mul<$t> for $bn {
            type Output = Self;
            #[inline]
            fn mul(mut self, rhs: $t) -> Self {
                self *= rhs;
                self
            }
        }

        impl MulAssign for $bn {
            #[inline]
            fn mul_assign(&mut self, rhs: Self) {
                self.xy *= rhs.xy;
                self.xz *= rhs.xz;
                self.yz *= rhs.yz;
            }
        }

        impl MulAssign<$t> for $bn {
            #[inline]
            fn mul_assign(&mut self, rhs: $t) {
                self.xy *= rhs;
                self.xz *= rhs;
                self.yz *= rhs;
            }
        }

        impl Div for $bn {
            type Output = Self;
            #[inline]
            fn div(mut self, rhs: $bn) -> Self {
                self /= rhs;
                self
            }
        }

        impl Div<$t> for $bn {
            type Output = $bn;
            #[inline]
            fn div(mut self, rhs: $t) -> $bn {
                self /= rhs;
                self
            }
        }

        impl DivAssign for $bn {
            #[inline]
            fn div_assign(&mut self, rhs: $bn) {
                self.xy /= rhs.xy;
                self.xz /= rhs.xz;
                self.yz /= rhs.yz;
            }
        }

        impl DivAssign<$t> for $bn {
            #[inline]
            fn div_assign(&mut self, rhs: $t) {
                self.xy /= rhs;
                self.xz /= rhs;
                self.yz /= rhs;
            }
        }

        impl Neg for $bn {
            type Output = Self;
            #[inline]
            fn neg(mut self) -> Self {
                self.xy = -self.xy;
                self.xz = -self.xz;
                self.yz = -self.yz;
                self
            }
        }
        )+
    }
}

bivec2s!(
    (Bivec2) => f32,
    (Bivec2x4) => f32x4,
    (Bivec2x8) => f32x8
);

#[cfg(feature = "f64")]
bivec2s!(
    (DBivec2) => f64,
    (DBivec2x2) => f64x2,
    (DBivec2x4) => f64x4
);

bivec3s!(
    Bivec3 => (Vec3, Rotor3, f32),
    Bivec3x4 => (Vec3x4, Rotor3x4, f32x4),
    Bivec3x8 => (Vec3x8, Rotor3x8, f32x8)
);

#[cfg(feature = "f64")]
bivec3s!(
    DBivec3 => (DVec3, DRotor3, f64),
    DBivec3x2 => (DVec3x2, DRotor3x2, f64x2),
    DBivec3x4 => (DVec3x4, DRotor3x4, f64x4)
);

#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck::TestResult;
    use quickcheck_macros::quickcheck;

    #[quickcheck]
    fn test_bivector_inverse(bv: (f32, f32, f32)) -> TestResult {
        let bivec = Bivec3::new(bv.0, bv.1, bv.2);

        if bivec.mag().is_nan() || bivec.mag().is_infinite() || bivec.mag() < 0.1e5 {
            return TestResult::discard();
        }

        let scale = bivec.gmul(&bivec.reverse());
        if scale.bv.mag() > 0.1e6 {
            return TestResult::from_bool(false);
        }
        let inverse = bivec.reverse() / scale.s;

        let unit = bivec.gmul(&inverse);

        TestResult::from_bool(unit.bv.mag().abs() < 0.1e6 && (unit.s - 1.0).abs() < 0.1e6)
    }
}
