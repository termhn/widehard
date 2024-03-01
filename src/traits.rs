use std::ops;

/// Arbitrary geometric product support for a type.
///
/// The geometric product is the sum of a grade-lowering (dot) and grade-raising
/// (wedge) product. These products are quite useful individually, as well as in
/// concert as the geometric product.
pub trait GeometricMul<Rhs>
where
    // NOTE: we'd ideally prefer to define add in terms of Self::Lower, as
    // the typical definition is "dot" + "wedge".  Unfortunately, the Lower type
    // can often be the base type, (e.g. f32), which would make this trait much
    // less ergonomic.
    Self::Upper: ops::Add<Self::Lower, Output = Self::Full>,
{
    type Upper;
    type Lower;
    type Full;

    fn dot(&self, v: &Rhs) -> Self::Lower;
    fn wedge(&self, v: &Rhs) -> Self::Upper;

    fn gmul(&self, other: &Rhs) -> Self::Full {
        self.wedge(other) + self.dot(other)
    }
}
