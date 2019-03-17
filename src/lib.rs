/// Rust affine transforms
///
///
/// Based in:
/// Python affine transforms
/// by Matthew T. Perry on Sun 13 September 2015
/// http://www.perrygeo.com/python-affine-transforms.html
///
/// | x' |   | a  b  c | | x |
/// | y' | = | d  e  f | | y |
/// | 1  |   | 0  0  1 | | 1 |
///
/// Affine Python library
/// a = width of a pixel
/// b = row rotation (typically zero)
/// c = x-coordinate of the upper-left corner of the upper-left pixel
/// d = column rotation (typically zero)
/// e = height of a pixel (typically negative)
/// f = y-coordinate of the of the upper-left corner of the upper-left pixel
///
/// ESRI World File
/// a = width of a pixel
/// d = column rotation (typically zero)
/// b = row rotation (typically zero)
/// e = height of a pixel (typically negative)
/// c = x-coordinate of the center of the upper-left pixel
/// f = y-coordinate of the center of the upper-left pixel
///
/// It's important to note that the c and f parameters refer to the center of the cell, not the origin!
///
/// GDAL with the "Geotransform" array
/// c = x-coordinate of the upper-left corner of the upper-left pixel
/// a = width of a pixel
/// b = row rotation (typically zero)
/// f = y-coordinate of the of the upper-left corner of the upper-left pixel
/// d = column rotation (typically zero)
/// e = height of a pixel (typically negative)
///
///
/// a = width of a pixel
/// b = row rotation (typically zero)
/// c = x-coordinate of the upper-left corner of the upper-left pixel
/// d = column rotation (typically zero)
/// e = height of a pixel (typically negative)
/// f = y-coordinate of the of the upper-left corner of the upper-left pixel
/// g = 0
/// h = 0
/// i = 1
/// ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i')
///
/// x' = ax + c (x' = ax + by + c with x-rotation)
/// y' = ey + f (y' = dx + ey + d with y-rotation)
///
/// x = (x' - c) / a
/// y = (y' - f) / e
///
/// raster (x: column, y: row)
/// geospatial coordinate (x': longitude, y': latitude)
/// pixel resolution (a: width, e: height)
///

use std::ops::{Mul, Add, Sub, Neg, Not};

///
///
///
#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Affine<T> {
    // pixel's width: x-resolution
    width: T,

    // row rotation (typically zero)
    row_rotation: T,

    // col origin: x-offset
    col_origin: T,

    // column rotation (typically zero)
    col_rotation: T,

    // pixel's height: y-resolution
    height: T,

    // row origin: y-offset
    row_origin: T
}

///
///
///
pub trait AffineTransform<T>
    where T: Mul<Output = T> + Add<Output = T> + Sub<Output = T> + Neg<Output = T> + Copy + PartialEq + PartialOrd + Default{

    ///
    ///
    ///
    fn new(width: T,
           row_rotation: T,
           col_origin: T,
           col_rotation: T,
           height: T,
           row_origin: T) -> Affine<T>;

    ///
    ///
    ///
    fn from_gdal(col_origin: T,
                 width: T,
                 row_rotation: T,
                 row_origin: T,
                 col_rotation: T,
                 height: T) -> Affine<T>;

    /// Transform determinant.
    fn determinant(&self) -> T;

    /// Checks if the transform is degenarate (cannot be inverted).
    fn is_degenerate(&self) -> bool;

    /// Checks if the transform is proper (positive define).
    fn is_proper(&self) -> bool;

    /// Transform adjugate.
    fn adjugate(&self) -> Affine<T>;
}

impl<T> AffineTransform<T> for Affine<T>
    where T: Mul<Output = T> + Add<Output = T> + Sub<Output = T> + Neg<Output = T> + Copy + PartialEq + PartialOrd + Default {
    fn new(width: T,
           row_rotation: T,
           col_origin: T,
           col_rotation: T,
           height: T,
           row_origin: T) -> Affine<T> {
        Affine {
            width: width,
            row_rotation: row_rotation,
            col_origin: col_origin,
            col_rotation: col_rotation,
            height: height,
            row_origin: row_origin
        }
    }

    fn from_gdal(col_origin: T,
                 width: T,
                 row_rotation: T,
                 row_origin: T,
                 col_rotation: T,
                 height: T) -> Affine<T> {
        Affine {
            width: width,
            row_rotation: row_rotation,
            col_origin: col_origin,
            col_rotation: col_rotation,
            height: height,
            row_origin: row_origin
        }
    }

    fn determinant(&self) -> T {
        self.width * self.height - self.row_rotation * self.col_rotation
    }

    fn is_degenerate(&self) -> bool {
        let zero = T::default();
        self.determinant() == zero
    }

    fn is_proper(&self) -> bool {
        let zero = T::default();
        self.determinant() > zero
    }

    fn adjugate(&self) -> Affine<T> {
        // Sign of elements is positive or negative if coeficients sum is pair or impair respectively.
        let width = self.height;
        let row_rotation = -self.row_rotation;
        let col_origin = self.row_rotation * self.row_origin - self.col_origin * self.height;
        let col_rotation = -self.col_rotation;
        let height = self.width;
        let row_origin = self.col_origin * self.col_rotation - self.width * self.row_origin;
        Affine::new(width, row_rotation, col_origin, col_rotation, height, row_origin)
    }
}

/// Scalar multiplication LHS (left-hand side)
macro_rules! impl_multiplication_lhs {
    ($type: ty) => {
        impl Mul<Affine<$type>> for $type {
            type Output = Affine<$type>;

            fn mul(self, matrix: Affine<$type>) -> Affine<$type> {

                Affine::new(self * matrix.width,
                            self * matrix.row_rotation,
                            self * matrix.col_origin,
                            self * matrix.col_rotation,
                            self * matrix.height,
                            self * matrix.row_origin)
            }
        }
    }
}

/// Scalar multiplication RHS (right-hand side)
macro_rules! impl_multiplication_rhs {
    ($type: ty) => {
        impl Mul<$type> for Affine<$type> {
            type Output = Affine<$type>;

            fn mul(self, scalar: $type) -> Affine<$type> {

                Affine::new(scalar * self.width,
                            scalar * self.row_rotation,
                            scalar * self.col_origin,
                            scalar * self.col_rotation,
                            scalar * self.height,
                            scalar * self.row_origin)
            }
        }
    }
}

impl_multiplication_lhs!(f32);
impl_multiplication_lhs!(f64);
impl_multiplication_rhs!(f32);
impl_multiplication_rhs!(f64);

/// Applies the transform
macro_rules! impl_transform {
    ($type: ty) => {
        impl Mul<($type, $type)> for Affine<$type>{
            type Output = ($type, $type);

            /// Transform the coordinates.
            fn mul(self, domain: ($type, $type)) -> ($type, $type) {
                let x = domain.0;
                let y = domain.1;
                let image_x = self.width * x + self.col_origin;
                let image_y = self.height * y + self.row_origin;

                if (self.row_rotation, self.col_rotation) == (0.0, 0.0) {
                    return (image_x, image_y);
                }

                let rotated_x = image_x + self.row_rotation * y;
                let rotated_y = image_y + self.col_rotation * x;
                (rotated_x, rotated_y)
            }
        }
    }
}

impl_transform!(f32);
impl_transform!(f64);

/// Inverse transform
macro_rules! impl_inverse_transform {
    ($type: ty) => {
        impl Not for Affine<$type>
        {
            type Output = Affine<$type>;

            fn not(self) -> Affine<$type> {
                if self.is_degenerate() {
                    // TransformNotInvertibleError
                    panic!("Cannot invert degenerate transform");
                }

                // The reciprocal (inverse) of determinant.
                let inverse_determinant = self.determinant().recip();
                let adjugate = self.adjugate();

                inverse_determinant * adjugate
            }
        }
    }
}

impl_inverse_transform!(f32);
impl_inverse_transform!(f64);

/// Applies the transform
macro_rules! impl_compostion {
    ($type: ty) => {
        impl Mul<Affine<$type>> for Affine<$type>{
            type Output = Affine<$type>;

            /// Transform the coordinates.
            fn mul(self, other: Affine<$type>) -> Affine<$type> {
                let width: $type = self.width * other.width + self.row_rotation * other.col_rotation;
                let row_rotation: $type = self.width * other.row_rotation + self.row_rotation * other.height;
                let col_origin: $type = self.width * other.col_origin + self.row_rotation * other.row_origin + self.col_origin;
                let col_rotation: $type = self.col_rotation * other.width + self.height * other.col_rotation;
                let height: $type = self.col_rotation * other.row_rotation + self.height * other.height;
                let row_origin: $type = self.col_rotation * other.col_origin + self.height * other.row_origin + self.row_origin;
                Affine::new(width, row_rotation, col_origin, col_rotation, height, row_origin)
            }
        }
    }
}

impl_compostion!(f32);
impl_compostion!(f64);

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_new() {
        let affine = Affine::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0);
        let result = Affine { width: 1.0, row_rotation: 2.0, col_origin: 3.0,
                              col_rotation: 4.0, height: 5.0, row_origin: 6.0 };
        assert_eq!(affine, result);
    }

    #[test]
    fn test_diff() {
        let affine = Affine::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0);
        let result = Affine { width: 2.0, row_rotation: 2.0, col_origin: 3.0,
                              col_rotation: 4.0, height: 5.0, row_origin: 6.0 };
        assert_ne!(affine, result);
    }

    #[test]
    fn test_new_f64_items() {
        let affine = Affine::new(1.0f64, 2.0f64, 3.0f64,
                                 4.0f64, 5.0f64, 6.0f64);
        let result = Affine::<f64> { width: 1.0, row_rotation: 2.0, col_origin: 3.0,
                              col_rotation: 4.0, height: 5.0, row_origin: 6.0 };
        assert_eq!(affine, result);
    }

    #[test]
    fn test_new_f64_struct() {
        let affine = Affine::<f64>::new(1.0, 2.0, 3.0,
                                        4.0, 5.0, 6.0);
        let result = Affine::<f64> { width: 1.0, row_rotation: 2.0, col_origin: 3.0,
                              col_rotation: 4.0, height: 5.0, row_origin: 6.0 };
        assert_eq!(affine, result);
    }

    #[test]
    fn test_from_gdal() {
        let affine = Affine::from_gdal(3.0, 1.0, 2.0,
                                       6.0, 4.0, 5.0);
        let result = Affine::<f64> { width: 1.0, row_rotation: 2.0, col_origin: 3.0,
                              col_rotation: 4.0, height: 5.0, row_origin: 6.0 };
        assert_eq!(affine, result);
    }

    #[test]
    fn test_determinant() {
        let affine = Affine::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0);
        let determinant = -3.0;
        let result = affine.determinant();
        assert_eq!(result, determinant);
    }

    #[test]
    fn test_degenerate() {
        let affine = Affine::new(1.0, 2.0, 3.0,
                                 2.0, 4.0, 6.0);
        let result = affine.is_degenerate();
        assert_eq!(result, true);
    }

    #[test]
    fn test_not_degenerate() {
        let affine = Affine::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0);
        let result = affine.is_degenerate();
        assert_eq!(result, false);
    }

    #[test]
    fn test_proper() {
        let affine = Affine::new(1.0, 0.0, 3.0,
                                 0.0, 4.0, 6.0);
        let result = affine.is_proper();
        assert_eq!(result, true);
    }

    #[test]
    fn test_not_proper() {
        let affine = Affine::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0);
        let result = affine.is_proper();
        assert_eq!(result, false);
    }

    #[test]
    fn test_adjugate() {
        let affine = Affine::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0);
        let adjugate = Affine::new(5.0, -2.0, -3.0,
                                   -4.0, 1.0, 6.0);
        let result = affine.adjugate();
        assert_eq!(adjugate, result);
    }


    #[test]
    fn test_left_multiplication() {
        let scalar = 3.0;
        let affine = scalar * Affine::new(1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0);
        let result = Affine::new(3.0, 6.0, 9.0,
                                 12.0, 15.0, 18.0);
        assert_eq!(affine, result);
    }

    #[test]
    fn test_left_multiplication_f64() {
        let scalar = 3.0f64;
        let affine = scalar * Affine::<f64>::new(1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0);
        let result = Affine::<f64>::new(3.0, 6.0, 9.0,
                                 12.0, 15.0, 18.0);
        assert_eq!(affine, result);
    }

    #[test]
    fn test_right_multiplication() {
        let scalar = 4.0;
        let affine = Affine::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0) * scalar;
        let result = Affine::new(4.0, 8.0, 12.0,
                                 16.0, 20.0, 24.0);
        assert_eq!(affine, result);
    }

    #[test]
    fn test_right_multiplication_f64() {
        let scalar = 4.0f64;
        let affine = Affine::<f64>::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0) * scalar;
        let result = Affine::<f64>::new(4.0, 8.0, 12.0,
                                 16.0, 20.0, 24.0);
        assert_eq!(affine, result);
    }

    #[test]
    fn test_transform_direct() {
        let affine = Affine::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0);
        let (row, col) = (2.0, 3.0);
        let (lon, lat) = affine * (row, col);
        assert_eq!(lon, 11.0);
        assert_eq!(lat, 29.0);
    }

    #[test]
    fn test_transform_direct_f64() {
        let affine = Affine::<f64>::new(1.0, 2.0, 3.0,
                                        4.0, 5.0, 6.0);
        let (row, col) = (2.0f64, 3.0f64);
        let (lon, lat) = affine * (row, col);
        assert_eq!(lon, 11.0f64);
        assert_eq!(lat, 29.0f64);
    }

    #[test]
    fn test_transform_direct_ordinary() {
        let affine = Affine::new(0.5, 0.0, -80.0,
                                 0.0, -0.8, -40.0);
        let (row, col) = (2.0, 3.0);
        let (lon, lat) = affine * (row, col);
        assert_eq!(lon, -79.0);
        assert_eq!(lat, -42.4);
    }

    #[test]
    fn test_inverse() {
        let affine = Affine::new(1.0, 0.0, 3.0,
                                 0.0, 5.0, 5.0);
        let inverse = Affine::new(1.0, 0.0, -3.0,
                                  0.0, 0.2, -1.0);
        let result = !affine;
        assert_eq!(inverse, result);
    }

    #[test]
    fn test_transform_inverse() {
        let affine = Affine::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0);
        let (lon, lat) = (11.0, 29.0);
        let (row, col): (f64, f64) = !affine * (lon, lat);
        assert_eq!(row, 2.0);
        assert_eq!(col, 3.0);
    }

    #[test]
    fn test_transform_inverse_ordinary() {
        let affine = Affine::new(0.5, 0.0, -80.0,
                                 0.0, -0.8, -40.0);
        let (lon, lat) = (-78.5, -41.6);
        let (row, col): (f64, f64) = !affine * (lon, lat);
        assert_eq!(row, 3.0);
        assert_eq!(col, 2.0);
    }

    #[test]
    fn test_transform_inverse_center() {
        let affine = Affine::<f32>::new(0.5, 0.0, -80.0,
                                 0.0, -0.8, -40.0);
        let (lon, lat) = (-78.25, -42.0);
        let (row, col): (f32, f32) = !affine * (lon, lat);
        assert_eq!(row.floor(), 3.0);
        assert_eq!(col.floor(), 2.0);
    }

    #[test]
    fn test_transform_inverse_center_f64() {
        let affine = Affine::<f64>::new(0.5, 0.0, -80.0,
                                 0.0, -0.8, -40.0);
        let (lon, lat) = (-78.25, -42.0);
        let (row, col): (f64, f64) = !affine * (lon, lat);
        assert_eq!(row.floor(), 3.0);
        assert_eq!(col.floor(), 2.0);
    }

    #[test]
    #[should_panic]
    fn test_transform_inverse_degenerate() {
        let affine = Affine::new(0.0, 0.0, -80.0,
                                 0.0, -0.8, -40.0);
        let (lon, lat) = (-78.25, -42.0);
        let (row, col): (f64, f64) = !affine * (lon, lat);
        assert_eq!(row, 3.0);
        assert_eq!(col, 2.0);
    }

    #[test]
    fn test_transform_identity() {
        let affine = Affine::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0);
        let identity = Affine::new(1.0, 0.0, 0.0,
                                 0.0, 1.0, 0.0);
        let inverse = !affine;
        let result = affine * inverse;
        assert_eq!(result, identity);
    }
}

