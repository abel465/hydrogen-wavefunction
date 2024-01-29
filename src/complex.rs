use core::convert::From;
use core::ops::*;

#[derive(Copy, Clone, PartialEq)]
pub struct Complex {
    pub x: f64,
    pub y: f64,
}

impl Complex {
    pub fn new(x: f64, y: f64) -> Self {
        Complex { x, y }
    }
    pub const ZERO: Complex = Complex { x: 0.0, y: 0.0 };
    pub const ONE: Complex = Complex { x: 1.0, y: 0.0 };
    pub const I: Complex = Complex { x: 0.0, y: 1.0 };
}

impl Complex {
    pub fn conjugate(&self) -> Self {
        Self::new(self.x, -self.y)
    }

    pub fn powf(self, exp: f64) -> Self {
        let (r, theta) = self.to_polar();
        Self::from_polar(r.powf(exp), theta * exp)
    }

    pub fn norm(self) -> f64 {
        self.norm_squared().sqrt()
    }

    pub fn norm_squared(self) -> f64 {
        (self * self.conjugate()).x
    }

    pub fn arg(self) -> f64 {
        self.y.atan2(self.x)
    }

    pub fn to_polar(self) -> (f64, f64) {
        (self.norm(), self.arg())
    }

    pub fn from_polar(r: f64, theta: f64) -> Self {
        Self::new(r * theta.cos(), r * theta.sin())
    }

    pub fn sqrt(self) -> Self {
        Self::new(
            ((self.norm() + self.x) / 2.0).sqrt(),
            self.y.signum() * ((self.norm() - self.x) / 2.0).sqrt(),
        )
    }

    pub fn exp(self) -> Self {
        Self::from_polar(self.x.exp(), self.y)
    }
}

impl Add for Complex {
    type Output = Self;
    fn add(self, other: Self) -> Self::Output {
        Complex::new(self.x + other.x, self.y + other.y)
    }
}

impl Sub for Complex {
    type Output = Self;
    fn sub(self, other: Self) -> Self::Output {
        Complex::new(self.x - other.x, self.y - other.y)
    }
}

impl Mul for Complex {
    type Output = Self;
    fn mul(self, other: Self) -> Self::Output {
        Complex::new(
            self.x * other.x - self.y * other.y,
            self.x * other.y + self.y * other.x,
        )
    }
}

impl Mul<f64> for Complex {
    type Output = Self;
    fn mul(self, other: f64) -> Self::Output {
        Complex::new(self.x * other, self.y * other)
    }
}

impl Mul<Complex> for f64 {
    type Output = Complex;
    fn mul(self, other: Complex) -> Self::Output {
        Complex::new(self * other.x, self * other.y)
    }
}

impl Div<f64> for Complex {
    type Output = Self;
    fn div(self, other: f64) -> Self::Output {
        Complex::new(self.x / other, self.y / other)
    }
}

impl Div<Complex> for Complex {
    type Output = Self;
    fn div(self, other: Complex) -> Self::Output {
        let d = other.x * other.x + other.y * other.y;
        Complex::new(
            (self.x * other.x + self.y * other.y) / d,
            (self.y * other.x - self.x * other.y) / d,
        )
    }
}

impl Neg for Complex {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Complex::new(-self.x, -self.y)
    }
}

impl PartialEq<f64> for Complex {
    fn eq(&self, other: &f64) -> bool {
        self.y == 0.0 && self.x == *other
    }
}

impl std::fmt::Debug for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.y == 0.0 {
            write!(f, "{}", self.x)
        } else if self.y < 0.0 {
            write!(f, "{} - {}i", self.x, -self.y)
        } else {
            write!(f, "{} + {}i", self.x, self.y)
        }
    }
}

impl From<f64> for Complex {
    fn from(value: f64) -> Complex {
        Complex::new(value, 0.0)
    }
}
