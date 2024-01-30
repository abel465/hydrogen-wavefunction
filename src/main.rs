mod complex;
pub mod integrate;

use complex::Complex;
use core::f64::consts::PI;

const A: f64 = 1.00054 * 5.29177210903e-11; // Bohr radius

fn factorialu(n: u32) -> f64 {
    if n > 12 {
        panic!("");
    }
    [
        1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600,
    ][n as usize] as f64
}

fn laguerre_polynomial(r: u32, s: u32, x: f64) -> f64 {
    let mut sum = 0.0;
    for q in 0..=s {
        sum += (-1.0_f64).powi(q as i32) * factorialu(s + r) * factorialu(s + r) * x.powi(q as i32)
            / (factorialu(s - q) * factorialu(r + q) * factorialu(q));
    }
    sum
}

fn binomial(n: u32, k: u32) -> f64 {
    let mut x = 1.0;
    for i in 1..=k {
        x *= (n + 1 - i) as f64 / i as f64;
    }
    x
}

fn general_binomial(n: f64, k: u32) -> f64 {
    let mut x = 1.0;
    for i in 0..k {
        x *= n - i as f64;
    }
    x / factorialu(k)
}

fn legendre_polynomial(m: i32, mut l: i32, x: f64) -> Complex {
    fn legendre_polynomial_positive(m: u32, l: u32, x: f64) -> Complex {
        let mut sm = 0.0;
        for k in m..=l {
            sm += factorialu(k) / factorialu(k - m)
                * x.powi((k - m) as i32)
                * binomial(l, k)
                * general_binomial(((l + k) as f64 - 1.0) / 2.0, l);
        }
        let bb = Complex::new(1.0 - x * x, 0.0).powf(m as f64 / 2.0);
        (-1.0_f64).powi(m as i32) * 2.0_f64.powi(l as i32) * bb * sm
    }
    if l < 0 {
        l = -l - 1;
    }
    if m < 0 {
        (-1.0_f64).powi(-m) * factorialu((l + m) as u32) / factorialu((l - m) as u32)
            * legendre_polynomial_positive((-m) as u32, l as u32, x)
    } else {
        legendre_polynomial_positive(m as u32, l as u32, x)
    }
}

fn spherical_harmonic(m: i32, l: i32, theta: f64, phi: f64) -> Complex {
    let normalization_constant = (((2 * l + 1) as f64 * factorialu((l - m) as u32))
        / (4.0 * PI * factorialu((l + m) as u32)))
    .sqrt();
    let angular = (Complex::I * phi * m as f64).exp();
    let lp = legendre_polynomial(m, l, theta.cos());
    normalization_constant * lp * angular
}

fn radial_wavefunction(n: u32, l: u32, r: f64) -> f64 {
    let p = (2.0 * r) / (n as f64 * A);
    let normalization_constant = ((2.0 / (n as f64 * A)).powi(3) * factorialu(n - l - 1)
        / (2.0 * n as f64 * factorialu(n + l).powi(3)))
    .sqrt();
    let asymptotic_forms = (-r / (n as f64 * A)).exp() * p.powi(l as i32);
    let lp = laguerre_polynomial(2 * l + 1, n - l - 1, p);
    normalization_constant * asymptotic_forms * lp
}

fn hydrogen_wavefunction(n: u32, l: u32, m: i32, theta: f64, phi: f64, r: f64) -> Complex {
    let radial = radial_wavefunction(n, l, r);
    let angular = spherical_harmonic(m, l as i32, theta, phi);
    radial * angular
}

fn main() {
    println!("{:?}", hydrogen_wavefunction(1, 0, 0, 0.0, 0.0, A));
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;
    use integrate::*;

    macro_rules! assert_similar {
        ($x:expr, $y:expr) => {
            let delta = ($x - $y).abs();
            let e = ($x * f64::EPSILON).abs().max(($y * f64::EPSILON).abs());
            if delta > e {
                panic!(
                    "assertion `left ~= right` failed
  left: {}
 right: {}",
                    $x, $y
                )
            }
        };
    }

    #[test]
    fn test_hydrogen_wavefunction_integral() {
        for n in 0..3 {
            for l in 0..n {
                for m in 0..=l {
                    let f = |x: f64, y: f64, z: f64| {
                        let r = (x * x + y * y + z * z).sqrt();
                        let theta = (z / r).acos();
                        let phi = y.signum() * (x / (x * x + y * y).sqrt()).acos();
                        let v = hydrogen_wavefunction(n, l, m as i32, theta, phi, r);
                        v.norm_squared()
                    };
                    let d = A * 50.0;
                    let total_probability = integrate3(&f, [-d, -d, -d], [d, d, d], 140);
                    assert_approx_eq!(total_probability, 1.0, 0.05);
                }
            }
        }
    }

    #[test]
    fn test_radial_wavefunction_integral() {
        for n in 1..6 {
            for l in 0..n {
                let f = |r: f64| {
                    let v = radial_wavefunction(n, l, r);
                    v * v * r * r
                };
                let total_probability = integrate1(&f, 0.0, A * 100.0, 10000);
                assert_approx_eq!(total_probability, 1.0, 1e-5);
            }
        }
    }

    #[test]
    fn test_radial_wavefunction() {
        for r in [
            0.0,
            A / 99.9,
            A / 3.0,
            A / 2.0,
            A,
            2.0 * A,
            4.0 * A,
            40.0 * A,
            999.9 * A,
        ] {
            assert_eq!(
                radial_wavefunction(1, 0, r),
                2.0 * A.powf(-3.0 / 2.0) * (-r / A).exp()
            );
            assert_similar!(
                radial_wavefunction(2, 0, r),
                1.0 / 2.0_f64.sqrt()
                    * A.powf(-3.0 / 2.0)
                    * (1.0 - r / (2.0 * A))
                    * (-r / (2.0 * A)).exp()
            );
            assert_similar!(
                radial_wavefunction(2, 1, r),
                1.0 / 24.0_f64.sqrt() * A.powf(-3.0 / 2.0) * (r / A) * (-r / (2.0 * A)).exp()
            );
        }
    }

    #[test]
    fn test_spherical_harmonic_integral() {
        for l in 1..4 {
            for m in -l..=l {
                let f = |theta: f64, phi: f64| {
                    let v = spherical_harmonic(m, l, theta, phi);
                    v.norm_squared() * theta.sin()
                };
                let total_probability = integrate2(&f, [0.0, 0.0], [PI, 2.0 * PI], 1000);
                assert_approx_eq!(total_probability, 1.0, 0.002);
            }
        }
    }

    #[test]
    fn test_spherical_harmonic() {
        let i = Complex::I;
        for theta in [
            0.0,
            PI / 8.0,
            2.0 * PI / 8.0,
            3.0 * PI / 8.0,
            4.0 * PI / 8.0,
            5.0 * PI / 8.0,
            6.0 * PI / 8.0,
            7.0 * PI / 8.0,
            PI,
        ] {
            for phi in [
                -PI,
                -7.0 * PI / 8.0,
                -6.0 * PI / 8.0,
                -5.0 * PI / 8.0,
                -4.0 * PI / 8.0,
                -3.0 * PI / 8.0,
                -2.0 * PI / 8.0,
                -PI / 8.0,
                0.0,
                PI / 8.0,
                2.0 * PI / 8.0,
                3.0 * PI / 8.0,
                4.0 * PI / 8.0,
                5.0 * PI / 8.0,
                6.0 * PI / 8.0,
                7.0 * PI / 8.0,
                PI,
            ] {
                assert_eq!(
                    spherical_harmonic(0, 0, theta, phi),
                    0.5 * (1.0 / PI).sqrt()
                );
                assert_eq!(
                    spherical_harmonic(0, 1, theta, phi),
                    0.5 * (3.0 / PI).sqrt() * theta.cos()
                );
                {
                    let result = spherical_harmonic(0, 2, theta, phi);
                    let truth = Complex::from(
                        0.25 * (5.0 / PI).sqrt() * (3.0 * theta.cos() * theta.cos() - 1.0),
                    );
                    assert_approx_eq!(result.x, truth.x, f64::EPSILON);
                    assert_approx_eq!(result.y, truth.y, f64::EPSILON);
                }
                {
                    let result = spherical_harmonic(1, 1, theta, phi);
                    let truth = -0.5 * (3.0 / (2.0 * PI)).sqrt() * theta.sin() * (i * phi).exp();
                    assert_approx_eq!(result.x, truth.x, f64::EPSILON);
                    assert_approx_eq!(result.y, truth.y, f64::EPSILON);
                }
                {
                    let result = spherical_harmonic(1, 2, theta, phi);
                    let truth = -0.5
                        * (15.0 / (2.0 * PI)).sqrt()
                        * (theta).sin()
                        * (theta).cos()
                        * (i * phi).exp();
                    assert_approx_eq!(result.x, truth.x, f64::EPSILON);
                    assert_approx_eq!(result.y, truth.y, f64::EPSILON);
                }
                {
                    let result = spherical_harmonic(2, 2, theta, phi);
                    let truth = 0.25
                        * (15.0 / (2.0 * PI)).sqrt()
                        * (theta).sin()
                        * (theta).sin()
                        * (2.0 * i * phi).exp();
                    assert_approx_eq!(result.x, truth.x, f64::EPSILON);
                    assert_approx_eq!(result.y, truth.y, f64::EPSILON);
                }
                {
                    let result = spherical_harmonic(-1, 1, theta, phi);
                    let truth = 0.5 * (3.0 / (2.0 * PI)).sqrt() * theta.sin() * (-i * phi).exp();
                    assert_approx_eq!(result.x, truth.x, f64::EPSILON);
                    assert_approx_eq!(result.y, truth.y, f64::EPSILON);
                }
                {
                    let result = spherical_harmonic(-1, 2, theta, phi);
                    let truth = 0.5
                        * (15.0 / (2.0 * PI)).sqrt()
                        * theta.sin()
                        * theta.cos()
                        * (-i * phi).exp();
                    assert_approx_eq!(result.x, truth.x, f64::EPSILON);
                    assert_approx_eq!(result.y, truth.y, f64::EPSILON);
                }
                {
                    let result = spherical_harmonic(-2, 2, theta, phi);
                    let truth = 0.25
                        * (15.0 / (2.0 * PI)).sqrt()
                        * theta.sin()
                        * theta.sin()
                        * (-2.0 * i * phi).exp();
                    assert_approx_eq!(result.x, truth.x, f64::EPSILON);
                    assert_approx_eq!(result.y, truth.y, f64::EPSILON);
                }
            }
        }
    }

    #[test]
    fn test_legendre_polynomial() {
        for x in [-1.4, -1.0, 0.0, 1.0, 1.4] {
            {
                let result = legendre_polynomial(0, 0, x);
                let truth = Complex::ONE;
                assert_approx_eq!(result.x, truth.x, f64::EPSILON);
                assert_approx_eq!(result.y, truth.y, f64::EPSILON);
            }
            {
                let result = legendre_polynomial(0, 1, x);
                let truth = x * Complex::ONE;
                assert_approx_eq!(result.x, truth.x, f64::EPSILON);
                assert_approx_eq!(result.y, truth.y, f64::EPSILON);
            }
            {
                let result = legendre_polynomial(1, 1, x);
                let truth = -Complex::new(1.0 - x * x, 0.0).sqrt();
                assert_approx_eq!(result.x, truth.x, f64::EPSILON);
                assert_approx_eq!(result.y, truth.y, f64::EPSILON);
            }
            {
                let result = legendre_polynomial(-1, 1, x);
                let truth = -0.5 * legendre_polynomial(1, 1, x);
                assert_approx_eq!(result.x, truth.x, f64::EPSILON);
                assert_approx_eq!(result.y, truth.y, f64::EPSILON);
            }
            {
                let result = legendre_polynomial(-1, 2, x);
                let truth = -1.0 / 6.0 * legendre_polynomial(1, 2, x);
                assert_approx_eq!(result.x, truth.x, f64::EPSILON);
                assert_approx_eq!(result.y, truth.y, f64::EPSILON);
            }
            {
                let result = legendre_polynomial(-2, 2, x);
                let truth = 1.0 / 24.0 * legendre_polynomial(2, 2, x);
                assert_approx_eq!(result.x, truth.x, f64::EPSILON);
                assert_approx_eq!(result.y, truth.y, f64::EPSILON);
            }
            {
                let result = legendre_polynomial(0, 2, x);
                let truth = Complex::ONE * 0.5 * (3.0 * x * x - 1.0);
                assert_approx_eq!(result.x, truth.x, f64::EPSILON);
                assert_approx_eq!(result.y, truth.y, f64::EPSILON);
            }
            {
                let result = legendre_polynomial(1, 2, x);
                let truth = -3.0 * x * Complex::from(1.0 - x * x).sqrt();
                assert_approx_eq!(result.x, truth.x, f64::EPSILON * 2.0);
                assert_approx_eq!(result.y, truth.y, f64::EPSILON);
            }
            {
                let result = legendre_polynomial(2, 2, x);
                let truth = 3.0 * Complex::from(1.0 - x * x);
                assert_approx_eq!(result.x, truth.x, f64::EPSILON);
                assert_approx_eq!(result.y, truth.y, f64::EPSILON * 3.0);
            }
        }
    }

    #[test]
    fn test_laguerre_polynomial() {
        for x in [-1.4, -1.0, 0.0, 1.0, 1.4] {
            assert_eq!(laguerre_polynomial(0, 0, x), 1.0);
            assert_eq!(laguerre_polynomial(0, 1, x), 1.0 - x);
            assert_eq!(laguerre_polynomial(1, 0, x), 1.0);
            assert_eq!(laguerre_polynomial(1, 1, x), -2.0 * x + 4.0);
            assert_approx_eq!(
                laguerre_polynomial(0, 2, x),
                x * x - 4.0 * x + 2.0,
                f64::EPSILON * 2.0
            );
            assert_approx_eq!(
                laguerre_polynomial(1, 2, x),
                3.0 * x * x - 18.0 * x + 18.0,
                f64::EPSILON * 5.0
            );
            assert_approx_eq!(laguerre_polynomial(2, 1, x), -6.0 * x + 18.0);
        }
    }

    #[test]
    fn test_binomial() {
        assert_eq!(binomial(0, 1), 0.0);
        assert_eq!(binomial(0, 0), 1.0);
        assert_eq!(binomial(1, 0), 1.0);
        assert_eq!(binomial(2, 0), 1.0);
        assert_eq!(binomial(1, 1), 1.0);
        assert_eq!(binomial(1, 1), 1.0);
        assert_eq!(binomial(2, 1), 2.0);
        assert_eq!(binomial(2, 2), 1.0);
        assert_eq!(binomial(4, 2), 6.0);
    }

    #[test]
    fn test_general_binomial() {
        for (n, k) in [
            (0, 1),
            (0, 0),
            (1, 0),
            (2, 0),
            (1, 1),
            (1, 1),
            (2, 1),
            (2, 2),
            (4, 2),
        ] {
            assert_eq!(general_binomial(n as f64, k), binomial(n, k));
        }
    }

    #[test]
    fn test_factorial() {
        assert_eq!(factorialu(0), 1.0);
        assert_eq!(factorialu(1), 1.0);
        assert_eq!(factorialu(2), 2.0);
        assert_eq!(factorialu(3), 6.0);
        assert_eq!(factorialu(4), 24.0);
        assert_eq!(factorialu(5), 120.0);
        assert_eq!(factorialu(6), 720.0);
    }
}
