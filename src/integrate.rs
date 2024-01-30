pub fn integrate1(f: &dyn Fn(f64) -> f64, a: f64, b: f64, p: u32) -> f64 {
    if p == 0 || (a - b).abs() == 0.0 {
        return 0.0;
    }

    let delta = (b - a) / p as f64;

    let mut integral = f(a) + f(b);
    let mut pos = a;
    for _ in 1..p {
        pos += delta;
        integral += f(pos);
    }

    integral * delta
}

pub fn integrate2(f: &dyn Fn(f64, f64) -> f64, a: [f64; 2], b: [f64; 2], p: u32) -> f64 {
    if p == 0 || a.iter().zip(b).any(|(a, b)| (a - b).abs() == 0.0) {
        return 0.0;
    }

    let delta = [(b[0] - a[0]) / p as f64, (b[1] - a[1]) / p as f64];

    let mut integral = f(a[0], a[1]) + f(b[0], b[1]);
    let mut pos = [a[0], a[1]];
    for _ in 1..p {
        pos[0] += delta[0];
        for _ in 1..p {
            pos[1] += delta[1];
            integral += f(pos[0], pos[1]);
        }
        pos[1] = a[1];
    }

    integral * delta[0] * delta[1]
}

pub fn integrate3(f: &dyn Fn(f64, f64, f64) -> f64, a: [f64; 3], b: [f64; 3], p: u32) -> f64 {
    if p == 0 || a.iter().zip(b).any(|(a, b)| (a - b).abs() == 0.0) {
        return 0.0;
    }

    let delta = [
        (b[0] - a[0]) / p as f64,
        (b[1] - a[1]) / p as f64,
        (b[2] - a[2]) / p as f64,
    ];

    let mut integral = f(a[0], a[1], a[2]) + f(b[0], b[1], b[2]);
    let mut pos = [a[0], a[1], a[2]];
    for _ in 1..p {
        pos[0] += delta[0];
        for _ in 1..p {
            pos[1] += delta[1];
            for _ in 1..p {
                pos[2] += delta[2];
                integral += f(pos[0], pos[1], pos[2]);
            }
            pos[2] = a[2];
        }
        pos[1] = a[1];
    }

    integral * delta[0] * delta[1] * delta[2]
}
