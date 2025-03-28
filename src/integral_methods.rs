use core::f64;
use std::fmt;

/// Versión modificada del crate integrate (https://docs.rs/integrate/0.1.4/integrate/index.html) para permitir aceptar closures (lambda functions en python). 
/// también permite calcular integrales indefinidas con límites infinitos con un cambio de variable.
//---------------------------------------------------------
#[derive(Clone, Debug)]
struct SubInterval<F> {
    upper_limit: F,
    lower_limit: F,
    function: [F; 5],
    interval: Option<Box<SubInterval<F>>>,
}

type Result<T> = std::result::Result<T, AdaptiveSimpsonError>;

#[derive(Debug, Clone)]
pub struct AdaptiveSimpsonError;

impl fmt::Display for AdaptiveSimpsonError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let msg = "No subinterval of length > min_h was found for which the estimated error was less that the pro-rated tolerance";
        write!(f, "{}", msg)
    }
}


pub fn adaptive_simpson_method(
    f: impl Fn(f64) -> f64,
    a: f64,
    b: f64,
    min_h: f64,
    tolerance: f64,
) -> Result<f64> {
    // Si alguno de los límites es infinito arctan(inf) = pi/2
    // y la integral ya es definida
    if a.is_infinite() || b.is_infinite(){
        return helper_func(|x| f(x.tan())/x.cos().powi(2), a.atan(), b.atan(), min_h, tolerance)
    }
    return helper_func(f, a, b, min_h, tolerance);
    fn helper_func(f: impl Fn(f64) -> f64,
    a: f64,
    b: f64,
    min_h: f64,
    tolerance: f64,) -> Result<f64>{
    let two = 2.;

    let mut integral: f64 = 0.0;
    let epsilon_density = two * tolerance / (b - a);

    // Create the initial level, with lower_limit = a, upper_limit = b,
    // and f(x) evaluated at a, b, and (a + b) / 2.

    let interval: SubInterval<f64> = SubInterval {
        upper_limit: b,
        lower_limit: a,
        function: [f(a), f64::NAN, f((a + b) / two), f64::NAN, f(b)],
        interval: None,
    };

    let mut pinterval = Box::new(interval);

    // Calculate the tolerance for the current interval.
    // calculate the single subinterval Simpson rule,
    // and the two subintervals composite Simpson rule.

    let mut epsilon = epsilon_density * (b - a);
    let (mut s1, mut s2) = simpson_rule_update(&f, &mut pinterval);

    let mut qinterval: SubInterval<f64>;

    while pinterval.upper_limit - pinterval.lower_limit > min_h {
        if (s1 - s2).abs() < epsilon {
            // If the two estimates are close, then increment the
            // integral and if we are not at the right end, set the
            // left end of the new interval to the right end of the
            // old interval and the right end of the new interval
            // remains the same (as the previous right end for this
            // interval.

            integral += s2;

            if pinterval.interval.is_none() {
                return Ok(integral);
            }

            // Move to the next interval
            qinterval = *pinterval.interval.take().unwrap();
            qinterval.lower_limit = pinterval.upper_limit;
            qinterval.function[0] = qinterval.function[2];
            qinterval.function[2] = qinterval.function[3];

            pinterval = Box::new(qinterval);
        } else {
            // If the two estimates are not close, then create a new
            // interval with same left end point and right end point
            // at the midpoint of the current interval.

            let limit1 = pinterval.lower_limit;
            let limit2 = (pinterval.upper_limit + pinterval.lower_limit) / two;

            let upper_limit = if limit1 > limit2 { limit1 } else { limit2 };
            let lower_limit = if limit1 > limit2 { limit2 } else { limit1 };

            qinterval = SubInterval {
                lower_limit,
                upper_limit,
                function: [f64::NAN; 5],
                interval: None,
            };

            qinterval.function[0] = pinterval.function[0];
            qinterval.function[2] = pinterval.function[1];
            qinterval.function[4] = pinterval.function[2];

            qinterval.interval = Some(pinterval);

            pinterval = Box::new(qinterval);
        }

        // Update Simpson's rule for the new interval
        (s1, s2) = simpson_rule_update(&f, &mut pinterval);
        epsilon = epsilon_density * (pinterval.upper_limit - pinterval.lower_limit);
    }
    Err(AdaptiveSimpsonError)
}
}

fn simpson_rule_update(
    f: &dyn Fn(f64) -> f64,
    pinterval: &mut SubInterval<f64>,
) -> (f64, f64) {
    let two = 2.;
    let four = two + two;
    let six = four + two;

    let h = pinterval.upper_limit - pinterval.lower_limit;
    let h4 = h / four;

    pinterval.function[1] = f(pinterval.lower_limit + h4);
    pinterval.function[3] = f(pinterval.upper_limit - h4);

    let mut s1 = pinterval.function[0] + four * pinterval.function[2] + pinterval.function[4];
    s1 *= h / six;

    let mut s2 = pinterval.function[0]
        + four * pinterval.function[1]
        + two * pinterval.function[2]
        + four * pinterval.function[3]
        + pinterval.function[4];
    s2 *= h / (six * two);

    (s1, s2)
}

//---------------------------------------

pub fn runge_kutta<T>(f: impl Fn(f64, f64, &mut T) -> f64, delta_x: f64, initial_x: f64, initial_condition: f64, steps: u64, params: &mut T) -> f64{
    if initial_condition.is_nan(){
        return f64::NAN;
    }
    let h = delta_x/steps as f64;
    let half_h = h/2.0;
    let mut y = initial_condition;
    let mut k1;
    let mut k2;
    let mut k3;
    let mut k4;
    let mut x = initial_x;
    for i in 0..steps{
        k1 = f(x, y, params);
        k2 = f(x+half_h, y+half_h*k1, params);
        k3 = f(x+half_h, y+half_h*k2, params);
        k4 = f(x+h, y+k3*h, params);
        
        y = y + h*(k1 +2.*k2+2.*k3+k4)/6.;
        x = initial_x + i as f64*h;
    }
    y
}


pub fn runge_kutta_at_points<T>(f: impl Fn(f64, f64, &mut T) -> f64, xx: &Vec<f64>, initial_condition: f64, sub_steps: u64, mut params: T) -> Vec<f64>{
    let mut yy = Vec::with_capacity(xx.len());
    yy.push(initial_condition);
    for i in 1..xx.len(){
        yy.push(runge_kutta(&f, xx[i]-xx[i-1], xx[i-1], yy[i-1], sub_steps, &mut params));
    }
    yy
}

pub fn trapezoid_method(xx: &[f64], yy: &[f64]) -> f64{
    let dx: Vec<_> = xx.windows(2).map(|sl| sl[1]-sl[0]).collect();
    let parl: Vec<_> = yy.windows(2).map(|sl| (sl[0]+sl[1])/2.).collect();
    dx.into_iter().zip(parl.into_iter()).map(|(a, b)| a*b).sum()
}
