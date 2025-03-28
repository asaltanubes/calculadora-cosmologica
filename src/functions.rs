use core::f64;

use crate::input::ParametrosCosmologicos;
use crate::integral_methods::adaptive_simpson_method;
use crate::{runge_kutta_at_points, trapezoid_method};
use crate::plotting::linspace;


pub fn s_k(xi: f64, k: f64) -> f64{
    if k > 0.{
        return xi.sin()
    }
    if k <0.{
        return xi.sinh()
    }
    
    return xi
}

/// Calcula el valor de E(x) donde x puede ser 1+z o a'/a al pasarse como input
/// si el valor resultaría en un e^2 <0 se devuleve un número negativo
pub fn e(x: f64, parametros: ParametrosCosmologicos) -> f64{
    let arg = parametros.densidad_energia_oscura + parametros.densidad_curvatura*x.powi(2)
   + parametros.densidad_materia*x.powi(3) + parametros.densidad_radiacion*x.powi(4);
    let ret = if arg <0. {-(-arg).sqrt()} else {arg.sqrt()};
    // if ret <=0. {println!("ret. {ret}")}
    return ret
}


pub fn luminosity_distance(z: f64, parametros: ParametrosCosmologicos) -> f64{
    parametros.factor_escala*(1.0+z) * s_k(parametros.c/(parametros.factor_escala * parametros.constante_hubble) * adaptive_simpson_method(|zp: f64| 1./e(1.0+zp, parametros), 0.0, z, 1e-6, 1e-7).unwrap_or(0.), parametros.signo_curvatura)
}

pub fn angular_distance(z: f64, parametros: ParametrosCosmologicos) -> f64{
    (1.0+z).powi(-2) * luminosity_distance(z, parametros)
}

pub fn lookback_time(z:f64, parametros: ParametrosCosmologicos) -> f64{
    adaptive_simpson_method(|zp: f64| 1./((1.+zp)*e(1.+zp, parametros)), z, f64::INFINITY, 1e-20, 1e-7).unwrap_or(f64::INFINITY)/parametros.constante_hubble
    
}

/// Edo que representa el factor de escala. s representa si el universo está en expansión o en contracción o si ya ha acabado
/// esto es importante en universos con curvatura positiva y constante cosmológica suficientemente pequeña o universos
/// con constante cosmológica negativa. Cuando se ha llegado al tamaño máximo del universo lo cual significa que E(a) = 0
/// o en nuestro caso si a es mayor de lo que debería (porque estamos usando un método numérico) obteniendo e(a)<0
/// se cambia el signo de s.
/// además si a = 0 entonces a' = 0 y una solución posible es a=0 para siempre (porque E(t, a) no es lipchitz entorno a a=0 hay varias soluciones).
pub fn edo_para_a(_: f64, a: f64, sig_and_params: &mut (i8, bool, ParametrosCosmologicos)) -> f64{
    let s = &mut sig_and_params.0;
    let last_iteration_changed_sign = &mut sig_and_params.1;
    let parametros = sig_and_params.2;
    let res = e(parametros.factor_escala/a, parametros);
    if  res <= 0. {
        if !*last_iteration_changed_sign{*s*=-1;}
        *last_iteration_changed_sign=true;
    }
    else {*last_iteration_changed_sign=false}
        if a < 1e-8{*s = 0}
        *s as f64*a*parametros.constante_hubble*res.abs()
}


// pub fn z(tp: f64, eps: f64, parametros: ParametrosCosmologicos) -> f64{
//     return newton_method(|z| t(z, parametros), |zp| (1.+zp) * e(1.+zp, parametros), tp, 0., eps);
// }

// pub fn t(z: f64, parametros: ParametrosCosmologicos) -> f64{
//     adaptive_simpson_method(|zp| 1./((1.+zp) * e(1.+zp, parametros)), z, 0., 1e-10, 1e-7).unwrap_or(0.)/parametros.constante_hubble
// }

pub fn horizonte_de_partículas(t: f64, parametros: ParametrosCosmologicos) -> (f64, Vec<f64>, Vec<f64>){
    
    let universe_lifetime = lookback_time(0., parametros);
    if universe_lifetime.is_finite(){
        if t <=0.{
            let elems = 5000;
            let mut tt_at = linspace(0., t, elems);
            let mut aa_at = runge_kutta_at_points(edo_para_a, &tt_at, parametros.factor_escala, 10, (1, true, parametros));
            let mut tt_bt = linspace(t, -universe_lifetime, elems);
            let mut aa_bt = runge_kutta_at_points(edo_para_a, &tt_bt, aa_at[elems as usize-1], 10, (1, true, parametros));
            tt_at.reverse();
            aa_at.reverse(); 
            tt_bt.reverse();
            aa_bt.reverse();
            aa_bt.extend_from_slice(&aa_at[1..]);
            tt_bt.extend_from_slice(&tt_at[1..]);
            let func_bt: Vec<_> = aa_bt.iter().map(|a| 1./a).collect();
            let mut distance = vec![0.; 2*elems as usize-1];
            for i in 1..tt_bt.len(){
                distance[i] = aa_bt[i]*parametros.c*trapezoid_method(&tt_bt[i-1..=i], &func_bt[i-1..=i]) + distance[i-1];
            }
            return (distance[elems as usize-1], tt_bt, distance);
        }
        else {
            let elems = 5000;
            let tt_at = linspace(0., t, elems);
            let aa_at = runge_kutta_at_points(edo_para_a, &tt_at, parametros.factor_escala, 10, (1, true, parametros));
            let mut tt_bt = linspace(0., -universe_lifetime, elems);
            let mut aa_bt = runge_kutta_at_points(edo_para_a, &tt_bt, parametros.factor_escala, 10, (1, true, parametros));
            tt_bt.reverse();
            aa_bt.reverse();
            aa_bt.extend_from_slice(&aa_at[1..]);
            tt_bt.extend_from_slice(&tt_at[1..]);
            let func_bt: Vec<_> = aa_bt.iter().map(|a| 1./a).collect();
            let mut distance = vec![0.; 2*elems as usize-1];
            for i in 1..tt_bt.len(){
                distance[i] = aa_bt[i]*parametros.c*trapezoid_method(&tt_bt[0..=i], &func_bt[0..=i]);
            }
            return (distance[2*elems as usize-2], tt_bt, distance);
       }


    }
    return (f64::INFINITY, Vec::new(), Vec::new());
}