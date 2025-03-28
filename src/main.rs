
use core::f64;
use std::fs;

use anyhow::Result;
use toml;
// use gnuplot::{AxesCommon, Caption, Color, Figure, LabelOption};
// use plotters::prelude::*;
// use full_palette::ORANGE;
use std::process::Command;


mod input;
mod functions;
mod integral_methods;
mod consts;
mod plotting;


use functions::*;
use input::*;
use integral_methods::*;
use consts::*;
use plotting::*;



// const C: f64 = 1.0;
// const FONT: &str = "New Computer Modern";

fn main() -> Result<()>{
    let parametros: ParametrosCosmologicos = toml::from_str::<InputFile>(fs::read_to_string("input_data.toml")?.as_str())?.into(); 
    println!("{:?}", parametros);
    println!("1lc =  {:e} m", C_GY/parametros.c);
    println!("c: {:e} leguas cósmicas/giga año", parametros.c);
    println!("k: {:e}", parametros.signo_curvatura);
    println!("edad actual del universo: {} giga años", lookback_time(0., parametros));


    let zz: Vec<f64> = linspace(0., parametros.z_max_dist, 1000);
    let d1: Vec<_> = zz.iter().map(|&z| luminosity_distance(z, parametros)*parametros.legua_cosmica_a_mpc).collect();
    let d2: Vec<_> = zz.iter().map(|&z| angular_distance(z ,parametros)*parametros.legua_cosmica_a_mpc).collect();


    let py_program = format!(r#"
import numpy as np
import matplotlib.pyplot as plt
plt.plot({xx}, {y1}, c="tab:orange", label="Distancia luminosidad")
plt.plot({xx}, {y2}, c="tab:blue", label="Distancia angular")
plt.xlabel("z", fontsize=14)
plt.ylabel("D[Mpc]", fontsize=14)
plt.legend()
plt.yscale("log")
plt.savefig("plots/distancias.svg")
    "#, xx=format_as_list(&zz), y1 = format_as_list(&d1), y2=format_as_list(&d2));
    std::fs::write("plots/distancias.py", py_program)?;
    let mut distancias = Command::new("py").arg("plots/distancias.py").spawn()?;

    let mut tt: Vec<f64>;
    let mut aa: Vec<f64>;
    if parametros.rango_factor_escala.0 *parametros.rango_factor_escala.1 <0.{
        tt = linspace(0., parametros.rango_factor_escala.0, 1000);
        aa = runge_kutta_at_points(edo_para_a, &tt, parametros.factor_escala, 10, (1, false, parametros));
        tt.reverse();
        aa.reverse();
        let tt2: Vec<f64> = linspace(0., parametros.rango_factor_escala.1, 1000);
        let aa2: Vec<_> = runge_kutta_at_points(edo_para_a, &tt2, parametros.factor_escala, 10, (1, false, parametros));
        tt.extend_from_slice(&tt2[1..]);
        aa.extend_from_slice(&aa2[1..]);
    }
    else{
        let extr = parametros.rango_factor_escala.0.signum() *  parametros.rango_factor_escala.0.abs().max(parametros.rango_factor_escala.1.abs());
        tt = linspace(0., extr, 1000);
        aa = runge_kutta_at_points(edo_para_a, &tt, parametros.factor_escala, 10, (1, false, parametros));
    }

    let py_program = format!(r#"
import numpy as np
import matplotlib.pyplot as plt
plt.plot({xx}, {yy})
plt.ticklabel_format(style='plain')
plt.ylabel("a/lc", fontsize=14)
plt.xlabel("$t$/Giga años", fontsize=14)
plt.savefig("plots/factor_escala.svg")
    "#, xx=format_as_list(&tt), yy = format_as_list(&aa));
    std::fs::write("plots/factor_escala.py", py_program)?;
    let mut escala = Command::new("py").arg("plots/factor_escala.py").spawn()?;

    let hh: Vec<_> = aa.as_slice().windows(2).map(|sl|  if sl[1]>sl[0] {1.} else {-1.}*parametros.constante_hubble *e(parametros.factor_escala/sl[0], parametros)*S_TO_GY*3.08567758e19).collect();
    let py_program = format!(r#"
import numpy as np
import matplotlib.pyplot as plt
plt.plot({xx}, {yy})
plt.ylabel("H/km/s/Mpc", fontsize=14)
plt.xlabel("$t$/Giga años", fontsize=14)
# plt.ticklabel_format(style='sci', axis='y', scilimits=(5,-5))
plt.yscale("log")
plt.savefig("plots/constante_hubble.svg")
    "#, xx=format_as_list(&tt[0..tt.len()-1]), yy = format_as_list(&hh));
    std::fs::write("plots/constante_hubble.py", py_program)?;
    let mut hubble = Command::new("py").arg("plots/constante_hubble.py").spawn()?;
    
    // let xx: Vec<f64> = linspace(parametros.rango_radio_hubble.0, parametros.rango_radio_hubble.1, 10000);
    let yy: Vec<_> = hh.iter().map(|a| parametros.c/a*parametros.legua_cosmica_a_mpc).collect();

    let py_program = format!(r#"
import numpy as np
import matplotlib.pyplot as plt
plt.plot({xx}, {yy})
plt.ticklabel_format(style='sci', axis='y', scilimits=(-3, 3))
plt.ylabel(r"$D_{{\rm H}}$/Mpc", fontsize=14)
plt.xlabel("$t$/Giga años", fontsize=14)
plt.savefig("plots/radio de hubble.svg")
    "#, xx=format_as_list(&tt[0..tt.len()-1]), yy = format_as_list(&yy));
    std::fs::write("plots/radio_hubble.py", py_program)?;
    let mut hubble_rad = Command::new("py").arg("plots/radio_hubble.py").spawn()?;
 

    // println!("root 2: {}", newton_method(|x| x.sqrt(), |x| 2.*x.sqrt(), 2., 1., 1e-6));
    let (dist, tt, dd) = horizonte_de_partículas(parametros.t_horizonte_particulas, parametros);

    println!("horizonte de particulas en t={} Giga años es: {} Mpc", parametros.t_horizonte_particulas, dist*parametros.legua_cosmica_a_mpc);

    let py_program = format!(r#"
import numpy as np
import matplotlib.pyplot as plt
plt.plot({xx}, {yy})
plt.ticklabel_format(style='sci', axis='y', scilimits=(-3, 3))
plt.ylabel(r"$D_{{\rm HP}}$/Mpc", fontsize=14)
plt.xlabel("$t$/Giga años", fontsize=14)
plt.savefig("plots/horizonte particulas.svg")
    "#, xx=format_as_list(&tt), yy = format_as_list(&dd.iter().map(|x| x*parametros.legua_cosmica_a_mpc).collect::<Vec<_>>()));
    std::fs::write("plots/horizonte_de_partículas.py", py_program)?;
    let mut horiz = Command::new("py").arg("plots/horizonte_de_partículas.py").spawn()?;

    let mut xx = linspace(0., parametros.z_max_edad, 1000);
    let age: Vec<_>;
    if lookback_time(0., parametros).is_finite(){
        age = xx.iter().map(|z| lookback_time(*z, parametros)).collect();
        println!("Edad del universo en z={} es: {} giga años", parametros.z_max_edad, age[age.len()-1])
    }
    else{
        xx = vec![];
        age = vec![];
    }
    let py_program = format!(r#"
import numpy as np
import matplotlib.pyplot as plt
plt.plot({xx}, {yy})
plt.ticklabel_format(style='plain')
plt.ylabel("edad del universo en giga años", fontsize=14)
plt.xlabel("$z$", fontsize=14)
plt.savefig("plots/edad del universo.svg")
    "#, xx=format_as_list(&xx), yy = format_as_list(&age));
    std::fs::write("plots/edad_universo.py", py_program)?;
    let mut edad = Command::new("py").arg("plots/edad_universo.py").spawn()?;



    distancias.wait()?;
    hubble.wait()?;
    escala.wait()?;
    hubble_rad.wait()?;
    horiz.wait()?;
    edad.wait()?;

    println!("Pulse enter para salir...");
    std::io::stdin().read_line(&mut String::new()).unwrap();
    
    Ok(())
}
