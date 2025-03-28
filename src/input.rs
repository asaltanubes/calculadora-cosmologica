
use serde::{Serialize, Deserialize};
use crate::{C_MPC_GY, S_TO_GY};

#[derive(Clone, Copy, Serialize, Deserialize, Debug, Default)]
pub struct InputFile{
    densidad_materia: f64,
    densidad_radiacion: f64,
    densidad_energia_oscura: f64,
    constante_hubble: f64,

    rango_factor_escala: (f64, f64),
    t_horizonte_particulas: f64,
    z_max_dist: f64,
    z_max_edad: f64
}

#[derive(Clone, Copy, Debug, Default)]
pub struct ParametrosCosmologicos{
    pub densidad_materia: f64,
    pub densidad_radiacion: f64,
    pub densidad_energia_oscura: f64,
    pub densidad_curvatura: f64,
    pub constante_hubble: f64,
    pub factor_escala: f64,
    pub signo_curvatura: f64,
    pub c: f64,
    pub legua_cosmica_a_mpc: f64,

    pub rango_factor_escala: (f64, f64),
    pub t_horizonte_particulas: f64,
    pub z_max_dist: f64,
    pub z_max_edad: f64
}

impl std::convert::From<InputFile> for ParametrosCosmologicos{
    fn from(value: InputFile) -> Self{
        let densidad_materia = value.densidad_materia;
        let densidad_radiacion = value.densidad_radiacion;
        let densidad_energia_oscura = value.densidad_energia_oscura;
        let constante_hubble = value.constante_hubble/3.08567758e19/S_TO_GY; // km/s/Mpc a 1/GY
        let densidad_curvatura = 1. - densidad_energia_oscura - densidad_materia - densidad_radiacion;
        let signo_curvatura = -if densidad_curvatura == 0. {0.} else {densidad_curvatura.signum()};

        let c =  if densidad_curvatura != 0. {densidad_curvatura.abs().sqrt() * 1. * constante_hubble} //  leguas cósmica / GY
                     else{C_MPC_GY}; // legua cósmica = Mpc
        
        let legua_cosmica_a_mpc = C_MPC_GY/c;
        ParametrosCosmologicos{densidad_materia, densidad_radiacion, densidad_energia_oscura, 
                                densidad_curvatura, constante_hubble, signo_curvatura, c, factor_escala: 1., legua_cosmica_a_mpc, 
                                rango_factor_escala: value.rango_factor_escala, t_horizonte_particulas: value.t_horizonte_particulas,
                            z_max_dist: value.z_max_dist, z_max_edad: value.z_max_edad}
    } 
}