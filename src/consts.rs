pub const C: f64 = 299792458.;// m/s
pub const M_TO_MPC: f64 = 1./3.0857e22;
pub const C_MPC: f64 = C*M_TO_MPC;
pub const C_GY: f64 = C/S_TO_GY;
pub const C_MPC_GY: f64 = C_MPC/S_TO_GY;
pub const S_TO_GY: f64 = 1./(3600.*24.*365.25*1e9);