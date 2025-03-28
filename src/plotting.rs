
pub fn linspace(min: f64, max: f64, steps: i64) -> Vec<f64>{
    (0..steps).map(|x| x as f64/steps as f64 * (max-min) + min).collect()
}


pub fn format_as_list<'a, T: IntoIterator<Item = I>, I: Into<&'a f64>>(list: T) -> String{
    format!("[{}]", list.into_iter().map(|x| x.into().to_string()).collect::<Vec<String>>().join(", ")).replace("NaN", "np.nan")
}