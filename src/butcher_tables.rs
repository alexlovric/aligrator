// Coefficients for RK4
#[rustfmt::skip]
pub fn rk4_classical_table() -> ([f64; 4], [f64; 6], [f64; 4]) {
    let a = [0., 1./2., 1./2., 1.];
    let b = [
        1./2.,
        0.0, 1./2.,
        0., 0., 1.,
    ];
    let c = [1./6., 2./6., 2./6., 1./6.];

    (a, b, c)
}

// Coefficients for RK45 (Fehlberg)
#[rustfmt::skip]
pub fn rk45_fehlberg_alfa_one_third_table() -> ([f64; 6], [f64; 15], [f64; 6], [f64; 6]) {
    let a = [0., 2./9., 1./3., 3./4., 1., 5./6.];
    let b = [
        2./9.,
        1./12., 1./4., 
        69./128., -243./128., 135./64., 
        -17./12., 27./4., -27./5., 16./15.,
        65./432., -5./16., 13./16., 4./27., 5./144.,
    ];
    let ch = [47./450., 0., 12./25., 32./225., 1./30., 6./25.];
    let ct = [1./150., 0., -3./100., 16./75., 1./20., -6./25.];

    (a, b, ch, ct)
}

#[rustfmt::skip]
pub fn _rk45_fehlberg_alfa_three_eigths_table() -> ([f64; 6], [f64; 15], [f64; 6], [f64; 6]) {
    let a = [0., 1./4., 3./8., 12./13., 1., 1./2.];
    let b = [
        1./4.,
        3./32., 9./32.,
        1932./2197., -7200./2197., 7296./2197., 
        439./216., -8., 3680./513., -845./4104., 
        -8./27., 2., -3544./2565., 1859./4104., -11./40.
    ];
    let ch = [16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55.];
    let ct = [-1./360., 0., 128./4275., 2197./75240., -1./50., -2./55.];

    (a, b, ch, ct)
}

// Explicit Runge-Kutta Methods with Estimates of the Local Truncation Error
// Author(s): J. H. Verner
#[rustfmt::skip]
pub fn rk89_verner_table() -> ([f64; 16], [f64; 120], [f64; 16], [f64; 16]) {
    let rt6 = 6.0_f64.sqrt();

    let a = [
        0., 1./12., 1./9., 1./6., (2. + 2. * rt6)/15., (6. + rt6)/15., (6. - rt6)/15., 
        2./3., 1./2., 1./3., 1./4., 4./3., 5./6., 1., 1./6., 1.
    ];

    let b = [
        1./12.,
        1./27., 2./27.,
        1./24., 0., 1./8.,
        (4.+94.*rt6)/375., 0., (-94. - 84.*rt6)/125., (328.+208.*rt6)/375.,
        (9.-rt6)/150., 0., 0., (312.+32.*rt6)/1425., (69.+29.*rt6)/570.,
        (927.-347.*rt6)/1250., 0., 0., (-16248.+7328.*rt6)/9375., (-489.+179.*rt6)/3750., (14268.-5798.*rt6)/9375.,
        2./27., 0., 0., 0., 0., (16.-rt6)/54., (16.+rt6)/54.,
        19./256., 0., 0., 0., 0., (118.-23.*rt6)/512., (118.+23.*rt6)/512., -9./256.,
        11./144., 0., 0., 0., 0., (266.-rt6)/864., (266.+rt6)/864., -1./16., -8./27.,
        (5034.-271.*rt6)/61440., 0., 0., 0., 0., 0., (7859.-1626.*rt6)/10240., (-2232.+813.*rt6)/20480., (-594.+271.*rt6)/960., (657.-813.*rt6)/5120.,
        (5996.-3794.*rt6)/405., 0., 0., 0., 0., (-4342.-338.*rt6)/9., (154922.-40458.*rt6)/135., (-4176.+3794.*rt6)/45., (-340864.+242816.*rt6)/405., (26304.-15176.*rt6)/45., -26624./81.,
        (3793.+2168.*rt6)/103680., 0., 0., 0., 0., (4042.+2263.*rt6)/13824., (-231278.+40717.*rt6)/69120., (7947.-2168.*rt6)/11520., (1048.-542.*rt6)/405., (-1383.+542.*rt6)/720., 2624./1053., 3./1664.,
        -137./1296., 0., 0., 0., 0., (5642.-337.*rt6)/864., (5642.+337.*rt6)/864., -299./48., 184./81., -44./9., -5120./1053., -11./468., 16./9.,
        (33617.-2168.*rt6)/518400., 0., 0., 0., 0., (-3846.+31.*rt6)/13824., (155338.-52807.*rt6)/345600., (-12537.+2168.*rt6)/57600., (92.+542.*rt6)/2025., (-1797.-542.*rt6)/3600., 320./567., -1./1920., 4./105., 0.,
        (-36487.-30352.*rt6)/279600., 0., 0., 0., 0., (-29666.-4499.*rt6)/7456., (2779182.-615973.*rt6)/186400., (-94329.+91056.*rt6)/93200., (-232192.+121408.*rt6)/17475., (101226.-22764.*rt6)/5825., -169984./9087., -87./30290., 492./1165., 0., 1260./233.,
    ];

    let c = [
        23./525., 0., 0., 0., 0., 0., 0., 171./1400., 86./525., 93./280., -2048./6825., 
        -3./18200., 39./175., 0., 9./25., 233./4200.
    ];

    let e = [
        -7./400., 0., 0., 0., 0., 0., 0., 63./200., -14./25., 21./20., -1024./975., 
        -21./36400., -3./25., -9./280., 9./25., 233./4200.
    ];

    (a, b, c, e)
}

// #[rustfmt::skip]
// fn get_classical_rk4_butcher_table<const O: usize, const S: usize>(order: usize) -> ([f64; O], [f64; S], [f64; O]) {
//     let result = match order {
//         4 => {
//             let a: [f64; O] = [0., 1./2., 1./2., 1.];
//             let b: [f64; 6] = [
//                 1./2.,
//                 0.0, 1./2.,
//                 0., 0., 1.,
//             ];
//             let c: [f64; 4] = [1./6., 2./6., 2./6., 1./6.];

//             return (a,b,c)
//         },
//         // Add more cases here for different orders and stages
//         _ => todo!("🐒 Not implemented unsupported order and stages combination"),
//     };
// }
