pub fn compute_l2_error(numerical: &[f64], analytical: &[f64]) -> f64 {
    let sum_of_squares = numerical
        .iter()
        .zip(analytical.iter())
        .map(|(num, anal)| (num - anal).powi(2))
        .sum::<f64>();
    (sum_of_squares / numerical.len() as f64).sqrt()
}

pub fn compute_order_of_accuracy(dts: &[f64], l2_errors: &[f64]) -> Option<f64> {
    if dts.len() != l2_errors.len() || dts.len() < 2 {
        return None;
    }

    let mut orders_of_accuracy = vec![];
    for i in 1..dts.len() {
        let h1 = dts[i - 1];
        let e1 = l2_errors[i - 1];
        let h2 = dts[i];
        let e2 = l2_errors[i];

        // Compute the order of accuracy
        let order = (e1.log10() - e2.log10()) / (h1 / h2).log10();
        orders_of_accuracy.push(order);
    }

    // Calculate the average order of accuracy
    let average_order: f64 = orders_of_accuracy.iter().sum::<f64>() / orders_of_accuracy.len() as f64;

    Some(average_order)
}
