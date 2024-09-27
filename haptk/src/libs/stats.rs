use statrs::distribution::{ContinuousCDF, Normal, StudentsT};
use statrs::statistics::Statistics;

// #[cfg(feature = "linear-regression")]
// pub fn linear_regression(x: Vec<f64>, y: Vec<f64>) -> f64 {
//     use ndarray_glm::{Linear, ModelBuilder};
//     let y = ndarray::Array1::from_vec(y);
//     let x = ndarray::Array2::from_shape_vec((y.len(), 1), x).unwrap();

//     let model = ModelBuilder::<Linear>::data(&y, &x).build().unwrap();
//     let fit = model.fit().unwrap();
//     fit.lr_test()
// }

pub fn two_tail_welch_t_test2(x: &[f64], y: &[f64]) -> f64 {
    let n1 = x.len() as f64;
    let n2 = y.len() as f64;
    let mean1 = x.mean();
    let mean2 = y.mean();
    let variance1 = f64::powi(x.std_dev(), 2);
    let variance2 = f64::powi(y.std_dev(), 2);
    let t = (mean1 - mean2) / (variance1 / n1 + variance2 / n2).sqrt();
    let var_by_n1 = variance1 / n1;
    let var_by_n2 = variance2 / n2;
    let upper_part = f64::powi(var_by_n1 + var_by_n2, 2);
    let lower_part = f64::powi(var_by_n1, 2) / (n1 - 1.0) + f64::powi(var_by_n2, 2) / (n2 - 1.0);
    let df = upper_part / lower_part;
    let tdist = StudentsT::new(0.0, 1.0, df).unwrap();
    tdist.cdf(t) * 2.0
}

pub fn two_tail_welch_t_test(x: &[f64], y: &[f64]) -> f64 {
    let x: Summary = x.iter().collect();
    let y: Summary = y.iter().collect();
    let difference = x.compare(&y, 95.0);
    difference.p_value
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Group {
    X,
    Y,
}

//https://docs.rs/rustats/0.1.1/src/rustats/hypothesis_testings.rs.html
pub fn two_tail_mann_whitney_u_test_no_corr(xs: &[i64], ys: &[i64]) -> f64 {
    let mut vs = xs
        .iter()
        .map(|x| (x, Group::X))
        .chain(ys.iter().map(|y| (y, Group::Y)))
        .collect::<Vec<_>>();
    vs.sort_by(|a, b| a.0.partial_cmp(b.0).unwrap());

    let n = vs.len();
    let xn = vs.iter().filter(|t| t.1 == Group::X).count();
    let yn = n - xn;

    let mut counts: Vec<(usize, usize)> = Vec::with_capacity(vs.len());
    let mut prev = None;
    for (v, group) in vs {
        if prev.as_ref() != Some(&v) {
            counts.push((0, 0));
        }
        match group {
            Group::X => counts.last_mut().unwrap_or_else(|| unreachable!()).0 += 1,
            Group::Y => counts.last_mut().unwrap_or_else(|| unreachable!()).1 += 1,
        }
        prev = Some(v);
    }
    let u = calculate_u(n, xn, yn, counts.clone());
    let au = calculate_au(n as f64, xn as f64, yn as f64, counts);
    let mu = (xn * yn / 2) as f64;
    let z = (u - mu) / au;
    let normal_distribution = Normal::new(0.0, 1.0).unwrap();
    (1.0 - normal_distribution.cdf(z.abs())) * 2.0
}

fn calculate_au(n: f64, xn: f64, yn: f64, counts: Vec<(usize, usize)>) -> f64 {
    let t = counts
        .iter()
        .map(|&(x, y)| x + y)
        .map(|t| t * t * t - t)
        .sum::<usize>() as f64;
    ((xn * yn * ((n + 1.0) - t / (n * (n - 1.0)))) / 12.0).sqrt()
}

fn calculate_u(n: usize, xn: usize, yn: usize, counts: Vec<(usize, usize)>) -> f64 {
    let mut xr = 0.0;
    let mut rank = 1;
    for (x, y) in counts.iter().cloned() {
        let temp = (rank..).take(x + y).map(|x| x as f64);
        xr += temp.mean() * x as f64;
        rank += x + y;
    }
    let yr = (n * (n + 1) / 2) as f64 - xr;

    let xu = xr - (xn * (xn + 1) / 2) as f64;
    let yu = yr - (yn * (yn + 1) / 2) as f64;
    xu.min(yu)
}
use std::iter::FromIterator;

/// The statistical difference between two [Summary] instances.
#[derive(Copy, Clone, Debug)]
pub struct Difference {
    /// The absolute difference between the samples' means.
    pub effect: f64,

    /// The difference in means between the two samples, normalized for variance. Technically, this
    /// is Cohen's d.
    pub effect_size: f64,

    /// The minimum allowed effect at the given confidence level.
    pub critical_value: f64,

    /// The p-value for the test: the probability that accepting the results of this test will be a
    /// Type 1 error, in which the null hypothesis (i.e. there is no difference between the means of
    /// the two samples) will be rejected when it is in fact true.
    pub p_value: f64,

    /// The significance level of the test. It is the maximum allowed value of the p-value.
    pub alpha: f64,

    /// The probability of a Type 2 error: the probability that the null hypothesis will be retained
    /// despite it not being true.
    pub beta: f64,
}

impl Difference {
    /// Whether or not the difference is statistically significant.
    pub fn is_significant(&self) -> bool {
        self.effect > self.critical_value
    }
}

/// A statistical summary of a normally distributed data set.
///
/// Created from an iterable of `f64`s:
///
/// ```
/// let summary: nanostat::Summary = vec![0.1, 0.45, 0.42].iter().collect();
/// ```
#[derive(Copy, Clone, Debug)]
pub struct Summary {
    /// The number of measurements in the set.
    pub n: f64,
    /// The arithmetic mean of the measurements.
    pub mean: f64,
    /// The sample variance of the data set.
    pub variance: f64,
}

impl<'a> FromIterator<&'a f64> for Summary {
    fn from_iter<T: IntoIterator<Item = &'a f64>>(iter: T) -> Self {
        // Welford's one-pass algorithm for corrected variance
        let (mut mean, mut s, mut n) = (0.0, 0.0, 0.0);
        for x in iter {
            n += 1.0;
            let delta = x - mean;
            mean += delta / n;
            s += delta * (x - mean);
        }
        let variance = s / (n - 1.0); // Bessel's correction
        Summary { n, mean, variance }
    }
}

impl Summary {
    /// The standard deviation of the sample.
    pub fn std_dev(&self) -> f64 {
        self.variance.sqrt()
    }

    /// The standard error of the sample.
    pub fn std_err(&self) -> f64 {
        self.std_dev() / self.n.sqrt()
    }

    /// Calculate the statistical difference between the two summaries using a two-tailed Welch's
    /// t-test. The confidence level must be in the range `(0, 100)`.
    pub fn compare(&self, other: &Summary, confidence: f64) -> Difference {
        assert!(
            0.0 < confidence && confidence < 100.0,
            "confidence must be (0,100)"
        );

        let (a, b) = (self, other);

        // Calculate the significance level.
        let alpha = 1.0 - (confidence / 100.0);

        // Calculate the degrees of freedom.
        let nu = (a.variance / a.n + b.variance / b.n).powf(2.0)
            / ((a.variance).powf(2.0) / ((a.n).powf(2.0) * (a.n - 1.0))
                + (b.variance).powf(2.0) / ((b.n).powf(2.0) * (b.n - 1.0)));

        // Create a Student's T distribution with location of 0, a scale of 1, and the same number
        // of degrees of freedom as in the test.
        let dist_st = StudentsT::new(0.0, 1.0, nu).unwrap();

        // Calculate the hypothetical two-tailed t-value for the given significance level.
        let t_hyp = dist_st.inverse_cdf(1.0 - (alpha / TAILS));

        // Calculate the absolute difference between the means of the two samples.
        let effect = (a.mean - b.mean).abs();

        // Calculate the standard error.
        let std_err = (a.variance / a.n + b.variance / b.n).sqrt();

        // Calculate the experimental t-value.
        let t_exp = effect / std_err;

        // Calculate the p-value given the experimental t-value.
        let p_value = dist_st.cdf(-t_exp) * TAILS;

        // Calculate the critical value.
        let critical_value = t_hyp * std_err;

        // Calculate the standard deviation using mean variance.
        let std_dev = ((a.variance + b.variance) / 2.0).sqrt();

        // Calculate Cohen's d for the effect size.
        let effect_size = effect / std_dev;

        // Calculate the statistical power.
        let z = effect / (std_dev * (1.0 / a.n + 1.0 / b.n).sqrt());
        let dist_norm = Normal::new(0.0, 1.0).unwrap();
        let za = dist_norm.inverse_cdf(1.0 - alpha / TAILS);
        let beta = dist_norm.cdf(z - za) - dist_norm.cdf(-z - za);

        Difference {
            effect,
            effect_size,
            critical_value,
            p_value,
            alpha,
            beta,
        }
    }
}

/// The number of distribution tails used to determine significance. In this case, we always use a
/// two-tailed test because our null hypothesis is that the samples are not different.
const TAILS: f64 = 2.0;

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;

    #[test]
    fn test_two_tail_welch_t_test() {
        let vec1 = &[1.0, 2.0, 1.0, 3.0];
        let vec2 = &[1.0, 2.0, 1.0, 5.0];
        let result = two_tail_welch_t_test(vec1, vec2);
        assert_eq!("0.6596", &format!("{result:.4}"));

        let vec1 = &[1.0, -2.0, 1.0, 3.0];
        let vec2 = &[-1.0, 2.0, 1.0, 5.0];
        let result = two_tail_welch_t_test(vec1, vec2);
        assert_eq!("0.5606", &format!("{result:.4}"));

        let vec1 = &[1.0, -2.0, 1.0, 3.0];
        let vec2 = &[-1.0, 2.0, 1.0, 5.0, 5.0, 5.0];
        let result = two_tail_welch_t_test(vec1, vec2);
        assert_eq!("0.1959", &format!("{result:.4}"));
    }

    #[test]
    fn two_tail_mwu_test() {
        let vec1 = &[1, 2, 1, 3, 2, 1];
        let vec2 = &[1, 2, 1, 100, 50, 30];
        let result = two_tail_mann_whitney_u_test_no_corr(vec1, vec2);
        assert_eq!("0.2416", &format!("{result:.4}"));

        let vec1 = &[-1, 2, 1, 3, 2, 1];
        let vec2 = &[1, -2, 1, 100, 50, 30];
        let result = two_tail_mann_whitney_u_test_no_corr(vec1, vec2);
        assert_eq!("0.5136", &format!("{result:.4}"));

        let vec1 = &[1, -2, 1, 3, 4, 4];
        let vec2 = &[-1, 2, 1, 5, 5, 5, 5, 5, 5, 5];
        let result = two_tail_mann_whitney_u_test_no_corr(vec1, vec2);
        assert_eq!("0.0407", &format!("{result:.4}"));
    }
}
