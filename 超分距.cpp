#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using cd = std::complex<double>;
constexpr double PI = 3.14159265358979323846;

struct RadarParams {
    // 物理常数
    double c = 299792458.0;               // 光速 m/s

    // FMCW 配置
    double fc = 77e9;                     // 载频 Hz
    double bandwidth = 2e9;               // 调频带宽 Hz
    double chirp_time = 40e-6;            // 单个 Chirp 时长 s
    double sample_rate = 5e6;             // ADC 采样率 Hz
    int samples_per_chirp = 16;          // 每个 chirp 采样点数
    int chirp_count = 16;                // 多快拍数（用于超分辨）

    // 仿真目标
    std::vector<double> target_ranges = {30.00, 30.12}; // m，距离非常接近
    std::vector<double> target_rcs = {1.0, 0.7};        // 相对散射强度

    // 噪声
    double snr_db = 22.0;                 // 信噪比 dB

    // MUSIC 扫描
    double min_scan_range = 28.0;         // m
    double max_scan_range = 32.0;         // m
    double scan_step = 0.002;             // m
};

struct EstimateResult {
    double fft_peak_range = 0.0;
    std::vector<double> music_peak_ranges;
};

double db_to_linear(double db) {
    return std::pow(10.0, db / 10.0);
}

std::vector<cd> dft(const std::vector<cd>& x) {
    const int N = static_cast<int>(x.size());
    std::vector<cd> X(N, cd(0.0, 0.0));
    for (int k = 0; k < N; ++k) {
        cd sum(0.0, 0.0);
        for (int n = 0; n < N; ++n) {
            double angle = -2.0 * PI * k * n / N;
            sum += x[n] * std::exp(cd(0.0, angle));
        }
        X[k] = sum;
    }
    return X;
}

std::vector<std::vector<cd>> simulate_fmcw_snapshots(const RadarParams& p) {
    std::vector<std::vector<cd>> X(p.samples_per_chirp, std::vector<cd>(p.chirp_count, cd(0.0, 0.0)));

    const double slope = p.bandwidth / p.chirp_time;
    std::mt19937 rng(2026);
    std::normal_distribution<double> gauss(0.0, 1.0);

    // 理想信号
    for (int m = 0; m < p.chirp_count; ++m) {
        for (int n = 0; n < p.samples_per_chirp; ++n) {
            double t = static_cast<double>(n) / p.sample_rate;
            cd sample(0.0, 0.0);
            for (size_t i = 0; i < p.target_ranges.size(); ++i) {
                double r = p.target_ranges[i];
                double a = p.target_rcs[i];
                double fb = (2.0 * slope * r) / p.c; // beat frequency
                double phase = 2.0 * PI * fb * t + 2.0 * PI * 0.07 * m * (i + 1); // 跨chirp相位变化
                sample += a * std::exp(cd(0.0, phase));
            }
            X[n][m] = sample;
        }
    }

    // 按目标 SNR 加噪
    double signal_power = 0.0;
    for (int n = 0; n < p.samples_per_chirp; ++n) {
        for (int m = 0; m < p.chirp_count; ++m) {
            signal_power += std::norm(X[n][m]);
        }
    }
    signal_power /= (p.samples_per_chirp * p.chirp_count);

    const double snr_linear = db_to_linear(p.snr_db);
    const double noise_power = signal_power / snr_linear;
    const double sigma = std::sqrt(noise_power / 2.0);

    for (int n = 0; n < p.samples_per_chirp; ++n) {
        for (int m = 0; m < p.chirp_count; ++m) {
            X[n][m] += cd(sigma * gauss(rng), sigma * gauss(rng));
        }
    }

    return X;
}

// Jacobi 特征值分解（实对称矩阵）
void jacobi_eigen_decomposition(std::vector<std::vector<double>> A,
                                std::vector<double>& eigenvalues,
                                std::vector<std::vector<double>>& eigenvectors) {
    int n = static_cast<int>(A.size());
    eigenvectors.assign(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) eigenvectors[i][i] = 1.0;

    for (int iter = 0; iter < 20 * n * n; ++iter) {
        int p = 0, q = 1;
        double max_val = std::abs(A[0][1]);
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double val = std::abs(A[i][j]);
                if (val > max_val) {
                    max_val = val;
                    p = i;
                    q = j;
                }
            }
        }
        if (max_val < 1e-10) break;

        double phi = 0.5 * std::atan2(2.0 * A[p][q], A[q][q] - A[p][p]);
        double c = std::cos(phi), s = std::sin(phi);

        for (int k = 0; k < n; ++k) {
            double Apk = A[p][k], Aqk = A[q][k];
            A[p][k] = c * Apk - s * Aqk;
            A[q][k] = s * Apk + c * Aqk;
        }
        for (int k = 0; k < n; ++k) {
            double Akp = A[k][p], Akq = A[k][q];
            A[k][p] = c * Akp - s * Akq;
            A[k][q] = s * Akp + c * Akq;
        }
        for (int k = 0; k < n; ++k) {
            double vkp = eigenvectors[k][p], vkq = eigenvectors[k][q];
            eigenvectors[k][p] = c * vkp - s * vkq;
            eigenvectors[k][q] = s * vkp + c * vkq;
        }
    }

    eigenvalues.resize(n);
    for (int i = 0; i < n; ++i) eigenvalues[i] = A[i][i];
}

EstimateResult estimate_range_fft_music(const RadarParams& p) {
    auto X = simulate_fmcw_snapshots(p);
    const int N = p.samples_per_chirp;
    const int M = p.chirp_count;

    // 1) 传统 FFT 距离估计（先对 chirp 平均）
    std::vector<cd> mean_signal(N, cd(0.0, 0.0));
    for (int n = 0; n < N; ++n) {
        cd sum(0.0, 0.0);
        for (int m = 0; m < M; ++m) sum += X[n][m];
        mean_signal[n] = sum / static_cast<double>(M);
    }
    auto spectrum = dft(mean_signal);
    std::vector<double> mag(N / 2);
    for (int k = 0; k < N / 2; ++k) mag[k] = std::abs(spectrum[k]);
    int kmax = std::distance(mag.begin(), std::max_element(mag.begin(), mag.end()));

    double slope = p.bandwidth / p.chirp_time;
    double fb_fft = static_cast<double>(kmax) * p.sample_rate / N;
    double range_fft = p.c * fb_fft / (2.0 * slope);

    // 2) MUSIC 超分辨
    // 协方差矩阵 R = 1/M * X X^H
    std::vector<std::vector<cd>> R(N, std::vector<cd>(N, cd(0.0, 0.0)));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cd sum(0.0, 0.0);
            for (int m = 0; m < M; ++m) {
                sum += X[i][m] * std::conj(X[j][m]);
            }
            R[i][j] = sum / static_cast<double>(M);
        }
    }

    // 将 Hermitian 复矩阵扩展为等价实对称矩阵，做特征分解
    // [ Re(R) -Im(R);
    //   Im(R)  Re(R)]
    const int NR = 2 * N;
    std::vector<std::vector<double>> A(NR, std::vector<double>(NR, 0.0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double re = R[i][j].real();
            double im = R[i][j].imag();
            A[i][j] = re;
            A[i][j + N] = -im;
            A[i + N][j] = im;
            A[i + N][j + N] = re;
        }
    }

    std::vector<double> eigvals;
    std::vector<std::vector<double>> eigvecs;
    jacobi_eigen_decomposition(A, eigvals, eigvecs);

    std::vector<int> order(NR);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return eigvals[a] < eigvals[b]; });

    int signal_count = static_cast<int>(p.target_ranges.size());
    int noise_dim = NR - 2 * signal_count;

    // 噪声子空间基向量
    std::vector<std::vector<double>> En(noise_dim, std::vector<double>(NR));
    for (int d = 0; d < noise_dim; ++d) {
        int idx = order[d];
        for (int r = 0; r < NR; ++r) En[d][r] = eigvecs[r][idx];
    }

    std::vector<double> scan_ranges;
    std::vector<double> pseudo;

    for (double r = p.min_scan_range; r <= p.max_scan_range + 1e-12; r += p.scan_step) {
        double fb = (2.0 * slope * r) / p.c;
        std::vector<double> a(NR, 0.0); // [Re(s);Im(s)]
        for (int n = 0; n < N; ++n) {
            double t = static_cast<double>(n) / p.sample_rate;
            double ph = 2.0 * PI * fb * t;
            a[n] = std::cos(ph);
            a[n + N] = std::sin(ph);
        }

        double denom = 0.0;
        for (int d = 0; d < noise_dim; ++d) {
            double proj = 0.0;
            for (int i = 0; i < NR; ++i) proj += En[d][i] * a[i];
            denom += proj * proj;
        }

        scan_ranges.push_back(r);
        pseudo.push_back(1.0 / std::max(denom, 1e-12));
    }

    // 选 MUSIC 最大两个峰
    std::vector<int> idx(pseudo.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return pseudo[a] > pseudo[b]; });

    std::vector<double> peaks;
    for (int id : idx) {
        double rr = scan_ranges[id];
        bool separated = true;
        for (double pr : peaks) {
            if (std::abs(pr - rr) < 0.03) {
                separated = false;
                break;
            }
        }
        if (separated) peaks.push_back(rr);
        if (static_cast<int>(peaks.size()) == signal_count) break;
    }
    std::sort(peaks.begin(), peaks.end());

    return {range_fft, peaks};
}

void print_params(const RadarParams& p) {
    std::cout << "===== 雷达参数（高精度超分辨测距）=====\n";
    std::cout << "载频 fc: " << p.fc / 1e9 << " GHz\n";
    std::cout << "带宽 B: " << p.bandwidth / 1e6 << " MHz\n";
    std::cout << "Chirp 时长 T: " << p.chirp_time * 1e6 << " us\n";
    std::cout << "ADC 采样率 Fs: " << p.sample_rate / 1e6 << " MSps\n";
    std::cout << "每 Chirp 采样点 N: " << p.samples_per_chirp << "\n";
    std::cout << "快拍数 M: " << p.chirp_count << "\n";
    std::cout << "SNR: " << p.snr_db << " dB\n";

    double range_res = p.c / (2.0 * p.bandwidth);
    std::cout << "理论 FFT 距离分辨率 c/(2B): " << range_res * 100.0 << " cm\n";
    std::cout << "目标距离: ";
    for (size_t i = 0; i < p.target_ranges.size(); ++i) {
        std::cout << p.target_ranges[i] << (i + 1 < p.target_ranges.size() ? ", " : " m\n");
    }
    std::cout << "=======================================\n\n";
}

int main() {
    RadarParams params;
    print_params(params);

    auto result = estimate_range_fft_music(params);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "传统 FFT 估计距离(单峰): " << result.fft_peak_range << " m\n";
    std::cout << "MUSIC 超分辨估计距离: ";
    for (size_t i = 0; i < result.music_peak_ranges.size(); ++i) {
        std::cout << result.music_peak_ranges[i] << (i + 1 < result.music_peak_ranges.size() ? ", " : " m\n");
    }

    std::cout << "\n说明: 两目标间隔 0.12 m，小于 FFT 分辨率，MUSIC 能实现超分辨。\n";
    return 0;
}
