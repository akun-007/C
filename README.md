# C

FMCW 雷达超分辨测距（MUSIC）示例代码（C++17，单文件可运行）。

## 功能
- 生成带噪声的 FMCW 拍频信号（多 chirp 快拍）。
- 给出**传统 FFT 距离估计**结果（基线）。
- 给出 **MUSIC 超分辨距离估计**结果（可分离近距离目标）。
- 在程序中集中定义雷达关键参数，便于按项目需求调参。

## 雷达参数（代码内 `RadarParams`）
- 载频 `fc = 77e9` Hz
- 带宽 `bandwidth = 2e9` Hz
- Chirp 时长 `chirp_time = 40e-6` s
- ADC 采样率 `sample_rate = 5e6` Hz
- 每 Chirp 采样点 `samples_per_chirp = 16`
- 快拍数 `chirp_count = 16`
- 仿真目标距离 `target_ranges = {30.00, 30.12}` m
- 目标强度 `target_rcs = {1.0, 0.7}`
- 信噪比 `snr_db = 22`
- MUSIC 扫描区间 `[28, 32]` m，步长 `0.002` m

> 说明：为了让示例在受限环境下快速运行，默认采样点和快拍数设置较小；工程中请增大 `samples_per_chirp/chirp_count` 并采用 FFTW/Eigen/Armadillo 等高性能库。

## 编译运行
```bash
g++ -std=c++17 -O2 test.cpp -o radar_demo
./radar_demo
```

## 输出解读
- `传统 FFT 估计距离(单峰)`：低复杂度基线；近距双目标下可能只能给出单峰或偏差较大。
- `MUSIC 超分辨估计距离`：可在同一距离单元中分离更近的目标。

## 工程落地建议（高精度）
1. 增大带宽（提升理论距离分辨率 `c/(2B)`）。
2. 增加快拍数与相干积累时间（提升子空间稳定性）。
3. 使用更高性能线代与 FFT 库，替换示例中的教学版 DFT/Jacobi。
4. 加入硬件标定（时钟、IQ 不平衡、相位噪声、斜率非线性）。
5. 做 CFAR + 跟踪滤波（卡尔曼/IMM）提升稳态测距精度。
