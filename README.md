# Colocation Mining with Fast Maximal Clique & Rare Feature

Cài đặt C++17 của thuật toán khai phá **mẫu đồng vị trí phổ biến (Prevalent Co-location Patterns – PCP)** trên dữ liệu không gian, kết hợp hai ý tưởng chính:

1. **Fast Maximal Clique Enumeration (MCE)** — liệt kê clique cực đại bằng thuật toán lai *Degeneracy ordering + Bron–Kerbosch Pivot + RCD (Recursive Core Decomposition)*, tự động chuyển đổi giữa hai nhánh dựa trên cấu trúc đồ thị con.
2. **Cauchy kernel weighting cho đặc trưng hiếm (rare feature)** — dùng chỉ số *Rare Intensity* để tính **Weighted Participation Index (WPI)**, tránh việc các đặc trưng hiếm bị các đặc trưng phổ biến lấn át.

---

## Luồng xử lý

```
CSV dataset
   │
   ├─ (1) DataLoader::load_csv          → danh sách SpatialInstance (+ lấy mẫu theo %)
   ├─ (2) countFeatures                 → số lượng instance của từng feature
   ├─ (3) calculateDispersion           → delta (độ phân tán tần suất giữa các feature)
   ├─ (4) NeighborGraph::buildNeighborGraph
   │        plane sweep theo trục X + khoảng cách Euclid ≤ neighbor_distance
   │        (chỉ nối các instance KHÁC feature type)
   ├─ (5) MaximalCliqueHashmap::executeBK
   │        degeneracy ordering → với mỗi đỉnh chọn BK-RCD hoặc BK-Pivot
   │        kết quả: hashmap  Colocation → {FeatureType → tập instance}
   ├─ (6) extractInitialCandidates      → hàng đợi ưu tiên các ứng viên (size lớn trước)
   └─ (7) Miner::minePCPs               → lọc theo WPI ≥ min_prevalence
                                           + cắt tỉa bằng downward closure (Lemma 2)
   →  results.txt
```

### Cơ chế chuyển đổi thuật toán (hybrid switch)

Với mỗi đồ thị con `P`, chương trình tính bộ ba `(s, k, div)` trong [maximal_clique_hashmap.cpp](src/maximal_clique_hashmap.cpp#L263-L300):

| Ký hiệu | Ý nghĩa |
|---|---|
| `s` | số đỉnh nối với **tất cả** các đỉnh còn lại (kernel) |
| `k` | số đỉnh còn lại (shell) |
| `div` | độ phân kỳ `max(deg) − s` trong tập shell |

Ngưỡng chọn nhánh (theo paper, có điều chỉnh theo `div`):

```
div == 0 → threshold = 2.8k − 4.5
div == 1 → threshold = 2.8k − 8.0
div >= 2 → threshold = 2.8k − 11.0

s >= threshold  →  BK-RCD   (vùng đặc)
s <  threshold  →  BK-Pivot (vùng thưa)
```

### Trọng số Cauchy cho đặc trưng hiếm

Cho colocation `C`, gọi `N(f)` là số instance của feature `f` và `N(f_min)` là giá trị nhỏ nhất trong `C`:

```
delta = 2/(m(m−1)) · Σ_{i<j} num(f_j)/num(f_i)      (m = số feature toàn cục, tần suất sắp tăng dần)

v(f)  = N(f) / N(f_min)
w(f)  = 1 + ((v(f) − 1) / delta)²                    ← Cauchy kernel
RI(f) = 1 / w(f)

WPR(f) = PR(f) · (1 / RI(f))       với PR(f) = |instance của f tham gia| / N(f)
WPI(C) = min_{f ∈ C} WPR(f)
```

`C` là mẫu phổ biến khi `WPI(C) ≥ min_prevalence`.

---

## Cấu trúc thư mục

```
├── CMakeLists.txt
├── config/
│   └── config.txt                   # tham số chạy
├── data/                            # các bộ dữ liệu CSV
├── include/
│   ├── config.h                     # AppConfig + ConfigLoader
│   ├── csv.hpp                      # thư viện đọc CSV (header-only, vincentlaucsb/csv-parser)
│   ├── data_loader.h
│   ├── maximal_clique_hashmap.h     # liệt kê clique cực đại
│   ├── miner.h                      # khai phá PCP
│   ├── neighbor_graph.h
│   ├── types.h                      # SpatialInstance, NeighborSet, Colocation…
│   └── utils.h                      # countFeatures, calculateDispersion, calcRareIntensity
└── src/
    ├── config.cpp
    ├── data_loader.cpp
    ├── main.cpp                     # entry point + báo cáo kết quả
    ├── maximal_clique_hashmap.cpp
    ├── miner.cpp
    ├── neighbor_graph.cpp
    └── utils.cpp
```

---

## Yêu cầu

- CMake ≥ 3.10
- Trình biên dịch hỗ trợ C++17 (MSVC / Visual Studio 2019+)
- **Windows**: [main.cpp](src/main.cpp#L20-L23) dùng `windows.h`, `psapi.h` để đo peak memory nên hiện tại chỉ build được trên Windows. Muốn chạy Linux/macOS cần bỏ hoặc thay thế phần đo bộ nhớ này.

## Build & chạy

Mở thư mục bằng Visual Studio (hỗ trợ CMake sẵn), hoặc dùng dòng lệnh:

```powershell
cmake -S . -B out/build/x64-Release -DCMAKE_BUILD_TYPE=Release
cmake --build out/build/x64-Release --config Release
```

Target `copy_resources` tự động copy `config/` và `data/` sang thư mục build, nên chỉ cần chạy:

```powershell
cd out/build/x64-Release
./main                     # dùng ./config/config.txt
./main path/to/config.txt  # hoặc chỉ định file config khác
```

> Nên chạy ở cấu hình **Release**: các bộ dữ liệu lớn (ví dụ `plant_big_100_ex_checkin.csv` ~409k dòng) chậm hơn nhiều ở Debug.

---

## Cấu hình ([config/config.txt](config/config.txt))

Định dạng `key=value`, dòng bắt đầu bằng `#` là chú thích. Khóa không nhận diện được sẽ bị bỏ qua, thiếu file thì dùng giá trị mặc định.

| Khóa | Kiểu | Mặc định | Ý nghĩa |
|---|---|---|---|
| `dataset_path` | string | `data/sample_data.csv` | đường dẫn CSV đầu vào |
| `neighbor_distance` | double | `5.0` | ngưỡng khoảng cách Euclid để coi hai instance là láng giềng |
| `min_prevalence` | double | `0.6` | ngưỡng WPI tối thiểu để một colocation là phổ biến |
| `min_cond_prob` | double | `0.5` | ngưỡng xác suất có điều kiện (đã đọc, chưa dùng trong luồng hiện tại) |
| `percentage_instances` | double | `1.0` | tỉ lệ lấy mẫu ngẫu nhiên theo từng feature (`1.0` = dùng toàn bộ) |
| `debug_mode` | bool | `false` | bật log gỡ lỗi |

`output_path` có trong file config nhưng **chưa** được `ConfigLoader` đọc — đường dẫn kết quả hiện đang cố định trong `main.cpp`.

---

## Định dạng dữ liệu

CSV có header, bắt buộc các cột `Feature`, `Instance`, và cặp toạ độ `LocX`/`LocY` **hoặc** `X`/`Y` (loader tự nhận diện). Các cột thừa như `Checkin` được bỏ qua.

```csv
Feature,Instance,LocX,LocY
A,1,9,8
A,2,3,4
B,1,7,4
```

ID instance được sinh bằng `Feature + Instance` (ví dụ `A1`, `B2`).

### Các bộ dữ liệu kèm theo

| File | Số dòng | Ghi chú |
|---|---|---|
| `sample_data.csv` | 15 | dữ liệu đồ chơi để kiểm thử nhanh |
| `gau_mountain.csv` | ~11k | toạ độ lớn (~1e7), cần `neighbor_distance` tương ứng |
| `LasVegas_x_y_alphabet_version_03_2.csv` | ~22.7k | dùng mặc định trong config (`neighbor_distance=160`) |
| `5k_15f_40k.csv` | 40k | 15 feature, dữ liệu tổng hợp |
| `Shanghai_POI_sample_035_checkin_exp_dis.csv` | ~41k | POI Thượng Hải kèm số check-in |
| `plant_big_100_ex_checkin.csv` | ~409k | bộ lớn nhất, dùng `percentage_instances` để giảm tải |

> `neighbor_distance` phụ thuộc hoàn toàn vào thang toạ độ của từng bộ dữ liệu — đổi dataset thì phải chỉnh lại tham số này.

---

## Kết quả

Chương trình ghi báo cáo ra `results.txt` ở **thư mục cha của thư mục chạy** (`../results.txt` so với working directory), gồm:

- thông tin dataset và tham số đã dùng
- tổng thời gian thực thi (giây)
- peak memory (MB, qua `GetProcessMemoryInfo`)
- số mẫu tìm được và danh sách đầy đủ các mẫu

```
=== FINAL REPORT ===
Dataset Path:      data/LasVegas_x_y_alphabet_version_03_2.csv
Total Instances:   22724
Neighbor Distance: 160
Min Prevalence:    0.1
Percentage Data:    100%
----------------------------------------
Execution Time: 12.345 s
Peak Memory Usage: 512 MB
Patterns Found: 37
----------------------------------------
[1] {A, B}
[2] {A, B, C}
...
```

---

## Ghi chú về cài đặt

- Đồ thị láng giềng **chỉ nối các instance khác feature type** — đúng với ngữ nghĩa colocation, nhưng nghĩa là đồ thị luôn là đa phần (multipartite).
- Clique có kích thước < 2 không được ghi vào hashmap ([maximal_clique_hashmap.cpp:59](src/maximal_clique_hashmap.cpp#L59)).
- `Miner::generateSubsets` chỉ sinh tập con bớt đúng một feature và dừng khi `|C| <= 2`, nên mẫu nhỏ nhất được xét là size-2.
- Lemma 2 (`deducePrevalentSubsets`): nếu `C` phổ biến thì mọi tập con của `C` chứa `f_min` cũng phổ biến — dùng để bỏ qua việc tính lại WPI.
- Khi `percentage_instances < 1.0`, việc lấy mẫu dùng `std::random_device` nên **kết quả giữa các lần chạy sẽ khác nhau**.
