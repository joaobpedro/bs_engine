#pragma once

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <type_traits>

#define DISCRETIZATION 50
#define NUMBER_OF_STATES 4

// Bend stiffener geometry and engineering design parameters
struct BendStiffner {
    float length;              // Total length L [m]
    float root_dia;            // Outer diameter at root D_root [m]
    float tip_dia;             // Outer diameter at tip D_tip [m]
    float inner_dia;           // Inner bore diameter D_i [m]
    float root_cyl_length = 0.0f;  // Cylindrical base section length L_cyl [m] (0 <= L_cyl < L)
    float pipe_EI_kNm2 = 15.0f;    // Internal umbilical/pipe bending stiffness EI_pipe [kN.m^2]
    float allowable_strain = 0.07f; // Max allowable polyurethane strain [-] (e.g. 0.07 = 7%)
    float allowable_mbr_m = 2.5f;   // Min allowable pipe bend radius MBR [m]
    bool use_arc_length = true;     // Integrate along physical arc length s in [0, L]
    bool use_far_field_bc = false;  // false: Prescribed Tip Angle; true: Far-Field Tension & Angle (OrcaFlex BC)
    float pipe_od = 0.110f;                 // Flexible pipe outer diameter D_pipe [m] (<= inner_dia)
    float insert_length = 0.0f;             // Embedded steel root insert sleeve length L_insert [m]
    float insert_dia = 0.220f;              // Embedded steel root insert outer diameter D_insert [m]
    float insert_EI_kNm2 = 450.0f;          // Steel root insert base flexural stiffness [kN.m^2]
    float allowable_bond_stress_mpa = 1.50f; // Allowable PU-to-steel insert interface bond stress [MPa]
};

inline bool is_valid_geometry(const BendStiffner &bs) {
    return std::isfinite(bs.length) &&
           std::isfinite(bs.root_dia) &&
           std::isfinite(bs.tip_dia) &&
           std::isfinite(bs.inner_dia) &&
           (bs.length > 1e-4f) &&
           (bs.inner_dia >= 0.0f) &&
           (bs.root_dia > bs.inner_dia + 1e-5f) &&
           (bs.tip_dia > bs.inner_dia + 1e-5f) &&
           (bs.root_dia >= bs.tip_dia - 1e-5f) &&
           (bs.root_cyl_length >= 0.0f) &&
           (bs.root_cyl_length < bs.length) &&
           (bs.allowable_strain > 1e-5f) &&
           (bs.allowable_mbr_m > 1e-4f);
}

// Effective inner diameter of the polyurethane sleeve at station s (accounting for embedded steel root insert)
inline float get_pu_inner_dia(const BendStiffner &bs, float s) {
    if (bs.insert_length > 1e-4f && s >= 0.0f && s < bs.insert_length) {
        float hi_ins = std::max(bs.inner_dia, bs.root_dia - 0.02f);
        float d_ins = std::clamp(bs.insert_dia, bs.inner_dia, hi_ins);
        // Smooth transition over the last 15% of insert length
        float t = s / bs.insert_length;
        if (t > 0.85f) {
            float alpha = (t - 0.85f) / 0.15f;
            return d_ins + (bs.inner_dia - d_ins) * alpha;
        }
        return d_ins;
    }
    return bs.inner_dia;
}

// Effective flexural stiffness [N.m^2] of the embedded steel root insert at station s
inline float get_insert_EI_Pa(const BendStiffner &bs, float s) {
    if (bs.insert_length <= 1e-4f || s < 0.0f || s >= bs.insert_length) return 0.0f;
    float t = std::clamp(s / bs.insert_length, 0.0f, 1.0f);
    // Slotted/fingered steel insert cone tapers its flexural contribution smoothly to 0 at s = L_insert
    return std::max(0.0f, bs.insert_EI_kNm2) * 1000.0f * (1.0f - t) * (1.0f - t);
}

// Calculates peak PU-to-Steel interface bond shear stress [MPa] along the embedded steel root insert
inline float get_insert_bond_stress_mpa(const BendStiffner &bs, float root_moment_Nm, float root_shear_N) {
    if (bs.insert_length <= 1e-4f) return 0.0f;
    constexpr float pi = 3.14159265358979323846f;
    float lo_ins = bs.inner_dia + 0.01f;
    float hi_ins = std::max(lo_ins, bs.root_dia - 0.02f);
    float d_ins = std::clamp(bs.insert_dia, lo_ins, hi_ins);
    float L_ins = bs.insert_length;
    float area_bond = pi * d_ins * L_ins;
    float z_bond = (pi / 4.0f) * d_ins * d_ins * L_ins;
    float tau_pa = (std::fabs(root_moment_Nm) / std::max(1e-6f, z_bond)) +
                   (std::fabs(root_shear_N) / std::max(1e-6f, area_bond));
    return tau_pa / 1e6f;
}

// Calculates annular clearance gap [mm], effective contact span [m], and free-span contact angle [deg]
inline float get_annular_gap_mm(const BendStiffner &bs) {
    if (bs.pipe_od <= 0.0f || bs.inner_dia <= bs.pipe_od) return 0.0f;
    return (bs.inner_dia - bs.pipe_od) * 1000.0f;
}

inline float get_gap_contact_span_m(const BendStiffner &bs) {
    return std::max(0.25f, (bs.root_cyl_length > 0.05f) ? bs.root_cyl_length : (0.15f * bs.length));
}

inline float get_gap_contact_angle_deg(const BendStiffner &bs) {
    if (bs.pipe_od <= 0.0f || bs.inner_dia <= bs.pipe_od) return 0.0f;
    float gap_m = bs.inner_dia - bs.pipe_od;
    float L_span = get_gap_contact_span_m(bs);
    return std::atan(gap_m / L_span) * (180.0f / 3.14159265358979323846f);
}

struct Vec2 {
    float x;
    float y;
};

// Stack-allocated state vector [y, theta, M, V] plus axial position x_pos for arc-length tracking
struct State4 {
    float items[NUMBER_OF_STATES]; // [0]=y, [1]=theta, [2]=M, [3]=V
    float x_pos = 0.0f;            // Deformed X coordinate [m]
};

template<typename T>
struct Vec1 {
    T* items;
    size_t count;
    size_t capacity;

    Vec1() : items(nullptr), count(0), capacity(0) {}

    explicit Vec1(size_t size) : items(nullptr), count(0), capacity(size) {
        if (capacity > 0) {
            items = new T[capacity]();
        }
    }

    Vec1(const Vec1& other) : items(nullptr), count(other.count), capacity(other.capacity) {
        if (capacity > 0) {
            items = new T[capacity]();
            for (size_t i = 0; i < count; ++i) {
                items[i] = other.items[i];
            }
        }
    }

    Vec1(Vec1&& other) noexcept : items(other.items), count(other.count), capacity(other.capacity) {
        other.items = nullptr;
        other.count = 0;
        other.capacity = 0;
    }

    Vec1& operator=(const Vec1& other) {
        if (this != &other) {
            if (capacity < other.capacity) {
                delete[] items;
                capacity = other.capacity;
                items = (capacity > 0) ? new T[capacity]() : nullptr;
            }
            count = other.count;
            for (size_t i = 0; i < count; ++i) {
                items[i] = other.items[i];
            }
        }
        return *this;
    }

    Vec1& operator=(Vec1&& other) noexcept {
        if (this != &other) {
            delete[] items;
            items = other.items;
            count = other.count;
            capacity = other.capacity;
            other.items = nullptr;
            other.count = 0;
            other.capacity = 0;
        }
        return *this;
    }

    ~Vec1() {
        delete[] items;
        items = nullptr;
        count = 0;
        capacity = 0;
    }

    T& operator[](size_t index) {
        assert(index < capacity && capacity > 0);
        return items[(index < capacity) ? index : (capacity > 0 ? capacity - 1 : 0)];
    }

    const T& operator[](size_t index) const {
        assert(index < capacity && capacity > 0);
        return items[(index < capacity) ? index : (capacity > 0 ? capacity - 1 : 0)];
    }

    void reserve(size_t new_capacity) {
        if (new_capacity <= capacity) return;
        T* new_items = new T[new_capacity]();
        for (size_t i = 0; i < count; ++i) {
            new_items[i] = std::move(items[i]);
        }
        delete[] items;
        items = new_items;
        capacity = new_capacity;
    }

    void append(const T& item) {
        if (count >= capacity) {
            size_t new_cap = (capacity == 0) ? 8 : (capacity * 2);
            reserve(new_cap);
        }
        items[count] = item;
        count++;
    }

    void clear() {
        if constexpr (std::is_arithmetic<T>::value) {
            if (items && count > 0) {
                std::memset(items, 0, count * sizeof(T));
            }
        } else {
            for (size_t i = 0; i < count; ++i) {
                items[i] = T();
            }
        }
        count = 0;
    }

    float max_pos() const {
        if constexpr (std::is_arithmetic<T>::value) {
            float max_element = 0.0f;
            if (count == 0) return 0.0f;
            for (size_t i = 0; i < count; i++) {
                if ((float)items[i] > max_element) {
                    max_element = (float)items[i];
                }
            }
            return max_element;
        } else {
            return 0.0f;
        }
    }

    float max_abs_signed(float fallback = 0.0f) const {
        if constexpr (std::is_arithmetic<T>::value) {
            if (count == 0) return fallback;
            float best = (float)items[0];
            for (size_t i = 1; i < count; ++i) {
                float v = (float)items[i];
                if (std::fabs(v) > std::fabs(best)) best = v;
            }
            return best;
        } else {
            return fallback;
        }
    }

    float min_val(float fallback = 0.0f) const {
        if constexpr (std::is_arithmetic<T>::value) {
            if (count == 0) return fallback;
            float m = (float)items[0];
            for (size_t i = 1; i < count; ++i) {
                float v = (float)items[i];
                if (v < m) m = v;
            }
            return m;
        } else {
            return fallback;
        }
    }

    float max_val(float fallback = 0.0f) const {
        if constexpr (std::is_arithmetic<T>::value) {
            if (count == 0) return fallback;
            float m = (float)items[0];
            for (size_t i = 1; i < count; ++i) {
                float v = (float)items[i];
                if (v > m) m = v;
            }
            return m;
        } else {
            return fallback;
        }
    }

    float min_pos(float fallback = 0.0f) const {
        if constexpr (std::is_arithmetic<T>::value) {
            if (count == 0) return fallback;
            float min_element = 1e30f;
            bool found = false;
            for (size_t i = 0; i < count; i++) {
                float v = (float)items[i];
                if (v > 0.0f && v < min_element) {
                    min_element = v;
                    found = true;
                }
            }
            return found ? min_element : fallback;
        } else {
            return fallback;
        }
    }
};

// Default polyurethane strain vs tangent modulus table (in kPa)
static constexpr int DEFAULT_MATERIAL_DATA_SIZE = 21;
static constexpr float DEFAULT_STRAIN_VECTOR[DEFAULT_MATERIAL_DATA_SIZE] = {
    0.00f, 0.01f, 0.02f, 0.03f, 0.04f, 0.05f, 0.06f, 0.07f, 0.08f, 0.09f,
    0.10f, 0.11f, 0.12f, 0.13f, 0.14f, 0.15f, 0.16f, 0.17f, 0.18f, 0.19f, 0.20f
};
static constexpr float DEFAULT_EMOD_VECTOR_KPA[DEFAULT_MATERIAL_DATA_SIZE] = {
    215800.0f, 215800.0f, 179500.0f, 138100.0f, 94300.0f, 59500.0f, 39700.0f,
    28200.0f,  21100.0f,  17000.0f,  14300.0f,  12700.0f, 11000.0f, 10600.0f,
    10000.0f,  9700.0f,   9600.0f,   8600.0f,   9400.0f,  8600.0f,  8500.0f
};

enum MaterialTempMode {
    MAT_TEMP_NOMINAL = 0,
    MAT_TEMP_COLD    = 1,
    MAT_TEMP_HOT     = 2
};

struct Material {
    bool is_nonlinear;
    bool table_is_tangent_modulus;       // true: table is d(sigma)/d(eps), integrate for stress & secant E
    bool use_cross_section_integration;  // true: integrate stress across annular depth y in [0, R_o]
    int active_temp_mode;                // MAT_TEMP_NOMINAL, MAT_TEMP_COLD, or MAT_TEMP_HOT

    float temp_nominal_c;
    float temp_cold_c;
    float temp_hot_c;

    float linear_E_nominal_kpa;
    float linear_E_cold_kpa;
    float linear_E_hot_kpa;
    float linear_E_kpa; // Active linear E [kPa]

    Vec1<float> strain_table;
    Vec1<float> emod_nominal_kpa;
    Vec1<float> emod_cold_kpa;
    Vec1<float> emod_hot_kpa;
    Vec1<float> emod_table_kpa;          // Active table E [kPa]
    Vec1<float> stress_table_kpa;        // Precomputed cumulative stress integral S(strain_table[i]) [kPa]
    bool is_uniform_grid = false;        // O(1) direct index lookup when strain_table is uniformly spaced
    float inv_grid_step = 100.0f;

    Material()
        : is_nonlinear(true),
          table_is_tangent_modulus(true),
          use_cross_section_integration(false),
          active_temp_mode(MAT_TEMP_NOMINAL),
          temp_nominal_c(20.0f),
          temp_cold_c(4.0f),
          temp_hot_c(35.0f),
          linear_E_nominal_kpa(215800.0f),
          linear_E_cold_kpa(215800.0f * 1.35f),
          linear_E_hot_kpa(215800.0f * 0.75f),
          linear_E_kpa(215800.0f),
          strain_table(DEFAULT_MATERIAL_DATA_SIZE),
          emod_nominal_kpa(DEFAULT_MATERIAL_DATA_SIZE),
          emod_cold_kpa(DEFAULT_MATERIAL_DATA_SIZE),
          emod_hot_kpa(DEFAULT_MATERIAL_DATA_SIZE),
          emod_table_kpa(DEFAULT_MATERIAL_DATA_SIZE),
          stress_table_kpa(DEFAULT_MATERIAL_DATA_SIZE) {
        reset_to_default();
    }

    void sync_active_table() {
        if (active_temp_mode == MAT_TEMP_COLD) {
            emod_table_kpa = emod_cold_kpa;
            linear_E_kpa = linear_E_cold_kpa;
        } else if (active_temp_mode == MAT_TEMP_HOT) {
            emod_table_kpa = emod_hot_kpa;
            linear_E_kpa = linear_E_hot_kpa;
        } else {
            emod_table_kpa = emod_nominal_kpa;
            linear_E_kpa = linear_E_nominal_kpa;
        }

        // Precompute cumulative trapezoidal stress table S[i] [kPa] and check uniform grid spacing for O(1) lookup
        size_t n = std::min(strain_table.count, emod_table_kpa.count);
        stress_table_kpa.clear();
        if (stress_table_kpa.capacity < n) {
            stress_table_kpa.reserve(n);
        }
        is_uniform_grid = false;
        inv_grid_step = 100.0f;

        if (n >= 2) {
            float dx_avg = (strain_table.items[n - 1] - strain_table.items[0]) / (float)(n - 1);
            if (dx_avg > 1e-6f) {
                bool uniform = true;
                for (size_t i = 0; i < n - 1; ++i) {
                    float dxi = strain_table.items[i + 1] - strain_table.items[i];
                    if (std::fabs(dxi - dx_avg) > 1e-5f) {
                        uniform = false;
                        break;
                    }
                }
                if (uniform) {
                    is_uniform_grid = true;
                    inv_grid_step = 1.0f / dx_avg;
                }
            }

            double acc = 0.0;
            if (strain_table.items[0] > 0.0f) {
                acc = (double)strain_table.items[0] * (double)std::max(1.0f, emod_table_kpa.items[0]);
            }
            stress_table_kpa.append((float)acc);
            for (size_t i = 0; i < n - 1; ++i) {
                double dx = std::max(0.0, (double)strain_table.items[i + 1] - (double)strain_table.items[i]);
                double y0 = std::max(1.0, (double)emod_table_kpa.items[i]);
                double y1 = std::max(1.0, (double)emod_table_kpa.items[i + 1]);
                acc += 0.5 * (y0 + y1) * dx;
                stress_table_kpa.append((float)acc);
            }
        }
    }

    void set_active_temp_mode(int mode) {
        active_temp_mode = std::clamp(mode, 0, 2);
        sync_active_table();
    }

    Material with_temp_mode(int mode) const {
        Material m = *this;
        m.set_active_temp_mode(mode);
        return m;
    }

    void reset_to_default(float nominal_scale = 1.0f, float cold_scale = 1.35f, float hot_scale = 0.75f) {
        strain_table.clear();
        emod_nominal_kpa.clear();
        emod_cold_kpa.clear();
        emod_hot_kpa.clear();
        for (int i = 0; i < DEFAULT_MATERIAL_DATA_SIZE; ++i) {
            strain_table.append(DEFAULT_STRAIN_VECTOR[i]);
            emod_nominal_kpa.append(DEFAULT_EMOD_VECTOR_KPA[i] * nominal_scale);
            emod_cold_kpa.append(DEFAULT_EMOD_VECTOR_KPA[i] * cold_scale);
            emod_hot_kpa.append(DEFAULT_EMOD_VECTOR_KPA[i] * hot_scale);
        }
        linear_E_nominal_kpa = DEFAULT_EMOD_VECTOR_KPA[0] * nominal_scale;
        linear_E_cold_kpa = DEFAULT_EMOD_VECTOR_KPA[0] * cold_scale;
        linear_E_hot_kpa = DEFAULT_EMOD_VECTOR_KPA[0] * hot_scale;
        sync_active_table();
    }

    const Vec1<float>& get_active_emod_table() const {
        if (active_temp_mode == MAT_TEMP_COLD && emod_cold_kpa.count == strain_table.count) return emod_cold_kpa;
        if (active_temp_mode == MAT_TEMP_HOT && emod_hot_kpa.count == strain_table.count) return emod_hot_kpa;
        if (emod_nominal_kpa.count == strain_table.count) return emod_nominal_kpa;
        return emod_table_kpa;
    }

    float get_active_linear_E_kpa() const {
        if (active_temp_mode == MAT_TEMP_COLD) return linear_E_cold_kpa;
        if (active_temp_mode == MAT_TEMP_HOT) return linear_E_hot_kpa;
        return linear_E_nominal_kpa;
    }

    bool has_valid_table() const {
        const Vec1<float>& active_t = get_active_emod_table();
        return strain_table.count >= 2 &&
               strain_table.count == active_t.count &&
               active_t.items[0] > 0.0f;
    }
};

// Piecewise linear interpolation of default non-linear E table (returns Pa)
inline float get_non_linear_E(float strain) {
    float s = std::fabs(strain);
    if (s <= DEFAULT_STRAIN_VECTOR[0]) {
        return DEFAULT_EMOD_VECTOR_KPA[0] * 1000.0f;
    }
    if (s >= DEFAULT_STRAIN_VECTOR[DEFAULT_MATERIAL_DATA_SIZE - 1]) {
        return DEFAULT_EMOD_VECTOR_KPA[DEFAULT_MATERIAL_DATA_SIZE - 1] * 1000.0f;
    }

    int idx = std::clamp((int)(s * 100.0f), 0, DEFAULT_MATERIAL_DATA_SIZE - 2);
    float x0 = DEFAULT_STRAIN_VECTOR[idx];
    float x1 = DEFAULT_STRAIN_VECTOR[idx + 1];
    float y0 = DEFAULT_EMOD_VECTOR_KPA[idx];
    float y1 = DEFAULT_EMOD_VECTOR_KPA[idx + 1];
    float dx = x1 - x0;
    float result = (dx > 1e-12f) ? (y0 + (s - x0) * (y1 - y0) / dx) : y0;
    return result * 1000.0f;
}

// Piecewise linear interpolation of custom non-linear E table (returns Pa)
inline float get_non_linear_E_custom(float strain, const Vec1<float> &strain_input_vector, const Vec1<float> &E_mod_input_vector) {
    assert(strain_input_vector.count == E_mod_input_vector.count);
    size_t n = std::min(strain_input_vector.count, E_mod_input_vector.count);
    if (n == 0) return DEFAULT_EMOD_VECTOR_KPA[0] * 1000.0f;
    if (n == 1) return E_mod_input_vector.items[0] * 1000.0f;

    float s = std::fabs(strain);
    const float* xs = strain_input_vector.items;
    const float* ys = E_mod_input_vector.items;
    if (s <= xs[0]) return ys[0] * 1000.0f;
    if (s >= xs[n - 1]) return ys[n - 1] * 1000.0f;

    const float* it = std::upper_bound(xs, xs + n, s);
    size_t i = (it > xs) ? (size_t)(it - xs - 1) : 0;
    if (i >= n - 1) i = n - 2;
    float x0 = xs[i];
    float x1 = xs[i + 1];
    float y0 = ys[i];
    float y1 = ys[i + 1];
    float dx = x1 - x0;
    float result = (dx > 1e-12f) ? (y0 + (s - x0) * (y1 - y0) / dx) : y0;
    return result * 1000.0f;
}

// Evaluates point-wise fiber Secant Young's Modulus E_sec(strain) in Pa (O(1) on uniform grids, O(log N) on general grids)
inline float get_fiber_E_secant(float strain, const Material &mat) {
    if (!mat.is_nonlinear || !mat.has_valid_table()) {
        return std::max(1.0f, mat.get_active_linear_E_kpa()) * 1000.0f;
    }

    float s = std::fabs(strain);
    const Vec1<float>& active_emod = mat.get_active_emod_table();
    size_t n = mat.strain_table.count;
    const float* xs = mat.strain_table.items;
    const float* ys = active_emod.items;

    if (!mat.table_is_tangent_modulus) {
        return get_non_linear_E_custom(s, mat.strain_table, active_emod);
    }

    if (s <= 1e-7f || s <= xs[0]) {
        return std::max(1.0f, ys[0]) * 1000.0f;
    }

    // Fast path using precomputed prefix stress table S[0..n-1]
    if (mat.stress_table_kpa.count == n) {
        const float* ss = mat.stress_table_kpa.items;
        if (s >= xs[n - 1]) {
            double y_last = std::max(1.0f, ys[n - 1]);
            double sigma_kpa = (double)ss[n - 1] + y_last * ((double)s - (double)xs[n - 1]);
            return (float)((sigma_kpa / (double)s) * 1000.0);
        }

        size_t i = 0;
        if (mat.is_uniform_grid) {
            int idx = (int)((s - xs[0]) * mat.inv_grid_step);
            i = (size_t)std::clamp(idx, 0, (int)n - 2);
        } else {
            const float* it = std::upper_bound(xs, xs + n, s);
            i = (it > xs) ? (size_t)(it - xs - 1) : 0;
            if (i >= n - 1) i = n - 2;
        }

        double x0 = xs[i];
        double x1 = xs[i + 1];
        double dx = x1 - x0;
        double y0 = std::max(1.0f, ys[i]);
        double y1 = std::max(1.0f, ys[i + 1]);
        double ds = (double)s - x0;
        double ys_interp = (dx > 1e-12) ? (y0 + (ds / dx) * (y1 - y0)) : y0;
        double sigma_kpa = (double)ss[i] + 0.5 * (y0 + ys_interp) * ds;
        return (float)((sigma_kpa / (double)s) * 1000.0);
    }

    // Fallback prefix sum if table was mutated without sync_active_table()
    double sigma_kpa = 0.0;
    if (xs[0] > 0.0f) {
        double clamp_s = std::min((double)s, (double)xs[0]);
        sigma_kpa += clamp_s * (double)ys[0];
    }
    for (size_t i = 0; i < n - 1; ++i) {
        double x0 = xs[i];
        double x1 = xs[i + 1];
        double dx = x1 - x0;
        if (dx <= 1e-12) continue;
        double y0 = std::max(1.0f, ys[i]);
        double y1 = std::max(1.0f, ys[i + 1]);
        if (s >= x1) {
            sigma_kpa += 0.5 * (y0 + y1) * dx;
        } else if (s > x0) {
            double ds = s - x0;
            double ys_interp = y0 + (ds / dx) * (y1 - y0);
            sigma_kpa += 0.5 * (y0 + ys_interp) * ds;
            break;
        }
    }
    if (s > xs[n - 1]) {
        double y_last = std::max(1.0f, ys[n - 1]);
        sigma_kpa += y_last * ((double)s - (double)xs[n - 1]);
    }
    return (float)((sigma_kpa / (double)s) * 1000.0);
}

// Evaluates effective cross-sectional Young's Modulus E(strain) in Pa (either outer-fiber secant or annular depth integration)
inline float get_material_E(float strain, const Material &mat) {
    if (!mat.is_nonlinear || !mat.has_valid_table() || !mat.use_cross_section_integration) {
        return get_fiber_E_secant(strain, mat);
    }
    // 4-point moment-of-area weighted quadrature across normalized depth xi in (0, 1]
    static constexpr float xi_q[4] = { 0.32425342f, 0.61337143f, 0.83603111f, 0.96816024f };
    static constexpr float w_q[4]  = { 0.10122854f, 0.31234708f, 0.43392346f, 0.15250092f };
    float s = std::fabs(strain);
    float e_acc = 0.0f;
    for (int q = 0; q < 4; ++q) {
        e_acc += w_q[q] * get_fiber_E_secant(s * xi_q[q], mat);
    }
    return e_acc;
}

// Analytically inverts the monotonic stress-strain law sigma(eps) = eps * E_eff(eps) = sigma_req_pa
inline float solve_strain_from_stress(double sigma_req_pa, const Material &mat) {
    if (sigma_req_pa <= 0.0) return 0.0f;

    if (!mat.is_nonlinear || !mat.has_valid_table()) {
        double E_pa = std::max(1.0, (double)mat.get_active_linear_E_kpa()) * 1000.0;
        return (float)(sigma_req_pa / E_pa);
    }

    const Vec1<float>& active_emod = mat.get_active_emod_table();
    size_t n = mat.strain_table.count;
    const float* xs = mat.strain_table.items;
    const float* ys = active_emod.items;
    double sigma_kpa = sigma_req_pa / 1000.0;

    float eps_est = 0.0f;
    if (!mat.table_is_tangent_modulus) {
        double prev_x = 0.0;
        double prev_sig = 0.0;
        double min_slope = std::max(1.0, (double)ys[n - 1] * 0.1);
        bool found = false;
        for (size_t i = 0; i < n; ++i) {
            double xi = xs[i];
            double yi = std::max(1.0, (double)ys[i]);
            if (xi <= prev_x) continue;
            double sig_i = std::max(xi * yi, prev_sig + min_slope * (xi - prev_x));
            if (sigma_kpa <= sig_i) {
                double t = (sigma_kpa - prev_sig) / (sig_i - prev_sig);
                eps_est = (float)(prev_x + t * (xi - prev_x));
                found = true;
                break;
            }
            prev_x = xi;
            prev_sig = sig_i;
        }
        if (!found) {
            double last_slope = std::max(1.0, (double)ys[n - 1]);
            eps_est = (float)(prev_x + (sigma_kpa - prev_sig) / last_slope);
        }
    } else if (mat.stress_table_kpa.count == n) {
        // Fast O(log N) inversion using precomputed cumulative stress table S[0..n-1]
        const float* ss = mat.stress_table_kpa.items;
        if (xs[0] > 0.0f && sigma_kpa <= (double)ss[0]) {
            double y0 = std::max(1.0, (double)ys[0]);
            eps_est = (float)(sigma_kpa / y0);
        } else if (sigma_kpa >= (double)ss[n - 1]) {
            double y_last = std::max(1.0, (double)ys[n - 1]);
            eps_est = (float)((double)xs[n - 1] + (sigma_kpa - (double)ss[n - 1]) / y_last);
        } else {
            const float* it = std::upper_bound(ss, ss + n, (float)sigma_kpa);
            size_t i = (it > ss) ? (size_t)(it - ss - 1) : 0;
            if (i >= n - 1) i = n - 2;
            double x0 = xs[i];
            double x1 = xs[i + 1];
            double dx = x1 - x0;
            if (dx <= 1e-12) {
                eps_est = (float)x0;
            } else {
                double y0 = std::max(1.0, (double)ys[i]);
                double y1 = std::max(1.0, (double)ys[i + 1]);
                double rem = std::max(0.0, sigma_kpa - (double)ss[i]);
                double slope = (y1 - y0) / dx;
                double ds = 0.0;
                if (std::fabs(slope) < 1e-9) {
                    ds = rem / y0;
                } else {
                    double disc = y0 * y0 + 2.0 * slope * rem;
                    ds = (-y0 + std::sqrt(std::max(0.0, disc))) / slope;
                }
                eps_est = (float)(x0 + std::clamp(ds, 0.0, dx));
            }
        }
    } else {
        double acc_sigma = 0.0;
        bool found = false;
        if (xs[0] > 0.0f) {
            double y0 = std::max(1.0, (double)ys[0]);
            double d_sig0 = (double)xs[0] * y0;
            if (sigma_kpa <= d_sig0) {
                eps_est = (float)(sigma_kpa / y0);
                found = true;
            }
            acc_sigma += d_sig0;
        }

        if (!found) {
            for (size_t i = 0; i < n - 1; ++i) {
                double x0 = xs[i];
                double x1 = xs[i + 1];
                double dx = x1 - x0;
                if (dx <= 1e-12) continue;
                double y0 = std::max(1.0, (double)ys[i]);
                double y1 = std::max(1.0, (double)ys[i + 1]);
                double d_sigma = 0.5 * (y0 + y1) * dx;

                if (sigma_kpa <= acc_sigma + d_sigma) {
                    double rem = sigma_kpa - acc_sigma;
                    double slope = (y1 - y0) / dx;
                    double ds = 0.0;
                    if (std::fabs(slope) < 1e-9) {
                        ds = rem / y0;
                    } else {
                        double disc = y0 * y0 + 2.0 * slope * rem;
                        ds = (-y0 + std::sqrt(std::max(0.0, disc))) / slope;
                    }
                    eps_est = (float)(x0 + std::clamp(ds, 0.0, dx));
                    found = true;
                    break;
                }
                acc_sigma += d_sigma;
            }
        }
        if (!found) {
            double y_last = std::max(1.0, (double)ys[n - 1]);
            eps_est = (float)((double)xs[n - 1] + (sigma_kpa - acc_sigma) / y_last);
        }
    }

    // If cross-sectional annular integration is active, refine eps_est with 3 fast secant/Newton steps
    if (mat.use_cross_section_integration) {
        double e = std::max(1e-7, (double)eps_est);
        for (int k = 0; k < 3; ++k) {
            double E_curr = (double)get_material_E((float)e, mat);
            double f = e * E_curr - sigma_req_pa;
            if (std::fabs(f) < 1.0) break;
            double de = std::max(1e-6, e * 0.01);
            double f_plus = (e + de) * (double)get_material_E((float)(e + de), mat) - sigma_req_pa;
            double df = (f_plus - f) / de;
            if (df > 1.0) {
                e = std::max(0.0, e - f / df);
            }
        }
        return (float)e;
    }

    return eps_est;
}

// Outer diameter at station x (supports cylindrical root section L_cyl followed by conical taper)
inline float get_dia(const BendStiffner &dimensions, float x) {
    if (dimensions.length <= 1e-6f) return dimensions.root_dia;
    float s = std::clamp(x, 0.0f, dimensions.length);
    float max_cyl = std::max(0.0f, dimensions.length - 1e-5f);
    float cyl_len = std::clamp(dimensions.root_cyl_length, 0.0f, max_cyl);
    if (s <= cyl_len) {
        return dimensions.root_dia;
    }
    float cone_len = std::max(1e-6f, dimensions.length - cyl_len);
    float t = std::clamp((s - cyl_len) / cone_len, 0.0f, 1.0f);
    return dimensions.root_dia - (dimensions.root_dia - dimensions.tip_dia) * t;
}

inline float get_Inertia(const BendStiffner& dimensions, float x) {
    constexpr float pi_over_64 = 3.14159265358979323846f / 64.0f;
    float outer_dia = get_dia(dimensions, x);
    float inner_dia = std::min(outer_dia - 0.01f, get_pu_inner_dia(dimensions, x));
    float do2 = outer_dia * outer_dia;
    float di2 = inner_dia * inner_dia;
    float inertia = pi_over_64 * (do2 * do2 - di2 * di2);
    return std::max(1e-12f, inertia);
}

inline float get_Area(const BendStiffner& dimensions, float x) {
    constexpr float pi_over_4 = 3.14159265358979323846f / 4.0f;
    float outer_dia = get_dia(dimensions, x);
    float inner_dia = std::min(outer_dia - 0.01f, get_pu_inner_dia(dimensions, x));
    float area = pi_over_4 * (outer_dia * outer_dia - inner_dia * inner_dia);
    return std::max(1e-10f, area);
}

inline float get_EI(const BendStiffner& dimensions, float strain, float x) {
    float inertia = get_Inertia(dimensions, x);
    float passive_EI_Pa = std::max(0.0f, dimensions.pipe_EI_kNm2) * 1000.0f + get_insert_EI_Pa(dimensions, x);
    return get_non_linear_E(strain) * inertia + passive_EI_Pa;
}

inline float get_EI(const BendStiffner& dimensions, float strain, float x, const Material &mat) {
    float inertia = get_Inertia(dimensions, x);
    float passive_EI_Pa = std::max(0.0f, dimensions.pipe_EI_kNm2) * 1000.0f + get_insert_EI_Pa(dimensions, x);
    return get_material_E(strain, mat) * inertia + passive_EI_Pa;
}

// Solves the coupled composite cross-section (Polyurethane Cone + Internal Pipe EI_pipe + Steel Root Insert) at station s
// Given total bending moment M_total [N.m] and axial tension N_axial [N], returns:
//   out_curvature [1/m], out_total_strain [-], out_EI_total [N.m^2]
inline void solve_composite_section(const BendStiffner &bs,
                                    const Material &mat,
                                    float s_pos,
                                    float m_total_Nm,
                                    float axial_N,
                                    float &out_curvature,
                                    float &out_total_strain,
                                    float &out_EI_total) {
    float d_o = get_dia(bs, s_pos);
    float r_o = 0.5f * d_o;
    float I_bs = get_Inertia(bs, s_pos);
    float A_bs = get_Area(bs, s_pos);
    float passive_EI = std::max(0.0f, bs.pipe_EI_kNm2) * 1000.0f + get_insert_EI_Pa(bs, s_pos); // N.m^2

    double abs_M = std::fabs((double)m_total_Nm);
    double sigma_axial = (double)std::fabs(axial_N) / (double)A_bs;

    if (abs_M <= 1e-6 && sigma_axial <= 1e-4) {
        out_EI_total = get_material_E(0.0f, mat) * I_bs + passive_EI;
        out_curvature = 0.0f;
        out_total_strain = 0.0f;
        return;
    }

    if (passive_EI <= 1e-3f) {
        // Pure bend stiffener section
        double sigma_req = (abs_M * (double)r_o) / (double)I_bs + sigma_axial;
        float eps = solve_strain_from_stress(sigma_req, mat);
        float E_bs = get_material_E(eps, mat);
        out_EI_total = E_bs * I_bs;
        out_curvature = m_total_Nm / out_EI_total;
        out_total_strain = eps;
        return;
    }

    // Composite section: M_total = (E_bs(eps) * I_bs + EI_passive) * kappa, with eps = |kappa| * r_o + eps_axial
    // Initialize from uncoupled estimate and refine in up to 3 fast iterations
    double E0 = (double)get_material_E(0.01f, mat);
    double frac_bs = (E0 * I_bs) / (E0 * I_bs + passive_EI);
    double sigma_req = (abs_M * frac_bs * (double)r_o) / (double)I_bs + sigma_axial;
    float eps = solve_strain_from_stress(sigma_req, mat);

    for (int k = 0; k < 3; ++k) {
        double E_bs = (double)get_material_E(eps, mat);
        double EI_tot = E_bs * (double)I_bs + (double)passive_EI;
        double kappa_abs = abs_M / EI_tot;
        double sigma_b = E_bs * kappa_abs * (double)r_o;
        float next_eps = solve_strain_from_stress(sigma_b + sigma_axial, mat);
        if (std::fabs(next_eps - eps) <= 1e-6f) {
            eps = next_eps;
            break;
        }
        eps = next_eps;
    }

    float E_final = get_material_E(eps, mat);
    out_EI_total = E_final * I_bs + passive_EI;
    out_curvature = m_total_Nm / out_EI_total;
    out_total_strain = eps;
}

// Fast zero-allocation ODE right-hand side supporting both Arc-Length (s) and Cartesian (x) formulations
inline State4 equations_fast(float s_or_x,
                             const State4 &y,
                             const BendStiffner &bs,
                             const Material &mat,
                             float tension_N = 0.0f,
                             float target_angle_rad = 0.0f) {
    float T_abs = std::fabs(tension_N);
    float angle_diff = target_angle_rad - y.items[1];
    float cos_diff = std::cos(angle_diff);
    float axial_N = T_abs * std::max(0.0f, cos_diff);
    float kappa = 0.0f, local_strain = 0.0f, EI_total = 1.0f;
    solve_composite_section(bs, mat, s_or_x, y.items[2], axial_N, kappa, local_strain, EI_total);

    bool far_field = bs.use_far_field_bc && (T_abs > 1e-2f);
    float shear_v = far_field
        ? (-T_abs * std::sin(angle_diff))
        : y.items[3];

    State4 ds;
    if (bs.use_arc_length) {
        // Exact Elastica in arc-length s in [0, L]: singularity-free at all angles!
        ds.x_pos = std::cos(y.items[1]);
        ds.items[0] = std::sin(y.items[1]);          // dy/ds = sin(theta)
        ds.items[1] = kappa;                         // dtheta/ds = M / EI_total
        ds.items[2] = shear_v;                       // dM/ds = V(s)
        ds.items[3] = far_field ? (T_abs * cos_diff * kappa) : 0.0f; // dV/ds = T * cos(phi_inf - theta) * dtheta/ds
    } else {
        // Legacy Cartesian x-projection formulation
        constexpr float max_theta = 1.553343f;
        float theta_clamped = std::clamp(y.items[1], -max_theta, max_theta);
        float cos_t = std::cos(theta_clamped);
        float tan_t = std::tan(theta_clamped);
        ds.x_pos = 1.0f;
        ds.items[0] = tan_t;
        ds.items[1] = kappa / cos_t;
        ds.items[2] = shear_v;
        ds.items[3] = far_field ? (T_abs * cos_diff * ds.items[1]) : 0.0f;
    }
    return ds;
}

// Fast zero-allocation 4th-order Runge-Kutta step
inline State4 RK4_fast(float s_or_x,
                       const State4 &y,
                       float h,
                       const BendStiffner &bs,
                       const Material &mat,
                       float tension_N = 0.0f,
                       float target_angle_rad = 0.0f) {
    State4 k1 = equations_fast(s_or_x, y, bs, mat, tension_N, target_angle_rad);

    State4 y_temp;
    for (int i = 0; i < NUMBER_OF_STATES; ++i) {
        y_temp.items[i] = y.items[i] + 0.5f * h * k1.items[i];
    }
    y_temp.x_pos = y.x_pos + 0.5f * h * k1.x_pos;
    State4 k2 = equations_fast(s_or_x + 0.5f * h, y_temp, bs, mat, tension_N, target_angle_rad);

    for (int i = 0; i < NUMBER_OF_STATES; ++i) {
        y_temp.items[i] = y.items[i] + 0.5f * h * k2.items[i];
    }
    y_temp.x_pos = y.x_pos + 0.5f * h * k2.x_pos;
    State4 k3 = equations_fast(s_or_x + 0.5f * h, y_temp, bs, mat, tension_N, target_angle_rad);

    for (int i = 0; i < NUMBER_OF_STATES; ++i) {
        y_temp.items[i] = y.items[i] + h * k3.items[i];
    }
    y_temp.x_pos = y.x_pos + h * k3.x_pos;
    State4 k4 = equations_fast(s_or_x + h, y_temp, bs, mat, tension_N, target_angle_rad);

    State4 y_next;
    for (int i = 0; i < NUMBER_OF_STATES; ++i) {
        y_next.items[i] = y.items[i] + (h / 6.0f) * (k1.items[i] + 2.0f * k2.items[i] + 2.0f * k3.items[i] + k4.items[i]);
    }
    y_next.x_pos = y.x_pos + (h / 6.0f) * (k1.x_pos + 2.0f * k2.x_pos + 2.0f * k3.x_pos + k4.x_pos);
    return y_next;
}

inline State4 shoot_fast(const BendStiffner &bs,
                         size_t steps,
                         float y0,
                         float theta0,
                         float guessed_m0,
                         float guessed_v0,
                         const Material &mat,
                         float tension_N = 0.0f,
                         float target_angle_rad = 0.0f) {
    float h = bs.length / (float)steps;
    float s = 0.0f;
    State4 y = { { y0, theta0, guessed_m0, guessed_v0 }, 0.0f };
    for (size_t i = 0; i < steps; ++i) {
        y = RK4_fast(s, y, h, bs, mat, tension_N, target_angle_rad);
        s += h;
    }
    return y;
}

// Legacy Vec1-based ODE functions preserved for backwards compatibility
inline Vec1<float> equations(float x, const Vec1<float> &y, const BendStiffner bs, float strain) {
    Vec1<float> dydx(NUMBER_OF_STATES);
    float EI = get_EI(bs, strain, x);
    constexpr float max_theta = 1.553343f;
    float theta_clamped = std::clamp(y.items[1], -max_theta, max_theta);
    dydx.append(std::tan(theta_clamped));
    dydx.append(y.items[2] / (EI * std::cos(theta_clamped)));
    dydx.append(y.items[3]);
    dydx.append(0.0f);
    return dydx;
}

inline Vec1<float> RK4(float x, const Vec1<float> &y, float h, const BendStiffner &bs, float strain) {
    (void)strain;
    Material default_mat;
    State4 s = { { y.items[0], y.items[1], y.items[2], y.items[3] }, x };
    State4 s_next = RK4_fast(x, s, h, bs, default_mat, 0.0f, 0.0f);
    Vec1<float> y_next(NUMBER_OF_STATES);
    for (int i = 0; i < NUMBER_OF_STATES; ++i) {
        y_next.append(s_next.items[i]);
    }
    return y_next;
}

inline Vec1<float> shoot(const BendStiffner &bs,
                         size_t steps,
                         float y0,
                         float theta0,
                         float guessed_m0,
                         float guessed_v0,
                         const Vec1<float> &strain) {
    (void)strain;
    Material default_mat;
    State4 s_end = shoot_fast(bs, steps, y0, theta0, guessed_m0, guessed_v0, default_mat, 0.0f, 0.0f);
    Vec1<float> y(NUMBER_OF_STATES);
    for (int i = 0; i < NUMBER_OF_STATES; ++i) {
        y.append(s_end.items[i]);
    }
    return y;
}

// Solves the BVP with Ridders' bracketed root-finding method (supports Prescribed Tip Angle & Far-Field Tension Elastica BCs)
inline Vec2 solve_bvp_material(const BendStiffner &bs,
                               size_t steps,
                               float y0,
                               float theta0,
                               float targetML,
                               float thetaL,
                               const Material &mat,
                               float tension_N = 0.0f) {
    if (!is_valid_geometry(bs) || steps == 0) {
        return { 0.0f, 0.0f };
    }

    double dtheta = (double)thetaL - (double)theta0;
    if (std::fabs(dtheta) < 1e-8 && std::fabs(targetML) < 1e-8) {
        return { 0.0f, 0.0f };
    }

    double sgn = (dtheta >= 0.0) ? 1.0 : -1.0;
    double target_abs = std::fabs(dtheta);
    bool far_field = bs.use_far_field_bc && (std::fabs(tension_N) > 1e-2f);
    double pipe_EI_Pa = std::max(0.0, (double)bs.pipe_EI_kNm2) * 1000.0;
    double T_abs = std::fabs((double)tension_N);

    auto v0_for_m0 = [&](double m0) -> double {
        if (far_field) {
            return -T_abs * std::sin((double)thetaL - (double)theta0);
        }
        return ((double)targetML - m0) / (double)bs.length;
    };

    auto eval_residual = [&](double m_abs) -> double {
        double m0 = sgn * m_abs;
        double v0 = v0_for_m0(m0);
        State4 end_state = shoot_fast(bs, steps, y0, theta0, (float)m0, (float)v0, mat, tension_N, thetaL);

        if (far_field) {
            // Exact tensioned-elastica tip moment continuity at s = L:
            // M(L) = 2 * sqrt(EI_pipe * T) * sin((phi_inf - theta(L)) / 2)
            constexpr double half_pi = 1.5707963267948966;
            double half_angle_diff = std::clamp(0.5 * ((double)thetaL - (double)end_state.items[1]), -half_pi, half_pi);
            double m_tip_target = 2.0 * std::sqrt(pipe_EI_Pa * T_abs) * std::sin(half_angle_diff);
            return sgn * ((double)end_state.items[2] - m_tip_target);
        } else {
            return sgn * ((double)end_state.items[1] - (double)theta0) - target_abs;
        }
    };

    double a = 0.0;
    double fa = eval_residual(a);
    if (std::fabs(fa) < 1e-7) {
        return { 0.0f, (float)v0_for_m0(0.0) };
    }

    double EI_mid = (double)get_EI(bs, 0.02f, 0.5f * bs.length, mat);
    double b = far_field
        ? std::max(100.0, 0.5 * T_abs * (double)bs.length * std::fabs(std::sin(dtheta)))
        : std::max(100.0, 0.8 * EI_mid * target_abs / (double)bs.length);
    double fb = eval_residual(b);

    for (int expand = 0; expand < 30 && fb < 0.0; ++expand) {
        a = b;
        fa = fb;
        b *= 1.6;
        fb = eval_residual(b);
    }

    // Ridders' bracketed root-finding method on [x_lo, x_hi] where f_lo <= 0 <= f_hi
    double x_lo = a, f_lo = fa;
    double x_hi = b, f_hi = fb;
    double m_best = (std::fabs(f_lo) < std::fabs(f_hi)) ? x_lo : x_hi;
    double f_best = (std::fabs(f_lo) < std::fabs(f_hi)) ? std::fabs(f_lo) : std::fabs(f_hi);
    double tol = far_field ? 1e-3 : 1e-6;

    for (int k = 0; k < 28; ++k) {
        if (f_best < tol || std::fabs(x_hi - x_lo) <= 1e-6 * (1.0 + std::fabs(m_best))) {
            break;
        }

        double x_mid = 0.5 * (x_lo + x_hi);
        double f_mid = eval_residual(x_mid);
        if (std::fabs(f_mid) < f_best) {
            f_best = std::fabs(f_mid);
            m_best = x_mid;
        }
        if (f_best < tol) break;

        double denom = std::sqrt(std::max(0.0, f_mid * f_mid - f_lo * f_hi));
        if (denom <= 1e-14) {
            if (f_mid > 0.0) { x_hi = x_mid; f_hi = f_mid; }
            else { x_lo = x_mid; f_lo = f_mid; }
            continue;
        }

        double x_new = x_mid - (x_mid - x_lo) * (f_mid / denom);
        double f_new = eval_residual(x_new);
        if (std::fabs(f_new) < f_best) {
            f_best = std::fabs(f_new);
            m_best = x_new;
        }
        if (f_best < tol) break;

        if (f_mid * f_new < 0.0) {
            if (x_mid < x_new) {
                x_lo = x_mid; f_lo = f_mid;
                x_hi = x_new; f_hi = f_new;
            } else {
                x_lo = x_new; f_lo = f_new;
                x_hi = x_mid; f_hi = f_mid;
            }
        } else {
            if (f_mid > 0.0) { x_hi = x_mid; f_hi = f_mid; }
            else { x_lo = x_mid; f_lo = f_mid; }

            if (f_new > 0.0 && x_new < x_hi) { x_hi = x_new; f_hi = f_new; }
            else if (f_new < 0.0 && x_new > x_lo) { x_lo = x_new; f_lo = f_new; }
        }
    }

    double final_m0 = sgn * m_best;
    double final_v0 = v0_for_m0(final_m0);
    return { (float)final_m0, (float)final_v0 };
}

inline float solve_V0(const BendStiffner &bs,
                      size_t steps,
                      float y0,
                      float theta0,
                      float curr_guessed_m0,
                      float thetaL,
                      float guessed_v0,
                      const Vec1<float> &strain) {
    (void)steps; (void)y0; (void)theta0; (void)thetaL; (void)guessed_v0; (void)strain;
    if (bs.length <= 1e-6f) return 0.0f;
    return -curr_guessed_m0 / bs.length;
}

inline Vec2 solve_bvp(const BendStiffner &bs,
                      size_t steps,
                      float y0,
                      float theta0,
                      float targetML,
                      float thetaL,
                      float guessed_m0,
                      float guessed_v0,
                      const Vec1<float> &strain) {
    (void)guessed_m0; (void)guessed_v0; (void)strain;
    Material default_mat;
    return solve_bvp_material(bs, steps, y0, theta0, targetML, thetaL, default_mat, 0.0f);
}

// Local pipe curvature [1/m] accounting for annular clearance gap contact transition inside root entry span
inline float get_local_pipe_curvature(const BendStiffner &bs, float s, float abs_kappa_beam, float abs_angle_rad) {
    if (bs.pipe_od <= 0.0f || bs.inner_dia <= bs.pipe_od + 1e-4f || abs_angle_rad <= 1e-4f) {
        return abs_kappa_beam;
    }
    float gap_m = bs.inner_dia - bs.pipe_od;
    float L_span = get_gap_contact_span_m(bs);
    if (s >= L_span) return abs_kappa_beam;
    float theta_gap_rad = std::atan(gap_m / L_span);
    float activation = std::clamp(abs_angle_rad / std::max(1e-4f, theta_gap_rad), 0.0f, 1.0f);
    float shape = 1.0f - (s / L_span);
    float delta_kappa = 0.25f * (gap_m / (L_span * L_span)) * activation * shape;
    return abs_kappa_beam + delta_kappa;
}

// Calculates total outer-fiber strain, curvature [1/m], and bend radius [m] at all discretization nodes
inline void calculate_strain_material(size_t steps,
                                      const BendStiffner &bs,
                                      const Vec1<Vec1<float>> &results,
                                      Vec1<float> &strain,
                                      const Material &mat,
                                      float tension_N = 0.0f,
                                      float tip_angle_rad = 0.0f,
                                      Vec1<float>* out_curvature = nullptr,
                                      Vec1<float>* out_bend_radius = nullptr) {
    size_t n = std::min(steps, results.count);
    if (n == 0) return;
    if (strain.capacity < n) strain.reserve(n);
    if (out_curvature) {
        out_curvature->clear();
        if (out_curvature->capacity < n) out_curvature->reserve(n);
    }
    if (out_bend_radius) {
        out_bend_radius->clear();
        if (out_bend_radius->capacity < n) out_bend_radius->reserve(n);
    }

    float h = (n > 1) ? (bs.length / (float)(n - 1)) : 0.0f;
    float s = 0.0f;
    float abs_tip_angle = std::fabs(tip_angle_rad);

    for (size_t i = 0; i < n; ++i) {
        float theta_s = results.items[i].items[1];
        float m_s = results.items[i].items[2];
        float axial_N = std::fabs(tension_N) * std::max(0.0f, std::cos(tip_angle_rad - theta_s));

        float kappa = 0.0f, total_eps = 0.0f, EI_tot = 1.0f;
        solve_composite_section(bs, mat, s, m_s, axial_N, kappa, total_eps, EI_tot);

        strain.items[i] = total_eps;
        float abs_kappa = get_local_pipe_curvature(bs, s, std::fabs(kappa), abs_tip_angle);
        // In Far-Field Tension mode, also check bare pipe curvature right outside the tip exit (s = L+)
        if (i + 1 == n && bs.use_far_field_bc && std::fabs(tension_N) > 1e-2f && bs.pipe_EI_kNm2 > 1e-4f) {
            float kappa_pipe_exit = std::fabs(m_s) / (bs.pipe_EI_kNm2 * 1000.0f);
            abs_kappa = std::max(abs_kappa, kappa_pipe_exit);
        }
        if (out_curvature) {
            out_curvature->append(abs_kappa);
        }
        if (out_bend_radius) {
            float mbr = (abs_kappa > 1e-6f) ? (1.0f / abs_kappa) : 1000.0f;
            out_bend_radius->append(std::min(1000.0f, mbr));
        }
        s += h;
    }
    strain.count = n;
}

inline void calculate_strain(size_t steps,
                             const BendStiffner &bs,
                             const Vec1<Vec1<float>> &results,
                             Vec1<float> &strain) {
    Material default_mat;
    calculate_strain_material(steps, bs, results, strain, default_mat, 0.0f, 0.0f, nullptr, nullptr);
}

// Calculates Polyurethane volume [m^3] (cylindrical root + conical frustum - inner bore - embedded steel insert)
inline float get_pu_volume_m3(const BendStiffner &bs) {
    if (!is_valid_geometry(bs)) return 0.0f;
    constexpr float pi = 3.14159265358979323846f;
    float l_cyl = std::clamp(bs.root_cyl_length, 0.0f, bs.length);
    float l_cone = std::max(0.0f, bs.length - l_cyl);
    float dr = bs.root_dia;
    float dt = bs.tip_dia;
    float di = bs.inner_dia;

    float v_cyl_out = (pi / 4.0f) * dr * dr * l_cyl;
    float v_cone_out = (pi / 12.0f) * l_cone * (dr * dr + dr * dt + dt * dt);
    float v_bore = (pi / 4.0f) * di * di * bs.length;
    float v_insert = 0.0f;
    if (bs.insert_length > 1e-4f) {
        float l_ins = std::min(bs.insert_length, bs.length);
        float hi_ins = std::max(di, dr - 0.02f);
        float d_ins = std::clamp(bs.insert_dia, di, hi_ins);
        v_insert = (pi / 4.0f) * std::max(0.0f, d_ins * d_ins - di * di) * l_ins * 0.925f;
    }
    return std::max(0.0f, v_cyl_out + v_cone_out - v_bore - v_insert);
}

inline float get_pu_mass_kg(const BendStiffner &bs, float density_kg_m3 = 1150.0f) {
    return get_pu_volume_m3(bs) * density_kg_m3;
}

struct OptimizationConfig {
    float min_length = 1.0f;
    float max_length = 6.0f;
    float min_root_dia = 0.25f;
    float max_root_dia = 1.20f;
    float cyl_length_ratio = 0.12f;        // L_cyl / L
    float safety_margin = 0.95f;           // Target utilization <= 95% of allowable limits
    float min_root_stiffness_ratio = 15.0f; // Min EI_bs(0) / EI_pipe when T == 0
    bool optimize_length = true;
    bool optimize_root_dia = true;
    bool use_all_pairs = false;
    bool force_far_field_when_tension = true;
    bool include_temp_envelope = false;    // Check Cold (1.35x), Nominal (1.0x), and Warm (0.75x) during optimization
};

inline Material scale_material(const Material &base_mat, float factor) {
    Material m = base_mat;
    m.linear_E_nominal_kpa = base_mat.linear_E_nominal_kpa * factor;
    m.linear_E_cold_kpa = base_mat.linear_E_cold_kpa * factor;
    m.linear_E_hot_kpa = base_mat.linear_E_hot_kpa * factor;
    for (size_t i = 0; i < m.emod_nominal_kpa.count; ++i) m.emod_nominal_kpa.items[i] = base_mat.emod_nominal_kpa.items[i] * factor;
    for (size_t i = 0; i < m.emod_cold_kpa.count; ++i)    m.emod_cold_kpa.items[i]    = base_mat.emod_cold_kpa.items[i] * factor;
    for (size_t i = 0; i < m.emod_hot_kpa.count; ++i)     m.emod_hot_kpa.items[i]     = base_mat.emod_hot_kpa.items[i] * factor;
    m.sync_active_table();
    return m;
}

#define OPT_GRID_N 11

struct OptimizationResult {
    bool converged = false;
    bool feasible = false;
    BendStiffner best_bs = {};
    float volume_m3 = 0.0f;
    float mass_kg = 0.0f;
    float max_strain = 0.0f;
    float min_mbr_m = 1000.0f;
    float max_root_moment_kNm = 0.0f;
    float max_bond_stress_mpa = 0.0f;
    float strain_utilization = 0.0f;
    float mbr_utilization = 0.0f;
    float bond_utilization = 0.0f;
    int evaluations = 0;

    // 2D Design Space (L vs D_root) grid from Stage-1 global exploration
    int grid_nL = 0;
    int grid_nDr = 0;
    float grid_L[OPT_GRID_N] = { 0.0f };
    float grid_Dr[OPT_GRID_N] = { 0.0f };
    bool grid_feasible[OPT_GRID_N][OPT_GRID_N] = { { false } };
    float grid_max_strain[OPT_GRID_N][OPT_GRID_N] = { { 0.0f } };
    float grid_min_mbr[OPT_GRID_N][OPT_GRID_N] = { { 0.0f } };
    float grid_mass_kg[OPT_GRID_N][OPT_GRID_N] = { { 0.0f } };

    // Interpolated Feasibility Limit Boundary Curves: Dr_min(L) for Strain Limit and MBR Limit
    int strain_bound_count = 0;
    float strain_bound_L[OPT_GRID_N] = { 0.0f };
    float strain_bound_Dr[OPT_GRID_N] = { 0.0f };

    int mbr_bound_count = 0;
    float mbr_bound_L[OPT_GRID_N] = { 0.0f };
    float mbr_bound_Dr[OPT_GRID_N] = { 0.0f };
};

// Fast stack-only evaluation of a candidate BendStiffner design for a single (angle_rad, tension_N) load case
inline void evaluate_candidate_single(const BendStiffner &cand,
                                      const Material &mat,
                                      size_t steps,
                                      float angle_rad,
                                      float tension_N,
                                      float &out_max_strain,
                                      float &out_avg_strain,
                                      float &out_tip_half_strain,
                                      float &out_root_half_strain,
                                      float &out_min_mbr,
                                      float &out_root_m0_kNm,
                                      float* out_bond_stress_mpa = nullptr) {
    Vec2 m0_v0 = solve_bvp_material(cand, steps, 0.0f, 0.0f, 0.0f, angle_rad, mat, tension_N);
    State4 state = { { 0.0f, 0.0f, m0_v0.x, m0_v0.y }, 0.0f };
    float h = cand.length / (float)steps;
    float s = 0.0f;

    out_max_strain = 0.0f;
    out_avg_strain = 0.0f;
    out_tip_half_strain = 0.0f;
    out_root_half_strain = 0.0f;
    float max_kappa = 0.0f;
    out_root_m0_kNm = std::fabs(m0_v0.x) / 1000.0f;
    if (out_bond_stress_mpa) {
        *out_bond_stress_mpa = get_insert_bond_stress_mpa(cand, m0_v0.x, m0_v0.y);
    }

    float abs_angle = std::fabs(angle_rad);
    size_t n_nodes = steps + 1;
    for (size_t i = 0; i < n_nodes; ++i) {
        float axial_N = std::fabs(tension_N) * std::max(0.0f, std::cos(angle_rad - state.items[1]));
        float kappa = 0.0f, eps = 0.0f, EI_tot = 1.0f;
        solve_composite_section(cand, mat, s, state.items[2], axial_N, kappa, eps, EI_tot);

        if (eps > out_max_strain) out_max_strain = eps;
        if (i < n_nodes / 2) {
            if (eps > out_root_half_strain) out_root_half_strain = eps;
        } else {
            if (eps > out_tip_half_strain) out_tip_half_strain = eps;
        }
        out_avg_strain += eps;

        float abs_k = get_local_pipe_curvature(cand, s, std::fabs(kappa), abs_angle);
        if (i + 1 == n_nodes && cand.use_far_field_bc && std::fabs(tension_N) > 1e-2f && cand.pipe_EI_kNm2 > 1e-4f) {
            float k_exit = std::fabs(state.items[2]) / (cand.pipe_EI_kNm2 * 1000.0f);
            abs_k = std::max(abs_k, k_exit);
        }
        if (abs_k > max_kappa) max_kappa = abs_k;

        if (i < steps) {
            state = RK4_fast(s, state, h, cand, mat, tension_N, angle_rad);
            s += h;
        }
    }
    out_avg_strain /= (float)n_nodes;
    out_min_mbr = (max_kappa > 1e-6f) ? (1.0f / max_kappa) : 1000.0f;
}

// Multi-stage global + local pattern search optimizer for Bend Stiffener sizing
inline OptimizationResult optimize_bend_stiffener(const BendStiffner &initial_bs,
                                                  const Material &mat,
                                                  const OptimizationConfig &cfg,
                                                  float single_angle_deg,
                                                  float single_tension_kN,
                                                  const Vec1<float>* pairs_angles_deg = nullptr,
                                                  const Vec1<float>* pairs_tensions_kN = nullptr) {
    OptimizationResult res;
    res.best_bs = initial_bs;

    // Build effective load case list (keep up to 6 governing envelope pairs for ultra-fast evaluation)
    struct LoadCase { float angle_rad; float tension_N; };
    LoadCase cases[8];
    int num_cases = 0;

    if (cfg.use_all_pairs && pairs_angles_deg && pairs_tensions_kN && pairs_angles_deg->count > 0) {
        size_t n_pairs = std::min(pairs_angles_deg->count, pairs_tensions_kN->count);
        if (n_pairs <= 6) {
            for (size_t i = 0; i < n_pairs; ++i) {
                cases[num_cases++] = {
                    pairs_angles_deg->items[i] * (3.14159265358979323846f / 180.0f),
                    pairs_tensions_kN->items[i] * 1000.0f
                };
            }
        } else {
            // Select governing envelope cases: max angle, max tension, max(angle*tension), max(angle/sqrt(T)), etc.
            size_t idx_max_a = 0, idx_max_t = 0, idx_max_at = 0, idx_min_t_high_a = 0;
            for (size_t i = 0; i < n_pairs; ++i) {
                float a = std::fabs(pairs_angles_deg->items[i]);
                float t = std::fabs(pairs_tensions_kN->items[i]);
                if (a > std::fabs(pairs_angles_deg->items[idx_max_a])) idx_max_a = i;
                if (t > std::fabs(pairs_tensions_kN->items[idx_max_t])) idx_max_t = i;
                if (a * t > std::fabs(pairs_angles_deg->items[idx_max_at] * pairs_tensions_kN->items[idx_max_at])) idx_max_at = i;
                if (a / std::sqrt(t + 1.0f) > std::fabs(pairs_angles_deg->items[idx_min_t_high_a]) / std::sqrt(std::fabs(pairs_tensions_kN->items[idx_min_t_high_a]) + 1.0f)) {
                    idx_min_t_high_a = i;
                }
            }
            size_t picked[4] = { idx_max_a, idx_max_t, idx_max_at, idx_min_t_high_a };
            for (int p = 0; p < 4; ++p) {
                bool dup = false;
                for (int q = 0; q < p; ++q) if (picked[p] == picked[q]) dup = true;
                if (!dup) {
                    cases[num_cases++] = {
                        pairs_angles_deg->items[picked[p]] * (3.14159265358979323846f / 180.0f),
                        pairs_tensions_kN->items[picked[p]] * 1000.0f
                    };
                }
            }
        }
    } else {
        cases[num_cases++] = {
            single_angle_deg * (3.14159265358979323846f / 180.0f),
            single_tension_kN * 1000.0f
        };
    }

    // Inner diameter and Tip diameter are fixed by the flexible pipe and tip wall thickness
    float fixed_Dt = std::max(initial_bs.inner_dia + 0.01f, initial_bs.tip_dia);
    float L_lo = cfg.optimize_length ? std::max(0.3f, cfg.min_length) : initial_bs.length;
    float L_hi = cfg.optimize_length ? std::max(L_lo + 0.1f, cfg.max_length) : initial_bs.length;
    float Dr_lo = cfg.optimize_root_dia ? std::max(fixed_Dt + 0.02f, cfg.min_root_dia) : initial_bs.root_dia;
    float Dr_hi = cfg.optimize_root_dia ? std::max(Dr_lo + 0.02f, cfg.max_root_dia) : initial_bs.root_dia;

    float allow_eps_eff = std::max(1e-4f, initial_bs.allowable_strain * std::clamp(cfg.safety_margin, 0.5f, 1.0f));
    float allow_mbr_eff = std::max(0.1f, initial_bs.allowable_mbr_m / std::clamp(cfg.safety_margin, 0.5f, 1.0f));
    float allow_bond_eff = std::max(0.05f, initial_bs.allowable_bond_stress_mpa * std::clamp(cfg.safety_margin, 0.5f, 1.0f));
    float pipe_EI_Pa = std::max(0.0f, initial_bs.pipe_EI_kNm2) * 1000.0f;

    // Pre-build temperature material instances ONCE outside eval_objective (zero heap allocations during optimization!)
    int n_temps = cfg.include_temp_envelope ? 3 : 1;
    Material pre_mats[3] = {
        cfg.include_temp_envelope ? mat.with_temp_mode(MAT_TEMP_NOMINAL) : mat,
        mat.with_temp_mode(MAT_TEMP_COLD),
        mat.with_temp_mode(MAT_TEMP_HOT)
    };
    double e0_root_pa = (double)get_material_E(0.01f, pre_mats[0]);

    auto eval_objective = [&](float L, float Dr, size_t steps,
                              bool &out_feas, float &out_max_eps, float &out_min_mbr,
                              float &out_max_m0, float &out_max_bond) -> double {
        res.evaluations++;
        Dr = std::max(Dr, fixed_Dt + 0.02f);

        BendStiffner cand = initial_bs;
        cand.length = L;
        cand.root_dia = Dr;
        cand.tip_dia = fixed_Dt; // Strictly preserve user's tip_dia
        cand.root_cyl_length = std::clamp(cfg.cyl_length_ratio * L, 0.0f, 0.5f * L);

        float worst_eps = 0.0f;
        float worst_mbr = 1000.0f;
        float worst_m0 = 0.0f;
        float worst_bond = 0.0f;
        float avg_eps_accum = 0.0f;

        for (int t_idx = 0; t_idx < n_temps; ++t_idx) {
            const Material &case_mat = pre_mats[t_idx];
            for (int c = 0; c < num_cases; ++c) {
                BendStiffner case_bs = cand;
                if (cfg.force_far_field_when_tension && std::fabs(cases[c].tension_N) > 10.0f) {
                    case_bs.use_far_field_bc = true;
                }
                float c_max_eps = 0.0f, c_avg_eps = 0.0f, c_tip_eps = 0.0f, c_root_eps = 0.0f, c_mbr = 1000.0f, c_m0 = 0.0f, c_bond = 0.0f;
                evaluate_candidate_single(case_bs, case_mat, steps, cases[c].angle_rad, cases[c].tension_N,
                                          c_max_eps, c_avg_eps, c_tip_eps, c_root_eps, c_mbr, c_m0, &c_bond);
                if (c_max_eps > worst_eps) worst_eps = c_max_eps;
                if (c_mbr < worst_mbr) worst_mbr = c_mbr;
                if (c_m0 > worst_m0) worst_m0 = c_m0;
                if (c_bond > worst_bond) worst_bond = c_bond;
                avg_eps_accum += c_avg_eps;
            }
        }
        avg_eps_accum /= (float)(num_cases * n_temps);

        out_max_eps = worst_eps;
        out_min_mbr = worst_mbr;
        out_max_m0 = worst_m0;
        out_max_bond = worst_bond;

        double u_eps = (double)worst_eps / (double)allow_eps_eff;
        double u_mbr = (double)allow_mbr_eff / std::max(1e-4, (double)worst_mbr);
        double u_bond = (cand.insert_length > 1e-4f) ? ((double)worst_bond / (double)allow_bond_eff) : 0.0;

        // Minimum root stiffness requirement relative to internal pipe when tension is zero
        double EI_bs_root = e0_root_pa * (double)get_Inertia(cand, 0.0f) + (double)get_insert_EI_Pa(cand, 0.0f);
        double u_stiff = 0.0;
        if (pipe_EI_Pa > 1.0f) {
            double req_root_EI = (double)cfg.min_root_stiffness_ratio * (double)pipe_EI_Pa;
            u_stiff = req_root_EI / std::max(1.0, EI_bs_root);
        }

        out_feas = (u_eps <= 1.0) && (u_mbr <= 1.0) && (u_stiff <= 1.0) && (u_bond <= 1.0);

        double vol = (double)get_pu_volume_m3(cand);
        double non_uniformity = (worst_eps > 1e-6f) ? ((double)(worst_eps - avg_eps_accum) / (double)worst_eps) : 0.0;
        double base_cost = vol * (1.0 + 0.20 * non_uniformity);

        if (out_feas) {
            return base_cost;
        }

        double pen = 0.0;
        if (u_eps > 1.0)   pen += (u_eps - 1.0) * (u_eps - 1.0) * 100.0;
        if (u_mbr > 1.0)   pen += (u_mbr - 1.0) * (u_mbr - 1.0) * 100.0;
        if (u_bond > 1.0)  pen += (u_bond - 1.0) * (u_bond - 1.0) * 50.0;
        if (u_stiff > 1.0) pen += (u_stiff - 1.0) * (u_stiff - 1.0) * 25.0;
        return 1000.0 + pen + 0.1 * base_cost;
    };

    // Stage 1: 2D Global Grid Search over (Length L, Root Diameter D_root)
    const size_t coarse_steps = 20;
    int nL  = cfg.optimize_length ? OPT_GRID_N : 1;
    int nDr = cfg.optimize_root_dia ? OPT_GRID_N : 1;
    res.grid_nL = nL;
    res.grid_nDr = nDr;

    float best_L = std::clamp(initial_bs.length, L_lo, L_hi);
    float best_Dr = std::clamp(initial_bs.root_dia, Dr_lo, Dr_hi);
    double best_cost = 1e30;

    for (int iDr = 0; iDr < nDr; ++iDr) {
        res.grid_Dr[iDr] = (nDr > 1) ? (Dr_lo + (Dr_hi - Dr_lo) * (float)iDr / (float)(nDr - 1)) : Dr_lo;
    }

    for (int iL = 0; iL < nL; ++iL) {
        float L = (nL > 1) ? (L_lo + (L_hi - L_lo) * (float)iL / (float)(nL - 1)) : L_lo;
        res.grid_L[iL] = L;
        for (int iDr = 0; iDr < nDr; ++iDr) {
            float Dr = res.grid_Dr[iDr];
            bool feas = false;
            float me = 0, mm = 0, m0 = 0, mb = 0;
            double cost = eval_objective(L, Dr, coarse_steps, feas, me, mm, m0, mb);

            BendStiffner g_bs = initial_bs;
            g_bs.length = L;
            g_bs.root_dia = std::max(Dr, fixed_Dt + 0.02f);
            g_bs.tip_dia = fixed_Dt;
            g_bs.root_cyl_length = std::clamp(cfg.cyl_length_ratio * L, 0.0f, 0.5f * L);

            res.grid_feasible[iL][iDr] = feas;
            res.grid_max_strain[iL][iDr] = me;
            res.grid_min_mbr[iL][iDr] = mm;
            res.grid_mass_kg[iL][iDr] = get_pu_mass_kg(g_bs);

            if (cost < best_cost) {
                best_cost = cost;
                best_L = L;
                best_Dr = Dr;
            }
        }
    }

    // Interpolate 1D Feasibility Boundary Curves Dr_min(L) for Strain Limit and MBR Limit
    if (nL > 1 && nDr > 1) {
        res.strain_bound_count = nL;
        res.mbr_bound_count = nL;
        for (int iL = 0; iL < nL; ++iL) {
            res.strain_bound_L[iL] = res.grid_L[iL];
            res.mbr_bound_L[iL] = res.grid_L[iL];

            // Strain Limit Boundary: find lowest Dr where max_strain <= allow_eps_eff
            float dr_strain = res.grid_Dr[nDr - 1];
            if (res.grid_max_strain[iL][0] <= allow_eps_eff) {
                dr_strain = res.grid_Dr[0];
            } else {
                for (int iDr = 0; iDr + 1 < nDr; ++iDr) {
                    float e0 = res.grid_max_strain[iL][iDr];
                    float e1 = res.grid_max_strain[iL][iDr + 1];
                    if ((e0 - allow_eps_eff) * (e1 - allow_eps_eff) <= 0.0f && std::fabs(e1 - e0) > 1e-6f) {
                        float t = std::clamp((allow_eps_eff - e0) / (e1 - e0), 0.0f, 1.0f);
                        dr_strain = res.grid_Dr[iDr] + t * (res.grid_Dr[iDr + 1] - res.grid_Dr[iDr]);
                        break;
                    }
                }
            }
            res.strain_bound_Dr[iL] = dr_strain;

            // MBR Limit Boundary: find lowest Dr where min_mbr >= allow_mbr_eff
            float dr_mbr = res.grid_Dr[nDr - 1];
            if (res.grid_min_mbr[iL][0] >= allow_mbr_eff) {
                dr_mbr = res.grid_Dr[0];
            } else {
                for (int iDr = 0; iDr + 1 < nDr; ++iDr) {
                    float m0 = res.grid_min_mbr[iL][iDr];
                    float m1 = res.grid_min_mbr[iL][iDr + 1];
                    if ((m0 - allow_mbr_eff) * (m1 - allow_mbr_eff) <= 0.0f && std::fabs(m1 - m0) > 1e-6f) {
                        float t = std::clamp((allow_mbr_eff - m0) / (m1 - m0), 0.0f, 1.0f);
                        dr_mbr = res.grid_Dr[iDr] + t * (res.grid_Dr[iDr + 1] - res.grid_Dr[iDr]);
                        break;
                    }
                }
            }
            res.mbr_bound_Dr[iL] = dr_mbr;
        }
    }

    // Stage 2: Multi-scale 2D coordinate pattern search refinement at full resolution
    const size_t fine_steps = (DISCRETIZATION > 1) ? (DISCRETIZATION - 1) : 24;
    bool feas_now = false;
    float eps_now = 0, mbr_now = 0, m0_now = 0, bond_now = 0;
    best_cost = eval_objective(best_L, best_Dr, fine_steps, feas_now, eps_now, mbr_now, m0_now, bond_now);

    float step_scales[4] = { 0.10f, 0.035f, 0.010f, 0.002f };
    for (int sc = 0; sc < 4; ++sc) {
        float dL  = cfg.optimize_length ? ((L_hi - L_lo) * step_scales[sc]) : 0.0f;
        float dDr = cfg.optimize_root_dia ? ((Dr_hi - Dr_lo) * step_scales[sc]) : 0.0f;

        for (int pass = 0; pass < 3; ++pass) {
            bool improved = false;
            if (cfg.optimize_length && dL > 1e-4f) {
                for (int dir = -1; dir <= 1; dir += 2) {
                    float cand_L = std::clamp(best_L + (float)dir * dL, L_lo, L_hi);
                    bool f = false; float me = 0, mm = 0, m0 = 0, mb = 0;
                    double c = eval_objective(cand_L, best_Dr, fine_steps, f, me, mm, m0, mb);
                    if (c < best_cost) {
                        best_cost = c; best_L = cand_L; improved = true;
                    }
                }
            }
            if (cfg.optimize_root_dia && dDr > 1e-4f) {
                for (int dir = -1; dir <= 1; dir += 2) {
                    float cand_Dr = std::clamp(best_Dr + (float)dir * dDr, Dr_lo, Dr_hi);
                    bool f = false; float me = 0, mm = 0, m0 = 0, mb = 0;
                    double c = eval_objective(best_L, cand_Dr, fine_steps, f, me, mm, m0, mb);
                    if (c < best_cost) {
                        best_cost = c; best_Dr = cand_Dr; improved = true;
                    }
                }
            }
            if (!improved) break;
        }
    }

    // Final evaluation of optimal design
    eval_objective(best_L, best_Dr, fine_steps, res.feasible, res.max_strain, res.min_mbr_m, res.max_root_moment_kNm, res.max_bond_stress_mpa);
    res.converged = true;
    res.best_bs = initial_bs;
    res.best_bs.length = std::round(best_L * 1000.0f) / 1000.0f;
    res.best_bs.root_dia = std::round(best_Dr * 1000.0f) / 1000.0f;
    res.best_bs.tip_dia = initial_bs.tip_dia; // Unchanged
    res.best_bs.inner_dia = initial_bs.inner_dia; // Unchanged
    res.best_bs.root_cyl_length = std::round(std::clamp(cfg.cyl_length_ratio * res.best_bs.length, 0.0f, 0.5f * res.best_bs.length) * 1000.0f) / 1000.0f;
    res.volume_m3 = get_pu_volume_m3(res.best_bs);
    res.mass_kg = get_pu_mass_kg(res.best_bs);
    res.strain_utilization = (initial_bs.allowable_strain > 1e-5f) ? (res.max_strain / initial_bs.allowable_strain) : 0.0f;
    res.mbr_utilization = (res.min_mbr_m > 1e-5f) ? (initial_bs.allowable_mbr_m / res.min_mbr_m) : 0.0f;
    res.bond_utilization = (initial_bs.allowable_bond_stress_mpa > 1e-5f) ? (res.max_bond_stress_mpa / initial_bs.allowable_bond_stress_mpa) : 0.0f;
    return res;
}

#define TAIL_NODES 25

// Integrates the bare flexible pipe tail beyond the bend stiffener tip (s in [L, L + tail_length])
inline void simulate_bare_pipe_tail(const BendStiffner &bs,
                                    float tip_x,
                                    float tip_y,
                                    float tip_theta,
                                    float tip_moment_Nm,
                                    float tension_N,
                                    float target_angle_rad,
                                    float tail_length,
                                    float out_s[TAIL_NODES],
                                    float out_x[TAIL_NODES],
                                    float out_y[TAIL_NODES],
                                    float out_theta[TAIL_NODES],
                                    float out_moment_kNm[TAIL_NODES],
                                    float out_curvature[TAIL_NODES],
                                    float out_mbr[TAIL_NODES]) {
    float pipe_EI_Pa = std::max(1.0f, bs.pipe_EI_kNm2 * 1000.0f);
    float T_abs = std::fabs(tension_N);
    bool far_field = bs.use_far_field_bc && (T_abs > 1e-2f);
    float h = tail_length / (float)(TAIL_NODES - 1);

    float s = bs.length;
    float x = tip_x;
    float y = tip_y;
    float theta = tip_theta;
    float m = far_field ? tip_moment_Nm : 0.0f;

    auto deriv = [&](float th, float mom, float &dx, float &dy, float &dth, float &dm) {
        dx = std::cos(th);
        dy = std::sin(th);
        dth = far_field ? (mom / pipe_EI_Pa) : 0.0f;
        dm = far_field ? (-T_abs * std::sin(target_angle_rad - th)) : 0.0f;
    };

    for (int i = 0; i < TAIL_NODES; ++i) {
        out_s[i] = s;
        out_x[i] = x;
        out_y[i] = y;
        out_theta[i] = theta;
        out_moment_kNm[i] = m / 1000.0f;
        float k_abs = far_field ? (std::fabs(m) / pipe_EI_Pa) : 0.0f;
        out_curvature[i] = k_abs;
        out_mbr[i] = (k_abs > 1e-6f) ? std::min(1000.0f, 1.0f / k_abs) : 1000.0f;

        if (i + 1 < TAIL_NODES) {
            float k1x, k1y, k1th, k1m;
            deriv(theta, m, k1x, k1y, k1th, k1m);
            float k2x, k2y, k2th, k2m;
            deriv(theta + 0.5f * h * k1th, m + 0.5f * h * k1m, k2x, k2y, k2th, k2m);
            float k3x, k3y, k3th, k3m;
            deriv(theta + 0.5f * h * k2th, m + 0.5f * h * k2m, k3x, k3y, k3th, k3m);
            float k4x, k4y, k4th, k4m;
            deriv(theta + h * k3th, m + h * k3m, k4x, k4y, k4th, k4m);

            x += (h / 6.0f) * (k1x + 2.0f * k2x + 2.0f * k3x + k4x);
            y += (h / 6.0f) * (k1y + 2.0f * k2y + 2.0f * k3y + k4y);
            theta += (h / 6.0f) * (k1th + 2.0f * k2th + 2.0f * k3th + k4th);
            m += (h / 6.0f) * (k1m + 2.0f * k2m + 2.0f * k3m + k4m);
            s += h;
        }
    }
}

struct TemperatureCheckRow {
    const char* label;
    float temp_c;
    float initial_E_mpa;
    float max_strain;
    float min_mbr_m;
    float root_moment_kNm;
    float root_shear_kN;
    float root_axial_kN;
    bool pass;
};

struct TemperatureEnvelopeSummary {
    bool evaluated = false;
    bool all_pass = false;
    TemperatureCheckRow rows[3]; // 0: Cold, 1: Nominal, 2: Hot
    float governing_strain = 0.0f;
    float governing_min_mbr_m = 1000.0f;
    float governing_root_moment_kNm = 0.0f;
    float governing_root_shear_kN = 0.0f;
    float governing_root_axial_kN = 0.0f;
};

inline TemperatureEnvelopeSummary evaluate_temperature_envelope(const BendStiffner &bs,
                                                                const Material &base_mat,
                                                                float angle_deg,
                                                                float tension_kN) {
    TemperatureEnvelopeSummary env;
    if (!is_valid_geometry(bs)) return env;

    static const char* labels[3] = {
        "Cold Temperature Curve",
        "Nominal Temperature Curve",
        "Hot / Warm Temperature Curve"
    };
    static const int modes[3] = { MAT_TEMP_COLD, MAT_TEMP_NOMINAL, MAT_TEMP_HOT };
    const float temps_c[3] = { base_mat.temp_cold_c, base_mat.temp_nominal_c, base_mat.temp_hot_c };

    float angle_rad = angle_deg * (3.14159265358979323846f / 180.0f);
    float tension_N = tension_kN * 1000.0f;
    const size_t steps = (DISCRETIZATION > 1) ? (DISCRETIZATION - 1) : 24;

    env.evaluated = true;
    env.all_pass = true;
    env.governing_strain = 0.0f;
    env.governing_min_mbr_m = 1000.0f;
    env.governing_root_moment_kNm = 0.0f;
    env.governing_root_shear_kN = 0.0f;
    env.governing_root_axial_kN = 0.0f;

    for (int i = 0; i < 3; ++i) {
        Material m_temp = base_mat.with_temp_mode(modes[i]);
        Vec2 m0_v0 = solve_bvp_material(bs, steps, 0.0f, 0.0f, 0.0f, angle_rad, m_temp, tension_N);

        float max_eps = 0.0f, avg_eps = 0.0f, tip_eps = 0.0f, root_eps = 0.0f, min_mbr = 1000.0f, m0_kNm = 0.0f;
        evaluate_candidate_single(bs, m_temp, steps, angle_rad, tension_N,
                                  max_eps, avg_eps, tip_eps, root_eps, min_mbr, m0_kNm);

        float v0_kN = std::fabs(m0_v0.y) / 1000.0f;
        float n0_kN = std::fabs(tension_kN) * std::max(0.0f, std::cos(angle_rad));
        bool row_pass = (max_eps <= bs.allowable_strain) && (min_mbr >= bs.allowable_mbr_m);
        float e0_mpa = get_material_E(0.0f, m_temp) / 1e6f;

        env.rows[i] = {
            labels[i],
            temps_c[i],
            e0_mpa,
            max_eps,
            min_mbr,
            m0_kNm,
            v0_kN,
            n0_kN,
            row_pass
        };

        if (!row_pass) env.all_pass = false;
        if (max_eps > env.governing_strain) env.governing_strain = max_eps;
        if (min_mbr < env.governing_min_mbr_m) env.governing_min_mbr_m = min_mbr;
        if (m0_kNm > env.governing_root_moment_kNm) env.governing_root_moment_kNm = m0_kNm;
        if (v0_kN > env.governing_root_shear_kN) env.governing_root_shear_kN = v0_kN;
        if (n0_kN > env.governing_root_axial_kN) env.governing_root_axial_kN = n0_kN;
    }
    return env;
}

// Save and Load project state (.bsproj text format)
inline bool save_project_file(const char* filename,
                              const BendStiffner &bs,
                              const Material &mat,
                              const OptimizationConfig &opt_cfg,
                              float target_angle_deg,
                              float target_tension_kN,
                              const Vec1<float> &pairs_angles_deg,
                              const Vec1<float> &pairs_tensions_kN) {
    FILE* f = std::fopen(filename, "w");
    if (!f) return false;

    std::fprintf(f, "BS_PROJECT_V2\n");
    std::fprintf(f, "geometry %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %d %d\n",
                 bs.length, bs.root_dia, bs.tip_dia, bs.inner_dia,
                 bs.root_cyl_length, bs.pipe_EI_kNm2, bs.allowable_strain, bs.allowable_mbr_m,
                 (int)bs.use_arc_length, (int)bs.use_far_field_bc);
    std::fprintf(f, "geometry_ext %.6f %.6f %.6f %.6f %.6f\n",
                 bs.pipe_od, bs.insert_length, bs.insert_dia, bs.insert_EI_kNm2, bs.allowable_bond_stress_mpa);
    std::fprintf(f, "loading %.6f %.6f\n", target_angle_deg, target_tension_kN);
    std::fprintf(f, "optimizer %.6f %.6f %.6f %.6f %.6f %.6f %.6f %d %d %d %d %d\n",
                 opt_cfg.min_length, opt_cfg.max_length, opt_cfg.min_root_dia, opt_cfg.max_root_dia,
                 opt_cfg.cyl_length_ratio, opt_cfg.safety_margin, opt_cfg.min_root_stiffness_ratio,
                 (int)opt_cfg.optimize_length, (int)opt_cfg.optimize_root_dia,
                 (int)opt_cfg.use_all_pairs, (int)opt_cfg.force_far_field_when_tension,
                 (int)opt_cfg.include_temp_envelope);

    size_t n_mat = mat.strain_table.count;
    std::fprintf(f, "material_v2 %d %d %d %d %.3f %.3f %.3f %.3f %.3f %.3f %zu\n",
                 (int)mat.is_nonlinear, (int)mat.table_is_tangent_modulus,
                 (int)mat.use_cross_section_integration, mat.active_temp_mode,
                 mat.temp_nominal_c, mat.temp_cold_c, mat.temp_hot_c,
                 mat.linear_E_nominal_kpa, mat.linear_E_cold_kpa, mat.linear_E_hot_kpa, n_mat);
    for (size_t i = 0; i < n_mat; ++i) {
        float en = (i < mat.emod_nominal_kpa.count) ? mat.emod_nominal_kpa.items[i] : 215800.0f;
        float ec = (i < mat.emod_cold_kpa.count) ? mat.emod_cold_kpa.items[i] : en;
        float eh = (i < mat.emod_hot_kpa.count) ? mat.emod_hot_kpa.items[i] : en;
        std::fprintf(f, "mat_row_v2 %.6f %.6f %.6f %.6f\n", mat.strain_table.items[i], en, ec, eh);
    }

    size_t n_pairs = std::min(pairs_angles_deg.count, pairs_tensions_kN.count);
    std::fprintf(f, "pairs %zu\n", n_pairs);
    for (size_t i = 0; i < n_pairs; ++i) {
        std::fprintf(f, "pair_row %.6f %.6f\n", pairs_angles_deg.items[i], pairs_tensions_kN.items[i]);
    }
    std::fclose(f);
    return true;
}

inline bool load_project_file(const char* filename,
                              BendStiffner &bs,
                              Material &mat,
                              OptimizationConfig &opt_cfg,
                              float &target_angle_deg,
                              float &target_tension_kN,
                              Vec1<float> &pairs_angles_deg,
                              Vec1<float> &pairs_tensions_kN) {
    FILE* f = std::fopen(filename, "r");
    if (!f) return false;

    char header[64] = { 0 };
    if (std::fscanf(f, "%63s", header) != 1 ||
        (std::strcmp(header, "BS_PROJECT_V1") != 0 && std::strcmp(header, "BS_PROJECT_V2") != 0)) {
        std::fclose(f);
        return false;
    }

    char tag[64] = { 0 };
    while (std::fscanf(f, "%63s", tag) == 1) {
        if (std::strcmp(tag, "geometry") == 0) {
            int arc = 1, ff = 0;
            if (std::fscanf(f, "%f %f %f %f %f %f %f %f %d %d",
                            &bs.length, &bs.root_dia, &bs.tip_dia, &bs.inner_dia,
                            &bs.root_cyl_length, &bs.pipe_EI_kNm2, &bs.allowable_strain, &bs.allowable_mbr_m,
                            &arc, &ff) == 10) {
                bs.use_arc_length = (arc != 0);
                bs.use_far_field_bc = (ff != 0);
            }
        } else if (std::strcmp(tag, "geometry_ext") == 0) {
            std::fscanf(f, "%f %f %f %f %f",
                        &bs.pipe_od, &bs.insert_length, &bs.insert_dia, &bs.insert_EI_kNm2, &bs.allowable_bond_stress_mpa);
        } else if (std::strcmp(tag, "loading") == 0) {
            std::fscanf(f, "%f %f", &target_angle_deg, &target_tension_kN);
        } else if (std::strcmp(tag, "optimizer") == 0) {
            int ol = 1, odr = 1, uap = 0, fff = 1, ite = 0;
            if (std::fscanf(f, "%f %f %f %f %f %f %f %d %d %d %d %d",
                            &opt_cfg.min_length, &opt_cfg.max_length, &opt_cfg.min_root_dia, &opt_cfg.max_root_dia,
                            &opt_cfg.cyl_length_ratio, &opt_cfg.safety_margin, &opt_cfg.min_root_stiffness_ratio,
                            &ol, &odr, &uap, &fff, &ite) == 12) {
                opt_cfg.optimize_length = (ol != 0);
                opt_cfg.optimize_root_dia = (odr != 0);
                opt_cfg.use_all_pairs = (uap != 0);
                opt_cfg.force_far_field_when_tension = (fff != 0);
                opt_cfg.include_temp_envelope = (ite != 0);
            }
        } else if (std::strcmp(tag, "material_v2") == 0) {
            int nl = 1, tan_m = 1, csi = 0, act_mode = 0;
            size_t n_mat = 0;
            if (std::fscanf(f, "%d %d %d %d %f %f %f %f %f %f %zu",
                            &nl, &tan_m, &csi, &act_mode,
                            &mat.temp_nominal_c, &mat.temp_cold_c, &mat.temp_hot_c,
                            &mat.linear_E_nominal_kpa, &mat.linear_E_cold_kpa, &mat.linear_E_hot_kpa,
                            &n_mat) == 11) {
                mat.is_nonlinear = (nl != 0);
                mat.table_is_tangent_modulus = (tan_m != 0);
                mat.use_cross_section_integration = (csi != 0);
                mat.strain_table.clear();
                mat.emod_nominal_kpa.clear();
                mat.emod_cold_kpa.clear();
                mat.emod_hot_kpa.clear();
                for (size_t i = 0; i < n_mat; ++i) {
                    char rtag[32] = { 0 };
                    float s_v = 0.0f, en = 0.0f, ec = 0.0f, eh = 0.0f;
                    if (std::fscanf(f, "%31s %f %f %f %f", rtag, &s_v, &en, &ec, &eh) == 5) {
                        mat.strain_table.append(s_v);
                        mat.emod_nominal_kpa.append(en);
                        mat.emod_cold_kpa.append(ec);
                        mat.emod_hot_kpa.append(eh);
                    }
                }
                mat.set_active_temp_mode(act_mode);
            }
        } else if (std::strcmp(tag, "material") == 0) {
            int nl = 1, tan_m = 1, csi = 0;
            size_t n_mat = 0;
            if (std::fscanf(f, "%d %d %d %f %zu", &nl, &tan_m, &csi, &mat.linear_E_nominal_kpa, &n_mat) == 5) {
                mat.is_nonlinear = (nl != 0);
                mat.table_is_tangent_modulus = (tan_m != 0);
                mat.use_cross_section_integration = (csi != 0);
                mat.strain_table.clear();
                mat.emod_nominal_kpa.clear();
                mat.emod_cold_kpa.clear();
                mat.emod_hot_kpa.clear();
                for (size_t i = 0; i < n_mat; ++i) {
                    char rtag[32] = { 0 };
                    float s_v = 0.0f, e_v = 0.0f;
                    if (std::fscanf(f, "%31s %f %f", rtag, &s_v, &e_v) == 3) {
                        mat.strain_table.append(s_v);
                        mat.emod_nominal_kpa.append(e_v);
                        mat.emod_cold_kpa.append(e_v * 1.35f);
                        mat.emod_hot_kpa.append(e_v * 0.75f);
                    }
                }
                mat.set_active_temp_mode(MAT_TEMP_NOMINAL);
            }
        } else if (std::strcmp(tag, "pairs") == 0) {
            size_t n_pairs = 0;
            if (std::fscanf(f, "%zu", &n_pairs) == 1) {
                pairs_angles_deg.clear();
                pairs_tensions_kN.clear();
                for (size_t i = 0; i < n_pairs; ++i) {
                    char rtag[32] = { 0 };
                    float a_v = 0.0f, t_v = 0.0f;
                    if (std::fscanf(f, "%31s %f %f", rtag, &a_v, &t_v) == 3) {
                        pairs_angles_deg.append(a_v);
                        pairs_tensions_kN.append(t_v);
                    }
                }
            }
        }
    }
    std::fclose(f);
    return true;
}

// Generates a self-contained HTML5 Engineering Calculation Report with embedded Vector SVG charts
inline bool generate_html_report(const char* filename,
                                 const BendStiffner &bs,
                                 const Material &mat,
                                 float target_angle_deg,
                                 float target_tension_kN,
                                 float root_moment_kNm,
                                 float root_shear_kN,
                                 float root_axial_kN,
                                 const float* x_stations,
                                 const float* def_x,
                                 const float* def_y,
                                 const float* angle_rad,
                                 const float* moment_Nm,
                                 const float* shear_N,
                                 const float* strain_arr,
                                 const float* curvature_arr,
                                 const float* mbr_arr,
                                 size_t n_nodes,
                                 const float* tail_x = nullptr,
                                 const float* tail_y = nullptr,
                                 size_t n_tail = 0) {
    if (!filename || n_nodes < 2) return false;
    FILE* f = std::fopen(filename, "w");
    if (!f) return false;

    float max_strain = 0.0f;
    float max_strain_s = 0.0f;
    float min_mbr = 1000.0f;
    float min_mbr_s = 0.0f;
    for (size_t i = 0; i < n_nodes; ++i) {
        if (strain_arr[i] > max_strain) {
            max_strain = strain_arr[i];
            max_strain_s = x_stations[i];
        }
        if (mbr_arr[i] < min_mbr) {
            min_mbr = mbr_arr[i];
            min_mbr_s = x_stations[i];
        }
    }

    float pu_vol_m3 = get_pu_volume_m3(bs);
    float pu_mass_kg = get_pu_mass_kg(bs);
    float bond_stress_mpa = get_insert_bond_stress_mpa(bs, root_moment_kNm * 1000.0f, root_shear_kN * 1000.0f);
    float gap_mm = get_annular_gap_mm(bs);
    float gap_angle_deg = get_gap_contact_angle_deg(bs);
    float resultant_kN = std::sqrt(root_shear_kN * root_shear_kN + root_axial_kN * root_axial_kN);

    float u_strain = (bs.allowable_strain > 1e-5f) ? (max_strain / bs.allowable_strain) : 0.0f;
    float u_mbr = (min_mbr > 1e-5f) ? (bs.allowable_mbr_m / min_mbr) : 0.0f;
    float u_bond = (bs.insert_length > 1e-4f && bs.allowable_bond_stress_mpa > 1e-5f)
                   ? (bond_stress_mpa / bs.allowable_bond_stress_mpa) : 0.0f;

    bool pass_strain = (max_strain <= bs.allowable_strain);
    bool pass_mbr = (min_mbr >= bs.allowable_mbr_m);
    bool pass_bond = (bs.insert_length <= 1e-4f) || (bond_stress_mpa <= bs.allowable_bond_stress_mpa);
    bool overall_pass = pass_strain && pass_mbr && pass_bond;

    TemperatureEnvelopeSummary temp_env = evaluate_temperature_envelope(bs, mat, target_angle_deg, target_tension_kN);

    std::fprintf(f, "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"UTF-8\">\n");
    std::fprintf(f, "<title>Bend Stiffener Engineering Calculation Report</title>\n");
    std::fprintf(f, "<style>\n"
                    ":root { --bg: #0f172a; --card: #1e293b; --border: #334155; --text: #f1f5f9; --muted: #94a3b8; "
                    "--accent: #38bdf8; --pass: #22c55e; --fail: #ef4444; --warn: #f59e0b; }\n"
                    "@media print { :root { --bg: #ffffff; --card: #f8fafc; --border: #cbd5e1; --text: #0f172a; --muted: #475569; } }\n"
                    "* { box-sizing: border-box; }\n"
                    "body { font-family: 'Segoe UI', system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--text); "
                    "margin: 0; padding: 28px 36px; line-height: 1.5; }\n"
                    ".container { max-width: 1240px; margin: 0 auto; }\n"
                    ".header { display: flex; justify-content: space-between; align-items: center; border-bottom: 2px solid var(--border); padding-bottom: 16px; margin-bottom: 24px; }\n"
                    "h1 { margin: 0; font-size: 24px; color: var(--accent); letter-spacing: 0.3px; }\n"
                    "h2 { font-size: 17px; color: var(--accent); border-left: 4px solid var(--accent); padding-left: 10px; margin: 26px 0 12px 0; }\n"
                    ".badge { padding: 6px 14px; border-radius: 6px; font-weight: 700; font-size: 14px; display: inline-block; }\n"
                    ".badge-pass { background: rgba(34, 197, 94, 0.18); color: var(--pass); border: 1px solid var(--pass); }\n"
                    ".badge-fail { background: rgba(239, 68, 68, 0.18); color: var(--fail); border: 1px solid var(--fail); }\n"
                    ".grid-4 { display: grid; grid-template-columns: repeat(4, 1fr); gap: 14px; margin-bottom: 20px; }\n"
                    ".grid-2 { display: grid; grid-template-columns: repeat(2, 1fr); gap: 18px; margin-bottom: 20px; }\n"
                    ".card { background: var(--card); border: 1px solid var(--border); border-radius: 8px; padding: 16px; }\n"
                    ".kpi-label { font-size: 12px; color: var(--muted); text-transform: uppercase; letter-spacing: 0.6px; }\n"
                    ".kpi-val { font-size: 22px; font-weight: 700; margin: 4px 0; }\n"
                    ".kpi-sub { font-size: 12px; color: var(--muted); }\n"
                    "table { width: 100%%; border-collapse: collapse; font-size: 13px; background: var(--card); border-radius: 8px; overflow: hidden; }\n"
                    "th, td { padding: 8px 12px; text-align: left; border-bottom: 1px solid var(--border); }\n"
                    "th { background: rgba(56, 189, 248, 0.1); color: var(--accent); font-weight: 600; }\n"
                    "svg { width: 100%%; height: auto; background: #0b1120; border: 1px solid var(--border); border-radius: 6px; }\n"
                    "</style>\n</head>\n<body>\n<div class=\"container\">\n");

    std::fprintf(f, "<div class=\"header\">\n"
                    "<div><h1>Bend Stiffener Nonlinear Elastica Engineering Report</h1>\n"
                    "<div style=\"color:var(--muted);font-size:13px;\">BS Engine &bull; Coupled Composite Polyurethane + Pipe + Root Insert Solver &bull; %s</div></div>\n"
                    "<div class=\"badge %s\">%s</div>\n</div>\n",
                 bs.use_far_field_bc ? "OrcaFlex Far-Field Tension Elastica BC" : "Prescribed Tip Angle BC",
                 overall_pass ? "badge-pass" : "badge-fail",
                 overall_pass ? "OVERALL DESIGN: PASS" : "OVERALL DESIGN: LIMIT EXCEEDED");

    // KPI Cards
    std::fprintf(f, "<div class=\"grid-4\">\n");
    std::fprintf(f, "<div class=\"card\"><div class=\"kpi-label\">Peak Polyurethane Strain</div>"
                    "<div class=\"kpi-val\" style=\"color:%s;\">%.2f%%</div>"
                    "<div class=\"kpi-sub\">Allowable: %.2f%% (Util: %.1f%% @ s=%.2fm)</div></div>\n",
                 pass_strain ? "var(--pass)" : "var(--fail)",
                 max_strain * 100.0f, bs.allowable_strain * 100.0f, u_strain * 100.0f, max_strain_s);

    std::fprintf(f, "<div class=\"card\"><div class=\"kpi-label\">Minimum Pipe Bend Radius</div>"
                    "<div class=\"kpi-val\" style=\"color:%s;\">%.2f m</div>"
                    "<div class=\"kpi-sub\">Allowable MBR: %.2f m (Util: %.1f%% @ s=%.2fm)</div></div>\n",
                 pass_mbr ? "var(--pass)" : "var(--fail)",
                 min_mbr, bs.allowable_mbr_m, u_mbr * 100.0f, min_mbr_s);

    std::fprintf(f, "<div class=\"card\"><div class=\"kpi-label\">Hang-Off Flange Reactions</div>"
                    "<div class=\"kpi-val\">M<sub>0</sub> = %.2f kN&middot;m</div>"
                    "<div class=\"kpi-sub\">V<sub>0</sub> = %.2f kN &bull; N<sub>0</sub> = %.2f kN (|R| = %.2f kN)</div></div>\n",
                 std::fabs(root_moment_kNm), std::fabs(root_shear_kN), std::fabs(root_axial_kN), resultant_kN);

    std::fprintf(f, "<div class=\"card\"><div class=\"kpi-label\">PU Volume, Mass &amp; Root Bond</div>"
                    "<div class=\"kpi-val\">%.1f kg <span style=\"font-size:14px;color:var(--muted);\">(%.3f m&sup3;)</span></div>"
                    "<div class=\"kpi-sub\">Insert Bond: %.2f / %.2f MPa &bull; Gap: %.1f mm</div></div>\n",
                 pu_mass_kg, pu_vol_m3, bond_stress_mpa, bs.allowable_bond_stress_mpa, gap_mm);
    std::fprintf(f, "</div>\n");

    // Input Summary Tables (2 columns)
    std::fprintf(f, "<div class=\"grid-2\">\n");
    std::fprintf(f, "<div class=\"card\"><h2 style=\"margin-top:0;\">1. Geometry &amp; Mechanical Parameters</h2><table>\n"
                    "<tr><th>Parameter</th><th>Symbol</th><th>Value</th><th>Units</th></tr>\n"
                    "<tr><td>Total Stiffener Length</td><td>L</td><td>%.3f</td><td>m</td></tr>\n"
                    "<tr><td>Cylindrical Base Length</td><td>L<sub>cyl</sub></td><td>%.3f</td><td>m</td></tr>\n"
                    "<tr><td>Outer Root Diameter</td><td>D<sub>root</sub></td><td>%.3f</td><td>m</td></tr>\n"
                    "<tr><td>Outer Tip Diameter</td><td>D<sub>tip</sub></td><td>%.3f</td><td>m</td></tr>\n"
                    "<tr><td>Inner Bore Diameter</td><td>D<sub>i</sub></td><td>%.3f</td><td>m</td></tr>\n"
                    "<tr><td>Flexible Pipe Outer Diameter</td><td>D<sub>pipe</sub></td><td>%.3f (Gap: %.1f mm, &Delta;&theta;<sub>gap</sub>: %.2f&deg;)</td><td>m</td></tr>\n"
                    "<tr><td>Embedded Steel Insert Sleeve</td><td>L<sub>ins</sub> &times; D<sub>ins</sub></td><td>%.3f m &times; %.3f m (EI<sub>ins</sub> = %.0f kN&middot;m&sup2;)</td><td>-</td></tr>\n"
                    "<tr><td>Internal Pipe Bending Stiffness</td><td>EI<sub>pipe</sub></td><td>%.2f</td><td>kN&middot;m&sup2;</td></tr>\n"
                    "</table></div>\n",
                 bs.length, bs.root_cyl_length, bs.root_dia, bs.tip_dia, bs.inner_dia,
                 bs.pipe_od, gap_mm, gap_angle_deg, bs.insert_length, bs.insert_dia, bs.insert_EI_kNm2, bs.pipe_EI_kNm2);

    const char* active_temp_label = (mat.active_temp_mode == MAT_TEMP_COLD) ? "Cold" :
                                    (mat.active_temp_mode == MAT_TEMP_HOT) ? "Hot" : "Nominal";
    float active_temp_c = (mat.active_temp_mode == MAT_TEMP_COLD) ? mat.temp_cold_c :
                          (mat.active_temp_mode == MAT_TEMP_HOT) ? mat.temp_hot_c : mat.temp_nominal_c;

    std::fprintf(f, "<div class=\"card\"><h2 style=\"margin-top:0;\">2. Loading, Boundary Conditions &amp; Criteria</h2><table>\n"
                    "<tr><th>Parameter</th><th>Value</th><th>Verification Status</th></tr>\n"
                    "<tr><td>Applied Bend Angle</td><td>%.2f deg</td><td>%s</td></tr>\n"
                    "<tr><td>Applied Far-Field Tension</td><td>%.2f kN</td><td>%s</td></tr>\n"
                    "<tr><td>Active Material Temperature Curve</td><td>%s (%.1f &deg;C)</td><td>Initial E<sub>0</sub> = %.1f MPa</td></tr>\n"
                    "<tr><td>Cross-Sectional Stress Integration</td><td>%s</td><td>%s</td></tr>\n"
                    "<tr><td>Polyurethane Outer-Fiber Strain Check</td><td>%.2f%% &le; %.2f%%</td><td><b style=\"color:%s;\">%s (%.1f%%)</b></td></tr>\n"
                    "<tr><td>Flexible Pipe Minimum Bend Radius Check</td><td>%.2f m &ge; %.2f m</td><td><b style=\"color:%s;\">%s (%.1f%%)</b></td></tr>\n"
                    "<tr><td>Steel Insert Interface Bond Shear Check</td><td>%.3f MPa &le; %.2f MPa</td><td><b style=\"color:%s;\">%s (%.1f%%)</b></td></tr>\n"
                    "</table></div>\n",
                 target_angle_deg, bs.use_arc_length ? "Exact Arc-Length Elastica (s)" : "Cartesian Projection",
                 target_tension_kN, bs.use_far_field_bc ? "OrcaFlex Far-Field Elastica BC" : "Prescribed Tip Angle BC",
                 active_temp_label, active_temp_c, get_material_E(0.0f, mat) / 1e6f,
                 mat.use_cross_section_integration ? "4-Point Gauss Annular Integration" : "Outer-Fiber Direct Evaluation",
                 mat.table_is_tangent_modulus ? "Tangent-to-Secant Integrated" : "Direct Secant Table",
                 max_strain * 100.0f, bs.allowable_strain * 100.0f, pass_strain ? "var(--pass)" : "var(--fail)", pass_strain ? "PASS" : "FAIL", u_strain * 100.0f,
                 min_mbr, bs.allowable_mbr_m, pass_mbr ? "var(--pass)" : "var(--fail)", pass_mbr ? "PASS" : "FAIL", u_mbr * 100.0f,
                 bond_stress_mpa, bs.allowable_bond_stress_mpa, pass_bond ? "var(--pass)" : "var(--fail)", pass_bond ? "PASS" : "FAIL", u_bond * 100.0f);
    std::fprintf(f, "</div>\n");

    // Vector SVG Charts (2x2 Grid)
    std::fprintf(f, "<h2>3. Engineering Verification Charts (Vector SVG)</h2>\n<div class=\"grid-2\">\n");

    // SVG 1: 2D Geometry, Steel Insert & Deformed Elastica + Pipe Tail
    {
        float x_max = bs.length;
        float y_min = -0.6f * bs.root_dia;
        float y_max = 0.6f * bs.root_dia;
        for (size_t i = 0; i < n_nodes; ++i) {
            x_max = std::max(x_max, def_x[i]);
            y_min = std::min(y_min, def_y[i]);
            y_max = std::max(y_max, def_y[i]);
        }
        if (tail_x && tail_y && n_tail > 0) {
            for (size_t i = 0; i < n_tail; ++i) {
                x_max = std::max(x_max, tail_x[i]);
                y_min = std::min(y_min, tail_y[i]);
                y_max = std::max(y_max, tail_y[i]);
            }
        }
        float pad_x = 0.12f * std::max(1.0f, x_max);
        float pad_y = 0.18f * std::max(0.5f, y_max - y_min);
        float vx0 = -0.1f, vx1 = x_max + pad_x;
        float vy0 = y_min - pad_y, vy1 = y_max + pad_y;
        const float W = 560.0f, H = 300.0f;
        const float ML = 52.0f, MR = 20.0f, MT = 26.0f, MB = 36.0f;
        auto sx = [&](float xv) { return ML + (xv - vx0) / std::max(1e-5f, vx1 - vx0) * (W - ML - MR); };
        auto sy = [&](float yv) { return MT + (vy1 - yv) / std::max(1e-5f, vy1 - vy0) * (H - MT - MB); };

        std::fprintf(f, "<div class=\"card\"><div style=\"font-weight:600;margin-bottom:8px;\">3.1 Bend Stiffener Profile, Steel Insert &amp; 2D Deformed Elastica</div>\n");
        std::fprintf(f, "<svg viewBox=\"0 0 %.0f %.0f\">\n", W, H);
        std::fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#334155\" stroke-dasharray=\"4,4\"/>\n",
                     sx(0.0f), sy(0.0f), sx(x_max), sy(0.0f));
        // Outer cone profile
        std::fprintf(f, "<polyline fill=\"none\" stroke=\"#38bdf8\" stroke-width=\"2\" points=\"");
        for (size_t i = 0; i < n_nodes; ++i) std::fprintf(f, "%.1f,%.1f ", sx(x_stations[i]), sy(0.5f * get_dia(bs, x_stations[i])));
        std::fprintf(f, "\"/>\n<polyline fill=\"none\" stroke=\"#38bdf8\" stroke-width=\"2\" points=\"");
        for (size_t i = 0; i < n_nodes; ++i) std::fprintf(f, "%.1f,%.1f ", sx(x_stations[i]), sy(-0.5f * get_dia(bs, x_stations[i])));
        std::fprintf(f, "\"/>\n");
        // Inner bore
        std::fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#64748b\" stroke-dasharray=\"3,3\"/>\n",
                     sx(0.0f), sy(0.5f * bs.inner_dia), sx(bs.length), sy(0.5f * bs.inner_dia));
        std::fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#64748b\" stroke-dasharray=\"3,3\"/>\n",
                     sx(0.0f), sy(-0.5f * bs.inner_dia), sx(bs.length), sy(-0.5f * bs.inner_dia));
        // Embedded steel insert
        if (bs.insert_length > 1e-4f) {
            float hi_ins = std::max(bs.inner_dia, bs.root_dia - 0.02f);
            float r_ins = 0.5f * std::clamp(bs.insert_dia, bs.inner_dia, hi_ins);
            float r_in = 0.5f * bs.inner_dia;
            std::fprintf(f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"rgba(245,158,11,0.35)\" stroke=\"#f59e0b\" stroke-width=\"1.5\"/>\n",
                         sx(0.0f), sy(r_ins), sx(bs.insert_length) - sx(0.0f), sy(r_in) - sy(r_ins));
            std::fprintf(f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"rgba(245,158,11,0.35)\" stroke=\"#f59e0b\" stroke-width=\"1.5\"/>\n",
                         sx(0.0f), sy(-r_in), sx(bs.insert_length) - sx(0.0f), sy(-r_ins) - sy(-r_in));
        }
        // Deformed elastica centerline
        std::fprintf(f, "<polyline fill=\"none\" stroke=\"#22c55e\" stroke-width=\"2.6\" points=\"");
        for (size_t i = 0; i < n_nodes; ++i) std::fprintf(f, "%.1f,%.1f ", sx(def_x[i]), sy(def_y[i]));
        std::fprintf(f, "\"/>\n");
        // Bare pipe tail
        if (tail_x && tail_y && n_tail > 0) {
            std::fprintf(f, "<polyline fill=\"none\" stroke=\"#f97316\" stroke-width=\"2.2\" stroke-dasharray=\"5,3\" points=\"");
            for (size_t i = 0; i < n_tail; ++i) std::fprintf(f, "%.1f,%.1f ", sx(tail_x[i]), sy(tail_y[i]));
            std::fprintf(f, "\"/>\n");
        }
        std::fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" fill=\"#94a3b8\" font-size=\"11\">X [m] (0 to %.2f m) | Blue: Cone | Amber: Steel Insert | Green: Deformed | Orange: Pipe Tail</text>\n",
                     ML, H - 8.0f, x_max);
        std::fprintf(f, "</svg></div>\n");
    }

    // SVG 2: Outer-Fiber Strain Distribution vs Allowable Limit
    {
        float y_top = std::max(bs.allowable_strain * 1.25f, max_strain * 1.15f) * 100.0f;
        if (y_top < 1.0f) y_top = 8.0f;
        const float W = 560.0f, H = 300.0f;
        const float ML = 52.0f, MR = 20.0f, MT = 26.0f, MB = 36.0f;
        auto sx = [&](float sv) { return ML + (sv / std::max(1e-5f, bs.length)) * (W - ML - MR); };
        auto sy = [&](float ev_pct) { return MT + (1.0f - (ev_pct / y_top)) * (H - MT - MB); };

        std::fprintf(f, "<div class=\"card\"><div style=\"font-weight:600;margin-bottom:8px;\">3.2 Outer-Fiber Polyurethane Strain Distribution &epsilon;(s) [%%]</div>\n");
        std::fprintf(f, "<svg viewBox=\"0 0 %.0f %.0f\">\n", W, H);
        // Allowable strain limit line
        float y_lim = sy(bs.allowable_strain * 100.0f);
        std::fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#ef4444\" stroke-width=\"1.8\" stroke-dasharray=\"6,4\"/>\n",
                     ML, y_lim, W - MR, y_lim);
        std::fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" fill=\"#ef4444\" font-size=\"11\">Allowable Limit (%.2f%%)</text>\n",
                     W - MR - 135.0f, y_lim - 6.0f, bs.allowable_strain * 100.0f);
        // Steel insert tip marker
        if (bs.insert_length > 1e-4f && bs.insert_length < bs.length) {
            float x_ins = sx(bs.insert_length);
            std::fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#f59e0b\" stroke-dasharray=\"3,3\"/>\n",
                         x_ins, MT, x_ins, H - MB);
            std::fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" fill=\"#f59e0b\" font-size=\"10\">Insert Tip (%.2fm)</text>\n",
                         x_ins + 4.0f, MT + 12.0f, bs.insert_length);
        }
        // Strain curve
        std::fprintf(f, "<polyline fill=\"none\" stroke=\"#38bdf8\" stroke-width=\"2.5\" points=\"");
        for (size_t i = 0; i < n_nodes; ++i) {
            std::fprintf(f, "%.1f,%.1f ", sx(x_stations[i]), sy(strain_arr[i] * 100.0f));
        }
        std::fprintf(f, "\"/>\n");
        std::fprintf(f, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"4.5\" fill=\"#f59e0b\"/>\n",
                     sx(max_strain_s), sy(max_strain * 100.0f));
        std::fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" fill=\"#94a3b8\" font-size=\"11\">Arc Length s [m] | Peak Strain = %.2f%% at s = %.2f m</text>\n",
                     ML, H - 8.0f, max_strain * 100.0f, max_strain_s);
        std::fprintf(f, "</svg></div>\n");
    }

    // SVG 3: Bending Moment M(s) & Curvature kappa(s)
    {
        float max_m_kNm = 0.1f;
        float max_k = 0.05f;
        for (size_t i = 0; i < n_nodes; ++i) {
            max_m_kNm = std::max(max_m_kNm, std::fabs(moment_Nm[i]) / 1000.0f);
            max_k = std::max(max_k, curvature_arr[i]);
        }
        float k_allow = (bs.allowable_mbr_m > 1e-4f) ? (1.0f / bs.allowable_mbr_m) : 0.4f;
        max_k = std::max(max_k * 1.15f, k_allow * 1.25f);

        const float W = 560.0f, H = 300.0f;
        const float ML = 52.0f, MR = 20.0f, MT = 26.0f, MB = 36.0f;
        auto sx = [&](float sv) { return ML + (sv / std::max(1e-5f, bs.length)) * (W - ML - MR); };
        auto sy_k = [&](float kv) { return MT + (1.0f - (kv / max_k)) * (H - MT - MB); };
        auto sy_m = [&](float mv) { return MT + (1.0f - (mv / (max_m_kNm * 1.15f))) * (H - MT - MB); };

        std::fprintf(f, "<div class=\"card\"><div style=\"font-weight:600;margin-bottom:8px;\">3.3 Pipe Curvature &kappa;(s) [1/m] &amp; Bending Moment |M(s)| [kN&middot;m]</div>\n");
        std::fprintf(f, "<svg viewBox=\"0 0 %.0f %.0f\">\n", W, H);
        float y_klim = sy_k(k_allow);
        std::fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#ef4444\" stroke-width=\"1.8\" stroke-dasharray=\"6,4\"/>\n",
                     ML, y_klim, W - MR, y_klim);
        std::fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" fill=\"#ef4444\" font-size=\"11\">Max Curvature 1/MBR (%.3f m⁻¹)</text>\n",
                     W - MR - 175.0f, y_klim - 6.0f, k_allow);
        // Curvature curve
        std::fprintf(f, "<polyline fill=\"none\" stroke=\"#a855f7\" stroke-width=\"2.5\" points=\"");
        for (size_t i = 0; i < n_nodes; ++i) std::fprintf(f, "%.1f,%.1f ", sx(x_stations[i]), sy_k(curvature_arr[i]));
        std::fprintf(f, "\"/>\n");
        // Bending moment curve
        std::fprintf(f, "<polyline fill=\"none\" stroke=\"#22c55e\" stroke-width=\"2.0\" stroke-dasharray=\"4,2\" points=\"");
        for (size_t i = 0; i < n_nodes; ++i) std::fprintf(f, "%.1f,%.1f ", sx(x_stations[i]), sy_m(std::fabs(moment_Nm[i]) / 1000.0f));
        std::fprintf(f, "\"/>\n");
        std::fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" fill=\"#94a3b8\" font-size=\"11\">Purple: Curvature &kappa;(s) (Min MBR = %.2f m) | Green Dashed: |M(s)| (M₀ = %.2f kN&middot;m)</text>\n",
                     ML, H - 8.0f, min_mbr, std::fabs(root_moment_kNm));
        std::fprintf(f, "</svg></div>\n");
    }

    // SVG 4: Multi-Temperature Polyurethane Stress-Strain Curves (Cold, Nominal, Hot)
    {
        Material m_cold = mat.with_temp_mode(MAT_TEMP_COLD);
        Material m_nom  = mat.with_temp_mode(MAT_TEMP_NOMINAL);
        Material m_hot  = mat.with_temp_mode(MAT_TEMP_HOT);
        float max_eps_plot = 0.25f;
        float max_sig_mpa = std::max(1.0f, (max_eps_plot * get_material_E(max_eps_plot, m_cold)) / 1e6f * 1.1f);

        const float W = 560.0f, H = 300.0f;
        const float ML = 52.0f, MR = 20.0f, MT = 26.0f, MB = 36.0f;
        auto sx = [&](float ev) { return ML + (ev / max_eps_plot) * (W - ML - MR); };
        auto sy = [&](float smpa) { return MT + (1.0f - (smpa / max_sig_mpa)) * (H - MT - MB); };

        std::fprintf(f, "<div class=\"card\"><div style=\"font-weight:600;margin-bottom:8px;\">3.4 Multi-Temperature Polyurethane Stress-Strain Curves &sigma;(&epsilon;) [MPa]</div>\n");
        std::fprintf(f, "<svg viewBox=\"0 0 %.0f %.0f\">\n", W, H);
        const Material* mats[3] = { &m_cold, &m_nom, &m_hot };
        const char* colors[3] = { "#38bdf8", "#22c55e", "#f97316" };
        for (int m_i = 0; m_i < 3; ++m_i) {
            std::fprintf(f, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"2.3\" points=\"", colors[m_i]);
            for (int k = 0; k <= 40; ++k) {
                float e = max_eps_plot * (float)k / 40.0f;
                float s_mpa = (e * get_material_E(e, *mats[m_i])) / 1e6f;
                std::fprintf(f, "%.1f,%.1f ", sx(e), sy(s_mpa));
            }
            std::fprintf(f, "\"/>\n");
        }
        std::fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" fill=\"#94a3b8\" font-size=\"11\">Strain &epsilon; [0 to 25%%] | Blue: Cold (%.0f&deg;C) | Green: Nominal (%.0f&deg;C) | Orange: Hot (%.0f&deg;C)</text>\n",
                     ML, H - 8.0f, mat.temp_cold_c, mat.temp_nominal_c, mat.temp_hot_c);
        std::fprintf(f, "</svg></div>\n");
    }
    std::fprintf(f, "</div>\n");

    // Section 4: Multi-Temperature Envelope Table
    std::fprintf(f, "<h2>4. Multi-Temperature Envelope Verification (Cold / Nominal / Hot)</h2>\n"
                    "<div class=\"card\"><table>\n"
                    "<tr><th>Temperature Condition</th><th>Temp [&deg;C]</th><th>Initial E<sub>0</sub> [MPa]</th>"
                    "<th>Peak PU Strain [%%]</th><th>Min Pipe MBR [m]</th><th>Root Moment M<sub>0</sub> [kN&middot;m]</th>"
                    "<th>Root Shear V<sub>0</sub> [kN]</th><th>Status</th></tr>\n");
    for (int i = 0; i < 3; ++i) {
        const auto &r = temp_env.rows[i];
        std::fprintf(f, "<tr><td>%s</td><td>%.1f</td><td>%.1f</td><td>%.2f%%</td><td>%.2f</td><td>%.2f</td><td>%.2f</td>"
                        "<td><b style=\"color:%s;\">%s</b></td></tr>\n",
                     r.label, r.temp_c, r.initial_E_mpa, r.max_strain * 100.0f, r.min_mbr_m,
                     r.root_moment_kNm, r.root_shear_kN, r.pass ? "var(--pass)" : "var(--fail)", r.pass ? "PASS" : "FAIL");
    }
    std::fprintf(f, "</table></div>\n");

    // Section 5: Nodal Profile Table
    std::fprintf(f, "<h2>5. Full Nodal Solution Profile (s &isin; [0, L])</h2>\n"
                    "<div class=\"card\"><table>\n"
                    "<tr><th>Node</th><th>s [m]</th><th>X<sub>def</sub> [m]</th><th>Y<sub>def</sub> [m]</th>"
                    "<th>D<sub>outer</sub> [m]</th><th>D<sub>PU,in</sub> [m]</th><th>&theta; [deg]</th>"
                    "<th>M [kN&middot;m]</th><th>V [kN]</th><th>Strain [%%]</th><th>&kappa; [1/m]</th><th>MBR [m]</th></tr>\n");
    for (size_t i = 0; i < n_nodes; ++i) {
        float s = x_stations[i];
        std::fprintf(f, "<tr><td>%zu</td><td>%.3f</td><td>%.3f</td><td>%.3f</td><td>%.3f</td><td>%.3f</td>"
                        "<td>%.2f</td><td>%.2f</td><td>%.2f</td><td>%.2f%%</td><td>%.4f</td><td>%.2f</td></tr>\n",
                     i, s, def_x[i], def_y[i], get_dia(bs, s), get_pu_inner_dia(bs, s),
                     angle_rad[i] * (180.0f / 3.14159265358979323846f),
                     moment_Nm[i] / 1000.0f, shear_N[i] / 1000.0f,
                     strain_arr[i] * 100.0f, curvature_arr[i], mbr_arr[i]);
    }
    std::fprintf(f, "</table></div>\n");

    std::fprintf(f, "</div>\n</body>\n</html>\n");
    std::fclose(f);
    return true;
}



