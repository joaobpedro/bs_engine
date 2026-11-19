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
};

inline bool is_valid_geometry(const BendStiffner &bs) {
    return (bs.length > 1e-4f) &&
           (bs.inner_dia >= 0.0f) &&
           (bs.root_dia > bs.inner_dia + 1e-5f) &&
           (bs.tip_dia > bs.inner_dia + 1e-5f) &&
           (bs.root_cyl_length >= 0.0f) &&
           (bs.root_cyl_length < bs.length);
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
        assert(index < capacity);
        return items[index < capacity ? index : (capacity - 1)];
    }

    const T& operator[](size_t index) const {
        assert(index < capacity);
        return items[index < capacity ? index : (capacity - 1)];
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
            if (items && capacity > 0) {
                std::memset(items, 0, capacity * sizeof(T));
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

struct Material {
    bool is_nonlinear;
    bool table_is_tangent_modulus;       // true: table is d(sigma)/d(eps), integrate for stress & secant E
    bool use_cross_section_integration;  // true: integrate stress across annular depth y in [0, R_o]
    float linear_E_kpa;
    Vec1<float> strain_table;
    Vec1<float> emod_table_kpa;

    Material()
        : is_nonlinear(true),
          table_is_tangent_modulus(true),
          use_cross_section_integration(false),
          linear_E_kpa(215800.0f),
          strain_table(DEFAULT_MATERIAL_DATA_SIZE),
          emod_table_kpa(DEFAULT_MATERIAL_DATA_SIZE) {
        reset_to_default();
    }

    void reset_to_default(float scale_factor = 1.0f) {
        strain_table.clear();
        emod_table_kpa.clear();
        for (int i = 0; i < DEFAULT_MATERIAL_DATA_SIZE; ++i) {
            strain_table.append(DEFAULT_STRAIN_VECTOR[i]);
            emod_table_kpa.append(DEFAULT_EMOD_VECTOR_KPA[i] * scale_factor);
        }
        linear_E_kpa = DEFAULT_EMOD_VECTOR_KPA[0] * scale_factor;
    }

    bool has_valid_table() const {
        return strain_table.count >= 2 &&
               strain_table.count == emod_table_kpa.count &&
               emod_table_kpa.items[0] > 0.0f;
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

    for (int i = 0; i < DEFAULT_MATERIAL_DATA_SIZE - 1; ++i) {
        float x0 = DEFAULT_STRAIN_VECTOR[i];
        float x1 = DEFAULT_STRAIN_VECTOR[i + 1];
        if (s >= x0 && s <= x1) {
            float y0 = DEFAULT_EMOD_VECTOR_KPA[i];
            float y1 = DEFAULT_EMOD_VECTOR_KPA[i + 1];
            float dx = x1 - x0;
            float result = (dx > 1e-12f) ? (y0 + (s - x0) * (y1 - y0) / dx) : y0;
            return result * 1000.0f;
        }
    }
    return DEFAULT_EMOD_VECTOR_KPA[DEFAULT_MATERIAL_DATA_SIZE - 1] * 1000.0f;
}

// Piecewise linear interpolation of custom non-linear E table (returns Pa)
inline float get_non_linear_E_custom(float strain, const Vec1<float> &strain_input_vector, const Vec1<float> &E_mod_input_vector) {
    assert(strain_input_vector.count == E_mod_input_vector.count);
    size_t n = std::min(strain_input_vector.count, E_mod_input_vector.count);
    if (n == 0) return DEFAULT_EMOD_VECTOR_KPA[0] * 1000.0f;
    if (n == 1) return E_mod_input_vector.items[0] * 1000.0f;

    float s = std::fabs(strain);
    if (s <= strain_input_vector.items[0]) {
        return E_mod_input_vector.items[0] * 1000.0f;
    }
    if (s >= strain_input_vector.items[n - 1]) {
        return E_mod_input_vector.items[n - 1] * 1000.0f;
    }

    for (size_t i = 0; i < n - 1; ++i) {
        float x0 = strain_input_vector.items[i];
        float x1 = strain_input_vector.items[i + 1];
        if (s >= x0 && s <= x1) {
            float y0 = E_mod_input_vector.items[i];
            float y1 = E_mod_input_vector.items[i + 1];
            float dx = x1 - x0;
            float result = (dx > 1e-12f) ? (y0 + (s - x0) * (y1 - y0) / dx) : y0;
            return result * 1000.0f;
        }
    }
    return E_mod_input_vector.items[n - 1] * 1000.0f;
}

// Evaluates point-wise fiber Secant Young's Modulus E_sec(strain) in Pa
inline float get_fiber_E_secant(float strain, const Material &mat) {
    if (!mat.is_nonlinear || !mat.has_valid_table()) {
        return std::max(1.0f, mat.linear_E_kpa) * 1000.0f;
    }

    float s = std::fabs(strain);
    size_t n = mat.strain_table.count;
    const float* xs = mat.strain_table.items;
    const float* ys = mat.emod_table_kpa.items;

    if (!mat.table_is_tangent_modulus) {
        return get_non_linear_E_custom(s, mat.strain_table, mat.emod_table_kpa);
    }

    if (s <= 1e-7f || s <= xs[0]) {
        return std::max(1.0f, ys[0]) * 1000.0f;
    }

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
    // Weights w_q integrate (16/pi) * xi^2 * sqrt(1 - xi^2) d(xi) = 1
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
        double E_pa = std::max(1.0, (double)mat.linear_E_kpa) * 1000.0;
        return (float)(sigma_req_pa / E_pa);
    }

    size_t n = mat.strain_table.count;
    const float* xs = mat.strain_table.items;
    const float* ys = mat.emod_table_kpa.items;
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
    float cyl_len = std::clamp(dimensions.root_cyl_length, 0.0f, dimensions.length - 1e-5f);
    if (s <= cyl_len) {
        return dimensions.root_dia;
    }
    float cone_len = dimensions.length - cyl_len;
    float t = std::clamp((s - cyl_len) / cone_len, 0.0f, 1.0f);
    return dimensions.root_dia - (dimensions.root_dia - dimensions.tip_dia) * t;
}

inline float get_Inertia(const BendStiffner& dimensions, float x) {
    constexpr float pi_over_64 = 3.14159265358979323846f / 64.0f;
    float outer_dia = get_dia(dimensions, x);
    float inner_dia = dimensions.inner_dia;
    float do2 = outer_dia * outer_dia;
    float di2 = inner_dia * inner_dia;
    float inertia = pi_over_64 * (do2 * do2 - di2 * di2);
    return std::max(1e-12f, inertia);
}

inline float get_Area(const BendStiffner& dimensions, float x) {
    constexpr float pi_over_4 = 3.14159265358979323846f / 4.0f;
    float outer_dia = get_dia(dimensions, x);
    float inner_dia = dimensions.inner_dia;
    float area = pi_over_4 * (outer_dia * outer_dia - inner_dia * inner_dia);
    return std::max(1e-10f, area);
}

inline float get_EI(const BendStiffner& dimensions, float strain, float x) {
    float inertia = get_Inertia(dimensions, x);
    float pipe_EI_Pa = std::max(0.0f, dimensions.pipe_EI_kNm2) * 1000.0f;
    return get_non_linear_E(strain) * inertia + pipe_EI_Pa;
}

inline float get_EI(const BendStiffner& dimensions, float strain, float x, const Material &mat) {
    float inertia = get_Inertia(dimensions, x);
    float pipe_EI_Pa = std::max(0.0f, dimensions.pipe_EI_kNm2) * 1000.0f;
    return get_material_E(strain, mat) * inertia + pipe_EI_Pa;
}

// Solves the coupled composite cross-section (Polyurethane Cone + Internal Pipe EI_pipe) at station s
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
    float pipe_EI = std::max(0.0f, bs.pipe_EI_kNm2) * 1000.0f; // N.m^2

    double abs_M = std::fabs((double)m_total_Nm);
    double sigma_axial = (double)std::fabs(axial_N) / (double)A_bs;

    if (pipe_EI <= 1e-3f) {
        // Pure bend stiffener section
        double sigma_req = (abs_M * (double)r_o) / (double)I_bs + sigma_axial;
        float eps = solve_strain_from_stress(sigma_req, mat);
        float E_bs = get_material_E(eps, mat);
        out_EI_total = E_bs * I_bs;
        out_curvature = m_total_Nm / out_EI_total;
        out_total_strain = eps;
        return;
    }

    // Composite section: M_total = (E_bs(eps) * I_bs + EI_pipe) * kappa, with eps = |kappa| * r_o + eps_axial
    // Initialize from uncoupled estimate and refine in 3 fast iterations
    double E0 = (double)get_material_E(0.01f, mat);
    double frac_bs = (E0 * I_bs) / (E0 * I_bs + pipe_EI);
    double sigma_req = (abs_M * frac_bs * (double)r_o) / (double)I_bs + sigma_axial;
    float eps = solve_strain_from_stress(sigma_req, mat);

    for (int k = 0; k < 3; ++k) {
        double E_bs = (double)get_material_E(eps, mat);
        double EI_tot = E_bs * (double)I_bs + (double)pipe_EI;
        double kappa_abs = abs_M / EI_tot;
        double sigma_b = E_bs * kappa_abs * (double)r_o;
        eps = solve_strain_from_stress(sigma_b + sigma_axial, mat);
    }

    float E_final = get_material_E(eps, mat);
    out_EI_total = E_final * I_bs + pipe_EI;
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
    float axial_N = std::fabs(tension_N) * std::max(0.0f, std::cos(target_angle_rad - y.items[1]));
    float kappa = 0.0f, local_strain = 0.0f, EI_total = 1.0f;
    solve_composite_section(bs, mat, s_or_x, y.items[2], axial_N, kappa, local_strain, EI_total);

    bool far_field = bs.use_far_field_bc && (std::fabs(tension_N) > 1e-2f);
    float shear_v = far_field
        ? (-std::fabs(tension_N) * std::sin(target_angle_rad - y.items[1]))
        : y.items[3];

    State4 ds;
    if (bs.use_arc_length) {
        // Exact Elastica in arc-length s in [0, L]: singularity-free at all angles!
        ds.x_pos = std::cos(y.items[1]);
        ds.items[0] = std::sin(y.items[1]);          // dy/ds = sin(theta)
        ds.items[1] = kappa;                         // dtheta/ds = M / EI_total
        ds.items[2] = shear_v;                       // dM/ds = V(s)
        ds.items[3] = 0.0f;
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
        ds.items[3] = 0.0f;
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
            double angle_diff = (double)thetaL - (double)end_state.items[1];
            double m_tip_target = 2.0 * std::sqrt(pipe_EI_Pa * T_abs) * std::sin(0.5 * angle_diff);
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

    for (size_t i = 0; i < n; ++i) {
        float theta_s = results.items[i].items[1];
        float m_s = results.items[i].items[2];
        float axial_N = std::fabs(tension_N) * std::max(0.0f, std::cos(tip_angle_rad - theta_s));

        float kappa = 0.0f, total_eps = 0.0f, EI_tot = 1.0f;
        solve_composite_section(bs, mat, s, m_s, axial_N, kappa, total_eps, EI_tot);

        strain.items[i] = total_eps;
        float abs_kappa = std::fabs(kappa);
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

// Calculates Polyurethane volume [m^3] (cylindrical root + conical frustum - inner bore)
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
    return std::max(0.0f, v_cyl_out + v_cone_out - v_bore);
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
};

struct OptimizationResult {
    bool converged = false;
    bool feasible = false;
    BendStiffner best_bs = {};
    float volume_m3 = 0.0f;
    float mass_kg = 0.0f;
    float max_strain = 0.0f;
    float min_mbr_m = 1000.0f;
    float max_root_moment_kNm = 0.0f;
    float strain_utilization = 0.0f;
    float mbr_utilization = 0.0f;
    int evaluations = 0;
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
                                      float &out_root_m0_kNm) {
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

        float abs_k = std::fabs(kappa);
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
    float pipe_EI_Pa = std::max(0.0f, initial_bs.pipe_EI_kNm2) * 1000.0f;

    auto eval_objective = [&](float L, float Dr, size_t steps,
                              bool &out_feas, float &out_max_eps, float &out_min_mbr, float &out_max_m0) -> double {
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
        float avg_eps_accum = 0.0f;

        for (int c = 0; c < num_cases; ++c) {
            BendStiffner case_bs = cand;
            if (cfg.force_far_field_when_tension && std::fabs(cases[c].tension_N) > 10.0f) {
                case_bs.use_far_field_bc = true;
            }
            float c_max_eps = 0.0f, c_avg_eps = 0.0f, c_tip_eps = 0.0f, c_root_eps = 0.0f, c_mbr = 1000.0f, c_m0 = 0.0f;
            evaluate_candidate_single(case_bs, mat, steps, cases[c].angle_rad, cases[c].tension_N,
                                      c_max_eps, c_avg_eps, c_tip_eps, c_root_eps, c_mbr, c_m0);
            if (c_max_eps > worst_eps) worst_eps = c_max_eps;
            if (c_mbr < worst_mbr) worst_mbr = c_mbr;
            if (c_m0 > worst_m0) worst_m0 = c_m0;
            avg_eps_accum += c_avg_eps;
        }
        avg_eps_accum /= (float)num_cases;

        out_max_eps = worst_eps;
        out_min_mbr = worst_mbr;
        out_max_m0 = worst_m0;

        double u_eps = (double)worst_eps / (double)allow_eps_eff;
        double u_mbr = (double)allow_mbr_eff / std::max(1e-4, (double)worst_mbr);

        // Minimum root stiffness requirement relative to internal pipe when tension is zero
        double EI_bs_root = (double)get_material_E(0.01f, mat) * (double)get_Inertia(cand, 0.0f);
        double u_stiff = 0.0;
        if (pipe_EI_Pa > 1.0f) {
            double req_root_EI = (double)cfg.min_root_stiffness_ratio * (double)pipe_EI_Pa;
            u_stiff = req_root_EI / std::max(1.0, EI_bs_root);
        }

        out_feas = (u_eps <= 1.0) && (u_mbr <= 1.0) && (u_stiff <= 1.0);

        double vol = (double)get_pu_volume_m3(cand);
        double non_uniformity = (worst_eps > 1e-6f) ? ((double)(worst_eps - avg_eps_accum) / (double)worst_eps) : 0.0;
        double base_cost = vol * (1.0 + 0.20 * non_uniformity);

        if (out_feas) {
            return base_cost;
        }

        double pen = 0.0;
        if (u_eps > 1.0)   pen += (u_eps - 1.0) * (u_eps - 1.0) * 100.0;
        if (u_mbr > 1.0)   pen += (u_mbr - 1.0) * (u_mbr - 1.0) * 100.0;
        if (u_stiff > 1.0) pen += (u_stiff - 1.0) * (u_stiff - 1.0) * 25.0;
        return 1000.0 + pen + 0.1 * base_cost;
    };

    // Stage 1: 2D Global Grid Search over (Length L, Root Diameter D_root)
    const size_t coarse_steps = 20;
    int nL  = cfg.optimize_length ? 11 : 1;
    int nDr = cfg.optimize_root_dia ? 11 : 1;

    float best_L = std::clamp(initial_bs.length, L_lo, L_hi);
    float best_Dr = std::clamp(initial_bs.root_dia, Dr_lo, Dr_hi);
    double best_cost = 1e30;

    for (int iL = 0; iL < nL; ++iL) {
        float L = (nL > 1) ? (L_lo + (L_hi - L_lo) * (float)iL / (float)(nL - 1)) : L_lo;
        for (int iDr = 0; iDr < nDr; ++iDr) {
            float Dr = (nDr > 1) ? (Dr_lo + (Dr_hi - Dr_lo) * (float)iDr / (float)(nDr - 1)) : Dr_lo;
            bool feas = false;
            float me = 0, mm = 0, m0 = 0;
            double cost = eval_objective(L, Dr, coarse_steps, feas, me, mm, m0);
            if (cost < best_cost) {
                best_cost = cost;
                best_L = L;
                best_Dr = Dr;
            }
        }
    }

    // Stage 2: Multi-scale 2D coordinate pattern search refinement at full resolution
    const size_t fine_steps = (DISCRETIZATION > 1) ? (DISCRETIZATION - 1) : 24;
    bool feas_now = false;
    float eps_now = 0, mbr_now = 0, m0_now = 0;
    best_cost = eval_objective(best_L, best_Dr, fine_steps, feas_now, eps_now, mbr_now, m0_now);

    float step_scales[4] = { 0.10f, 0.035f, 0.010f, 0.002f };
    for (int sc = 0; sc < 4; ++sc) {
        float dL  = cfg.optimize_length ? ((L_hi - L_lo) * step_scales[sc]) : 0.0f;
        float dDr = cfg.optimize_root_dia ? ((Dr_hi - Dr_lo) * step_scales[sc]) : 0.0f;

        for (int pass = 0; pass < 3; ++pass) {
            bool improved = false;
            if (cfg.optimize_length && dL > 1e-4f) {
                for (int dir = -1; dir <= 1; dir += 2) {
                    float cand_L = std::clamp(best_L + (float)dir * dL, L_lo, L_hi);
                    bool f = false; float me = 0, mm = 0, m0 = 0;
                    double c = eval_objective(cand_L, best_Dr, fine_steps, f, me, mm, m0);
                    if (c < best_cost) {
                        best_cost = c; best_L = cand_L; improved = true;
                    }
                }
            }
            if (cfg.optimize_root_dia && dDr > 1e-4f) {
                for (int dir = -1; dir <= 1; dir += 2) {
                    float cand_Dr = std::clamp(best_Dr + (float)dir * dDr, Dr_lo, Dr_hi);
                    bool f = false; float me = 0, mm = 0, m0 = 0;
                    double c = eval_objective(best_L, cand_Dr, fine_steps, f, me, mm, m0);
                    if (c < best_cost) {
                        best_cost = c; best_Dr = cand_Dr; improved = true;
                    }
                }
            }
            if (!improved) break;
        }
    }

    // Final evaluation of optimal design
    eval_objective(best_L, best_Dr, fine_steps, res.feasible, res.max_strain, res.min_mbr_m, res.max_root_moment_kNm);
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
    return res;
}

