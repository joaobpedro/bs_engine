#pragma once

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <type_traits>

#define DISCRETIZATION 50
#define NUMBER_OF_STATES 4

// Bend stiffener dimensions (in meters)
struct BendStiffner {
    float length;
    float root_dia;
    float tip_dia;
    float inner_dia;
};

inline bool is_valid_geometry(const BendStiffner &bs) {
    return (bs.length > 1e-4f) &&
           (bs.inner_dia >= 0.0f) &&
           (bs.root_dia > bs.inner_dia + 1e-5f) &&
           (bs.tip_dia > bs.inner_dia + 1e-5f);
}

struct Vec2 {
    float x;
    float y;
};

// Stack-allocated 4-state vector [y, theta, M, V] for zero-allocation ODE integration
struct State4 {
    float items[NUMBER_OF_STATES];
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

    // Bounds-guarded element access
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

    // Returns the largest positive value or zero if all elements are negative/empty
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
    bool table_is_tangent_modulus; // true: table is d(sigma)/d(eps), integrate for stress & secant E
    float linear_E_kpa;
    Vec1<float> strain_table;
    Vec1<float> emod_table_kpa;

    Material()
        : is_nonlinear(true),
          table_is_tangent_modulus(true),
          linear_E_kpa(215800.0f),
          strain_table(DEFAULT_MATERIAL_DATA_SIZE),
          emod_table_kpa(DEFAULT_MATERIAL_DATA_SIZE) {
        reset_to_default();
    }

    void reset_to_default() {
        strain_table.clear();
        emod_table_kpa.clear();
        for (int i = 0; i < DEFAULT_MATERIAL_DATA_SIZE; ++i) {
            strain_table.append(DEFAULT_STRAIN_VECTOR[i]);
            emod_table_kpa.append(DEFAULT_EMOD_VECTOR_KPA[i]);
        }
        linear_E_kpa = DEFAULT_EMOD_VECTOR_KPA[0];
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
            return result * 1000.0f; // kPa -> Pa
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
            return result * 1000.0f; // kPa -> Pa
        }
    }
    return E_mod_input_vector.items[n - 1] * 1000.0f;
}

// Evaluates effective Secant Young's Modulus E_sec(strain) in Pa for the total moment-curvature formulation
inline float get_material_E(float strain, const Material &mat) {
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

    // Integrate tangent modulus E_tan(u) from 0 to s to obtain true stress sigma(s),
    // then return secant modulus E_sec(s) = sigma(s) / s (in Pa).
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

// Analytically inverts the monotonic stress-strain law sigma(eps) = eps * E_sec(eps) = sigma_req_pa
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

    if (!mat.table_is_tangent_modulus) {
        // Direct secant table: build enforced-monotonic stress nodes sigma_i = x_i * y_i
        double prev_x = 0.0;
        double prev_sig = 0.0;
        double min_slope = std::max(1.0, (double)ys[n - 1] * 0.1);
        for (size_t i = 0; i < n; ++i) {
            double xi = xs[i];
            double yi = std::max(1.0, (double)ys[i]);
            if (xi <= prev_x) continue;
            double sig_i = std::max(xi * yi, prev_sig + min_slope * (xi - prev_x));
            if (sigma_kpa <= sig_i) {
                double t = (sigma_kpa - prev_sig) / (sig_i - prev_sig);
                return (float)(prev_x + t * (xi - prev_x));
            }
            prev_x = xi;
            prev_sig = sig_i;
        }
        double last_slope = std::max(1.0, (double)ys[n - 1]);
        return (float)(prev_x + (sigma_kpa - prev_sig) / last_slope);
    }

    // Tangent modulus table: exact piecewise quadratic stress-strain inversion
    double acc_sigma = 0.0;
    if (xs[0] > 0.0f) {
        double y0 = std::max(1.0, (double)ys[0]);
        double d_sig0 = (double)xs[0] * y0;
        if (sigma_kpa <= d_sig0) {
            return (float)(sigma_kpa / y0);
        }
        acc_sigma += d_sig0;
    }

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
            return (float)(x0 + std::clamp(ds, 0.0, dx));
        }
        acc_sigma += d_sigma;
    }

    double y_last = std::max(1.0, (double)ys[n - 1]);
    return (float)((double)xs[n - 1] + (sigma_kpa - acc_sigma) / y_last);
}

inline float get_dia(const BendStiffner &dimensions, float x) {
    if (dimensions.length <= 1e-6f) return dimensions.root_dia;
    float t = std::clamp(x / dimensions.length, 0.0f, 1.0f);
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
    return get_non_linear_E(strain) * inertia;
}

inline float get_EI(const BendStiffner& dimensions, float strain, float x, const Material &mat) {
    float inertia = get_Inertia(dimensions, x);
    return get_material_E(strain, mat) * inertia;
}

// Fast zero-allocation ODE right-hand side with local non-linear constitutive coupling
inline State4 equations_fast(float x,
                             const State4 &y,
                             const BendStiffner &bs,
                             const Material &mat,
                             float tension_N = 0.0f,
                             float tip_angle_rad = 0.0f) {
    float d_o = get_dia(bs, x);
    float I = get_Inertia(bs, x);
    float A = get_Area(bs, x);

    // Axial force component along current beam tangent
    float axial_N = std::fabs(tension_N) * std::max(0.0f, std::cos(tip_angle_rad - y.items[1]));
    double sigma_req = (std::fabs((double)y.items[2]) * (0.5 * (double)d_o)) / (double)I +
                       ((double)axial_N / (double)A);
    float local_strain = solve_strain_from_stress(sigma_req, mat);
    float EI = get_material_E(local_strain, mat) * I;

    // Guard angle away from +/- pi/2 singularity during shooting iterations
    constexpr float max_theta = 1.553343f; // ~89 deg
    float theta_clamped = std::clamp(y.items[1], -max_theta, max_theta);
    float cos_t = std::cos(theta_clamped);
    float tan_t = std::tan(theta_clamped);

    State4 dydx;
    dydx.items[0] = tan_t;
    dydx.items[1] = y.items[2] / (EI * cos_t);
    dydx.items[2] = y.items[3];
    dydx.items[3] = 0.0f;
    return dydx;
}

// Fast zero-allocation 4th-order Runge-Kutta step
inline State4 RK4_fast(float x,
                       const State4 &y,
                       float h,
                       const BendStiffner &bs,
                       const Material &mat,
                       float tension_N = 0.0f,
                       float tip_angle_rad = 0.0f) {
    State4 k1 = equations_fast(x, y, bs, mat, tension_N, tip_angle_rad);

    State4 y_temp;
    for (int i = 0; i < NUMBER_OF_STATES; ++i) {
        y_temp.items[i] = y.items[i] + 0.5f * h * k1.items[i];
    }
    State4 k2 = equations_fast(x + 0.5f * h, y_temp, bs, mat, tension_N, tip_angle_rad);

    for (int i = 0; i < NUMBER_OF_STATES; ++i) {
        y_temp.items[i] = y.items[i] + 0.5f * h * k2.items[i];
    }
    State4 k3 = equations_fast(x + 0.5f * h, y_temp, bs, mat, tension_N, tip_angle_rad);

    for (int i = 0; i < NUMBER_OF_STATES; ++i) {
        y_temp.items[i] = y.items[i] + h * k3.items[i];
    }
    State4 k4 = equations_fast(x + h, y_temp, bs, mat, tension_N, tip_angle_rad);

    State4 y_next;
    for (int i = 0; i < NUMBER_OF_STATES; ++i) {
        y_next.items[i] = y.items[i] + (h / 6.0f) * (k1.items[i] + 2.0f * k2.items[i] + 2.0f * k3.items[i] + k4.items[i]);
    }
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
                         float tip_angle_rad = 0.0f) {
    float h = bs.length / (float)steps;
    float x = 0.0f;
    State4 y = { { y0, theta0, guessed_m0, guessed_v0 } };
    for (size_t i = 0; i < steps; ++i) {
        y = RK4_fast(x, y, h, bs, mat, tension_N, tip_angle_rad);
        x += h;
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
    Material default_mat;
    State4 s = { { y.items[0], y.items[1], y.items[2], y.items[3] } };
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

// Solves the BVP with guaranteed bracketed Illinois-Secant/Bisection convergence
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

    // Because dM/dx = V and dV/dx = 0, M(x) = M0 + V0 * x,
    // so M(L) = targetML uniquely determines V0(M0) = (targetML - M0) / L.
    auto v0_for_m0 = [&](double m0) -> double {
        return ((double)targetML - m0) / (double)bs.length;
    };

    double dtheta = (double)thetaL - (double)theta0;
    if (std::fabs(dtheta) < 1e-8 && std::fabs(targetML) < 1e-8) {
        return { 0.0f, 0.0f };
    }

    double sgn = (dtheta >= 0.0) ? 1.0 : -1.0;
    double target_abs = std::fabs(dtheta);

    auto eval_residual = [&](double m_abs) -> double {
        double m0 = sgn * m_abs;
        double v0 = v0_for_m0(m0);
        State4 end_state = shoot_fast(bs, steps, y0, theta0, (float)m0, (float)v0, mat, tension_N, thetaL);
        return sgn * ((double)end_state.items[1] - (double)theta0) - target_abs;
    };

    double a = 0.0;
    double fa = eval_residual(a);
    if (std::fabs(fa) < 1e-7) {
        return { 0.0f, (float)v0_for_m0(0.0) };
    }

    double EI_mid = (double)get_material_E(0.02f, mat) * (double)get_Inertia(bs, 0.5f * bs.length);
    double b = std::max(100.0, 0.8 * EI_mid * target_abs / (double)bs.length);
    double fb = eval_residual(b);

    // Geometric bracket expansion until fa <= 0 <= fb
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

    for (int k = 0; k < 25; ++k) {
        if (f_best < 1e-6 || std::fabs(x_hi - x_lo) <= 1e-6 * (1.0 + std::fabs(m_best))) {
            break;
        }

        double x_mid = 0.5 * (x_lo + x_hi);
        double f_mid = eval_residual(x_mid);
        if (std::fabs(f_mid) < f_best) {
            f_best = std::fabs(f_mid);
            m_best = x_mid;
        }
        if (f_best < 1e-6) break;

        double denom = std::sqrt(std::max(0.0, f_mid * f_mid - f_lo * f_hi));
        if (denom <= 1e-14) {
            if (f_mid > 0.0) { x_hi = x_mid; f_hi = f_mid; }
            else { x_lo = x_mid; f_lo = f_mid; }
            continue;
        }

        // Since f_lo <= 0 <= f_hi, sgn(f_lo - f_hi) = -1
        double x_new = x_mid - (x_mid - x_lo) * (f_mid / denom);
        double f_new = eval_residual(x_new);
        if (std::fabs(f_new) < f_best) {
            f_best = std::fabs(f_new);
            m_best = x_new;
        }
        if (f_best < 1e-6) break;

        // Tighten bracket using x_mid and x_new
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

// Calculates total outer-fiber strain (bending + axial tension) at all discretization nodes
inline void calculate_strain_material(size_t steps,
                                      const BendStiffner &bs,
                                      const Vec1<Vec1<float>> &results,
                                      Vec1<float> &strain,
                                      const Material &mat,
                                      float tension_N = 0.0f,
                                      float tip_angle_rad = 0.0f) {
    size_t n = std::min(steps, results.count);
    if (n == 0) return;
    if (strain.capacity < n) {
        strain.reserve(n);
    }

    float h = (n > 1) ? (bs.length / (float)(n - 1)) : 0.0f;
    float x = 0.0f;

    for (size_t i = 0; i < n; ++i) {
        float d_o = get_dia(bs, x);
        float I = get_Inertia(bs, x);
        float A = get_Area(bs, x);
        float theta_x = results.items[i].items[1];
        float m_x = results.items[i].items[2];

        float axial_N = std::fabs(tension_N) * std::max(0.0f, std::cos(tip_angle_rad - theta_x));
        double sigma_req = (std::fabs((double)m_x) * (0.5 * (double)d_o)) / (double)I +
                           ((double)axial_N / (double)A);
        strain.items[i] = solve_strain_from_stress(sigma_req, mat);
        x += h;
    }
    strain.count = n;
}

inline void calculate_strain(size_t steps,
                             const BendStiffner &bs,
                             const Vec1<Vec1<float>> &results,
                             Vec1<float> &strain) {
    Material default_mat;
    calculate_strain_material(steps, bs, results, strain, default_mat, 0.0f, 0.0f);
}
