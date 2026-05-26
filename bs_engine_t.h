#pragma once

#include <cmath>
#include <cassert>

#define DISCRETIZATION 50
#define NUMBER_OF_STATES 4


// TODOS: when looping throught the Vec1 object
// need to put a guard so that no loops can go over the capacity / count of each
// Vec1. This may be done in a if statement? does this add any significant overhead?

// bend stiffner struc
struct BendStiffner {
    float length;
    float root_dia;
    float tip_dia;
    float inner_dia;
};

struct Vec2 {
    float x;
    float y;
};

template<typename T>
struct Vec1 {
    T* items;
    size_t count;
    size_t capacity;
    
    Vec1(size_t size)
    {
        capacity = size;
        count = 0;
        items = (T *)malloc(capacity * sizeof(T));
    };
    
    void append(T item)
    {
        if (count >= capacity)
        {
            capacity *= 2;
            items = (T*) realloc(items, capacity * sizeof(T));
        };
        
        items[count] = item;
        count++;
    };
    
    void clear() 
    {
        memset(items, 0, capacity*sizeof(float));
        count = 0;
    };
    
    // returns the largest positive value or zero if all elements are negative
    float max_pos()
    {
        float max_element = 0.0;
        // need to make sure that the count is higher than zero
        if (count == 0.0)
        {
            return 0.0;
        } else
        {
            for (int i = 0.0; i < count; i++)
            {
                if (items[i] > max_element) max_element = items[i];
            };
            return max_element;
        };
    };
};


// get_non_linear
float get_non_linear_E(float strain) { // need to pass material object in this call
    
    const int material_data_size = 21;
    
    float strain_vector[material_data_size] = { 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20 };
    float Emod_vector[material_data_size] = { 215800, 215800, 179500, 138100, 94300, 59500, 39700, 28200, 21100, 17000, 14300, 12700, 11000, 10600, 10000, 9700, 9600, 8600, 9400, 8600, 8500 };
    
    float distance = 1000;
    int min_index = 0;
    int max_index = 0;
    
    // need to use the absolute values here
    for (int i = 0; i < sizeof(strain_vector); i++)
    {
        float prev_distance = distance;
        distance = strain_vector[i] - strain;
        if (std::fabs(distance < prev_distance))
        {
            if (distance < 0)
            {
                min_index = i;
                max_index = i + 1;
            };
            if (distance > 0)
            {
                min_index = i - 1;
                max_index = i;
            }
        };
    }
    
    // need to deal with the edge cases
    if (min_index < 0) {
        min_index = 1;
        max_index = min_index + 1;
    };
    
    if (max_index > sizeof(strain_vector))
    {
        max_index = sizeof(strain_vector);
        min_index = max_index - 1;
    }
    
    float x0 = strain_vector[min_index];
    float x1 = strain_vector[max_index];
    float y0 = Emod_vector[min_index];
    float y1 = Emod_vector[max_index];
    float result = y1 + (strain - x0) * (y0 - y1) / (x1 - x0);
    return result * 1000.0; // convert to MPa
};

float get_non_linear_E_custom(float strain, Vec1<float> &strain_input_vector, Vec1<float> &E_mod_input_vector) { // need to pass material object in this call
    float distance = 1000;
    int min_index = 0;
    int max_index = 0;
    
    assert(strain_input_vector.count == E_mod_input_vector.count);

    // need to use the absolute values here
    for (int i = 0; i < sizeof(strain_input_vector.count); i++)
    {
        float prev_distance = distance;
        distance = strain_input_vector.items[i] - strain;
        if (std::fabs(distance < prev_distance))
        {
            if (distance < 0)
            {
                min_index = i;
                max_index = i + 1;
            };
            if (distance > 0)
            {
                min_index = i - 1;
                max_index = i;
            }
        };
    }
    
    // need to deal with the edge cases
    if (min_index < 0) {
        min_index = 1;
        max_index = min_index + 1;
    };
    
    if (max_index > sizeof(strain_input_vector.count))
    {
        max_index = strain_input_vector.count;
        min_index = max_index - 1;
    }
    
    float x0 = strain_input_vector.items[min_index];
    float x1 = strain_input_vector.items[max_index];
    float y0 = E_mod_input_vector.items[min_index];
    float y1 = E_mod_input_vector.items[max_index];
    float result = y1 + (strain - x0) * (y0 - y1) / (x1 - x0);
    return result * 1000.0; // convert to MPa
};

float get_dia(const BendStiffner &dimensions, float x) {
    
    float outer_dia = dimensions.root_dia - (dimensions.root_dia - dimensions.tip_dia) *
        (x / dimensions.length);
    return outer_dia;
};


float get_Inertia(const BendStiffner& dimensions, float x) {
    
    const float pi = 3.141592653599;
    // calculate outer diameter
    // this assumes the conic shape
    float outer_dia = dimensions.root_dia - (dimensions.root_dia - dimensions.tip_dia) *
        (x / dimensions.length);
    
    // inertia is given by  I = pi/64 *(Do^4 - Di^4)
    
    float inertia = (pi / 64) * (pow(outer_dia, 4) - pow(dimensions.inner_dia, 4));
    return inertia;
};


float get_EI(const BendStiffner& dimensions, float strain, float x) {
    
    const float m_E = 215800;
    float inertia = get_Inertia(dimensions, x);
    if (strain == 0) {
        return m_E * inertia;
    } else {
        return get_non_linear_E(strain) * inertia;
    }
};


Vec1<float> equations(float x, const Vec1<float> &y, const BendStiffner bs, float strain)
{
    Vec1<float> dydx(NUMBER_OF_STATES);
    
    float EI = get_EI(bs, strain, x);
    
    dydx.append(std::tan(y.items[1]));
    dydx.append(y.items[2] / (EI * std::cos(y.items[1])));
    dydx.append(y.items[3]);
    dydx.append(0);
    
    return dydx;
};


Vec1<float> RK4(float x, const Vec1<float> &y, float h, const BendStiffner &bs, float strain)
{
    size_t n = y.count;
    
    Vec1<float> k1 = equations(x, y, bs, strain);
    Vec1<float> y_temp(n);
    
    for (int i = 0; i < n; i++)
    {
        y_temp.items[i] = y.items[i] + 0.5 * h * k1.items[i];
        y_temp.count++;
    };
    Vec1<float> k2 = equations(x + 0.5*h, y_temp, bs, strain);
    
    for (int i = 0; i < n; i++)
    {
        y_temp.items[i] = y.items[i] + 0.5 * h * k2.items[i];
        // no need for the y_temp.count++ because we are overwriting the previous data points
    };
    Vec1<float> k3 = equations(x + 0.5*h, y_temp, bs, strain);
    
    for (int i = 0; i < n; i++)
    {
        y_temp.items[i] = y.items[i] + h * k3.items[i];
    };
    Vec1<float> k4 = equations(x + h, y_temp, bs, strain);
    Vec1<float> y_next(n);
    
    for (int i = 0; i < n; i++)
    {
        y_next.items[i] = 
            y.items[i] + (h / 6.0) * (k1.items[i] + 2.0 * k2.items[i] + 2.0 * k3.items[i] + k4.items[i]);
        y_next.count++;
    };
    
    return y_next;
};


Vec1<float> shoot(const BendStiffner &bs,
                  size_t steps,
                  float y0,
                  float theta0,
                  float guessed_m0,
                  float guessed_v0,
                  const Vec1<float> &strain)
{
    float h = bs.length / (float) steps;
    float x = 0.0;
    
    Vec1<float> y(NUMBER_OF_STATES);
    y.append(y0);
    y.append(theta0);
    y.append(guessed_m0);
    y.append(guessed_v0);
    
    for (int i = 0; i < steps; i++)
    {
        float strain_f = strain.items[i];
        y = RK4(x, y, h, bs, strain_f);
        x += h;
    }
    
    return y;
};


float solve_V0 (const BendStiffner &bs, 
                size_t steps, 
                float y0,
                float theta0,
                float curr_guessed_m0,
                float thetaL,
                float guessed_v0,
                const Vec1<float> &strain)
{
    float v0 = -guessed_v0;
    float v1 = guessed_v0;
    
    Vec1<float> f0_vec(NUMBER_OF_STATES);
    f0_vec = shoot(bs, steps, y0, theta0, curr_guessed_m0, v0, strain);
    float f0 = f0_vec.items[1] - thetaL;
    
    Vec1<float> f1_vec(NUMBER_OF_STATES);
    f1_vec = shoot(bs, steps, y0, theta0, curr_guessed_m0, v1, strain);
    float f1 = f1_vec.items[1] - thetaL;
    
    Vec1<float> temp_state(NUMBER_OF_STATES);
    
    for (int i = 0; i < 25; i++)
    {
        if (std::abs(f1) < 1e-6)
        {
            return v1;
        };
        if (std::abs(f1 - f0) < 1e-12)
        {
            break;
        };
        
        float v2 = v1 - f1 * (v1 - v0) / (f1 - f0);
        v0 = v1;
        f0 = f1;
        v1 = v2;
        
        temp_state = shoot(bs, steps, y0, theta0, curr_guessed_m0, v1, strain);
        f1 = temp_state.items[1] - thetaL;
    };
    return v1;
};


Vec2 solve_bvp (const BendStiffner &bs, 
                size_t steps,
                float y0,
                float theta0,
                float targetML,
                float thetaL,
                float guessed_m0,
                float guessed_v0,
                const Vec1<float> &strain)
{
    float u0 = -guessed_m0;
    float u1 = guessed_m0;
    
    float correct_v0_for_u0 = solve_V0(bs, steps, y0, theta0, u0, thetaL, guessed_v0, strain);
    
    Vec1<float> f0_vec(NUMBER_OF_STATES);
    f0_vec = shoot(bs, steps, y0, theta0, u0, correct_v0_for_u0, strain);
    float f0 = f0_vec.items[2] - targetML;
    
    float correct_v0_for_u1 = solve_V0(bs, steps, y0, theta0, u1, thetaL, guessed_v0, strain);
    
    Vec1<float> f1_vec(NUMBER_OF_STATES);
    f1_vec = shoot(bs, steps, y0, theta0, u1, correct_v0_for_u1, strain);
    float f1 = f1_vec.items[2] - targetML;
    
    for (int i = 0; i < 25; i++)
    {
        if (std::abs(f1) < 1e-6)
        {
            return { u1, correct_v0_for_u1 };
        };
        if (std::abs(f1 - f0) < 1e-12)
        {
            break;
        };
        
        float u2 = u1 - f1 * (u1 - u0) / (f1 - f0);
        u0 = u1;
        f0 = f1;
        u1 = u2;
        guessed_v0 = correct_v0_for_u1;
        
        correct_v0_for_u1 = solve_V0(bs, steps, y0, theta0, u1, thetaL, guessed_v0, strain);
        f1_vec = shoot(bs, steps, y0, theta0, u1, correct_v0_for_u1, strain);
        f1 = f1_vec.items[2] - targetML;
    };
    
    return { u1, correct_v0_for_u1 };
};


void calculate_strain(size_t steps,
                      const BendStiffner &bs,
                      const Vec1<Vec1<float>> &results, 
                      Vec1<float> &strain)
{
    float h = bs.length / (double) steps;
    float x = 0.0;
    
    size_t count = 0.0;
    
    for (int i = 0; i < DISCRETIZATION; i++)
    {
        //if (i > strain.capacity || i > results.capacity) break;
        float EI = get_EI(bs, strain.items[i], x);
        strain.items[i] = (get_dia(bs, x) * results.items[i].items[2]) / EI;
        strain.count = count;
        count++;
        x += h;
    }
};
