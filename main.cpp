#include <stdio.h>
#include <chrono>
#include <iostream>
#include "bs_engine_t.h"

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    // make sure we have a bend stiffner defined
    BendStiffner bs;
    bs.length = 3.0;
    bs.root_dia = 2.0;
    bs.tip_dia = 1.0;
    bs.inner_dia = 0.5;
    
    float target_angle = 0.2; // 0.2 radians for the target for now
    
    // where I am storing the results
    Vec1<Vec1<float>> results(DISCRETIZATION);
    Vec1<float> strain(DISCRETIZATION);
    
    // boundary conditions
    float y0 = 0.0;
    float theta0 = 0.0;
    float ml = 0.0;
    
    // inital guesses of moment and shear
    float guessed_m0 = 1000;
    float guessed_v0 = 1000;
    
    float angle = 0.01; // initial applied angle
    Vec2 m0_v0_guessed;
    
    while (target_angle > angle)
    {
        Vec1<Vec1<float>> temp_results(DISCRETIZATION);
        
        m0_v0_guessed = solve_bvp(bs, DISCRETIZATION, y0, theta0, ml, angle, guessed_m0, guessed_v0, strain);
        float m0 = m0_v0_guessed.x;
        float v0 = m0_v0_guessed.y;
        
        Vec1<float> y(NUMBER_OF_STATES);
        y.append(y0);
        y.append(theta0);
        y.append(m0);
        y.append(v0);
        
        float x = 0.0;
        float h = bs.length / (float) DISCRETIZATION;
        
        for (size_t i = 0; i < DISCRETIZATION; i++)
        {
            y = RK4(x, y, h, bs, strain.items[i]);
            temp_results.items[i] = y;
            x += h;
        }
        
        calculate_strain(DISCRETIZATION, bs, temp_results, strain);
        results = temp_results;
         
        angle += 0.01; // increment the applied angle
    };
    
    auto end = std::chrono::high_resolution_clock::now();
    
    printf("Strain:\n");
    for (size_t i = 0; i < DISCRETIZATION; i++)
    {
        printf("%f\n", strain.items[i]);
    };
    
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "Execution time: " << elapsed.count() << " ms" << std::endl;
    
    return 0;
}