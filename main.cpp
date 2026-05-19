#include "imgui.h"
#include "imgui_impl_sdl3.h"
#include "imgui_impl_sdlrenderer3.h"
#include "implot.h"
#include <stdio.h>
#include <SDL3/SDL.h>



#include <stdio.h>
#include <chrono>
#include <iostream>
#include "bs_engine_t.h"

#ifdef __EMSCRIPTEN__
#include "../libs/emscripten/emscripten_mainloop_stub.h"
#endif


bool solve(float target_angle, BendStiffner &bs, Vec1<Vec1<float>> &results, Vec1<float> &strain)
{
    auto start = std::chrono::high_resolution_clock::now();
    // make sure we have a bend stiffner defined
    //BendStiffner bs;
    //bs.length = 3.0;
    //bs.root_dia = 2.0;
    //bs.tip_dia = 1.0;
    //bs.inner_dia = 0.5;
    
    //float target_angle = 0.2; // 0.2 radians for the target for now
    
    // where I am storing the results
    //Vec1<Vec1<float>> results(DISCRETIZATION);
    //Vec1<float> strain(DISCRETIZATION);
    
    // boundary conditions, always the same for this software
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
        //for (int i = 0; i < 20; i++)
        //{
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
        //}
        
        angle += 0.01; // increment the applied angle
    };
    
    auto end = std::chrono::high_resolution_clock::now();
    
    printf("Strain:\n");
    for (size_t i = 0; i < DISCRETIZATION; i++)
    {
        printf("%f\n", strain.items[i]);
    };
    
    // print the moment
    printf("Bending Moment:\n");
    for (size_t i = 0; i < DISCRETIZATION; i++)
    {
        printf("%f\n", results.items[i].items[2]);
    }
    
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "Execution time: " << elapsed.count() << " ms" << std::endl;
    
    return true;
}



int main ()
{
    // setup SDL
    if (!SDL_Init(SDL_INIT_VIDEO | SDL_INIT_GAMEPAD))
    {
        printf("Error: SDL_Init(): %s\n", SDL_GetError());
        return 1;
    }
    
    // Create window with SDL_Renderer graphics context
    float main_scale = SDL_GetDisplayContentScale(SDL_GetPrimaryDisplay());
    SDL_WindowFlags window_flags = SDL_WINDOW_RESIZABLE | SDL_WINDOW_HIDDEN | SDL_WINDOW_HIGH_PIXEL_DENSITY;
    float width = 800;
    float height = 640;
    SDL_Window* window = SDL_CreateWindow("Bend Stiffner Calculator", (int)(width * main_scale), (int)(height * main_scale), window_flags);
    if (window == nullptr)
    {
        printf("Error: SDL_CreateWindow(): %s\n", SDL_GetError());
        return 1;
    }
    SDL_Renderer* renderer = SDL_CreateRenderer(window, nullptr);
    SDL_SetRenderVSync(renderer, 1);
    if (renderer == nullptr)
    {
        SDL_Log("Error: SDL_CreateRenderer(): %s\n", SDL_GetError());
        return 1;
    }
    SDL_SetWindowPosition(window, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED);
    SDL_ShowWindow(window);
    
    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    
    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    
    // Setup scaling
    ImGuiStyle& style = ImGui::GetStyle();
    style.ScaleAllSizes(main_scale);              // Bake a fixed style scale. (until we have a solution for dynamic style scaling, changing this requires resetting Style + calling this again)
    style.FontScaleDpi = main_scale * 1.0;        // Set initial font scale. (in docking branch: using io.ConfigDpiScaleFonts=true automatically overrides this for every window depending on the current monitor)
    
    // Setup Platform/Renderer backends
    ImGui_ImplSDL3_InitForSDLRenderer(window, renderer);
    ImGui_ImplSDLRenderer3_Init(renderer);
    
    
    // ###########################################################################################
    // State Data
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
    // Bend stiffner problem declarations
    BendStiffner bs = { 0 };
    float target_angle_deg = 0.0;
    float target_angle = 0.0;
    
    
    // variable to store results, both initalized to zeros
    // TODOS: need to guarantee the zero initialization
    Vec1<Vec1<float>> results(DISCRETIZATION);
    Vec1<float> strain(DISCRETIZATION);
    // this is hugly code but for now is good enough
    for (int i = 0; i < DISCRETIZATION; ++i)
    {
        strain.append(0);
    }
    
    // length array, this is basically for plotting
    float length[DISCRETIZATION];
    
    
    // Flags about state
    bool solved = false;
    bool extract_results = true;
    bool show_input_window = false;
    bool data_pasted = false;
    bool plot_bs = false;
    
    //Vec1<float> data1(256);
    //Vec1<float> data2(256);
    
    // ###########################################################################################
    
    // Main loop
    bool done = false;
    #ifdef __EMSCRIPTEN__
    io.IniFilename = nullptr;
    EMSCRIPTEN_MAINLOOP_BEGIN
    #else
    while (!done)
        #endif
    {
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            ImGui_ImplSDL3_ProcessEvent(&event);
            if (event.type == SDL_EVENT_QUIT)
                done = true;
            if (event.type == SDL_EVENT_WINDOW_CLOSE_REQUESTED && event.window.windowID == SDL_GetWindowID(window))
                done = true;
        }
        
        if (SDL_GetWindowFlags(window) & SDL_WINDOW_MINIMIZED)
        {
            SDL_Delay(10);
            continue;
        }
        
        // Start the Dear ImGui frame
        ImGui_ImplSDLRenderer3_NewFrame();
        ImGui_ImplSDL3_NewFrame();
        ImGui::NewFrame();
        
        static float f = 0.0f;
        static int counter = 0;
        
        ImGuiIO& io = ImGui::GetIO();
        
        // Force the next window to cover the entire SDL3 window
        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(io.DisplaySize);
        
        // Use flags to remove the title bar, resizing handles, and movement
        ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoDecoration |
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoSavedSettings;
        
        ImGui::Begin("Bend Stiffner Calculator", nullptr, window_flags);   // Create a window called "Hello, world!" and append into it.
        // ##########################################################################
        // here it goes my window button stuff, and my back end
        ImGui::Text("Please insert the dimensions of Bend Stiffner to be studied");
        ImGui::InputFloat("Length, in meters", &bs.length);
        ImGui::InputFloat("Root Diameter, in meters", &bs.root_dia);
        ImGui::InputFloat("Tip Diameter, in meters", &bs.tip_dia);
        ImGui::InputFloat("Inner Diameter, in meters", &bs.inner_dia);
        
        ImGui::Text("Please insert here the desired BS tip angle");
        ImGui::InputFloat("Target Angle, deg", &target_angle_deg);
        
        //if (ImGui::Button("Input window", ImVec2 { 0, 0 })) show_input_window = true;
        //ImGui::SameLine();
        
        target_angle = (3.14159265359 / 180.0) * target_angle_deg;
        if (ImGui::Button("Calculate Strain", ImVec2 { 0, 0 }))
        {
            solved = solve(target_angle, bs, results, strain);
            extract_results = true;
            //  this is here because if the lenght changes we need to update the length array
            float step = bs.length / DISCRETIZATION;
            for (int i = 0; i < DISCRETIZATION; i++)
            {
                length[i] = (float) i * step;
            };
        }
        ImGui::SameLine();
        if (ImGui::Button("Hide/show Bend Stiffner Profile", ImVec2 { 0, 0 })) plot_bs = !plot_bs;
        
        //ImGui::SameLine();
        //ImGui::Button("Calculate all input pairs", ImVec2 { 0, 0 })
        
        // display bend stiffner dimensions
        float x_data[5] = { 0, 0, bs.length, bs.length, 0 };
        float y_data[5] = { -bs.root_dia / 2, bs.root_dia / 2, bs.tip_dia / 2, -bs.tip_dia / 2 , -bs.root_dia / 2 };
        
        if (plot_bs)
        {
            ImPlot::BeginPlot("Bend Stiffner Profile", ImVec2(-1, 0.4 * height*main_scale));
            ImPlot::SetupAxes("Length, m", "Diameter, m");
            ImPlot::SetupAxesLimits(-0.1*bs.length, bs.length * 1.1, -1.1*bs.root_dia / 2, 1.1 * bs.root_dia / 2);
            // need more configuration then this, some axis scale and some
            ImPlotSpec spec;
            spec.LineWeight = 5.f;
            spec.Flags = ImPlotItemFlags_NoLegend;
            ImPlot::PlotLine("Bend Stiffner Plot", x_data, y_data, 5, spec);
            ImPlot::EndPlot();
        };
        
        //if (solved && extract_results)
        //{
        //    length_counter = 0;
        //    for (int i = 0; i < DISCRETIZATION; i++)
        //    {
        //        fstrain[i] = (float) 100.0 * strain[i];
        //        //printf("%f\n", fstrain[i]);
        //        length[i] = length_counter;
        //        length_counter += bs.length / DISCRETIZATION;
        //    }
        //    extract_results = false;
        //}
        
        if (solved && extract_results)
        {
            ImPlot::BeginPlot("Bend Stiffner Strain", ImVec2(-1, 0.4 * height*main_scale));
            ImPlot::SetupAxes("Length, m", "Strain, -");
            //ImPlot::SetupAxesLimits(-0.1*bs.length, bs.length * 1.1);
            // need more configuration then this, some axis scale and some
            ImPlotSpec spec2;
            spec2.LineWeight = 5.f;
            spec2.Flags = ImPlotItemFlags_NoLegend;
            ImPlot::PlotScatter("Bend Stiffner Plot", length, strain.items, DISCRETIZATION, spec2);
            ImPlot::EndPlot();
        }
        
        //ImGui::SetNextWindowPos(ImVec2(0, 0));
        //if (show_input_window)
        //{
        //    ImGui::Begin("Input Window", &show_input_window);
        //    if (ImGui::IsWindowFocused() && ImGui::GetIO().KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_V)) {
        //        data1.clear();
        //        data2.clear();
        //        const char* clipboardData = ImGui::GetClipboardText();
        //        if (clipboardData) {
        
        //            ParseTwoColumns(data1, data2, clipboardData);
        //        }
        //        data_pasted = true;
        //    }
        
        
        //    ImGui::Text("Copy and past the angle, tension pairs from excel");
        //    if (ImGui::BeginTable("InputTable", 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg))
        //    {
        
        //        // Set up headers
        //        ImGui::TableSetupColumn("Angle, deg");
        //        ImGui::TableSetupColumn("Tension, kN");
        //        ImGui::TableHeadersRow();
        
        //        if (data_pasted)
        //        {
        //            for (int i = 0; i < data1.size(); ++i)
        //            {
        //                ImGui::PushID(i);
        //                ImGui::TableNextRow();
        //                ImGui::TableNextColumn();
        //                ImGui::InputFloat("##my_hidden_id", &data1[i]);
        //                ImGui::PopID();
        //                ImGui::PushID(i + 10000);
        //                ImGui::TableNextColumn();
        //                ImGui::InputFloat("##my_hidden_id", &data2[i]);
        //                ImGui::PopID();
        //            }
        //        }
        
        //        ImGui::EndTable();
        //    }
        //    if (ImGui::Button("Close Me"))
        //        show_input_window = false;
        //    ImGui::End();
        //};
        
        //if (data1.size() > 0)
        //{
        //    auto it = std::max_element(data1.begin(), data1.end());
        //    target_angle_deg = *it;
        //}
        
        // ##########################################################################
        ImGui::End();
        
        // Rendering
        ImGui::Render();
        SDL_SetRenderScale(renderer, io.DisplayFramebufferScale.x, io.DisplayFramebufferScale.y);
        SDL_SetRenderDrawColorFloat(renderer, clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        SDL_RenderClear(renderer);
        ImGui_ImplSDLRenderer3_RenderDrawData(ImGui::GetDrawData(), renderer);
        SDL_RenderPresent(renderer);
    }
    #ifdef __EMSCRIPTEN__
    EMSCRIPTEN_MAINLOOP_END;
    #endif
    printf("diamter of bend stiffner: %f", bs.root_dia);
    
    // Cleanup
    // [If using SDL_MAIN_USE_CALLBACKS: all code below would likely be your SDL_AppQuit() function]
    ImGui_ImplSDLRenderer3_Shutdown();
    ImGui_ImplSDL3_Shutdown();
    ImPlot::DestroyContext(ImPlot::GetCurrentContext());
    ImGui::DestroyContext();
    
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    
    return 0;
}