#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

#include "imgui.h"
#include "imgui_impl_sdl3.h"
#include "imgui_impl_sdlrenderer3.h"
#include "implot.h"
#include <stdio.h>
#include <SDL3/SDL.h>

#include <chrono>
#include <iostream>
#include <string>
#include <cstring>
#include "bs_engine_t.h"
#include "bs_icon.h"

#ifdef __EMSCRIPTEN__
#include "../libs/emscripten/emscripten_mainloop_stub.h"
#endif

// Two-column clipboard parser (tab-separated rows from Excel/spreadsheets)
void ParseTwoColumns(Vec1<float> &data1, Vec1<float> &data2, const char* clipboard_data) {
    if (!clipboard_data || std::strlen(clipboard_data) == 0) return;

    size_t length = std::strlen(clipboard_data);
    size_t last_pos = 0;

    int column_index = 0;
    float first_col_value = 0.0f;

    for (size_t i = 0; i <= length; ++i) {
        char c = (i < length) ? clipboard_data[i] : '\n';

        if (c == '\t' || c == '\n' || c == '\r') {
            if (c == '\t') {
                if (column_index == 0) {
                    char* endptr = nullptr;
                    float val = std::strtof(clipboard_data + last_pos, &endptr);
                    if (endptr != clipboard_data + last_pos) {
                        first_col_value = val;
                        column_index = 1;
                    }
                }
                last_pos = i + 1;
            } else {
                if (c == '\r' && (i + 1 < length) && clipboard_data[i + 1] == '\n') {
                    i++;
                }

                if (column_index == 1) {
                    char* endptr = nullptr;
                    float second_col_value = std::strtof(clipboard_data + last_pos, &endptr);
                    if (endptr != clipboard_data + last_pos) {
                        data1.append(first_col_value);
                        data2.append(second_col_value);
                    }
                } else if (column_index == 0 && (i > last_pos)) {
                    char* endptr = nullptr;
                    float fallback_value = std::strtof(clipboard_data + last_pos, &endptr);
                    if (endptr != clipboard_data + last_pos) {
                        data1.append(fallback_value);
                        data2.append(0.0f);
                    }
                }

                column_index = 0;
                last_pos = i + 1;
            }
        }
    }
}

bool solve(float target_angle,
           BendStiffner &bs,
           Vec1<Vec1<float>> &results,
           Vec1<float> &strain,
           const Material &mat = Material(),
           float tension_kN = 0.0f,
           double* out_elapsed_ms = nullptr,
           bool verbose = true)
{
    if (!is_valid_geometry(bs)) {
        return false;
    }

    auto start = std::chrono::high_resolution_clock::now();

    float y0 = 0.0f;
    float theta0 = 0.0f;
    float ml = 0.0f;
    float tension_N = tension_kN * 1000.0f;

    // Number of integration intervals between DISCRETIZATION nodes [0 .. L]
    const size_t intervals = (DISCRETIZATION > 1) ? (DISCRETIZATION - 1) : 1;
    Vec2 m0_v0 = solve_bvp_material(bs, intervals, y0, theta0, ml, target_angle, mat, tension_N);
    float m0 = m0_v0.x;
    float v0 = m0_v0.y;

    results.clear();
    if (results.capacity < DISCRETIZATION) {
        results.reserve(DISCRETIZATION);
    }

    State4 state = { { y0, theta0, m0, v0 } };
    float x = 0.0f;
    float h = bs.length / (float)intervals;

    for (size_t i = 0; i < DISCRETIZATION; ++i) {
        Vec1<float> node_state(NUMBER_OF_STATES);
        for (int s = 0; s < NUMBER_OF_STATES; ++s) {
            node_state.append(state.items[s]);
        }
        results.append(node_state);

        if (i + 1 < DISCRETIZATION) {
            state = RK4_fast(x, state, h, bs, mat, tension_N, target_angle);
            x += h;
        }
    }

    calculate_strain_material(DISCRETIZATION, bs, results, strain, mat, tension_N, target_angle);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    if (out_elapsed_ms) {
        *out_elapsed_ms = elapsed.count();
    }

    if (verbose) {
        printf("Solved (angle = %.2f deg, tension = %.1f kN) in %.3f ms | M0 = %.2f N.m | V0 = %.2f N | Max Strain = %.5f\n",
               target_angle * 180.0f / 3.14159265359f,
               tension_kN,
               elapsed.count(),
               m0,
               v0,
               strain.max_pos());
    }

    return true;
}

void solve_pairs(Vec1<float> &strain_maxs,
                 const Vec1<float> &angles_deg,
                 const Vec1<float> &tensions_kN,
                 BendStiffner &bs,
                 Vec1<Vec1<float>> &results,
                 Vec1<float> &strain,
                 const Material &mat = Material(),
                 double* out_total_ms = nullptr)
{
    auto start = std::chrono::high_resolution_clock::now();
    strain_maxs.clear();

    for (size_t i = 0; i < angles_deg.count; ++i) {
        float angle_rad = angles_deg.items[i] * (3.14159265358979323846f / 180.0f);
        float tension_kN = (i < tensions_kN.count) ? tensions_kN.items[i] : 0.0f;
        solve(angle_rad, bs, results, strain, mat, tension_kN, nullptr, false);
        strain_maxs.append(strain.max_pos());
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    if (out_total_ms) {
        *out_total_ms = elapsed.count();
    }
}

int main()
{
#ifdef _WIN32
    // Enable Per-Monitor V2 DPI awareness so Windows DWM never bitmap-stretches the window
    SetProcessDpiAwarenessContext(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2);
#endif

    if (!SDL_Init(SDL_INIT_VIDEO | SDL_INIT_GAMEPAD))
    {
        printf("Error: SDL_Init(): %s\n", SDL_GetError());
        return 1;
    }

    float main_scale = SDL_GetDisplayContentScale(SDL_GetPrimaryDisplay());
    if (main_scale <= 0.0f) {
        main_scale = 1.0f;
    }
    SDL_WindowFlags sdl_window_flags = SDL_WINDOW_RESIZABLE | SDL_WINDOW_HIDDEN | SDL_WINDOW_HIGH_PIXEL_DENSITY;
    float width = 1280.0f;
    float height = 760.0f;
    SDL_Window* window = SDL_CreateWindow("Bend Stiffener Calculator", (int)(width * main_scale), (int)(height * main_scale), sdl_window_flags);
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
    if (SDL_Surface* icon_surface = SDL_CreateSurfaceFrom(
            BS_ICON_WIDTH, BS_ICON_HEIGHT, SDL_PIXELFORMAT_RGBA32,
            (void*)BS_ICON_RGBA, BS_ICON_WIDTH * 4))
    {
        SDL_SetWindowIcon(window, icon_surface);
        SDL_DestroySurface(icon_surface);
    }
    SDL_ShowWindow(window);

    float window_scale = SDL_GetWindowDisplayScale(window);
    if (window_scale > 0.0f) {
        main_scale = window_scale;
    }

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;

    ImGui::StyleColorsDark();

    ImGuiStyle& style = ImGui::GetStyle();
    style.ScaleAllSizes(main_scale);
    style.AntiAliasedLines = true;
    style.AntiAliasedLinesUseTex = true;
    style.AntiAliasedFill = true;

    // Load a crisp TrueType font at an exact integer pixel size with oversampling and pixel snapping
    ImFontConfig font_cfg;
    font_cfg.OversampleH = 3;
    font_cfg.OversampleV = 2;
    font_cfg.PixelSnapH = true;
    font_cfg.Flags |= ImFontFlags_NoLoadError;
    float base_font_size = std::round(17.0f * main_scale);

    ImFont* loaded_font = io.Fonts->AddFontFromFileTTF("C:\\Windows\\Fonts\\segoeui.ttf", base_font_size, &font_cfg);
    if (!loaded_font) {
        loaded_font = io.Fonts->AddFontFromFileTTF("imgui/misc/fonts/Roboto-Medium.ttf", base_font_size, &font_cfg);
    }
    if (!loaded_font) {
        loaded_font = io.Fonts->AddFontFromFileTTF("..\\imgui\\misc\\fonts\\Roboto-Medium.ttf", base_font_size, &font_cfg);
    }
    if (!loaded_font) {
        loaded_font = io.Fonts->AddFontDefaultVector(&font_cfg);
    }
    style.FontSizeBase = base_font_size;
    style.FontScaleMain = 1.0f;
    style.FontScaleDpi = 1.0f;

    ImGui_ImplSDL3_InitForSDLRenderer(window, renderer);
    ImGui_ImplSDLRenderer3_Init(renderer);

    // ###########################################################################################
    // Problem State Data (initialized with sensible default Bend Stiffener dimensions)
    ImVec4 clear_color = ImVec4(0.10f, 0.12f, 0.15f, 1.0f);
    BendStiffner bs = { 3.0f, 0.60f, 0.20f, 0.12f };
    Material material;
    float target_angle_deg = 20.0f;
    float target_tension_kN = 0.0f;
    float target_angle = 0.0f;

    Vec1<Vec1<float>> results(DISCRETIZATION);
    Vec1<float> strain(DISCRETIZATION);
    for (int i = 0; i < DISCRETIZATION; ++i)
    {
        strain.append(0.0f);
    }

    Vec1<float> strain_maxs(256);

    float length[DISCRETIZATION] = { 0.0f };
    float deflection_y[DISCRETIZATION] = { 0.0f };
    float moment_kNm[DISCRETIZATION] = { 0.0f };

    bool solved = false;
    bool extract_results = true;
    bool show_input_window = false;
    bool show_material_window = false;
    bool data_pasted = false;
    bool plot_bs = true;
    bool plot_pairs = false;
    bool plot_moment = false;
    double last_solve_ms = 0.0;
    double last_pairs_ms = 0.0;
    std::string status_message = "Ready.";

    Vec1<float> data1(256); // Angle, deg
    Vec1<float> data2(256); // Tension, kN

    // ###########################################################################################
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
            if (event.type == SDL_EVENT_WINDOW_DISPLAY_SCALE_CHANGED && event.window.windowID == SDL_GetWindowID(window))
            {
                float new_scale = SDL_GetWindowDisplayScale(window);
                if (new_scale > 0.0f) {
                    main_scale = new_scale;
                    style.FontSizeBase = std::round(17.0f * main_scale);
                }
            }
        }

        if (SDL_GetWindowFlags(window) & SDL_WINDOW_MINIMIZED)
        {
            SDL_Delay(10);
            continue;
        }

        ImGui_ImplSDLRenderer3_NewFrame();
        ImGui_ImplSDL3_NewFrame();
        ImGui::NewFrame();

        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(io.DisplaySize);

        ImGuiWindowFlags main_win_flags = ImGuiWindowFlags_NoDecoration |
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoSavedSettings |
            ImGuiWindowFlags_NoBringToFrontOnFocus;

        ImGui::Begin("Bend Stiffener Calculator", nullptr, main_win_flags);

        // --- Section 1: Geometry & Loading Inputs ---
        ImGui::SeparatorText("Bend Stiffener Geometry & Loading");
        ImGui::PushItemWidth(220.0f * main_scale);
        ImGui::InputFloat("Length [m]", &bs.length, 0.1f, 0.5f, "%.3f");
        ImGui::SameLine(420.0f * main_scale);
        ImGui::InputFloat("Root Diameter [m]", &bs.root_dia, 0.02f, 0.1f, "%.3f");

        ImGui::InputFloat("Tip Diameter [m]", &bs.tip_dia, 0.02f, 0.1f, "%.3f");
        ImGui::SameLine(420.0f * main_scale);
        ImGui::InputFloat("Inner Diameter [m]", &bs.inner_dia, 0.01f, 0.05f, "%.3f");

        ImGui::InputFloat("Target Tip Angle [deg]", &target_angle_deg, 1.0f, 5.0f, "%.2f");
        ImGui::SameLine(420.0f * main_scale);
        ImGui::InputFloat("Target Tension [kN]", &target_tension_kN, 5.0f, 25.0f, "%.2f");
        ImGui::PopItemWidth();

        ImGui::Checkbox("Non-Linear Material (Polyurethane)", &material.is_nonlinear);
        ImGui::SameLine();
        if (ImGui::Button("Material Properties...")) {
            show_material_window = true;
        }
        ImGui::SameLine();
        if (ImGui::Button("Angle/Tension Pairs Table...")) {
            show_input_window = true;
        }

        // --- Section 2: Actions ---
        ImGui::SeparatorText("Solver & Visualization Controls");
        target_angle = (3.14159265358979323846f / 180.0f) * target_angle_deg;

        if (ImGui::Button("Calculate Strain"))
        {
            if (!is_valid_geometry(bs)) {
                status_message = "Error: Invalid geometry! Ensure Length > 0 and Root/Tip Diameter > Inner Diameter >= 0.";
                solved = false;
            } else if (std::fabs(target_angle_deg) >= 89.5f) {
                status_message = "Error: Target Angle must be between -89 and +89 degrees.";
                solved = false;
            } else {
                solved = solve(target_angle, bs, results, strain, material, target_tension_kN, &last_solve_ms, true);
                extract_results = solved;
                if (solved) {
                    float step = (DISCRETIZATION > 1) ? (bs.length / (float)(DISCRETIZATION - 1)) : 0.0f;
                    for (int i = 0; i < DISCRETIZATION; i++)
                    {
                        length[i] = (float)i * step;
                        deflection_y[i] = results.items[i].items[0];
                        moment_kNm[i] = results.items[i].items[2] / 1000.0f;
                    }
                    char buf[256];
                    snprintf(buf, sizeof(buf),
                             "Solved in %.3f ms | Max Strain: %.4f (%.2f%%) | Root Moment M0: %.2f kN.m | Tip Deflection: %.3f m",
                             last_solve_ms, strain.max_pos(), strain.max_pos() * 100.0f,
                             results.items[0].items[2] / 1000.0f,
                             results.items[DISCRETIZATION - 1].items[0]);
                    status_message = buf;
                }
            }
        }

        ImGui::SameLine();
        if (ImGui::Button("Calculate All Input Pairs"))
        {
            if (!is_valid_geometry(bs)) {
                status_message = "Error: Invalid geometry!";
            } else if (data1.count == 0) {
                status_message = "No Angle/Tension pairs loaded. Open 'Angle/Tension Pairs Table...' and paste from Excel (Ctrl+V).";
            } else {
                solve_pairs(strain_maxs, data1, data2, bs, results, strain, material, &last_pairs_ms);
                plot_pairs = true;
                char buf[256];
                snprintf(buf, sizeof(buf),
                         "Solved %zu Angle/Tension pairs in %.2f ms (%.3f ms/pair) | Peak Strain across pairs: %.4f (%.2f%%)",
                         data1.count, last_pairs_ms, last_pairs_ms / (double)data1.count,
                         strain_maxs.max_pos(), strain_maxs.max_pos() * 100.0f);
                status_message = buf;
            }
        }

        ImGui::SameLine();
        ImGui::Checkbox("Show BS Profile", &plot_bs);
        ImGui::SameLine();
        ImGui::Checkbox("Show Bending Moment", &plot_moment);
        ImGui::SameLine();
        ImGui::Checkbox("Show Pairs Strain Plot", &plot_pairs);

        ImGui::TextColored(ImVec4(0.4f, 0.9f, 0.6f, 1.0f), "Status: %s", status_message.c_str());

        // --- Section 3: Plots ---
        float x_data[5] = { 0.0f, 0.0f, bs.length, bs.length, 0.0f };
        float y_data[5] = { -bs.root_dia / 2.0f, bs.root_dia / 2.0f, bs.tip_dia / 2.0f, -bs.tip_dia / 2.0f, -bs.root_dia / 2.0f };
        float inner_x[2] = { 0.0f, bs.length };
        float inner_y_top[2] = { bs.inner_dia / 2.0f, bs.inner_dia / 2.0f };
        float inner_y_bot[2] = { -bs.inner_dia / 2.0f, -bs.inner_dia / 2.0f };

        float plot_height = std::max(180.0f, 0.28f * io.DisplaySize.y);

        if (plot_bs)
        {
            if (ImPlot::BeginPlot("Bend Stiffener Geometry Profile", ImVec2(-1, plot_height))) {
                ImPlot::SetupAxes("Length [m]", "Radius [m]");
                ImPlot::SetupAxesLimits(-0.05 * bs.length, bs.length * 1.05, -0.65 * bs.root_dia, 0.65 * bs.root_dia, ImGuiCond_Always);
                ImPlotSpec spec;
                spec.LineWeight = 2.0f * main_scale;
                ImPlot::PlotLine("Outer Cone", x_data, y_data, 5, spec);
                ImPlotSpec inner_spec;
                inner_spec.LineWeight = 1.5f * main_scale;
                ImPlot::PlotLine("Inner Bore (+)", inner_x, inner_y_top, 2, inner_spec);
                ImPlot::PlotLine("Inner Bore (-)", inner_x, inner_y_bot, 2, inner_spec);
                ImPlot::EndPlot();
            }
        }

        if (solved && extract_results)
        {
            if (ImPlot::BeginPlot("Bend Stiffener Strain Distribution", ImVec2(-1, plot_height))) {
                ImPlot::SetupAxes("Length from Root [m]", "Outer Fiber Strain [-]");
                ImPlotSpec spec_line;
                spec_line.LineWeight = 2.0f * main_scale;
                spec_line.Marker = ImPlotMarker_Circle;
                spec_line.MarkerSize = 3.5f * main_scale;
                ImPlot::PlotLine("Strain Profile", length, strain.items, DISCRETIZATION, spec_line);
                ImPlot::EndPlot();
            }

            if (plot_moment)
            {
                if (ImPlot::BeginPlot("Bending Moment & Deflection Profile", ImVec2(-1, plot_height))) {
                    ImPlot::SetupAxes("Length from Root [m]", "Bending Moment [kN.m]");
                    ImPlotSpec spec_m;
                    spec_m.LineWeight = 2.0f * main_scale;
                    ImPlot::PlotLine("Bending Moment [kN.m]", length, moment_kNm, DISCRETIZATION, spec_m);
                    ImPlot::EndPlot();
                }
            }
        }

        if (plot_pairs && strain_maxs.count > 0)
        {
            if (ImPlot::BeginPlot("Max Strain for Tension / Angle Pairs", ImVec2(-1, plot_height))) {
                ImPlot::SetupAxes("Tip Angle [deg]", "Max Strain [-]");
                ImPlotSpec spec_p;
                spec_p.LineWeight = 1.5f * main_scale;
                spec_p.Marker = ImPlotMarker_Circle;
                spec_p.MarkerSize = 4.0f * main_scale;
                size_t plot_count = std::min(data1.count, strain_maxs.count);
                ImPlot::PlotScatter("Max Strain per Pair", data1.items, strain_maxs.items, (int)plot_count, spec_p);
                ImPlot::EndPlot();
            }
        }

        // --- Modal / Floating Window 1: Angle & Tension Pairs Input ---
        if (show_input_window)
        {
            ImGui::SetNextWindowSize(ImVec2(480 * main_scale, 420 * main_scale), ImGuiCond_FirstUseEver);
            ImGui::Begin("Angle & Tension Pairs Input", &show_input_window);
            if (ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows) &&
                ImGui::GetIO().KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_V))
            {
                data1.clear();
                data2.clear();
                const char* clipboardData = ImGui::GetClipboardText();
                if (clipboardData) {
                    ParseTwoColumns(data1, data2, clipboardData);
                }
                data_pasted = (data1.count > 0);
                if (data_pasted) {
                    target_angle_deg = data1.max_pos();
                }
            }

            ImGui::TextWrapped("Focus this window and press Ctrl+V to paste 2 columns (Angle [deg], Tension [kN]) from Excel.");
            if (ImGui::Button("Add Row")) {
                data1.append(10.0f);
                data2.append(0.0f);
                data_pasted = true;
            }
            ImGui::SameLine();
            if (ImGui::Button("Clear Table")) {
                data1.clear();
                data2.clear();
                strain_maxs.clear();
                data_pasted = false;
            }
            ImGui::SameLine();
            if (ImGui::Button("Close")) {
                show_input_window = false;
            }

            if (ImGui::BeginTable("InputTable", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY, ImVec2(0, 300 * main_scale)))
            {
                ImGui::TableSetupColumn("Angle [deg]");
                ImGui::TableSetupColumn("Tension [kN]");
                ImGui::TableSetupColumn("Max Strain [-]");
                ImGui::TableHeadersRow();

                size_t lowest_count = std::min(data1.count, data2.count);
                for (size_t i = 0; i < lowest_count; ++i)
                {
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::PushID((int)i);
                    ImGui::InputFloat("##angle_col", &data1.items[i]);
                    ImGui::PopID();

                    ImGui::TableNextColumn();
                    ImGui::PushID((int)i + 10000);
                    ImGui::InputFloat("##tension_col", &data2.items[i]);
                    ImGui::PopID();

                    ImGui::TableNextColumn();
                    if (i < strain_maxs.count) {
                        ImGui::Text("%.5f (%.2f%%)", strain_maxs.items[i], strain_maxs.items[i] * 100.0f);
                    } else {
                        ImGui::TextDisabled("-");
                    }
                }
                ImGui::EndTable();
            }
            ImGui::End();
        }

        // --- Modal / Floating Window 2: Material Properties Configuration ---
        if (show_material_window)
        {
            ImGui::SetNextWindowSize(ImVec2(540 * main_scale, 480 * main_scale), ImGuiCond_FirstUseEver);
            ImGui::Begin("Material Properties", &show_material_window);

            if (ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows) &&
                ImGui::GetIO().KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_V))
            {
                const char* clipboardData = ImGui::GetClipboardText();
                if (clipboardData) {
                    Vec1<float> new_strain(32);
                    Vec1<float> new_emod(32);
                    ParseTwoColumns(new_strain, new_emod, clipboardData);
                    if (new_strain.count >= 2 && new_strain.count == new_emod.count) {
                        material.strain_table = new_strain;
                        material.emod_table_kpa = new_emod;
                        material.linear_E_kpa = new_emod.items[0];
                        status_message = "Pasted custom non-linear material table from clipboard.";
                    }
                }
            }

            ImGui::Checkbox("Enable Non-Linear Material Curve", &material.is_nonlinear);
            if (!material.is_nonlinear) {
                ImGui::InputFloat("Linear Young's Modulus E [kPa]", &material.linear_E_kpa, 1000.0f, 10000.0f, "%.1f");
            } else {
                ImGui::Checkbox("Table values are Tangent Modulus d(sigma)/d(strain) [integrated for Secant E]", &material.table_is_tangent_modulus);
                ImGui::TextWrapped("Focus this window and press Ctrl+V to paste a 2-column (Strain [-], E-Modulus [kPa]) table from Excel.");

                if (ImGui::Button("Reset to Default Polyurethane")) {
                    material.reset_to_default();
                }
                ImGui::SameLine();
                if (ImGui::Button("Add Point")) {
                    float next_s = (material.strain_table.count > 0)
                        ? (material.strain_table.items[material.strain_table.count - 1] + 0.01f)
                        : 0.0f;
                    float last_e = (material.emod_table_kpa.count > 0)
                        ? material.emod_table_kpa.items[material.emod_table_kpa.count - 1]
                        : 215800.0f;
                    material.strain_table.append(next_s);
                    material.emod_table_kpa.append(last_e);
                }
                ImGui::SameLine();
                if (ImGui::Button("Close")) {
                    show_material_window = false;
                }

                if (ImGui::BeginTable("MaterialTable", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY, ImVec2(0, 300 * main_scale)))
                {
                    ImGui::TableSetupColumn("Strain [-]");
                    ImGui::TableSetupColumn("Table E [kPa]");
                    ImGui::TableSetupColumn("Effective Secant E [MPa]");
                    ImGui::TableHeadersRow();

                    size_t n_mat = std::min(material.strain_table.count, material.emod_table_kpa.count);
                    for (size_t i = 0; i < n_mat; ++i)
                    {
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn();
                        ImGui::PushID((int)i + 20000);
                        ImGui::InputFloat("##mat_strain", &material.strain_table.items[i], 0.0f, 0.0f, "%.4f");
                        ImGui::PopID();

                        ImGui::TableNextColumn();
                        ImGui::PushID((int)i + 30000);
                        ImGui::InputFloat("##mat_emod", &material.emod_table_kpa.items[i], 0.0f, 0.0f, "%.1f");
                        ImGui::PopID();

                        ImGui::TableNextColumn();
                        float e_sec_mpa = get_material_E(material.strain_table.items[i], material) / 1e6f;
                        ImGui::Text("%.2f MPa", e_sec_mpa);
                    }
                    ImGui::EndTable();
                }
            }
            ImGui::End();
        }

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

    // Cleanup
    ImGui_ImplSDLRenderer3_Shutdown();
    ImGui_ImplSDL3_Shutdown();
    ImPlot::DestroyContext(ImPlot::GetCurrentContext());
    ImGui::DestroyContext();

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
