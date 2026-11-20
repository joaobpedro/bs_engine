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

// Parses either 4 columns (Strain, Nominal_E, Cold_E, Hot_E) or 2 columns (Strain, Active_Temp_E) from Excel clipboard
int ParseMaterialClipboard(Material &mat, const char* clipboard_data) {
    if (!clipboard_data || std::strlen(clipboard_data) == 0) return 0;

    Vec1<float> col_s(32), col_n(32), col_c(32), col_h(32);
    size_t length = std::strlen(clipboard_data);
    size_t line_start = 0;
    int max_cols_found = 0;

    for (size_t i = 0; i <= length; ++i) {
        char c = (i < length) ? clipboard_data[i] : '\n';
        if (c == '\n' || c == '\r') {
            if (c == '\r' && (i + 1 < length) && clipboard_data[i + 1] == '\n') {
                i++;
            }
            // Parse tab-separated tokens on this line
            float vals[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
            int col_count = 0;
            size_t tok_start = line_start;
            for (size_t j = line_start; j <= i && col_count < 4; ++j) {
                char tc = (j < i && clipboard_data[j] != '\r') ? clipboard_data[j] : '\t';
                if (tc == '\t') {
                    if (j > tok_start) {
                        char* endptr = nullptr;
                        float v = std::strtof(clipboard_data + tok_start, &endptr);
                        if (endptr != clipboard_data + tok_start) {
                            vals[col_count++] = v;
                        }
                    }
                    tok_start = j + 1;
                }
            }
            if (col_count >= 4) {
                col_s.append(vals[0]);
                col_n.append(vals[1]);
                col_c.append(vals[2]);
                col_h.append(vals[3]);
                if (4 > max_cols_found) max_cols_found = 4;
            } else if (col_count >= 2) {
                col_s.append(vals[0]);
                col_n.append(vals[1]);
                if (2 > max_cols_found) max_cols_found = 2;
            }
            line_start = i + 1;
        }
    }

    if (col_s.count < 2) return 0;

    // Auto-detect percentage strain input (e.g. 0..20% pasted from Excel) and convert to fractional strain (0..0.20)
    float max_s = col_s.max_val(0.0f);
    float strain_scale = (max_s > 2.0f && max_s <= 100.0f) ? 0.01f : 1.0f;
    for (size_t k = 0; k < col_s.count; ++k) {
        col_s.items[k] = std::max(0.0f, col_s.items[k] * strain_scale);
        col_n.items[k] = std::max(1.0f, col_n.items[k]);
        if (k < col_c.count) col_c.items[k] = std::max(1.0f, col_c.items[k]);
        if (k < col_h.count) col_h.items[k] = std::max(1.0f, col_h.items[k]);
    }

    if (max_cols_found >= 4 && col_c.count == col_s.count && col_h.count == col_s.count) {
        mat.strain_table = col_s;
        mat.emod_nominal_kpa = col_n;
        mat.emod_cold_kpa = col_c;
        mat.emod_hot_kpa = col_h;
        mat.linear_E_nominal_kpa = col_n.items[0];
        mat.linear_E_cold_kpa = col_c.items[0];
        mat.linear_E_hot_kpa = col_h.items[0];
        mat.sync_active_table();
        return 4;
    }

    // 2-column paste: update strain_table (and resize other temp tables if count changed) and set active temp curve
    bool same_count = (col_s.count == mat.strain_table.count);
    mat.strain_table = col_s;
    if (!same_count) {
        mat.emod_nominal_kpa = col_n;
        mat.emod_cold_kpa = col_n;
        mat.emod_hot_kpa = col_n;
        for (size_t k = 0; k < col_n.count; ++k) {
            mat.emod_cold_kpa.items[k] = col_n.items[k] * 1.35f;
            mat.emod_hot_kpa.items[k] = col_n.items[k] * 0.75f;
        }
    }
    if (mat.active_temp_mode == MAT_TEMP_COLD) {
        mat.emod_cold_kpa = col_n;
        mat.linear_E_cold_kpa = col_n.items[0];
    } else if (mat.active_temp_mode == MAT_TEMP_HOT) {
        mat.emod_hot_kpa = col_n;
        mat.linear_E_hot_kpa = col_n.items[0];
    } else {
        mat.emod_nominal_kpa = col_n;
        mat.linear_E_nominal_kpa = col_n.items[0];
    }
    mat.sync_active_table();
    return 2;
}

bool solve(float target_angle,
           BendStiffner &bs,
           Vec1<Vec1<float>> &results,
           Vec1<float> &strain,
           const Material &mat = Material(),
           float tension_kN = 0.0f,
           Vec1<float>* out_x_coords = nullptr,
           Vec1<float>* out_curvature = nullptr,
           Vec1<float>* out_bend_radius = nullptr,
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

    const size_t intervals = (DISCRETIZATION > 1) ? (DISCRETIZATION - 1) : 1;
    Vec2 m0_v0 = solve_bvp_material(bs, intervals, y0, theta0, ml, target_angle, mat, tension_N);
    float m0 = m0_v0.x;
    float v0 = m0_v0.y;

    results.clear();
    if (results.capacity < DISCRETIZATION) {
        results.reserve(DISCRETIZATION);
    }
    if (out_x_coords) {
        out_x_coords->clear();
        if (out_x_coords->capacity < DISCRETIZATION) out_x_coords->reserve(DISCRETIZATION);
    }

    State4 state = { { y0, theta0, m0, v0 }, 0.0f };
    float s = 0.0f;
    float h = bs.length / (float)intervals;

    for (size_t i = 0; i < DISCRETIZATION; ++i) {
        Vec1<float> node_state(NUMBER_OF_STATES);
        for (int k = 0; k < NUMBER_OF_STATES; ++k) {
            node_state.append(state.items[k]);
        }
        results.append(node_state);
        if (out_x_coords) {
            out_x_coords->append(state.x_pos);
        }

        if (i + 1 < DISCRETIZATION) {
            state = RK4_fast(s, state, h, bs, mat, tension_N, target_angle);
            s += h;
        }
    }

    calculate_strain_material(DISCRETIZATION, bs, results, strain, mat, tension_N, target_angle, out_curvature, out_bend_radius);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    if (out_elapsed_ms) {
        *out_elapsed_ms = elapsed.count();
    }

    if (verbose) {
        float min_mbr = out_bend_radius ? out_bend_radius->min_pos(1000.0f) : 0.0f;
        printf("Solved (angle = %.2f deg, T = %.1f kN) in %.3f ms | M0 = %.2f kN.m | Max Strain = %.4f | Min MBR = %.2f m\n",
               target_angle * 180.0f / 3.14159265359f,
               tension_kN,
               elapsed.count(),
               m0 / 1000.0f,
               strain.max_pos(),
               min_mbr);
    }

    return true;
}

// Fast zero-allocation batch evaluation of (Angle, Tension) pairs without overwriting the active single-solve buffers
void solve_pairs(Vec1<float> &strain_maxs,
                 Vec1<float> &pair_min_mbrs,
                 Vec1<float> &pair_root_moments_kNm,
                 const Vec1<float> &angles_deg,
                 const Vec1<float> &tensions_kN,
                 BendStiffner &bs,
                 Vec1<Vec1<float>> &results,
                 Vec1<float> &strain,
                 const Material &mat = Material(),
                 double* out_total_ms = nullptr)
{
    (void)results;
    (void)strain;
    auto start = std::chrono::high_resolution_clock::now();
    strain_maxs.clear();
    pair_min_mbrs.clear();
    pair_root_moments_kNm.clear();

    const size_t steps = (DISCRETIZATION > 1) ? (DISCRETIZATION - 1) : 1;
    for (size_t i = 0; i < angles_deg.count; ++i) {
        float angle_rad = angles_deg.items[i] * (3.14159265358979323846f / 180.0f);
        float tension_N = ((i < tensions_kN.count) ? tensions_kN.items[i] : 0.0f) * 1000.0f;
        float max_eps = 0.0f, avg_eps = 0.0f, tip_eps = 0.0f, root_eps = 0.0f, min_mbr = 1000.0f, m0_kNm = 0.0f;
        evaluate_candidate_single(bs, mat, steps, angle_rad, tension_N,
                                  max_eps, avg_eps, tip_eps, root_eps, min_mbr, m0_kNm);
        strain_maxs.append(max_eps);
        pair_min_mbrs.append(min_mbr);
        pair_root_moments_kNm.append(m0_kNm);
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
    float width = 1360.0f;
    float height = 840.0f;
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
    // Problem State Data
    ImVec4 clear_color = ImVec4(0.10f, 0.12f, 0.15f, 1.0f);
    BendStiffner bs = { 3.0f, 0.60f, 0.20f, 0.12f, 0.35f, 15.0f, 0.07f, 2.50f, true, false };
    Material material;
    float target_angle_deg = 25.0f;
    float target_tension_kN = 50.0f;
    float target_angle = 0.0f;
    float allowable_strain_pct = bs.allowable_strain * 100.0f;

    Vec1<Vec1<float>> results(DISCRETIZATION);
    Vec1<float> strain(DISCRETIZATION);
    Vec1<float> x_coords(DISCRETIZATION);
    Vec1<float> curvature(DISCRETIZATION);
    Vec1<float> bend_radius(DISCRETIZATION);
    for (int i = 0; i < DISCRETIZATION; ++i)
    {
        strain.append(0.0f);
        x_coords.append(0.0f);
        curvature.append(0.0f);
        bend_radius.append(1000.0f);
    }

    Vec1<float> strain_maxs(256);
    Vec1<float> pair_min_mbrs(256);
    Vec1<float> pair_root_moments_kNm(256);

    float length[DISCRETIZATION] = { 0.0f };
    float deflection_y[DISCRETIZATION] = { 0.0f };
    float moment_kNm[DISCRETIZATION] = { 0.0f };
    float def_top_x[DISCRETIZATION] = { 0.0f };
    float def_top_y[DISCRETIZATION] = { 0.0f };
    float def_bot_x[DISCRETIZATION] = { 0.0f };
    float def_bot_y[DISCRETIZATION] = { 0.0f };

    // Bare flexible pipe tail arrays (s in [L, L + L_tail])
    float tail_s[TAIL_NODES] = { 0.0f };
    float tail_x[TAIL_NODES] = { 0.0f };
    float tail_y[TAIL_NODES] = { 0.0f };
    float tail_theta[TAIL_NODES] = { 0.0f };
    float tail_moment_kNm[TAIL_NODES] = { 0.0f };
    float tail_curvature[TAIL_NODES] = { 0.0f };
    float tail_mbr[TAIL_NODES] = { 0.0f };
    float tail_top_x[TAIL_NODES] = { 0.0f };
    float tail_top_y[TAIL_NODES] = { 0.0f };
    float tail_bot_x[TAIL_NODES] = { 0.0f };
    float tail_bot_y[TAIL_NODES] = { 0.0f };

    bool solved = false;
    bool extract_results = true;
    bool reset_plot_limits = true;
    bool show_input_window = false;
    bool show_material_window = false;
    bool show_optimizer_window = false;
    bool show_temp_env_window = false;
    bool data_pasted = false;
    bool plot_bs = true;
    bool plot_deformed = true;
    bool plot_curvature = true;
    bool plot_moment = false;
    bool plot_pairs = false;
    double last_solve_ms = 0.0;
    double last_pairs_ms = 0.0;
    double last_opt_ms = 0.0;
    OptimizationConfig opt_cfg;
    OptimizationResult last_opt_res;
    TemperatureEnvelopeSummary temp_env;
    char project_filename[128] = "project.bsproj";
    std::string status_message = "Ready. Click 'Calculate Strain' or 'Optimize Size' to run.";

    Vec1<float> data1(256); // Angle, deg
    Vec1<float> data2(256); // Tension, kN

    auto update_post_solve_data = [&]() {
        float step = (DISCRETIZATION > 1) ? (bs.length / (float)(DISCRETIZATION - 1)) : 0.0f;
        for (int i = 0; i < DISCRETIZATION; i++)
        {
            float s_i = (float)i * step;
            length[i] = s_i;
            deflection_y[i] = results.items[i].items[0];
            float theta_i = results.items[i].items[1];
            moment_kNm[i] = results.items[i].items[2] / 1000.0f;
            float r_o = 0.5f * get_dia(bs, s_i);
            float nx = -std::sin(theta_i);
            float ny = std::cos(theta_i);
            def_top_x[i] = x_coords.items[i] + r_o * nx;
            def_top_y[i] = deflection_y[i] + r_o * ny;
            def_bot_x[i] = x_coords.items[i] - r_o * nx;
            def_bot_y[i] = deflection_y[i] - r_o * ny;
        }

        // Simulate bare flexible pipe tail beyond the tip (s in [L, L + L_tail])
        float tail_len = std::clamp(0.60f * bs.length, 0.8f, 2.5f);
        float tip_x = x_coords.items[DISCRETIZATION - 1];
        float tip_y = deflection_y[DISCRETIZATION - 1];
        float tip_theta = results.items[DISCRETIZATION - 1].items[1];
        float tip_m_Nm = results.items[DISCRETIZATION - 1].items[2];
        float rad = (3.14159265358979323846f / 180.0f) * target_angle_deg;
        simulate_bare_pipe_tail(bs, tip_x, tip_y, tip_theta, tip_m_Nm,
                                target_tension_kN * 1000.0f, rad, tail_len,
                                tail_s, tail_x, tail_y, tail_theta,
                                tail_moment_kNm, tail_curvature, tail_mbr);
        float r_pipe = 0.5f * bs.inner_dia;
        for (int i = 0; i < TAIL_NODES; ++i) {
            float nx = -std::sin(tail_theta[i]);
            float ny = std::cos(tail_theta[i]);
            tail_top_x[i] = tail_x[i] + r_pipe * nx;
            tail_top_y[i] = tail_y[i] + r_pipe * ny;
            tail_bot_x[i] = tail_x[i] - r_pipe * nx;
            tail_bot_y[i] = tail_y[i] - r_pipe * ny;
        }

        // Evaluate Cold (4C) / Nominal (20C) / Warm (35C) temperature envelope & hang-off reactions
        temp_env = evaluate_temperature_envelope(bs, material, target_angle_deg, target_tension_kN);
    };

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

        // --- Section 1: Geometry, Pipe, Root Insert & Design Criteria Inputs ---
        ImGui::SeparatorText("Bend Stiffener Geometry, Steel Root Insert, Annular Gap & Design Limits");
        float col2_x = 420.0f * main_scale;
        float col3_x = 840.0f * main_scale;
        ImGui::PushItemWidth(165.0f * main_scale);

        ImGui::InputFloat("Total Length L [m]", &bs.length, 0.1f, 0.5f, "%.3f");
        ImGui::SameLine(col2_x);
        ImGui::InputFloat("Root Diameter [m]", &bs.root_dia, 0.02f, 0.1f, "%.3f");
        ImGui::SameLine(col3_x);
        ImGui::InputFloat("Internal Pipe EI [kN.m2]", &bs.pipe_EI_kNm2, 1.0f, 5.0f, "%.2f");

        ImGui::InputFloat("Root Cylinder L_cyl [m]", &bs.root_cyl_length, 0.05f, 0.2f, "%.3f");
        ImGui::SameLine(col2_x);
        ImGui::InputFloat("Tip Diameter [m]", &bs.tip_dia, 0.02f, 0.1f, "%.3f");
        ImGui::SameLine(col3_x);
        if (ImGui::InputFloat("Allowable Strain [%]", &allowable_strain_pct, 0.5f, 1.0f, "%.2f")) {
            bs.allowable_strain = std::max(0.001f, allowable_strain_pct / 100.0f);
        }

        ImGui::InputFloat("Inner Bore Dia [m]", &bs.inner_dia, 0.01f, 0.05f, "%.3f");
        ImGui::SameLine(col2_x);
        ImGui::InputFloat("Pipe Outer Dia [m]", &bs.pipe_od, 0.005f, 0.02f, "%.3f");
        ImGui::SameLine(col3_x);
        ImGui::InputFloat("Allowable Pipe MBR [m]", &bs.allowable_mbr_m, 0.1f, 0.5f, "%.2f");

        ImGui::InputFloat("Steel Insert L_ins [m]", &bs.insert_length, 0.05f, 0.1f, "%.3f");
        ImGui::SameLine(col2_x);
        ImGui::InputFloat("Steel Insert Dia [m]", &bs.insert_dia, 0.01f, 0.05f, "%.3f");
        ImGui::SameLine(col3_x);
        ImGui::InputFloat("Allow Bond Stress [MPa]", &bs.allowable_bond_stress_mpa, 0.1f, 0.5f, "%.2f");

        ImGui::InputFloat("Target Angle [deg]", &target_angle_deg, 1.0f, 5.0f, "%.2f");
        ImGui::SameLine(col2_x);
        ImGui::InputFloat("Target Tension [kN]", &target_tension_kN, 5.0f, 25.0f, "%.2f");
        ImGui::PopItemWidth();

        ImGui::SameLine(col3_x);
        ImGui::Checkbox("Arc-Length Elastica (s)", &bs.use_arc_length);
        ImGui::SameLine();
        ImGui::Checkbox("Far-Field Tension BC (OrcaFlex)", &bs.use_far_field_bc);

        ImGui::Checkbox("Non-Linear Material (PU)", &material.is_nonlinear);
        ImGui::SameLine();
        ImGui::Text("Active Temp:");
        ImGui::SameLine();
        if (ImGui::RadioButton("Nominal", material.active_temp_mode == MAT_TEMP_NOMINAL)) {
            material.set_active_temp_mode(MAT_TEMP_NOMINAL);
        }
        ImGui::SameLine();
        if (ImGui::RadioButton("Cold", material.active_temp_mode == MAT_TEMP_COLD)) {
            material.set_active_temp_mode(MAT_TEMP_COLD);
        }
        ImGui::SameLine();
        if (ImGui::RadioButton("Hot", material.active_temp_mode == MAT_TEMP_HOT)) {
            material.set_active_temp_mode(MAT_TEMP_HOT);
        }
        ImGui::SameLine();
        if (ImGui::Button("Material Properties...")) {
            show_material_window = true;
        }
        ImGui::SameLine();
        if (ImGui::Button("Angle/Tension Pairs Table...")) {
            show_input_window = true;
        }
        ImGui::SameLine();
        if (ImGui::Button("Size Optimizer Settings...")) {
            show_optimizer_window = true;
        }
        ImGui::SameLine();
        if (ImGui::Button("Temp Envelope & Reactions...")) {
            temp_env = evaluate_temperature_envelope(bs, material, target_angle_deg, target_tension_kN);
            show_temp_env_window = true;
        }
        ImGui::SameLine();
        if (ImGui::Button("Copy OrcaFlex Geometry")) {
            char buf[256];
            std::string tsv = "ProfileLength_m\tOuterDia_m\tInnerDia_m\n";
            snprintf(buf, sizeof(buf), "%.4f\t%.4f\t%.4f\n", 0.0f, bs.root_dia, bs.inner_dia);
            tsv += buf;
            if (bs.root_cyl_length > 1e-4f && bs.root_cyl_length < bs.length) {
                snprintf(buf, sizeof(buf), "%.4f\t%.4f\t%.4f\n", bs.root_cyl_length, bs.root_dia, bs.inner_dia);
                tsv += buf;
            }
            snprintf(buf, sizeof(buf), "%.4f\t%.4f\t%.4f\n", bs.length, bs.tip_dia, bs.inner_dia);
            tsv += buf;
            ImGui::SetClipboardText(tsv.c_str());
            status_message = "Copied Bend Stiffener Geometry table for OrcaFlex to clipboard.";
        }
        ImGui::SameLine();
        if (ImGui::Button("Save Project")) {
            if (save_project_file(project_filename, bs, material, opt_cfg, target_angle_deg, target_tension_kN, data1, data2)) {
                status_message = std::string("Saved project to '") + project_filename + "'.";
            } else {
                status_message = "Error: Could not save project file.";
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Load Project")) {
            if (load_project_file(project_filename, bs, material, opt_cfg, target_angle_deg, target_tension_kN, data1, data2)) {
                allowable_strain_pct = bs.allowable_strain * 100.0f;
                data_pasted = (data1.count > 0);
                target_angle = (3.14159265358979323846f / 180.0f) * target_angle_deg;
                solved = solve(target_angle, bs, results, strain, material, target_tension_kN,
                               &x_coords, &curvature, &bend_radius, &last_solve_ms, false);
                extract_results = solved;
                if (solved) {
                    update_post_solve_data();
                    reset_plot_limits = true;
                }
                status_message = std::string("Loaded project from '") + project_filename + "'.";
            } else {
                status_message = std::string("No existing '") + project_filename + "' found to load.";
            }
        }
        ImGui::SameLine();
        float live_bond_mpa = (solved && results.count > 0)
            ? get_insert_bond_stress_mpa(bs, results.items[0].items[2], results.items[0].items[3]) : 0.0f;
        ImGui::TextDisabled("| PU Mass: %.1f kg | Gap: %.1f mm (%.2f deg) | Bond: %.2f MPa",
                            get_pu_mass_kg(bs), get_annular_gap_mm(bs), get_gap_contact_angle_deg(bs), live_bond_mpa);

        // --- Section 2: Solver Controls & Export ---
        ImGui::SeparatorText("Solver, Verification & Export Controls");
        target_angle = (3.14159265358979323846f / 180.0f) * target_angle_deg;

        if (ImGui::Button("Calculate Strain & Curvature"))
        {
            if (!is_valid_geometry(bs)) {
                status_message = "Error: Invalid geometry! Ensure 0 <= L_cyl < L and Root Dia >= Tip Dia > Inner Dia >= 0.";
                solved = false;
            } else if (!bs.use_arc_length && std::fabs(target_angle_deg) >= 89.5f) {
                status_message = "Error: Cartesian X-mode requires Target Angle in (-89, +89) deg. Enable Arc-Length Elastica (s) for large angles.";
                solved = false;
            } else if (bs.use_far_field_bc && std::fabs(target_tension_kN) < 1e-3f) {
                status_message = "Error: Far-Field Tension BC requires Target Tension > 0 kN.";
                solved = false;
            } else {
                solved = solve(target_angle, bs, results, strain, material, target_tension_kN,
                               &x_coords, &curvature, &bend_radius, &last_solve_ms, true);
                extract_results = solved;
                if (solved) {
                    update_post_solve_data();
                    reset_plot_limits = true;
                    float max_eps = strain.max_pos();
                    float min_mbr = bend_radius.min_pos(1000.0f);
                    float bond_mpa = get_insert_bond_stress_mpa(bs, results.items[0].items[2], results.items[0].items[3]);
                    bool pass_strain = (max_eps <= bs.allowable_strain);
                    bool pass_mbr = (min_mbr >= bs.allowable_mbr_m);
                    bool pass_bond = (bs.insert_length <= 1e-4f) || (bond_mpa <= bs.allowable_bond_stress_mpa);
                    float v0_kN = std::fabs(results.items[0].items[3]) / 1000.0f;
                    char buf[380];
                    snprintf(buf, sizeof(buf),
                             "[%s] Solved in %.3f ms | Max Strain: %.2f%% (Allow %.1f%%) | Min MBR: %.2f m (Allow %.2f m) | Root M0: %.2f kN.m, V0: %.1f kN | Bond: %.2f MPa",
                             (pass_strain && pass_mbr && pass_bond) ? "PASS" : "FAIL",
                             last_solve_ms,
                             max_eps * 100.0f, bs.allowable_strain * 100.0f,
                             min_mbr, bs.allowable_mbr_m,
                             results.items[0].items[2] / 1000.0f, v0_kN, bond_mpa);
                    status_message = buf;
                }
            }
        }

        ImGui::SameLine();
        if (ImGui::Button("Optimize Bend Stiffener Size"))
        {
            auto t_opt0 = std::chrono::high_resolution_clock::now();
            opt_cfg.use_all_pairs = (data1.count > 0 && opt_cfg.use_all_pairs);
            last_opt_res = optimize_bend_stiffener(bs, material, opt_cfg, target_angle_deg, target_tension_kN, &data1, &data2);
            auto t_opt1 = std::chrono::high_resolution_clock::now();
            last_opt_ms = std::chrono::duration<double, std::milli>(t_opt1 - t_opt0).count();

            bs = last_opt_res.best_bs;
            if (opt_cfg.force_far_field_when_tension && std::fabs(target_tension_kN) > 0.01f) {
                bs.use_far_field_bc = true;
            }
            target_angle = (3.14159265358979323846f / 180.0f) * target_angle_deg;
            solved = solve(target_angle, bs, results, strain, material, target_tension_kN,
                           &x_coords, &curvature, &bend_radius, &last_solve_ms, true);
            extract_results = solved;
            if (solved) {
                update_post_solve_data();
                reset_plot_limits = true;
            }
            if (opt_cfg.use_all_pairs && data1.count > 0) {
                solve_pairs(strain_maxs, pair_min_mbrs, pair_root_moments_kNm,
                            data1, data2, bs, results, strain, material, &last_pairs_ms);
            }
            char buf[340];
            snprintf(buf, sizeof(buf),
                     "[%s] Optimized (%d evals in %.1f ms): L=%.3f m (L_cyl=%.3f m), D_root=%.3f m, D_tip=%.3f m | Mass=%.1f kg | Max Strain=%.2f%% | Min MBR=%.2f m",
                     last_opt_res.feasible ? "OPTIMAL" : "BEST-EFFORT",
                     last_opt_res.evaluations, last_opt_ms,
                     bs.length, bs.root_cyl_length, bs.root_dia, bs.tip_dia,
                     last_opt_res.mass_kg,
                     last_opt_res.max_strain * 100.0f,
                     last_opt_res.min_mbr_m);
            status_message = buf;
        }

        ImGui::SameLine();
        if (ImGui::Button("Calculate All Input Pairs"))
        {
            if (!is_valid_geometry(bs)) {
                status_message = "Error: Invalid geometry!";
            } else if (data1.count == 0) {
                status_message = "No Angle/Tension pairs loaded. Open 'Angle/Tension Pairs Table...' and paste from Excel (Ctrl+V).";
            } else {
                solve_pairs(strain_maxs, pair_min_mbrs, pair_root_moments_kNm,
                            data1, data2, bs, results, strain, material, &last_pairs_ms);
                plot_pairs = true;
                char buf[256];
                snprintf(buf, sizeof(buf),
                         "Solved %zu pairs in %.2f ms (%.3f ms/pair) | Peak Strain: %.2f%% | Min MBR: %.2f m | Max M0: %.2f kN.m",
                         data1.count, last_pairs_ms, last_pairs_ms / (double)data1.count,
                         strain_maxs.max_pos() * 100.0f,
                         pair_min_mbrs.min_pos(1000.0f),
                         pair_root_moments_kNm.max_pos());
                status_message = buf;
            }
        }

        if (solved) {
            ImGui::SameLine();
            if (ImGui::Button("Copy Profile (Excel)")) {
                std::string tsv = "ArcLength_m\tDeformedX_m\tDeformedY_m\tAngle_deg\tOuterDia_m\tMoment_kNm\tCurvature_1_m\tBendRadius_m\tStrain\n";
                char row[256];
                for (int i = 0; i < DISCRETIZATION; ++i) {
                    snprintf(row, sizeof(row), "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.6f\t%.4f\t%.6f\n",
                             length[i], x_coords.items[i], deflection_y[i],
                             results.items[i].items[1] * 180.0f / 3.14159265f,
                             get_dia(bs, length[i]), moment_kNm[i],
                             curvature.items[i], bend_radius.items[i], strain.items[i]);
                    tsv += row;
                }
                ImGui::SetClipboardText(tsv.c_str());
                status_message = "Copied full 9-column spatial profile to clipboard (ready to paste into Excel).";
            }
            ImGui::SameLine();
            if (ImGui::Button("Export CSV")) {
                FILE* fp = std::fopen("bs_profile_results.csv", "w");
                if (fp) {
                    std::fprintf(fp, "ArcLength_m,DeformedX_m,DeformedY_m,Angle_deg,OuterDia_m,Moment_kNm,Curvature_1_m,BendRadius_m,Strain\n");
                    for (int i = 0; i < DISCRETIZATION; ++i) {
                        std::fprintf(fp, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.6f,%.4f,%.6f\n",
                                     length[i], x_coords.items[i], deflection_y[i],
                                     results.items[i].items[1] * 180.0f / 3.14159265f,
                                     get_dia(bs, length[i]), moment_kNm[i],
                                     curvature.items[i], bend_radius.items[i], strain.items[i]);
                    }
                    std::fclose(fp);
                    status_message = "Exported spatial profile to 'bs_profile_results.csv'.";
                }
            }
            ImGui::SameLine();
            if (ImGui::Button("Generate HTML Report")) {
                float theta_arr[DISCRETIZATION];
                float mom_Nm_arr[DISCRETIZATION];
                float shear_N_arr[DISCRETIZATION];
                for (int i = 0; i < DISCRETIZATION; ++i) {
                    theta_arr[i] = results.items[i].items[1];
                    mom_Nm_arr[i] = results.items[i].items[2];
                    shear_N_arr[i] = results.items[i].items[3];
                }
                float m0_kNm = results.items[0].items[2] / 1000.0f;
                float v0_kN = results.items[0].items[3] / 1000.0f;
                float n0_kN = std::fabs(target_tension_kN) * std::max(0.0f, std::cos(target_angle));
                if (generate_html_report("bs_engineering_report.html", bs, material,
                                         target_angle_deg, target_tension_kN,
                                         m0_kNm, v0_kN, n0_kN,
                                         length, x_coords.items, deflection_y,
                                         theta_arr, mom_Nm_arr, shear_N_arr,
                                         strain.items, curvature.items, bend_radius.items,
                                         DISCRETIZATION, tail_x, tail_y, TAIL_NODES)) {
                    status_message = "Generated HTML Engineering Calculation Report: 'bs_engineering_report.html'.";
                } else {
                    status_message = "Error: Failed to write 'bs_engineering_report.html'.";
                }
            }
        }

        ImGui::Checkbox("Geometry & Deformed Shape", &plot_bs);
        if (plot_bs && solved) {
            ImGui::SameLine();
            ImGui::Checkbox("Overlay 2D Deformed Cone", &plot_deformed);
        }
        ImGui::SameLine();
        ImGui::Checkbox("Curvature & MBR Check", &plot_curvature);
        ImGui::SameLine();
        ImGui::Checkbox("Bending Moment", &plot_moment);
        ImGui::SameLine();
        ImGui::Checkbox("Pairs Envelope Plot", &plot_pairs);

        bool is_fail = (status_message.rfind("[FAIL]", 0) == 0 || status_message.rfind("Error:", 0) == 0);
        ImVec4 status_col = is_fail ? ImVec4(1.0f, 0.42f, 0.42f, 1.0f) : ImVec4(0.4f, 0.92f, 0.62f, 1.0f);
        ImGui::TextColored(status_col, "Status: %s", status_message.c_str());

        // --- Section 3: Plots ---
        float plot_height = std::max(175.0f, 0.25f * io.DisplaySize.y);

        // Undeformed polygon with cylindrical root L_cyl + conical taper
        float c_len = std::clamp(bs.root_cyl_length, 0.0f, std::max(0.0f, bs.length));
        float x_data[7] = { 0.0f, 0.0f, c_len, bs.length, bs.length, c_len, 0.0f };
        float y_data[7] = {
            -bs.root_dia / 2.0f, bs.root_dia / 2.0f, bs.root_dia / 2.0f,
             bs.tip_dia / 2.0f, -bs.tip_dia / 2.0f, -bs.root_dia / 2.0f, -bs.root_dia / 2.0f
        };
        float inner_x[2] = { 0.0f, bs.length };
        float inner_y_top[2] = { bs.inner_dia / 2.0f, bs.inner_dia / 2.0f };
        float inner_y_bot[2] = { -bs.inner_dia / 2.0f, -bs.inner_dia / 2.0f };

        if (plot_bs)
        {
            if (ImPlot::BeginPlot("Bend Stiffener Geometry & 2D Deformed Shape (with Steel Insert & Exiting Pipe Tail)", ImVec2(-1, plot_height))) {
                ImPlot::SetupAxes("Axial Position X [m]", "Transverse Position Y [m]");
                float max_x_extent = solved ? std::max(bs.length * 1.15f, tail_x[TAIL_NODES - 1] * 1.05f) : (bs.length * 1.10f);
                float max_y_extent = std::max(0.65f * bs.root_dia, solved ? (std::fabs(tail_y[TAIL_NODES - 1]) + bs.root_dia) : (0.65f * bs.root_dia));
                ImPlot::SetupAxesLimits(-0.05 * bs.length, max_x_extent, -max_y_extent, max_y_extent,
                                        reset_plot_limits ? ImGuiCond_Always : ImGuiCond_Once);
                reset_plot_limits = false;
                ImPlotSpec spec;
                spec.LineWeight = 2.0f * main_scale;
                ImPlot::PlotLine("Undeformed Cone", x_data, y_data, 7, spec);
                ImPlotSpec inner_spec;
                inner_spec.LineWeight = 1.2f * main_scale;
                ImPlot::PlotLine("Inner Bore (+)", inner_x, inner_y_top, 2, inner_spec);
                ImPlot::PlotLine("Inner Bore (-)", inner_x, inner_y_bot, 2, inner_spec);

                if (bs.pipe_od > 1e-4f && bs.pipe_od < bs.inner_dia) {
                    float pipe_y_top[2] = { bs.pipe_od / 2.0f, bs.pipe_od / 2.0f };
                    float pipe_y_bot[2] = { -bs.pipe_od / 2.0f, -bs.pipe_od / 2.0f };
                    ImPlot::PlotLine("Pipe OD (+)", inner_x, pipe_y_top, 2, inner_spec);
                    ImPlot::PlotLine("Pipe OD (-)", inner_x, pipe_y_bot, 2, inner_spec);
                }

                if (bs.insert_length > 1e-4f) {
                    float l_ins = std::min(bs.insert_length, bs.length);
                    float hi_ins = std::max(bs.inner_dia, bs.root_dia - 0.02f);
                    float r_ins = 0.5f * std::clamp(bs.insert_dia, bs.inner_dia, hi_ins);
                    float r_in = 0.5f * bs.inner_dia;
                    float ins_x[5] = { 0.0f, l_ins, l_ins, 0.0f, 0.0f };
                    float ins_y_top[5] = { r_in, r_in, r_ins, r_ins, r_in };
                    float ins_y_bot[5] = { -r_in, -r_in, -r_ins, -r_ins, -r_in };
                    ImPlotSpec ins_spec;
                    ins_spec.LineWeight = 1.8f * main_scale;
                    ImPlot::PlotLine("Steel Root Insert (+)", ins_x, ins_y_top, 5, ins_spec);
                    ImPlot::PlotLine("Steel Root Insert (-)", ins_x, ins_y_bot, 5, ins_spec);
                }

                if (solved && plot_deformed) {
                    ImPlotSpec def_center_spec;
                    def_center_spec.LineWeight = 2.2f * main_scale;
                    ImPlot::PlotLine("Deformed BS Centerline", x_coords.items, deflection_y, DISCRETIZATION, def_center_spec);
                    ImPlotSpec def_wall_spec;
                    def_wall_spec.LineWeight = 1.6f * main_scale;
                    ImPlot::PlotLine("Deformed BS Upper Wall", def_top_x, def_top_y, DISCRETIZATION, def_wall_spec);
                    ImPlot::PlotLine("Deformed BS Lower Wall", def_bot_x, def_bot_y, DISCRETIZATION, def_wall_spec);

                    // Exiting bare flexible pipe tail (s in [L, L + L_tail])
                    ImPlotSpec tail_spec;
                    tail_spec.LineWeight = 2.0f * main_scale;
                    ImPlot::PlotLine("Exiting Pipe Centerline", tail_x, tail_y, TAIL_NODES, tail_spec);
                    ImPlotSpec tail_wall_spec;
                    tail_wall_spec.LineWeight = 1.2f * main_scale;
                    ImPlot::PlotLine("Exiting Pipe Wall (+)", tail_top_x, tail_top_y, TAIL_NODES, tail_wall_spec);
                    ImPlot::PlotLine("Exiting Pipe Wall (-)", tail_bot_x, tail_bot_y, TAIL_NODES, tail_wall_spec);
                }
                ImPlot::EndPlot();
            }
        }

        if (solved && extract_results)
        {
            if (ImPlot::BeginPlot("Bend Stiffener Outer-Fiber Strain vs. Allowable Limit", ImVec2(-1, plot_height))) {
                ImPlot::SetupAxes("Arc Length from Root s [m]", "Outer Fiber Strain [-]");
                ImPlotSpec spec_line;
                spec_line.LineWeight = 2.0f * main_scale;
                spec_line.Marker = ImPlotMarker_Circle;
                spec_line.MarkerSize = 3.2f * main_scale;
                ImPlot::PlotLine("Total Strain", length, strain.items, DISCRETIZATION, spec_line);

                float lim_x[2] = { 0.0f, bs.length };
                float lim_y[2] = { bs.allowable_strain, bs.allowable_strain };
                ImPlotSpec lim_spec;
                lim_spec.LineWeight = 1.8f * main_scale;
                ImPlot::PlotLine("Allowable PU Strain Limit", lim_x, lim_y, 2, lim_spec);
                ImPlot::EndPlot();
            }

            if (plot_curvature)
            {
                if (ImPlot::BeginPlot("Curvature & Minimum Bend Radius (MBR) Verification (Inside BS + Exiting Pipe Tail)", ImVec2(-1, plot_height))) {
                    ImPlot::SetupAxes("Arc Length from Root s [m]", "Curvature kappa [1/m]");
                    ImPlotSpec spec_k;
                    spec_k.LineWeight = 2.0f * main_scale;
                    ImPlot::PlotLine("Curvature Inside BS [0..L]", length, curvature.items, DISCRETIZATION, spec_k);
                    ImPlot::PlotLine("Exiting Bare Pipe Tail [L..L+tail]", tail_s, tail_curvature, TAIL_NODES, spec_k);

                    float max_allow_kappa = (bs.allowable_mbr_m > 1e-4f) ? (1.0f / bs.allowable_mbr_m) : 1.0f;
                    float lim_x[2] = { 0.0f, tail_s[TAIL_NODES - 1] };
                    float lim_k[2] = { max_allow_kappa, max_allow_kappa };
                    ImPlotSpec lim_spec;
                    lim_spec.LineWeight = 1.8f * main_scale;
                    ImPlot::PlotLine("Max Allowable Curvature (1 / MBR)", lim_x, lim_k, 2, lim_spec);
                    ImPlot::EndPlot();
                }
            }

            if (plot_moment)
            {
                if (ImPlot::BeginPlot("Bending Moment Distribution (Inside BS + Exiting Pipe Tail)", ImVec2(-1, plot_height))) {
                    ImPlot::SetupAxes("Arc Length from Root s [m]", "Bending Moment [kN.m]");
                    ImPlotSpec spec_m;
                    spec_m.LineWeight = 2.0f * main_scale;
                    ImPlot::PlotLine("Bending Moment Inside BS [kN.m]", length, moment_kNm, DISCRETIZATION, spec_m);
                    ImPlot::PlotLine("Exiting Bare Pipe Moment [kN.m]", tail_s, tail_moment_kNm, TAIL_NODES, spec_m);
                    ImPlot::EndPlot();
                }
            }
        }

        if (plot_pairs && strain_maxs.count > 0)
        {
            if (ImPlot::BeginPlot("Max Strain Envelope for Tension / Angle Pairs", ImVec2(-1, plot_height))) {
                ImPlot::SetupAxes("Angle [deg]", "Max Strain [-]");
                ImPlotSpec spec_p;
                spec_p.LineWeight = 1.5f * main_scale;
                spec_p.Marker = ImPlotMarker_Circle;
                spec_p.MarkerSize = 4.0f * main_scale;
                size_t plot_count = std::min(data1.count, strain_maxs.count);
                ImPlot::PlotScatter("Max Strain per Pair", data1.items, strain_maxs.items, (int)plot_count, spec_p);

                if (data1.count > 0) {
                    float min_x = std::min(0.0f, data1.min_val(0.0f) * 1.05f);
                    float max_x = std::max(0.0f, data1.max_val(0.0f) * 1.05f);
                    if (std::fabs(max_x - min_x) < 1e-3f) max_x = min_x + 1.0f;
                    float lim_x[2] = { min_x, max_x };
                    float lim_y[2] = { bs.allowable_strain, bs.allowable_strain };
                    ImPlotSpec lim_spec;
                    lim_spec.LineWeight = 1.8f * main_scale;
                    ImPlot::PlotLine("Allowable Strain Limit", lim_x, lim_y, 2, lim_spec);
                }
                ImPlot::EndPlot();
            }
        }

        // --- Modal / Floating Window 1: Angle & Tension Pairs Input & Verification Table ---
        if (show_input_window)
        {
            ImGui::SetNextWindowSize(ImVec2(720 * main_scale, 460 * main_scale), ImGuiCond_FirstUseEver);
            ImGui::Begin("Angle & Tension Pairs Input & Envelope Verification", &show_input_window);
            if (ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows) &&
                ImGui::GetIO().KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_V))
            {
                data1.clear();
                data2.clear();
                strain_maxs.clear();
                pair_min_mbrs.clear();
                pair_root_moments_kNm.clear();
                const char* clipboardData = ImGui::GetClipboardText();
                if (clipboardData) {
                    ParseTwoColumns(data1, data2, clipboardData);
                }
                data_pasted = (data1.count > 0);
                if (data_pasted) {
                    target_angle_deg = data1.max_abs_signed(target_angle_deg);
                }
            }

            ImGui::TextWrapped("Focus this window and press Ctrl+V to paste 2 columns (Angle [deg], Tension [kN]) from Excel.");
            if (ImGui::Button("Add Row")) {
                data1.append(15.0f);
                data2.append(50.0f);
                data_pasted = true;
            }
            ImGui::SameLine();
            if (ImGui::Button("Run All Pairs")) {
                if (is_valid_geometry(bs) && data1.count > 0) {
                    solve_pairs(strain_maxs, pair_min_mbrs, pair_root_moments_kNm,
                                data1, data2, bs, results, strain, material, &last_pairs_ms);
                    plot_pairs = true;
                }
            }
            ImGui::SameLine();
            if (strain_maxs.count > 0 && ImGui::Button("Copy Results (Excel)")) {
                std::string tsv = "Angle_deg\tTension_kN\tMaxStrain\tMinMBR_m\tRootMoment_kNm\tStatus\n";
                char row[256];
                size_t n_rows = std::min(data1.count, strain_maxs.count);
                for (size_t i = 0; i < n_rows; ++i) {
                    float eps_i = strain_maxs.items[i];
                    float mbr_i = (i < pair_min_mbrs.count) ? pair_min_mbrs.items[i] : 1000.0f;
                    float m0_i = (i < pair_root_moments_kNm.count) ? pair_root_moments_kNm.items[i] : 0.0f;
                    bool ok = (eps_i <= bs.allowable_strain) && (mbr_i >= bs.allowable_mbr_m);
                    snprintf(row, sizeof(row), "%.3f\t%.3f\t%.6f\t%.3f\t%.3f\t%s\n",
                             data1.items[i], data2.items[i], eps_i, mbr_i, m0_i, ok ? "PASS" : "FAIL");
                    tsv += row;
                }
                ImGui::SetClipboardText(tsv.c_str());
                status_message = "Copied pairs summary table to clipboard (ready to paste into Excel).";
            }
            ImGui::SameLine();
            if (strain_maxs.count > 0 && ImGui::Button("Export Pairs CSV")) {
                FILE* fp = std::fopen("bs_pairs_results.csv", "w");
                if (fp) {
                    std::fprintf(fp, "Angle_deg,Tension_kN,MaxStrain,MinMBR_m,RootMoment_kNm,Status\n");
                    size_t n_rows = std::min(data1.count, strain_maxs.count);
                    for (size_t i = 0; i < n_rows; ++i) {
                        float eps_i = strain_maxs.items[i];
                        float mbr_i = (i < pair_min_mbrs.count) ? pair_min_mbrs.items[i] : 1000.0f;
                        float m0_i = (i < pair_root_moments_kNm.count) ? pair_root_moments_kNm.items[i] : 0.0f;
                        bool ok = (eps_i <= bs.allowable_strain) && (mbr_i >= bs.allowable_mbr_m);
                        std::fprintf(fp, "%.3f,%.3f,%.6f,%.3f,%.3f,%s\n",
                                     data1.items[i], data2.items[i], eps_i, mbr_i, m0_i, ok ? "PASS" : "FAIL");
                    }
                    std::fclose(fp);
                    status_message = "Exported pairs summary table to 'bs_pairs_results.csv'.";
                }
            }
            ImGui::SameLine();
            if (ImGui::Button("Clear Table")) {
                data1.clear();
                data2.clear();
                strain_maxs.clear();
                pair_min_mbrs.clear();
                pair_root_moments_kNm.clear();
                data_pasted = false;
            }
            ImGui::SameLine();
            if (ImGui::Button("Close")) {
                show_input_window = false;
            }

            if (ImGui::BeginTable("InputTable", 6, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY, ImVec2(0, 330 * main_scale)))
            {
                ImGui::TableSetupColumn("Angle [deg]");
                ImGui::TableSetupColumn("Tension [kN]");
                ImGui::TableSetupColumn("Max Strain [%]");
                ImGui::TableSetupColumn("Min MBR [m]");
                ImGui::TableSetupColumn("Root M0 [kN.m]");
                ImGui::TableSetupColumn("Check");
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
                        ImGui::Text("%.2f%% (%.4f)", strain_maxs.items[i] * 100.0f, strain_maxs.items[i]);
                    } else {
                        ImGui::TextDisabled("-");
                    }

                    ImGui::TableNextColumn();
                    if (i < pair_min_mbrs.count) {
                        ImGui::Text("%.2f m", pair_min_mbrs.items[i]);
                    } else {
                        ImGui::TextDisabled("-");
                    }

                    ImGui::TableNextColumn();
                    if (i < pair_root_moments_kNm.count) {
                        ImGui::Text("%.2f", pair_root_moments_kNm.items[i]);
                    } else {
                        ImGui::TextDisabled("-");
                    }

                    ImGui::TableNextColumn();
                    if (i < strain_maxs.count && i < pair_min_mbrs.count) {
                        bool ok = (strain_maxs.items[i] <= bs.allowable_strain) &&
                                  (pair_min_mbrs.items[i] >= bs.allowable_mbr_m);
                        if (ok) {
                            ImGui::TextColored(ImVec4(0.35f, 0.95f, 0.55f, 1.0f), "PASS");
                        } else {
                            ImGui::TextColored(ImVec4(1.0f, 0.40f, 0.40f, 1.0f), "FAIL");
                        }
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
            ImGui::SetNextWindowSize(ImVec2(760 * main_scale, 680 * main_scale), ImGuiCond_FirstUseEver);
            ImGui::Begin("Temperature-Specific Material Properties (Nominal / Cold / Hot)", &show_material_window);

            if (ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows) &&
                ImGui::GetIO().KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_V))
            {
                const char* clipboardData = ImGui::GetClipboardText();
                if (clipboardData) {
                    int cols_pasted = ParseMaterialClipboard(material, clipboardData);
                    if (cols_pasted == 4) {
                        status_message = "Pasted 4-column (Strain, Nominal E, Cold E, Hot E) table from clipboard.";
                    } else if (cols_pasted == 2) {
                        status_message = "Pasted 2-column (Strain, E) table into the active temperature curve.";
                    }
                }
            }

            ImGui::SeparatorText("Active Operating Temperature & Reference Temperatures");
            if (ImGui::RadioButton("Nominal Temp", material.active_temp_mode == MAT_TEMP_NOMINAL)) {
                material.set_active_temp_mode(MAT_TEMP_NOMINAL);
            }
            ImGui::SameLine();
            if (ImGui::RadioButton("Cold Temp", material.active_temp_mode == MAT_TEMP_COLD)) {
                material.set_active_temp_mode(MAT_TEMP_COLD);
            }
            ImGui::SameLine();
            if (ImGui::RadioButton("Hot Temp", material.active_temp_mode == MAT_TEMP_HOT)) {
                material.set_active_temp_mode(MAT_TEMP_HOT);
            }

            ImGui::PushItemWidth(110.0f * main_scale);
            ImGui::InputFloat("Nominal [C]", &material.temp_nominal_c, 1.0f, 5.0f, "%.1f");
            ImGui::SameLine();
            ImGui::InputFloat("Cold [C]", &material.temp_cold_c, 1.0f, 5.0f, "%.1f");
            ImGui::SameLine();
            ImGui::InputFloat("Hot [C]", &material.temp_hot_c, 1.0f, 5.0f, "%.1f");
            ImGui::PopItemWidth();

            ImGui::SeparatorText("Constitutive Model & E-Modulus Input");
            ImGui::Checkbox("Enable Non-Linear Material Curve", &material.is_nonlinear);
            if (!material.is_nonlinear) {
                ImGui::PushItemWidth(180.0f * main_scale);
                if (ImGui::InputFloat("Nominal Young's Modulus E [kPa]", &material.linear_E_nominal_kpa, 1000.0f, 10000.0f, "%.1f")) {
                    material.sync_active_table();
                }
                if (ImGui::InputFloat("Cold Young's Modulus E [kPa]", &material.linear_E_cold_kpa, 1000.0f, 10000.0f, "%.1f")) {
                    material.sync_active_table();
                }
                if (ImGui::InputFloat("Hot Young's Modulus E [kPa]", &material.linear_E_hot_kpa, 1000.0f, 10000.0f, "%.1f")) {
                    material.sync_active_table();
                }
                ImGui::PopItemWidth();
            } else {
                ImGui::Checkbox("Table values are Tangent Modulus d(sigma)/d(strain) [integrated for Secant E]", &material.table_is_tangent_modulus);
                ImGui::Checkbox("Cross-Sectional Annular Stress Integration (4-pt Gauss quadrature)", &material.use_cross_section_integration);

                ImGui::TextWrapped("Tip: Press Ctrl+V to paste either 4 columns (Strain, Nominal E, Cold E, Hot E) or 2 columns (Strain, Active Temp E) from Excel.");
                if (ImGui::Button("Reset All to Default Polyurethane (20C / 4C / 35C)")) {
                    material.reset_to_default(1.0f, 1.35f, 0.75f);
                }
                ImGui::SameLine();
                if (ImGui::Button("Add Row")) {
                    size_t n = material.strain_table.count;
                    float next_s = (n > 0) ? (material.strain_table.items[n - 1] + 0.01f) : 0.0f;
                    float last_en = (material.emod_nominal_kpa.count > 0) ? material.emod_nominal_kpa.items[material.emod_nominal_kpa.count - 1] : 215800.0f;
                    float last_ec = (material.emod_cold_kpa.count > 0) ? material.emod_cold_kpa.items[material.emod_cold_kpa.count - 1] : (last_en * 1.35f);
                    float last_eh = (material.emod_hot_kpa.count > 0) ? material.emod_hot_kpa.items[material.emod_hot_kpa.count - 1] : (last_en * 0.75f);
                    material.strain_table.append(next_s);
                    material.emod_nominal_kpa.append(last_en);
                    material.emod_cold_kpa.append(last_ec);
                    material.emod_hot_kpa.append(last_eh);
                    material.sync_active_table();
                }
                ImGui::SameLine();
                if (ImGui::Button("Copy OrcaFlex Stress-Strain (Active Temp)")) {
                    std::string tsv = "Strain_pct\tStress_kPa\tSecantE_kPa\n";
                    char row[128];
                    size_t n_mat = material.strain_table.count;
                    for (size_t i = 0; i < n_mat; ++i) {
                        float eps_i = material.strain_table.items[i];
                        float esec_kpa = get_material_E(eps_i, material) / 1000.0f;
                        float stress_kpa = eps_i * esec_kpa;
                        snprintf(row, sizeof(row), "%.2f\t%.2f\t%.2f\n", eps_i * 100.0f, stress_kpa, esec_kpa);
                        tsv += row;
                    }
                    ImGui::SetClipboardText(tsv.c_str());
                    status_message = "Copied Active Temperature Stress-Strain curve for OrcaFlex to clipboard.";
                }
                ImGui::SameLine();
                if (ImGui::Button("Close")) {
                    show_material_window = false;
                }

                // Live interactive Stress-Strain Plot showing Nominal, Cold, and Hot curves together
                constexpr int N_CURVE = 41;
                float curve_strain_pct[N_CURVE];
                float curve_stress_nom[N_CURVE];
                float curve_stress_cold[N_CURVE];
                float curve_stress_hot[N_CURVE];
                Material mat_nom = material.with_temp_mode(MAT_TEMP_NOMINAL);
                Material mat_cold = material.with_temp_mode(MAT_TEMP_COLD);
                Material mat_hot = material.with_temp_mode(MAT_TEMP_HOT);

                for (int i = 0; i < N_CURVE; ++i) {
                    float eps = 0.20f * (float)i / (float)(N_CURVE - 1);
                    curve_strain_pct[i] = eps * 100.0f;
                    curve_stress_nom[i]  = (eps * get_material_E(eps, mat_nom)) / 1e6f;
                    curve_stress_cold[i] = (eps * get_material_E(eps, mat_cold)) / 1e6f;
                    curve_stress_hot[i]  = (eps * get_material_E(eps, mat_hot)) / 1e6f;
                }

                if (ImPlot::BeginPlot("Temperature-Dependent Polyurethane Stress-Strain Curves (Cold / Nominal / Hot)", ImVec2(-1, 185 * main_scale))) {
                    ImPlot::SetupAxes("Strain [%]", "Stress [MPa]");
                    ImPlotSpec sp_c;
                    sp_c.LineWeight = 2.0f * main_scale;
                    ImPlot::PlotLine("Cold Stress [MPa]", curve_strain_pct, curve_stress_cold, N_CURVE, sp_c);
                    ImPlot::PlotLine("Nominal Stress [MPa]", curve_strain_pct, curve_stress_nom, N_CURVE, sp_c);
                    ImPlot::PlotLine("Hot Stress [MPa]", curve_strain_pct, curve_stress_hot, N_CURVE, sp_c);

                    if (solved) {
                        float op_eps = strain.max_pos();
                        float op_pct = op_eps * 100.0f;
                        float op_sig = (op_eps * get_material_E(op_eps, material)) / 1e6f;
                        ImPlotSpec sp_op;
                        sp_op.Marker = ImPlotMarker_Circle;
                        sp_op.MarkerSize = 5.5f * main_scale;
                        sp_op.LineWeight = 2.0f * main_scale;
                        ImPlot::PlotScatter("Peak Operating Strain (Active Temp)", &op_pct, &op_sig, 1, sp_op);
                    }
                    ImPlot::EndPlot();
                }

                if (ImGui::BeginTable("MaterialTable", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY, ImVec2(0, 220 * main_scale)))
                {
                    char col_nom[64], col_cold[64], col_hot[64];
                    snprintf(col_nom, sizeof(col_nom), "Nominal E (%.0f C) [kPa]", material.temp_nominal_c);
                    snprintf(col_cold, sizeof(col_cold), "Cold E (%.0f C) [kPa]", material.temp_cold_c);
                    snprintf(col_hot, sizeof(col_hot), "Hot E (%.0f C) [kPa]", material.temp_hot_c);

                    ImGui::TableSetupColumn("Strain [-]");
                    ImGui::TableSetupColumn(col_nom);
                    ImGui::TableSetupColumn(col_cold);
                    ImGui::TableSetupColumn(col_hot);
                    ImGui::TableSetupColumn("Active Secant E / Stress");
                    ImGui::TableHeadersRow();

                    size_t n_mat = material.strain_table.count;
                    bool edited = false;
                    for (size_t i = 0; i < n_mat; ++i)
                    {
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn();
                        ImGui::PushID((int)i + 20000);
                        if (ImGui::InputFloat("##mat_strain", &material.strain_table.items[i], 0.0f, 0.0f, "%.4f")) edited = true;
                        ImGui::PopID();

                        ImGui::TableNextColumn();
                        ImGui::PushID((int)i + 30000);
                        if (ImGui::InputFloat("##mat_enom", &material.emod_nominal_kpa.items[i], 0.0f, 0.0f, "%.1f")) {
                            if (i == 0) material.linear_E_nominal_kpa = material.emod_nominal_kpa.items[0];
                            edited = true;
                        }
                        ImGui::PopID();

                        ImGui::TableNextColumn();
                        ImGui::PushID((int)i + 40000);
                        if (ImGui::InputFloat("##mat_ecold", &material.emod_cold_kpa.items[i], 0.0f, 0.0f, "%.1f")) {
                            if (i == 0) material.linear_E_cold_kpa = material.emod_cold_kpa.items[0];
                            edited = true;
                        }
                        ImGui::PopID();

                        ImGui::TableNextColumn();
                        ImGui::PushID((int)i + 50000);
                        if (ImGui::InputFloat("##mat_ehot", &material.emod_hot_kpa.items[i], 0.0f, 0.0f, "%.1f")) {
                            if (i == 0) material.linear_E_hot_kpa = material.emod_hot_kpa.items[0];
                            edited = true;
                        }
                        ImGui::PopID();

                        ImGui::TableNextColumn();
                        float s_i = material.strain_table.items[i];
                        float e_sec_mpa = get_material_E(s_i, material) / 1e6f;
                        ImGui::Text("%.1f MPa (%.2f MPa)", e_sec_mpa, s_i * e_sec_mpa);
                    }
                    if (edited) {
                        material.sync_active_table();
                    }
                    ImGui::EndTable();
                }
            }
            ImGui::End();
        }

        // --- Modal / Floating Window 3: Bend Stiffener Size Optimizer ---
        if (show_optimizer_window)
        {
            ImGui::SetNextWindowSize(ImVec2(720 * main_scale, 690 * main_scale), ImGuiCond_FirstUseEver);
            ImGui::Begin("Bend Stiffener Size Optimizer", &show_optimizer_window);

            ImGui::TextWrapped("Automatically sizes Length (L) and Root Diameter (D_root) to minimize Polyurethane Volume & Mass while keeping Inner Diameter and Tip Diameter fixed by the flexible pipe.");
            ImGui::SeparatorText("Design Variables & Search Bounds (Tip & Inner Dia Fixed)");

            ImGui::PushItemWidth(130.0f * main_scale);
            ImGui::Checkbox("Optimize Length L [m]", &opt_cfg.optimize_length);
            if (opt_cfg.optimize_length) {
                ImGui::SameLine(250.0f * main_scale);
                ImGui::InputFloat("Min L", &opt_cfg.min_length, 0.1f, 0.5f, "%.2f");
                ImGui::SameLine();
                ImGui::InputFloat("Max L", &opt_cfg.max_length, 0.1f, 0.5f, "%.2f");
            }

            ImGui::Checkbox("Optimize Root Dia [m]", &opt_cfg.optimize_root_dia);
            if (opt_cfg.optimize_root_dia) {
                ImGui::SameLine(250.0f * main_scale);
                ImGui::InputFloat("Min D_root", &opt_cfg.min_root_dia, 0.02f, 0.1f, "%.2f");
                ImGui::SameLine();
                ImGui::InputFloat("Max D_root", &opt_cfg.max_root_dia, 0.02f, 0.1f, "%.2f");
            }
            ImGui::TextDisabled("Fixed by flexible pipe: Inner Dia = %.3f m, Tip Dia = %.3f m", bs.inner_dia, bs.tip_dia);

            ImGui::SeparatorText("Optimization Criteria & Loading Scope");
            ImGui::InputFloat("Cylindrical Root Ratio (L_cyl / L)", &opt_cfg.cyl_length_ratio, 0.02f, 0.05f, "%.2f");
            ImGui::InputFloat("Target Utilization Factor (0.5 - 1.0)", &opt_cfg.safety_margin, 0.02f, 0.05f, "%.2f");
            ImGui::InputFloat("Min Root/Pipe Stiffness Ratio", &opt_cfg.min_root_stiffness_ratio, 1.0f, 5.0f, "%.1f");
            ImGui::PopItemWidth();

            ImGui::Checkbox("Use Far-Field Tension Elastica Loading when Tension > 0 kN", &opt_cfg.force_far_field_when_tension);
            ImGui::Checkbox("Size Across Cold (4 C) / Nominal (20 C) / Warm (35 C) Temperature Envelope", &opt_cfg.include_temp_envelope);
            if (data1.count > 0) {
                ImGui::Checkbox("Optimize Across All Pasted Angle/Tension Pairs Envelope", &opt_cfg.use_all_pairs);
            } else {
                ImGui::TextDisabled("(Paste Angle/Tension pairs in the Pairs Table to enable multi-pair envelope optimization)");
            }

            if (ImGui::Button("Run Optimization & Apply Dimensions")) {
                auto t_opt0 = std::chrono::high_resolution_clock::now();
                last_opt_res = optimize_bend_stiffener(bs, material, opt_cfg, target_angle_deg, target_tension_kN, &data1, &data2);
                auto t_opt1 = std::chrono::high_resolution_clock::now();
                last_opt_ms = std::chrono::duration<double, std::milli>(t_opt1 - t_opt0).count();

                bs = last_opt_res.best_bs;
                if (opt_cfg.force_far_field_when_tension && std::fabs(target_tension_kN) > 0.01f) {
                    bs.use_far_field_bc = true;
                }
                target_angle = (3.14159265358979323846f / 180.0f) * target_angle_deg;
                solved = solve(target_angle, bs, results, strain, material, target_tension_kN,
                               &x_coords, &curvature, &bend_radius, &last_solve_ms, true);
                extract_results = solved;
                if (solved) {
                    update_post_solve_data();
                }
                if (opt_cfg.use_all_pairs && data1.count > 0) {
                    solve_pairs(strain_maxs, pair_min_mbrs, pair_root_moments_kNm,
                                data1, data2, bs, results, strain, material, &last_pairs_ms);
                }
            }
            ImGui::SameLine();
            if (ImGui::Button("Close")) {
                show_optimizer_window = false;
            }

            if (last_opt_res.converged) {
                ImGui::SeparatorText("Latest Optimization Result & 2D Design Space (L vs D_root)");
                ImGui::TextColored(
                    last_opt_res.feasible ? ImVec4(0.35f, 0.95f, 0.55f, 1.0f) : ImVec4(1.0f, 0.55f, 0.35f, 1.0f),
                    "Status: %s (%d candidate evaluations in %.2f ms)",
                    last_opt_res.feasible ? "FEASIBLE OPTIMUM FOUND" : "CONSTRAINTS ACTIVE AT BOUNDS (Increase Max Bounds)",
                    last_opt_res.evaluations, last_opt_ms);
                ImGui::Text("Optimal Dimensions: L = %.3f m (L_cyl = %.3f m) | D_root = %.3f m | D_tip = %.3f m",
                            last_opt_res.best_bs.length, last_opt_res.best_bs.root_cyl_length,
                            last_opt_res.best_bs.root_dia, last_opt_res.best_bs.tip_dia);
                ImGui::Text("Polyurethane Volume: %.4f m3 | Estimated PU Mass: %.1f kg | Bond Stress: %.2f MPa (%.1f%%)",
                            last_opt_res.volume_m3, last_opt_res.mass_kg,
                            last_opt_res.max_bond_stress_mpa, last_opt_res.bond_utilization * 100.0f);
                ImGui::Text("Peak Strain: %.2f%% (Utilization: %.1f%%) | Min MBR: %.2f m (Utilization: %.1f%%)",
                            last_opt_res.max_strain * 100.0f, last_opt_res.strain_utilization * 100.0f,
                            last_opt_res.min_mbr_m, last_opt_res.mbr_utilization * 100.0f);

                // 2D Optimization Design-Space Feasibility Map (L vs D_root)
                if (last_opt_res.grid_nL > 1 && last_opt_res.grid_nDr > 1) {
                    float feas_L[OPT_GRID_N * OPT_GRID_N];
                    float feas_Dr[OPT_GRID_N * OPT_GRID_N];
                    int n_feas = 0;
                    float infeas_L[OPT_GRID_N * OPT_GRID_N];
                    float infeas_Dr[OPT_GRID_N * OPT_GRID_N];
                    int n_infeas = 0;

                    for (int iL = 0; iL < last_opt_res.grid_nL; ++iL) {
                        for (int iDr = 0; iDr < last_opt_res.grid_nDr; ++iDr) {
                            if (last_opt_res.grid_feasible[iL][iDr]) {
                                feas_L[n_feas] = last_opt_res.grid_L[iL];
                                feas_Dr[n_feas] = last_opt_res.grid_Dr[iDr];
                                n_feas++;
                            } else {
                                infeas_L[n_infeas] = last_opt_res.grid_L[iL];
                                infeas_Dr[n_infeas] = last_opt_res.grid_Dr[iDr];
                                n_infeas++;
                            }
                        }
                    }

                    if (ImPlot::BeginPlot("2D Design-Space Feasibility Map (Length L vs. Root Diameter D_root)", ImVec2(-1, 215 * main_scale))) {
                        ImPlot::SetupAxes("Length L [m]", "Root Diameter D_root [m]");
                        if (n_infeas > 0) {
                            ImPlotSpec sp_inf;
                            sp_inf.Marker = ImPlotMarker_Cross;
                            sp_inf.MarkerSize = 4.0f * main_scale;
                            sp_inf.LineWeight = 1.5f * main_scale;
                            ImPlot::PlotScatter("Infeasible Grid", infeas_L, infeas_Dr, n_infeas, sp_inf);
                        }
                        if (n_feas > 0) {
                            ImPlotSpec sp_f;
                            sp_f.Marker = ImPlotMarker_Circle;
                            sp_f.MarkerSize = 4.2f * main_scale;
                            sp_f.LineWeight = 1.5f * main_scale;
                            ImPlot::PlotScatter("Feasible Grid", feas_L, feas_Dr, n_feas, sp_f);
                        }
                        if (last_opt_res.strain_bound_count > 1) {
                            ImPlotSpec sp_sb;
                            sp_sb.LineWeight = 2.0f * main_scale;
                            ImPlot::PlotLine("Strain Limit Boundary", last_opt_res.strain_bound_L, last_opt_res.strain_bound_Dr,
                                             last_opt_res.strain_bound_count, sp_sb);
                        }
                        if (last_opt_res.mbr_bound_count > 1) {
                            ImPlotSpec sp_mb;
                            sp_mb.LineWeight = 2.0f * main_scale;
                            ImPlot::PlotLine("Pipe MBR Limit Boundary", last_opt_res.mbr_bound_L, last_opt_res.mbr_bound_Dr,
                                             last_opt_res.mbr_bound_count, sp_mb);
                        }
                        float opt_L = last_opt_res.best_bs.length;
                        float opt_Dr = last_opt_res.best_bs.root_dia;
                        ImPlotSpec sp_opt;
                        sp_opt.Marker = ImPlotMarker_Diamond;
                        sp_opt.MarkerSize = 7.0f * main_scale;
                        sp_opt.LineWeight = 2.4f * main_scale;
                        ImPlot::PlotScatter("Optimal Design (L*, D_root*)", &opt_L, &opt_Dr, 1, sp_opt);
                        ImPlot::EndPlot();
                    }
                }
            }

            ImGui::End();
        }

        // --- Modal / Floating Window 4: Temperature Envelope & Hang-Off Interface Reactions ---
        if (show_temp_env_window)
        {
            ImGui::SetNextWindowSize(ImVec2(720 * main_scale, 380 * main_scale), ImGuiCond_FirstUseEver);
            ImGui::Begin("Temperature Envelope & Hang-Off Flange Reaction Loads", &show_temp_env_window);

            if (!temp_env.evaluated) {
                temp_env = evaluate_temperature_envelope(bs, material, target_angle_deg, target_tension_kN);
            }

            ImGui::SeparatorText("Governing Hang-Off Flange / Steel Collar Interface Reactions");
            ImGui::Text("Governing Root Bending Moment (M0): %.2f kN.m  (governed by Cold 4 C Stiff PU)", temp_env.governing_root_moment_kNm);
            ImGui::Text("Governing Root Shear Force (V0):    %.2f kN    | Root Axial Tension (N0): %.2f kN",
                        temp_env.governing_root_shear_kN, temp_env.governing_root_axial_kN);
            ImGui::Text("Governing Peak PU Strain:           %.2f%% (Allow %.1f%%) | Governing Min MBR: %.2f m (Allow %.2f m)",
                        temp_env.governing_strain * 100.0f, bs.allowable_strain * 100.0f,
                        temp_env.governing_min_mbr_m, bs.allowable_mbr_m);

            ImGui::SeparatorText("Multi-Temperature Verification (Using Your Cold / Nominal / Hot E-Modulus Curves)");
            if (ImGui::BeginTable("TempEnvTable", 7, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                ImGui::TableSetupColumn("Temperature Condition");
                ImGui::TableSetupColumn("E0 [MPa]");
                ImGui::TableSetupColumn("Max Strain [%]");
                ImGui::TableSetupColumn("Min MBR [m]");
                ImGui::TableSetupColumn("Root M0 [kN.m]");
                ImGui::TableSetupColumn("Root V0 [kN]");
                ImGui::TableSetupColumn("Status");
                ImGui::TableHeadersRow();

                for (int i = 0; i < 3; ++i) {
                    const auto &r = temp_env.rows[i];
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn(); ImGui::Text("%s (%.1f C)", r.label, r.temp_c);
                    ImGui::TableNextColumn(); ImGui::Text("%.1f", r.initial_E_mpa);
                    ImGui::TableNextColumn(); ImGui::Text("%.2f%%", r.max_strain * 100.0f);
                    ImGui::TableNextColumn(); ImGui::Text("%.2f m", r.min_mbr_m);
                    ImGui::TableNextColumn(); ImGui::Text("%.2f", r.root_moment_kNm);
                    ImGui::TableNextColumn(); ImGui::Text("%.2f", r.root_shear_kN);
                    ImGui::TableNextColumn();
                    if (r.pass) {
                        ImGui::TextColored(ImVec4(0.35f, 0.95f, 0.55f, 1.0f), "PASS");
                    } else {
                        ImGui::TextColored(ImVec4(1.0f, 0.40f, 0.40f, 1.0f), "FAIL");
                    }
                }
                ImGui::EndTable();
            }

            if (ImGui::Button("Re-Evaluate Envelope")) {
                temp_env = evaluate_temperature_envelope(bs, material, target_angle_deg, target_tension_kN);
            }
            ImGui::SameLine();
            if (ImGui::Button("Copy Envelope Table (Excel)")) {
                std::string tsv = "Condition\tTemp_C\tInitialE_MPa\tMaxStrain_pct\tMinMBR_m\tRootMoment_kNm\tRootShear_kN\tStatus\n";
                char row[256];
                for (int i = 0; i < 3; ++i) {
                    const auto &r = temp_env.rows[i];
                    snprintf(row, sizeof(row), "%s\t%.1f\t%.2f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n",
                             r.label, r.temp_c, r.initial_E_mpa, r.max_strain * 100.0f,
                             r.min_mbr_m, r.root_moment_kNm, r.root_shear_kN, r.pass ? "PASS" : "FAIL");
                    tsv += row;
                }
                ImGui::SetClipboardText(tsv.c_str());
                status_message = "Copied Temperature Envelope & Reaction table to clipboard.";
            }
            ImGui::SameLine();
            if (ImGui::Button("Close")) {
                show_temp_env_window = false;
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
