@set OUT_DIR=Release
@set OUT_EXE=bs_engine
@set SDL3_DIR=C:\Users\joao\source\SDL3-3.4.2
@set INCLUDES=/I.\imgui /I.\imgui\backends /I%SDL3_DIR%\include
@set SOURCES=main.cpp .\imgui\backends\imgui_impl_sdl3.cpp .\imgui\backends\imgui_impl_sdlrenderer3.cpp .\imgui\imgui*.cpp .\imgui\implot*.cpp
@set LIBS=C:\Users\joao\source\SDL3-3.4.2\lib\x64\SDL3.lib
mkdir %OUT_DIR%
@REM cl /nologo /Zi /MD /utf-8 %INCLUDES% %SOURCES% /Fe%OUT_DIR%/%OUT_EXE%.exe /Fo%OUT_DIR%/ /link %LIBS% /subsystem:console
cl /O2 /GL /DNDEBUG /MD /Gy /Zi /EHscD /utf-8 %INCLUDES% %SOURCES% /Fe%OUT_DIR%/%OUT_EXE%.exe /Fo%OUT_DIR%/ /link /OPT:REF /OPT:ICF /LTCG  %LIBS% /subsystem:console
