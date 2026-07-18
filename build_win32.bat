<<<<<<< HEAD
call "C:/Program Files/Microsoft Visual Studio/18/Community/VC/Auxiliary/Build/vcvars64.bat"
@set OUT_DIR=Debug
@set OUT_EXE=bs_engine
@set SDL3_DIR=C:\Users\joao\source\SDL3-3.4.2
@set INCLUDES=/I.\imgui /I.\imgui\backends /I%SDL3_DIR%\include
@set SOURCES=main.cpp .\imgui\backends\imgui_impl_sdl3.cpp .\imgui\backends\imgui_impl_sdlrenderer3.cpp .\imgui\imgui*.cpp .\imgui\implot*.cpp
@set LIBS=C:\Users\joao\source\SDL3-3.4.2\lib\x64\SDL3.lib
mkdir %OUT_DIR%
cl /nologo /EHsc /Zi /MD /utf-8 %INCLUDES% %SOURCES% /Fe%OUT_DIR%/%OUT_EXE%.exe /Fo%OUT_DIR%/ /link %LIBS% /subsystem:console
=======
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
@set OUT_DIR=Debug
@set OUT_EXE=bs_engine
@set SDL3_DIR=C:\source\personal\SDL3-3.4.8
@set INCLUDES=/I.\imgui /I.\imgui\backends /I%SDL3_DIR%\include
@set SOURCES=main.cpp .\imgui\backends\imgui_impl_sdl3.cpp .\imgui\backends\imgui_impl_sdlrenderer3.cpp .\imgui\imgui*.cpp .\imgui\implot*.cpp
@set LIBS=C:\source\personal\SDL3-3.4.8\lib\x64\SDL3.lib
mkdir %OUT_DIR%
cl /nologo /EHsc /Zi /MD /utf-8 %INCLUDES% %SOURCES% /Fe%OUT_DIR%/%OUT_EXE%.exe /Fo%OUT_DIR%/ /link %LIBS% /subsystem:console
>>>>>>> 21fef49406b1087c76701275d2d46200e10f971e
