# CMakeList.txt: EMEsolver 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)
 

# 如果支持，请为 MSVC 编译器启用热重载。
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("EMEsolver")



include_directories(
	E:/VisualStudioFiles/vcpkg/vcpkg/installed/x64-windows/include
	E:/VisualStudioFiles/vcpkg/ku/SuperLU/SRC
	E:/VisualStudioFiles/vcpkg/ku/armadillo-14.4.1/include
)
 


# 将源代码添加到此项目的可执行文件。
add_executable (EMEsolver 
src/main.cpp
    "include/sconnectright.h" "src/sconnectright.cpp" "include/emeSmatrix.h" "src/emeSmatrix.cpp" "include/emeTmaxtrix.h" "src/emeTmaxtrix.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET EMEsolver PROPERTY CXX_STANDARD 20)
endif()


target_link_libraries(${PROJECT_NAME}  
E:/VisualStudioFiles/vcpkg/ku/armadillo-14.4.1/lib_win64/libopenblas.lib
E:/VisualStudioFiles/vcpkg/ku/SuperLU/x64/Release/SuperLU.lib

E:/VisualStudioFiles/vcpkg/vcpkg/installed/x64-windows/lib/aec.lib
E:/VisualStudioFiles/vcpkg/vcpkg/installed/x64-windows/lib/hdf5.lib
E:/VisualStudioFiles/vcpkg/vcpkg/installed/x64-windows/lib/hdf5_hl.lib
E:/VisualStudioFiles/vcpkg/vcpkg/installed/x64-windows/lib/szip.lib
E:/VisualStudioFiles/vcpkg/vcpkg/installed/x64-windows/lib/zlib.lib
)
 



# TODO: 如有需要，请添加测试并安装目标。
