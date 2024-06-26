set(NativeArch "")
set(DefaultArch "")
set(SandyBridgeArch "")
set(HaswellArch "")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	set(NativeArch "-march=native")
	set(DefaultArch "-mtune=generic")
	set(SandyBridgeArch "-march=corei7-avx")
	set(HaswellArch "-march=core-avx2 -mfma") # -ffast-math is required to enable FMA
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(NativeArch "-march=native")
	set(DefaultArch "-mtune=generic")
	if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9")
		set(SandyBridgeArch "-march=corei7-avx")
		set(HaswellArch "-march=core-avx2")
	else(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9")
		set(SandyBridgeArch "-march=sandybridge")
		set(HaswellArch "-march=haswell")
	endif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	set(NativeArch "-march=native")
	set(DefaultArch "-mkl=sequential -xSSE4.2")
	set(SandyBridgeArch "-mkl=sequential -xAVX ")
	set(HaswellArch "-mkl=sequential -xCORE-AVX2 -mtune=core-avx2")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	# using Visual Studio C++
	# not supported yet
endif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

macro(set_cxx_flags_from_arch SIMD_ARCH)
	if("${SIMD_ARCH}" STREQUAL "NATIVE")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NativeArch}")
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${NativeArch}")
	elseif("${SIMD_ARCH}" STREQUAL "AVX2")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${HaswellArch}")
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${HaswellArch}")
	elseif("${SIMD_ARCH}" STREQUAL "AVX")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SandyBridgeArch}")
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SandyBridgeArch}")
	elseif("${SIMD_ARCH}" STREQUAL "SSE4")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DefaultArch}")
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${DefaultArch}")
	elseif(SIMD_ARCH STREQUAL "ARM")
		set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-register -fPIC -pipe -march=armv8-a")
		set(CMAKE_C_FLAGS "-Wall -Wno-deprecated-register -fPIC -pipe -march=armv8-a")
	endif("${SIMD_ARCH}" STREQUAL "NATIVE")
endmacro(set_cxx_flags_from_arch)