set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)

add_executable(miho main.cpp)

add_subdirectory(models)

target_link_libraries(miho
	PRIVATE
	models
	fmt::fmt
  CLI11::CLI11
)

